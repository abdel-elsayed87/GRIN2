
#' GRIN Statistics Lesions Plot
#'
#' @description
#' Generates a plot showing all lesion types that span a specified gene or regulatory feature with GRIN stats added.
#'
#' @param grin.res GRIN results for genes or regulatory elements, as returned by the \code{\link{grin.stats}} function.
#' @param feature Ensembl ID of a feature of interest. This can be either a gene (e.g., Ensembl gene ID) or a regulatory region from Ensembl Regulatory Build or FANTOM5 project.
#' @param lsn.clrs Named vector of colors for each lesion type. If not provided, colors are automatically assigned using \code{\link{default.grin.colors}}.
#' @param expand Numeric value that controls the proportion of the flanking region (upstream/downstream) around the gene to be included. Default is \code{0.0005}. Set to \code{0} to restrict the plot strictly to the locus boundaries.
#'
#' @details
#' The plot consists of two panels:
#' \itemize{
#'   \item \strong{Top panel}: Displays all lesion types overlapping the selected gene or regulatory feature. Lesions are color-coded by type, as indicated in the legend.
#'   \item \strong{Bottom panel}: Summarizes GRIN statistics for the feature, including the number of subjects affected per lesion type, and the corresponding \eqn{-\log_{10}(p)} and \eqn{-\log_{10}(q)} values for significance.
#' }
#'
#' This plot is particularly useful for regulatory features, which typically lack transcript structure. Therefore, no transcript or exon structure is shown.
#'
#' @return
#' A two-panel plot showing lesion distribution and GRIN statistics for a given gene or regulatory feature, without a transcript annotation panel.
#'
#' @export
#'
#' @importFrom graphics rect legend text segments
#'
#' @references
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' @author
#' Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org} and Stanley Pounds \email{stanley.pounds@stjude.org}
#'
#' @seealso \code{\link{grin.stats}}
#'
#' @examples
#' data(lesion_data)
#' data(hg38_gene_annotation)
#' data(hg38_chrom_size)
#'
#' # Run GRIN analysis
#' grin.results <- grin.stats(lesion_data,
#'                            hg38_gene_annotation,
#'                            hg38_chrom_size)
#'
#' # Plot lesion and GRIN stats for a gene of interest (e.g., WT1)
#' grin.stats.lsn.plot(grin.results, feature = "ENSG00000184937")
#'
#' # Can also be also used for regulatory features without transcript panels

grin.stats.lsn.plot=function(grin.res,          # GRIN results for regulatory features (output of the grin.stats function)
                             feature=NULL,      # Feature ensembl ID from Ensembl regulatory build or FANTOM5 project. An ensembl ID of a gene can be provided as well.
                             lsn.clrs=NULL,     # Specified colors per lesion types. If not specified, colors will be automatically assigned using default.grin.colors function
                             expand=0.0005)     # Controls ratio of the gene locus (start and end position) to the whole plot with default value = 0.0005 (setting expand=0 will only plot the gene locus from the start to the end position without any of the upstream or downstream regions of the gene)

{
  # regional regulatory feature plot
  if (is.character(feature))
  {
    # Find the requested regulatory feature
    gene.data=grin.res[["gene.data"]]
    gene.mtch=which(gene.data[,"gene"]==feature)

    if (length(feature)!=1)
      stop("Exactly one feature must be specified!")
    if (length(gene.mtch)==0)
      stop(paste0(feature," not found in gene.data."))

    if (length(gene.mtch)>1)
      stop(paste0("Multiple matches of ",feature," found in gene data."))

    gene.chr=gene.data[gene.mtch,"chrom"]
    gene.start=gene.data[gene.mtch,"loc.start"]
    gene.end=gene.data[gene.mtch,"loc.end"]
    gene.size=(gene.end-gene.start)+1

    # Find lesions on the same chromosome as the feature
    lsn.dset=order.index.lsn.data(grin.res[["lsn.data"]])
    lsn.data=lsn.dset$lsn.data

    lsn.types=unique(lsn.data$lsn.type)
    if (is.null(lsn.clrs))
      lsn.clrs=default.grin.colors(lsn.types)

    lsn.index=lsn.dset$lsn.index
    lsn.ind.mtch=which(lsn.index$chrom==gene.chr)
    if (length(lsn.ind.mtch)==0)
      stop(paste0("No lesions overlap ",feature,"."))

    lsn.chr.rows=NULL
    for (i in lsn.ind.mtch)
    {
      blk.rows=(lsn.index$row.start[i]:lsn.index$row.end[i])
      lsn.chr.rows=c(lsn.chr.rows,blk.rows)
    }
    lsn.chr.rows=unlist(lsn.chr.rows)

    lsn.chr.data=lsn.data[lsn.chr.rows,]

    if (any(lsn.chr.data$chrom!=gene.chr))
      stop(paste0("Error in finding lesions on same chromosome as gene ",feature,"."))

    # Find lesions at overlap the regulatory feature
    ov.rows=which((lsn.chr.data$loc.start<=gene.end)&(lsn.chr.data$loc.end>=gene.start))
    if (length(ov.rows)==0)
      stop(paste0("No lesions overlap gene ",feature,"."))

    lsn.gene=lsn.chr.data[ov.rows,]

    lsn.gene$size=lsn.gene$loc.end-
      lsn.gene$loc.start+1
    lsn.gene=lsn.gene[order(lsn.gene$lsn.type, lsn.gene$size), ]

    # define plotting data
    x.start=gene.start-expand*gene.size
    x.end=gene.end+expand*gene.size

    lsn.gene$index=1:nrow(lsn.gene)
    lsn.gene$subj.num=as.numeric(as.factor(lsn.gene$index))
    lsn.gene$lsn.clr=lsn.clrs[lsn.gene$lsn.type]
    lsn.gene$type.num=as.numeric(as.factor(lsn.gene$lsn.type))


    n.type=max(lsn.gene$type.num)
    n.subj=max(lsn.gene$subj.num)

    lsn.gene$y0=-lsn.gene$subj.num+(lsn.gene$type.num-1)/n.type
    lsn.gene$y1=-lsn.gene$subj.num+lsn.gene$type.num/n.type

    lsn.gene$x0=pmax(lsn.gene$loc.start,x.start)
    lsn.gene$x1=pmin(lsn.gene$loc.end,x.end)

    plot(c(x.start-0.055*(x.end-x.start),x.end+0.13*(x.end-x.start)),
         c(+0.05,-1.75)*n.subj,type="n",
         main="",
         xlab="",
         ylab="",axes=F)

    graphics::rect(x.start,
                   -(1:n.subj),
                   x.end,
                   -(1:n.subj)+1,
                   col=c("snow","gainsboro")[1+(1:n.subj)%%2],
                   border=NA)

    graphics::segments(c(gene.start,gene.end),
                       rep(-n.subj,2),
                       c(gene.start,gene.end),
                       rep(0,2),
                       col="darkgray")

    graphics::rect(lsn.gene$x0,
                   lsn.gene$y0,
                   lsn.gene$x1,
                   lsn.gene$y1,
                   col=lsn.gene$lsn.clr,
                   border=lsn.gene$lsn.clr)

    graphics::segments(c(gene.start,gene.end),
                       rep(-n.subj,2),
                       c(gene.start,gene.end),
                       rep(0,2),
                       col="darkgray",lty=2)

    graphics::text(c(gene.start,gene.end),
                   -n.subj,
                   c(gene.start,gene.end),
                   pos=1, cex=0.75)

    graphics::text((gene.start+gene.end)/2,
                   0,feature,pos=3,cex=1)


    lgd=graphics::legend((x.start+x.end)/2,-1.10*n.subj,
                         fill=lsn.clrs,
                         legend=names(lsn.clrs),
                         ncol=length(lsn.clrs),
                         xjust=0.4,border=NA,
                         cex=0.7,bty="n")

    graphics::text(lgd$text$x[1]-0.085*diff(range(lgd$text$x)),
                   -c(1.25,1.35,1.45)*n.subj,
                   c("n","-log10p","-log10q"),pos=2, cex=0.7, font=2)

    gene.stats=grin.res[["gene.hits"]]
    stat.mtch=which(gene.stats$gene==feature)
    gene.stats=gene.stats[stat.mtch,]

    graphics::text(lgd$text$x,-1.25*n.subj,
                   gene.stats[,paste0("nsubj.",names(lsn.clrs))],
                   cex=0.7)
    graphics::text(lgd$text$x,-1.35*n.subj,
                   round(-log10(gene.stats[,paste0("p.nsubj.",names(lsn.clrs))]),2),
                   cex=0.7)
    graphics::text(lgd$text$x,-1.45*n.subj,
                   round(-log10(gene.stats[,paste0("q.nsubj.",names(lsn.clrs))]),2),
                   cex=0.7)

    lgd2=graphics::legend((x.start+x.end)/2,-1.5*n.subj,
                          legend=paste0("const", 1:length(lsn.types), ".typ"),
                          ncol=length(lsn.types),
                          xjust=0.42,border=NA,
                          cex=0.7,bty="n")

    graphics::text(lgd2$text$x[1]-0.0005*diff(range(lgd2$text$x)),
                   -c(1.65,1.75)*n.subj,
                   c("const.p", "const.q"),pos=2, cex=0.7, font=2)

    graphics::text(lgd2$text$x,-1.65*n.subj,
                   round(-log10(gene.stats[,paste0("p",1:length(lsn.types), ".nsubj")]),2),
                   cex=0.7, pos=4)
    graphics::text(lgd2$text$x,-1.75*n.subj,
                   round(-log10(gene.stats[,paste0("q",1:length(lsn.types), ".nsubj")]),2),
                   cex=0.7, pos=4)
  }
}

