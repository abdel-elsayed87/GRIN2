
#' Chromosome Lesion Plot
#'
#' @description
#' Function plot lesion data of a selected lesion type or all lesion groups located on a certain region such as a chromosome band or the whole specified chromosome .
#'
#' @param grin.res GRIN results (output of the grin.stats function).
#' @param genome Function support either "hg19" or "hg38" genome assemblies based on the genome assembly that has been used to prepare the lesion data.
#' @param lsn.clrs Lesion colors for the regional gene plot (If not provided by the user, colors will be automatically assigned using default.grin.colors function).
#' @param chrom chromosome number
#' @param plot.start Start position of the locus of interest.
#' @param plot.end End position of the locus of interest.
#' @param lesion.grp Lesion type of interest (if not specified the plot will return all lesion groups)
#' @param spec.lsn.clr Assigned color for the lesion type of interest (should be specified when lesion.grp is specified).
#' @param expand Controls ratio of the specified locus (start and end position) to the whole plot with default value = 0.0005 (setting expand=0 will only plot the specified locus from the start to the end position without any of the upstream or downstream regions.
#' @param hg38.transcripts transcripts data retrieved from annotation hub for hg38 version 110 (should be only specified if genome="hg38").
#' @param hg19.cytoband hg19 chromosome bands start and end data in base pair (should be only specified if genome="hg19").
#' @param hg38.cytoband hg38 chromosome bands start and end data in base pair (should be only specified if genome="hg38").
#'
#' @details
#' Function will return a plot with either lesions of a certain lesion type specified by the user or all lesions of different lesion groups allocated to either a specific part of the chromosome such as a chromosome band or the whole length of the chromosme.As compared the the grin.gene.plot function, chrom.lsn.plot will not return a transcripts track which will allow the user to specify a large region or even the whole chromosome. A chromosome ideogram is added on the top of the lesion plot with the specified part of the chromosome encircled in red.
#'
#' @return
#' Function will return a plot of a pre-specfied lesion type or all lesion affecting a specific locus or the whole chromosome.
#'
#' @export
#'
#' @importFrom graphics rect legend text segments
#' @importFrom gridGraphics grid.echo
#' @importFrom grid grid.grab editGrob grid.draw popViewport grid.newpage pushViewport viewport gpar
#' @importFrom Gviz plotTracks IdeogramTrack
#' @importFrom ensembldb getGeneRegionTrackForGviz
#' @importFrom EnsDb.Hsapiens.v75 EnsDb.Hsapiens.v75
#' @importFrom GenomeInfoDb seqlevelsStyle
#'
#' @references
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' @author {Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org} and Stanley Pounds \email{stanley.pounds@stjude.org}}
#'
#' @seealso [grin.stats()]
#'
#' @examples
#' \dontrun{
#' data(lesion.data)
#' data(hg19.gene.annotation)
#' data(hg19.chrom.size)
#' data(hg19_cytoband)
#' data(hg38_cytoband)
#'
#' # run GRIN analysis using grin.stats function
#' grin.results=grin.stats(lesion.data,
#'                         genome.version="Human_GRCh37")
#'
#' # Plots Showing Different Types of Lesions Affecting a region of Interest without plotting the
#' # transcripts track (this will allow plotting a larger locus of the chromosome such as a
#' # chromosome band using hg19 genome assembly):
#' locus.plot=chrom.lsn.plot(grin.results, genome="hg19", hg19.cytoband=hg19_cytoband,
#'                           chrom=9, plot.start=19900000, plot.end=25600000,
#'                           lesion.grp = "loss", spec.lsn.clr = "blue")
#'
#' # Plots Showing Different Types of Lesions Affecting the whole chromosome:
#' chrom.lsn=chrom.lsn.plot(GRIN.results, genome="hg19", hg19.cytoband=hg19_cytoband,
#'                           chrom=9, plot.start=1, plot.end=141000000)
#'
#' # for hg38 genome assembly:
#' ah <- AnnotationHub()
#' gtf.V110 <- ah[["AH113665"]]
#'
#' # run GRIN analysis using grin.stats function
#' grin.results=grin.stats(lesion.data,
#'                        genome.version="Human_GRCh38")
#' # to plot deletions on the p21.3 band of chromosome 9:
#' chrom.lsn.plot(grin.results, genome="hg38", hg38.transcripts="gtf.v110",
#'                hg38.cytoband=hg38_cytoband, chrom=9, plot.start=19900000, plot.end=25600000,
#'                lesion.grp = "loss", spec.lsn.clr = "blue")
#'
#' # to plot all lesions on chromosome 9:
#' chrom.lsn.plot(grin.results, genome="hg38", hg38.transcripts="gtf.v110",
#'                hg38.cytoband=hg38_cytoband, chrom=9, plot.start=1, plot.end=141000000)
#' }

chrom.lsn.plot=function(grin.res,          # GRIN results (output of the grin.stats function)
                        genome,            # genome assembly (hg19 or hg38) gene=NULL,         # gene name (should be only specified in the regional gene plots)
                        lsn.clrs=NULL,     # Specified colors per lesion types (gene plots when gene name is specified). If not specified, colors will be automatically assigned using default.grin.colors function
                        chrom=NULL,        # chromosome number (should be only specified in the locus plots where plot.start and plot.end for the locus of interest are specified)
                        plot.start=NULL,   # start position of the locus of interest
                        plot.end=NULL,     # end position of the locus of interest
                        lesion.grp=NULL,   # lesion group of interest (specified to plot just one type of lesions)
                        spec.lsn.clr=NULL, # color of the lesion of interest (should be specified when lesion.grp is specified)
                        expand=0.0005,     # Controls ratio of the gene locus (start and end position) to the whole plot with default value = 0.3 (setting expand=0 will only plot the gene locus from the start to the end position without any of the upstream or downstream regions of the gene)
                        hg38.transcripts=NULL, # transcripts data retrieved from annotation hub for hg38 version 110 (should be only specified if genome="hg38").
                        hg19.cytoband=NULL,    # hg19 chromosome bands start and end data in base pair (should be only specified if genome="hg19").
                        hg38.cytoband=NULL)    # hg38 chromosome bands start and end data in base pair (should be only specified if genome="hg38").

{
  if (is.character(lesion.grp)) {
    # To specify a locus of interest for the plotting purpose
    locus.chr=chrom
    locus.start=plot.start
    locus.end=plot.end
    locus.size=(locus.end-locus.start)+1

    # Find lesions on the specified region of the chromosome
    lesion.data=grin.res[["lsn.data"]]
    lsn.type=lesion.grp
    lsn.data=lesion.data[lesion.data$lsn.type==lsn.type,]
    lsn.data=lsn.data[lsn.data$chrom==locus.chr,]
    lsn.dset=order.index.lsn.data(lsn.data)
    lsn.clr=spec.lsn.clr

    lsn.index=lsn.dset$lsn.index
    lsn.ind.mtch=which(lsn.index$chrom==locus.chr)
    if (length(lsn.ind.mtch)==0)
      stop(paste0("No lesions overlap"))

    lsn.chr.rows=NULL
    for (i in lsn.ind.mtch)
    {
      blk.rows=(lsn.index$row.start[i]:lsn.index$row.end[i])
      lsn.chr.rows=c(lsn.chr.rows,blk.rows)
    }
    lsn.chr.rows=unlist(lsn.chr.rows)

    lsn.chr.data=lsn.data[lsn.chr.rows,]
    lsn.chr.data=lsn.chr.data[lsn.chr.data$chrom==locus.chr,]

    if (any(lsn.chr.data$chrom!=locus.chr))
      stop(paste0("Error in finding lesions on same chromosome"))

    # Find lesions that overlap the locus
    ov.rows=which((lsn.chr.data$loc.start<=locus.end)&(lsn.chr.data$loc.end>=locus.start))
    if (length(ov.rows)==0)
      stop(paste0("No lesions overlap locus"))

    lsn.locus=lsn.chr.data[ov.rows,]

    lsn.locus$size=lsn.locus$loc.end-lsn.locus$loc.start+1
    lsn.locus=lsn.locus[order(lsn.locus$lsn.type, lsn.locus$size), ]

    # define plotting data
    x.start=locus.start-expand*locus.size
    x.end=locus.end+expand*locus.size

    lsn.locus$index=1:nrow(lsn.locus)
    lsn.locus$subj.num=as.numeric(as.factor(lsn.locus$index))
    lsn.locus$lsn.clr=lsn.clr
    lsn.locus$type.num=as.numeric(as.factor(lsn.locus$lsn.type))
    n.type=max(lsn.locus$type.num)
    n.subj=max(lsn.locus$subj.num)

    lsn.locus$y0=-lsn.locus$subj.num+(lsn.locus$type.num-1)/n.type
    lsn.locus$y1=-lsn.locus$subj.num+lsn.locus$type.num/n.type

    lsn.locus$x0=pmax(lsn.locus$loc.start,x.start)
    lsn.locus$x1=pmin(lsn.locus$loc.end,x.end)

    plot(c(x.start-0.02*(x.end-x.start),x.end+0.1*(x.end-x.start)),
         c(+0.1,-1.1)*n.subj,type="n",
         main="",
         xlab="",
         ylab="",axes=F)

    graphics::rect(x.start,
         -(1:n.subj),
         x.end,
         -(1:n.subj)+1,
         col=c("snow","gainsboro")[1+(1:n.subj)%%2],
         border=NA)

    graphics::segments(c(locus.start,locus.end),
             rep(-n.subj,2),
             c(locus.start,locus.end),
             rep(0,2),
             col="darkgray")

    graphics::rect(lsn.locus$x0,
         lsn.locus$y0,
         lsn.locus$x1,
         lsn.locus$y1,
         col=lsn.locus$lsn.clr,
         border=lsn.locus$lsn.clr)

    graphics::segments(c(locus.start,locus.end),
             rep(-n.subj,2),
             c(locus.start,locus.end),
             rep(0,2),
             col="darkgray",lty=2)

    graphics::text(c(locus.start,locus.end),
         -n.subj,
         c(locus.start,locus.end),
         pos=1, cex=0.85)
    graphics::text((locus.start+locus.end)/2,
         0,paste0(lsn.type, " _ ", "chr",locus.chr, ": ", locus.start, " - ", locus.end), pos=3,cex=1)

    gridGraphics::grid.echo()
    lesion.plt <- grid::grid.grab()
    lesion.plt <- grid::editGrob(lesion.plt, gp=grid::gpar(fontsize=12))

    ## add the transcripts track for all genes on the specified region  based on the genome assembly
    if (genome=="hg19")
    {
      edb <- EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
      GenomeInfoDb::seqlevelsStyle(edb) <- "UCSC"
      txs <- ensembldb::getGeneRegionTrackForGviz(edb, chromosome=locus.chr, start=plot.start,
                                       end=plot.end)

      ## Define the individual tracks:
      ## - Ideogram
      hg19_cytoband=hg19.cytoband
      ideo_track <- Gviz::IdeogramTrack(genome = "hg19", chromosome = locus.chr,
                                        bands = hg19_cytoband, from=plot.start, to=plot.end)

      grid::grid.newpage()

      # specify plotting regions for the locus lesions and transcript tracks
      grid::pushViewport(grid::viewport(height=1, width=1.25, y=0, just="bottom"))
      grid::grid.draw(lesion.plt)
      grid::popViewport(1)

      grid::pushViewport(grid::viewport(height=0.15, width=0.88, y=1, just="top"))
      Gviz::plotTracks(ideo_track, from = plot.start, to = plot.end, add=TRUE, showId = FALSE,
                 showBandId = TRUE, cex.bands = 0.4)
      grid::popViewport(1)

    }

    if (genome=="hg38")
    {
      trans.hg38=hg38.transcripts
      edb <- trans.hg38
      GenomeInfoDb::seqlevelsStyle(edb) <- "UCSC"
      txs <- ensembldb::getGeneRegionTrackForGviz(edb, chromosome=locus.chr, start=plot.start,
                                                   end=plot.end)

      ## Define the individual tracks:
      ## - Ideogram
      hg38_cytoband=hg38.cytoband
      ideo_track <- Gviz::IdeogramTrack(genome = "hg38", chromosome = locus.chr,
                                         bands = hg38_cytoband, from=plot.start, to=plot.end)

      grid::grid.newpage()

      grid::pushViewport(grid::viewport(height=1, width=1.25, y=0, just="bottom"))
      grid::grid.draw(lesion.plt)
      grid::popViewport(1)

      grid::pushViewport(grid::viewport(height=0.15, width=0.88, y=1, just="top"))
      Gviz::plotTracks(ideo_track, from = plot.start, to = plot.end, add=TRUE, showId = FALSE,
                        showBandId = TRUE, cex.bands = 0.4)
      grid::popViewport(1)
    }
  }
  else if (is.null(lesion.grp)) {
    # To specify a locus of interest for the plotting purpose
    locus.chr=chrom
    locus.start=plot.start
    locus.end=plot.end
    locus.size=(locus.end-locus.start)+1

    # Find lesions on the same chromosome as the gene
    lsn.dset=order.index.lsn.data(grin.res[["lsn.data"]])
    lsn.data=lsn.dset$lsn.data

    lsn.types=unique(lsn.data$lsn.type)
    if (is.null(lsn.clrs))
      lsn.clrs=default.grin.colors(lsn.types)

    lsn.index=lsn.dset$lsn.index
    lsn.ind.mtch=which(lsn.index$chrom==locus.chr)
    if (length(lsn.ind.mtch)==0)
      stop(paste0("No lesions overlap ",chrom,"."))

    lsn.chr.rows=NULL
    for (i in lsn.ind.mtch)
    {
      blk.rows=(lsn.index$row.start[i]:lsn.index$row.end[i])
      lsn.chr.rows=c(lsn.chr.rows,blk.rows)
    }
    lsn.chr.rows=unlist(lsn.chr.rows)

    lsn.chr.data=lsn.data[lsn.chr.rows,]

    if (any(lsn.chr.data$chrom!=locus.chr))
      stop(paste0("Error in finding lesions on same chromosome."))

    # Find lesions at overlap the locus
    ov.rows=which((lsn.chr.data$loc.start<=locus.end)&(lsn.chr.data$loc.end>=locus.start))
    if (length(ov.rows)==0)
      stop(paste0("No lesions overlap the chromosome."))

    lsn.locus=lsn.chr.data[ov.rows,]

    lsn.locus$size=lsn.locus$loc.end-lsn.locus$loc.start+1
    lsn.locus=lsn.locus[order(lsn.locus$lsn.type, lsn.locus$size), ]

    # define plotting data
    x.start=locus.start-expand*locus.size
    x.end=locus.end+expand*locus.size

    lsn.locus$index=1:nrow(lsn.locus)
    lsn.locus$subj.num=as.numeric(as.factor(lsn.locus$index))
    lsn.locus$lsn.clr=lsn.clrs[lsn.locus$lsn.type]
    lsn.locus$type.num=as.numeric(as.factor(lsn.locus$lsn.type))
    n.type=max(lsn.locus$type.num)
    n.subj=max(lsn.locus$subj.num)

    lsn.locus$y0=-lsn.locus$subj.num+(lsn.locus$type.num-1)/n.type
    lsn.locus$y1=-lsn.locus$subj.num+lsn.locus$type.num/n.type

    lsn.locus$x0=pmax(lsn.locus$loc.start,x.start)
    lsn.locus$x1=pmin(lsn.locus$loc.end,x.end)

    plot(c(x.start-0.02*(x.end-x.start),x.end+0.1*(x.end-x.start)),
         c(+0.1,-1.15)*n.subj,type="n",
         main="",
         xlab="",
         ylab="",axes=F)

    graphics::rect(x.start,
         -(1:n.subj),
         x.end,
         -(1:n.subj)+1,
         col=c("snow","gainsboro")[1+(1:n.subj)%%2],
         border=NA)

    graphics::segments(c(locus.start,locus.end),
             rep(-n.subj,2),
             c(locus.start,locus.end),
             rep(0,2),
             col="darkgray")

    graphics::rect(lsn.locus$x0,
         lsn.locus$y0,
         lsn.locus$x1,
         lsn.locus$y1,
         col=lsn.locus$lsn.clr,
         border=lsn.locus$lsn.clr)

    graphics::segments(c(locus.start,locus.end),
             rep(-n.subj,2),
             c(locus.start,locus.end),
             rep(0,2),
             col="darkgray",lty=2)

    graphics::text(c(locus.start,locus.end),
         -n.subj,
         c(locus.start,locus.end),
         pos=1, cex=0.85)
    graphics::text((locus.start+locus.end)/2,
         0,paste0("chr",locus.chr, ": ", locus.start, " - ", locus.end), pos=3,cex=1)

    lgd=graphics::legend((x.start+x.end)/2,-1.1*n.subj,
               fill=lsn.clrs,
               legend=names(lsn.clrs),
               ncol=length(lsn.clrs),
               xjust=0.45,border=NA,
               cex=0.75,bty="n")

    gridGraphics::grid.echo()
    lesion.plt <- grid::grid.grab()
    lesion.plt <- grid::editGrob(lesion.plt, gp=grid::gpar(fontsize=12))

    ## add the transcripts track for all genes on the specified region  based on the genome assembly
    if (genome=="hg19")
    {
      edb <- EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
      GenomeInfoDb::seqlevelsStyle(edb) <- "UCSC"
      txs <- ensembldb::getGeneRegionTrackForGviz(edb, chromosome=locus.chr, start=plot.start,
                                                  end=plot.end)

      ## Define the individual tracks:
      ## - Ideogram
      hg19_cytoband=hg19.cytoband
      ideo_track <- Gviz::IdeogramTrack(genome = "hg19", chromosome = locus.chr,
                                          bands = hg19_cytoband, from=plot.start, to=plot.end)

      grid::grid.newpage()

      # specify plotting regions for the locus lesions and transcript tracks
      grid::pushViewport(grid::viewport(height=1, width=1.25, y=0, just="bottom"))
      grid::grid.draw(lesion.plt)
      grid::popViewport(1)

      grid::pushViewport(grid::viewport(height=0.15, width=0.88, y=1, just="top"))
      Gviz::plotTracks(ideo_track, from = plot.start, to = plot.end, add=TRUE, showId = FALSE,
                 showBandId = TRUE, cex.bands = 0.4)
      grid::popViewport(1)

    }

    if (genome=="hg38")
    {
      trans.hg38=hg38.transcripts
      edb <- trans.hg38
      GenomeInfoDb::seqlevelsStyle(edb) <- "UCSC"
      txs <- ensembldb::getGeneRegionTrackForGviz(edb, chromosome=locus.chr, start=plot.start,
                                                      end=plot.end)

      ## Define the individual tracks:
      ## - Ideogram
      hg38_cytoband=hg38.cytoband
      ideo_track <- Gviz::IdeogramTrack(genome = "hg38", chromosome = locus.chr,
                                         bands = hg38_cytoband, from=plot.start, to=plot.end)

      grid::grid.newpage()

      grid::pushViewport(grid::viewport(height=1, width=1.25, y=0, just="bottom"))
      grid::grid.draw(lesion.plt)
      grid::popViewport(1)

      grid::pushViewport(grid::viewport(height=0.15, width=0.88, y=1, just="top"))
      Gviz::plotTracks(ideo_track, from = plot.start, to = plot.end, add=TRUE, showId = FALSE,
                         showBandId = TRUE, cex.bands = 0.4)
      grid::popViewport(1)
    }
  }
}
