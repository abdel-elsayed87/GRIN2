% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prob.hits.R
\name{prob.hits}
\alias{prob.hits}
\title{Find Probablity of Locus Hit}
\usage{
prob.hits(hit.cnt, chr.size = NULL)
}
\arguments{
\item{hit.cnt}{output results of the count.hits function with number of subjects and number of hits affecting each locus.}

\item{chr.size}{data.frame with the size of the 22 autosomes, in addition to X and Y chromosomes in base pairs. The data.frame should has two columns "chrom" with the chromosome number and "size" for the size of the chromosome in base pairs.}
}
\value{
A list with the following components:
\item{gene.hits}{data table of GRIN results that include gene annotation, number of subjects affected by each lesion type for example gain, loss, mutation, etc.., and number of hits affecting each locus. The GRIN results table will also include P and FDR adjusted q-values showing the probability of each locus of being affected by one or a constellation of multiple types of lesions.}
\item{lsn.data}{input lesion data}
\item{gene.data}{input gene annotation data}
\item{gene.lsn.data}{each row represent a gene overlapped by a certain lesion. Column "gene" shows the overlapped gene ensembl ID and "ID"" column has the patient ID.}
\item{chr.size}{data table showing the size of the 22 autosomes, in addition to X and Y chromosomes in base pairs.}
\item{gene.index}{data.frame with overlapped gene-lesion data rows that belong to each chromosome in the gene.lsn.data table.}
\item{lsn.index}{data.frame that shows the overlapped gene-lesion data rows taht belong to each lesion in the gene.lsn.data table.}
}
\description{
The function evaluates the probability of a locus to be affected by one or a constellation of multiple types of lesions.
}
\details{
The function computes p-value for the probability of each locus (gene or regulatory feature) to be affected by different types of lesions based on a convolution of independent but non-identical Bernoulli distributions to determine whether a certain locus has an abundance of lesions that is statistically significant.In addition, FDR-adjusted q value is computed for each locus based on Pounds & Cheng (2006) estimator of the proportion of tests with a true null (pi.hat). The function also evaluates if a certain locus is affected by a constellation of multiple types of lesions and computes a p and adjusted q values for the locus to be affected by one type of lesions (p1), two types of lesions (p2), etc...
}
\examples{
data(lesion.data)
data(hg19.gene.annotation)
data(hg19.chrom.size)

# prepare gene and lesion data for later computations:
prep.gene.lsn=prep.gene.lsn.data(lesion.data,
                                 hg19.gene.annotation)

# determine lesions that overlap each gene (locus):
gene.lsn.overlap=find.gene.lsn.overlaps(prep.gene.lsn)

# count number of subjects affected by different types of lesions and number of hits that affect
# each locus:
count.subj.hits=count.hits(gene.lsn.overlap)

# compute the probability of each locus to be affected by one or a constellation of multiple
# types of lesion
hits.prob=prob.hits(count.subj.hits, hg19.chrom.size)
}
\references{
Pounds, Stan, et al. (2013) A genomic random interval model for statistical analysis of genomic lesion data.

Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
}
\seealso{
\code{\link[=prep.gene.lsn.data]{prep.gene.lsn.data()}}, \code{\link[=find.gene.lsn.overlaps]{find.gene.lsn.overlaps()}}, \code{\link[=count.hits]{count.hits()}}
}
\author{
{Stanley Pounds \email{stanley.pounds@stjude.org}}
}
