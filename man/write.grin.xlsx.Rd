% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write.grin.xlsx.R
\name{write.grin.xlsx}
\alias{write.grin.xlsx}
\title{Write GRIN Results}
\usage{
write.grin.xlsx(grin.result, output.file)
}
\arguments{
\item{grin.result}{output results of the grin.stats function.}

\item{output.file}{output file name ".xlsx".}
}
\value{
This function return an excel file with seven sheets that include:
\item{gene.hits}{data table of GRIN results that include gene annotation, number of subjects affected by each lesion type for example gain, loss, mutation, etc.., and number of hits affecting each locus. The GRIN results table will also include P and FDR adjusted q-values showing the probability of each locus of being affected by one or a constellation of multiple types of lesions.}
\item{gene.lsn.data}{each row represent a gene overlapped by a certain lesion. Column "gene" shows the overlapped gene ensembl ID and "ID"" column has the patient ID.}
\item{lsn.data}{input lesion data}
\item{gene.data}{input gene annotation data}
\item{chr.size}{data table showing the size of the 22 autosomes, in addition to X and Y chromosomes in base pairs.}
\item{interpretation}{provides some details about  the content of each sheet in the output excel file and interpretation of each column in the "gene.hits" GRIN results table.}
\item{method.paragraph}{include a paragraph that explains the GRIN model and cite some references.}
}
\description{
The function Write GRIN results to an excel file with multiple sheets that include GRIN results, lesion data, gene annotation data, chromosome size, gene-lesion overlap and methods paragraph.
}
\examples{
data(lesion.data)
data(hg19.gene.annotation)
data(hg19.chrom.size)

# to directly retreive gene annotation and chromosome size files from Ensembl BioMart database,
# UCSC genome browsers and run the GRIN analysis:
grin.results=grin.stats(lesion.data,
                         hg19.gene.annotation,
                         hg19.chrom.size)

# Write GRIN results in to an excel sheet ".xlsx" using write.grin.xlsx function.
}
\references{
Pounds, Stan, et al. (2013) A genomic random interval model for statistical analysis of genomic lesion data.

Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
}
\seealso{
\code{\link[=grin.stats]{grin.stats()}}
}
\author{
{Stanley Pounds \email{stanley.pounds@stjude.org}}
}
