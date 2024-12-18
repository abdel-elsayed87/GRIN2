% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.ensembl.annotation.R
\name{get.ensembl.annotation}
\alias{get.ensembl.annotation}
\title{Get Ensembl Gene and Regulatory Features Annotation Data}
\usage{
get.ensembl.annotation(genome.assembly)
}
\arguments{
\item{genome.assembly}{User can specify one of four genome assemblies that include "Human_GRCh38", "Human_GRCh37", "Mouse_HGCm39" and "Mouse_HGCm38".}
}
\value{
A list of three components:
\item{gene.annotation}{A 9 columns data.frame with gene annotation data that include enesembl ID, chromosome, gene start and gene end position, gene name, gene description, biotype, chromosome strand, and chromosome band.}
\item{reg.annotation.predicted}{A 5 columns data.frame with regulatory features annotation data directly retreived from Ensembl regulatory build that include enesembl ID, chromosome, description(promoter, enhancer, etc..), feature start and end positions.}
\item{reg.annotation.validated}{A 5 columns data.frame with regulatory features annotation data for experimentally verified features retreived from FANTOM5 project that include feature ID, chromosome, description(enhancer, transcription start site (TSS)), feature start and end positions.}
}
\description{
Function directly retrieve gene and regulatory features annotation data from Ensembl BioMart database based on the specified genome assembly.
}
\details{
Based on the genome assembly specified by the user, the function will directly retrieve gene and regulatory features annotation data from ensembl BioMart database. Annotation data include enesembl ID, the chromosome on which the gene is located, gene start and gene end position, gene name, gene description, biotype, chromosome strand and chromosome band. Gene classes (biotypes) include protein coding genes, long noncoding RNAs (lncRNAs), microRNAs (miRNAs), small nuclear RNAs (snRNA), small nucleolar RNAs (snoRNA), immunoglobulins (IGs), T-cell receptors (TCRs) and pseudogens. Regulatory features data retrieved from Ensembl regulatory build are categorized in 6 classes that include promoters, promoter flanking regions, predicted enhancers, CTCF binding sites, transcription factor (TF) binding sites and the open chromatin regions.Ensembl first imports publicly available data from different large epigenomic consortia such as ENCODE, Roadmap Epigenomics and Blueprint. All high-throughput sequencing data sets are then uniformly processed using the Ensembl Regulation Sequence Analysis (ERSA) pipeline to generate signal tracks for enriched regions also referred to as annotated features or peaks. Segmentation data provide information about promoter, promoter flanking regions, enhancers and CTCF binding sites. If TF binding probability is >0 in areas outside previously mentioned regions, it takes a label of TF binding site. If any open chromatin region did not overlap with the above features, it takes a label of unannotated open chromatin. users will also have the chance to use a list of experimentally verified enhancers/transcription start sites (TSS) using the CAGE (Cap Analysis of Gene Expression) experiment on a multitude of different primary cells and tissues from the Functional Annotation of the Mouse/Mammalian Genome (FANTOM5) project.
}
\examples{
\donttest{
# Retrieve annotation data for human hg19 genome assembly:
hg19.ann=get.ensembl.annotation("Human_GRCh37")
# gene annotation data:
hg19.gene.annotation=hg19.ann$gene.annotation
# regulatory features annotation data retrieved from Ensembl regulatory build:
hg19.reg.annotation=hg19.ann$reg.annotation.predicted
# regulatory features annotation data retrieved from FANTOM5 project:
hg19.fantom.annotation=hg19.ann$reg.annotation.validated

# "Human_GRCh38" can be used instead of "Human_GRCh37" to retrieve gene and regulatory features
# annotation data for human hg38 genome assembly.
}
}
\references{
Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.

Zerbino, Daniel R., et al. (2015). The ensembl regulatory build.

Kinsella, Rhoda J., et al. (2011). Ensembl BioMarts: a hub for data retrieval across taxonomic space.
}
\seealso{
\code{\link[biomaRt:useEnsembl]{biomaRt::useEnsembl()}}, \code{\link[biomaRt:getBM]{biomaRt::getBM()}}
}
\author{
{Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org} and Stanley Pounds \email{stanley.pounds@stjude.org}}
}
