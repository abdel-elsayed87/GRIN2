
#' Example T-ALL Lesion Dataset
#'
#' Genomic lesion dataset including copy number variations, single nucleotide variants, and structural rearrangements affecting 265 newly diagnosed T-cell Acute Lymphoblastic Leukemia (T-ALL) patients, as reported by Liu, Yu, et al. (2017). The original lesion coordinates were based on the GRCh37 (hg19) human genome assembly. We converted these coordinates to GRCh38 (hg38) using the UCSC LiftOver tool (<https://genome.ucsc.edu/cgi-bin/hgLiftOver>) prior to running GRIN2 analyses.
#'
#' @format ## `lesion_data`
#' A data frame with 6,861 rows and 5 columns:
#' \describe{
#'   \item{ID}{Patient identifier for the individual affected by the lesion}
#'   \item{chrom}{Chromosome on which the lesion is located}
#'   \item{loc.start}{Lesion start position (in base pairs, hg38)}
#'   \item{loc.end}{Lesion end position (in base pairs, hg38)}
#'   \item{lsn.type}{Type of lesion (e.g., gain, loss, mutation, fusion, etc.)}
#' }
#'
#' @source Adapted from the supplementary tables of Liu, Yu, et al. (2017), *Nature Genetics*
#' (<https://www.nature.com/articles/ng.3909#Sec27>)
"lesion_data"
