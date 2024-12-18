% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/onco.print.props.R
\name{onco.print.props}
\alias{onco.print.props}
\title{Oncoprint proportions}
\usage{
onco.print.props(lsn.data, clr = NULL, hgt = NULL)
}
\arguments{
\item{lsn.data}{data.frame with 5 columns including "ID" which is the subject identifier, "chrom" which is the chromosome on which the lesion is located, "loc.start" with lesion start position, "loc.end" which is the lesion end position), and "lsn.type" which is the lesion category for example gain, mutation, etc..).}

\item{clr}{Lesion colors (If not provided by the user, colors will be automatically assigned using default.grin.colors function).}

\item{hgt}{Manually assign the proportion of the oncoprint rectangle that should be color filled for each lesion group.}
}
\value{
Function return a list of three lists specifying the color assigned to each lesion type, the proportion of the rectangle that should be color filled in the oncoprint based on the average size of each lesion type, and the legend parameters.
}
\description{
The function order lesion types based on their average size and assign the proportion of the oncoprint rectangle that should be color filled based on the average size of each lesion type.
}
\details{
Some patients might be affected by two or more lesion types in the same gene for example gain AND mutations. To get all lesion types represented in the same rectangle in the oncoprint, this function order lesion types based on the average size of each type and assign the proportion of the oncoprint rectangle that should be color filled based on the average size of each lesion type. Color filled proportion of the oncoprint rectangles can be also specified by the user for  each lesion type based on the hgt parameter.
}
\examples{
data(lesion.data)
onco.props=onco.print.props(lesion.data, hgt = c("gain"=4, "loss"=3, "mutation"=2, "fusion"=1))
# if hgt argument is not specified, the lesion category "mutation" for single point mutations will
# be assigned size=1 because it has the smallest average lesion size and will have the smallest
# proportion of the filled oncoprint rectangles 1/4=0.25
}
\references{
Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
}
\author{
{Lakshmi Patibandla \email{LakshmiAnuhya.Patibandla@stjude.org}, Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org} and Stanley Pounds \email{stanley.pounds@stjude.org}}
}
