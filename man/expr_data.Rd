% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expr_data.R
\docType{data}
\name{expr_data}
\alias{expr_data}
\title{Example T-ALL Gene Expression Dataset}
\format{
\subsection{\code{expr_data}}{

A data frame with 417 rows and 265 columns:
\describe{
\item{gene}{Ensembl gene IDs of the 417 selected genes included in the dataset.}
\item{...}{Each remaining column represents a T-ALL patient with log2-normalized expression values.}
}
}
}
\source{
Data extracted from the supplementary materials of Liu, Yu, et al. (2017), \emph{Nature Genetics}. \url{https://www.nature.com/articles/ng.3909#Sec27}
}
\usage{
expr_data
}
\description{
Log2-normalized gene expression data for 417 genes across 265 newly diagnosed T-cell Acute Lymphoblastic Leukemia (T-ALL) patients, as reported by Liu, Yu, et al. (2017).
}
\keyword{datasets}
