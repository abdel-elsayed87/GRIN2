% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clin.data.R
\docType{data}
\name{clin.data}
\alias{clin.data}
\title{Example T-ALL Dataset Clinical Data}
\format{
\subsection{\code{clin.data}}{

A data frame with 265 rows and 11 columns:
\describe{
\item{ID}{Patient identifier}
\item{Sex}{Patient gender}
\item{Race}{Patient race}
\item{Age_Days}{Patient age in days}
\item{WBC}{White Blood Cell (WBC) count}
\item{MRD29}{Minimal Residual Disease (MRD) percentage}
\item{MRD.binary}{MRD as a categorical variable (0 if MRD<=0.1 or 1 if MRD>0.1)}
\item{os.time}{Overall survival time in years (time between diagnosis and either the last follow-up or death)}
\item{os.censor}{Survival status (0=alive at the last follow-up or, 1=dead)}
\item{efs.time}{Event-free survival time in years}
\item{efs.censor}{Event indicator (0=censored without event or, 1=event)}
}
}
}
\source{
Data was extracted from the supplementary material tables of the published Liu, Yu, et al. (2017) manuscript \url{https://www.nature.com/articles/ng.3909#Sec27} and the publicly available clinical data on TARGET database. The two files were merged and selected list of variables were kept in the final clinical data file.
}
\usage{
clin.data
}
\description{
Clinical data file showing demographic and clinical outcomes of 265 newly diagnosed T-cell Acute Lymphoblastic Leukemia (T-ALL) patients that was reported by Liu, Yu, et al. (2017).
}
\keyword{datasets}
