% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getters.R
\name{occurrences}
\alias{occurrences}
\title{Get occurrences of taxa in general or per condition}
\usage{
occurrences(ta, condition = NULL, pres_abs = F)
}
\arguments{
\item{ta}{a tidyamplicons object}

\item{condition}{a string denoting a categorical variable in the sample table}

\item{pres_abs}{wether to resort to presence/absense screening}
}
\description{
Returns a tidy table of occurrences: taxon presence counts in samples,
overall or per condition.
}
\details{
Condition should be a categorical variable present in the samples table.
Supply condition as a string.
}
