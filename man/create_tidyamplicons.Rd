% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/converters.R
\name{create_tidyamplicons}
\alias{create_tidyamplicons}
\title{Initiate tidyamplicons object}
\usage{
create_tidyamplicons(abundance_matrix, taxa_are_columns = TRUE)
}
\arguments{
\item{abundance_matrix}{Numerical matrix containing the abundance data.}

\item{taxa_are_columns}{A logical scalar. Are the taxa defined in columns?}
}
\description{
\code{tidyamplicons} returns a tidyamplicons object given a numeric matrix.
}
\details{
This function initiates a tidyamplicons object based on a numeric matrix. It
will automatically create a dummy taxa and sample table which will need to be
updated using the functions \code{\link{add_taxon_tibble}} and
\code{\link{add_sample_tibble}}.
}
\examples{
# Initiate abundance matrix
x <- matrix(
 c(1500, 1300, 280, 356),
 ncol = 2
)
rownames(x) <- c("taxon1", "taxon2")
colnames(x) <- c("sample1", "sample2")

# Convert to tidyamplicons object
data <- create_tidyamplicons(x,
                     taxa_are_columns = FALSE)


\dontrun{
tidyamplicons("a")
}

}
