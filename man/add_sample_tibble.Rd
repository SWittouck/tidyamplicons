% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adders_samples.R
\name{add_sample_tibble}
\alias{add_sample_tibble}
\title{Add sample table to the tidyamplicons object}
\usage{
add_sample_tibble(ta, sample_tibble)
}
\arguments{
\item{ta}{Tidyamplicons object.}

\item{sample_tibble}{A tibble containing sample data for each sample. samples
should be rows, while sample data should be columns. At least one column
name needs to be shared with the sample tibble of ta. The default shared
column name is 'sample'.}
}
\description{
\code{add_sample_tibble} adds a sample tibble to the tidyamplicons object.
}
\details{
This function adds a sample tibble containing metadata for each sample to the
tidyamplicons object. It is used after initiating a tidyamplicons object
using a numerical abundance matrix and the function
\code{\link{create_tidyamplicons}}. Also see \code{\link{add_taxon_tibble}}
to update the taxon data of the tidyamplicons object.
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
                     taxa_are_columns = FALSE
                     )

# Initiate sample tibble
sample <- c("sample1", "sample2")
environment <- c("food fermentation", "human stool")
sample_tibble <- tibble::tibble(sample, environment)

# Add sample tibble to tidyamplicons object
data <- data \%>\%
add_sample_tibble(sample_tibble)

}
