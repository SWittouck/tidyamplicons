# convert matrix with abundances to tidy data frame
as_abundances <- function(abundances_matrix, taxa_are_columns = TRUE, value = "abundance") {

  if (
    ! is.matrix(abundances_matrix) |
    ! is.numeric(abundances_matrix)
  ) stop("first argument should be an abundances matrix")

  if (! taxa_are_columns) abundances_matrix = t(abundances_matrix)

  abundances_matrix %>%
    as_tibble() %>%
    mutate(sample = row.names(abundances_matrix)) %>%
    gather(key = "taxon", value = !! value, - sample) %>%
    filter(!! value > 0)

}

# convert abundances tidy data frame to matrix
as_abundances_matrix <- function(abundances, value = abundance) {

  if (
    ! is.data.frame(abundances) |
    is.null(abundances$taxon) |
    is.null(abundances$sample)
  ) stop("first argument should be an abundances table (data frame)")

  value <- enquo(value)

  abundances_wide <- abundances %>%
    select(sample, taxon, !! value) %>%
    spread(key = taxon, value = !! value, fill = 0)

  abundances_wide %>%
    select(- sample) %>%
    as.matrix() %>%
    `row.names<-`(abundances_wide$sample)

}
