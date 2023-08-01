#' Apply a sample filtering to the taxon and abundance tables
#'
#' DEPRECATED, see \code{\link{filter_samples}}
#'
#' Should only be used internally.
#'
#' @export
process_sample_selection <- function(ta) {

  # filter abundance table
  selected_samples <- ta$samples$sample_id
  ta$abundances <- ta$abundances %>%
    filter(sample_id %in% selected_samples)

  # filter taxon table
  selected_taxa <- ta$abundances$taxon_id %>% unique()
  ta$taxa <- ta$taxa %>%
    filter(taxon_id %in% selected_taxa)

  # return ta object
  ta

}

#' Apply a taxon filtering to the abundance table
#'
#' DEPRECATED, see \code{\link{filter_taxa}}
#'
#' Should only be used internally.
#'
#' @export
process_taxon_selection <- function(ta) {

  # filter abundance table
  selected_taxa <- ta$taxa$taxon_id
  ta$abundances <- ta$abundances %>%
    filter(taxon_id %in% selected_taxa)

  # return ta object
  ta

}

#' Apply an abundance filtering to the taxon table
#'
#' DEPRECATED, see \code{\link{filter_abundances}}
#'
#' Should only be used internally.
#'
#' @export
process_abundance_selection <- function(ta) {

  # filter taxon table
  selected_taxa <- ta$abundances$taxon_id %>% unique()
  ta$taxa <- ta$taxa %>%
    filter(taxon_id %in% selected_taxa)

  # return ta object
  ta

}

#' Return rank names associated with a tidytacos object
#'
#' @export
rank_names <- function(ta) {

  if (! is.null(ta$rank_names)) {
    ta$rank_names
  } else {
    c("kingdom", "phylum", "class", "order", "family", "genus")
  }

}

#' Set rank names for a tidytacos object
#' @param ta a tidytacos object
#' @param rank_names a vector containing the names of the ranks used to describe the taxa
#' @export
set_rank_names <- function(ta, rank_names) {

  ta$rank_names <- rank_names

  ta

}

# for internal use, I think
merge_redundant_taxa <- function(ta) {

  # merge taxa in taxon table
  ta$taxa <- ta$taxa %>%
    group_by(taxon_id) %>%
    summarize_all(function(x) {
      x <- unique(x)
      x <- x[! is.na(x)]
      if (length(x) == 1) return(x)
      as.character(NA)
    })

  # merge taxa in abundances table
  ta$abundances <- ta$abundances %>%
    group_by(sample_id, taxon_id) %>%
    summarise(abundance = sum(abundance)) %>%
    ungroup()

  ta

}

retain_taxon_id <- function(ta) {
  if ((! "taxon_id" %in% names(ta$taxa)) ||
  (! "taxon_id" %in% names(ta$abundances)) || 
  (is.null(ta$abundances$taxon_id)) ||
  (is.null(ta$taxa$taxon_id))) {
    stop("You cannot delete the taxon_id column")
  }
}

retain_sample_id <- function(ta) {
  if ((! "sample_id" %in% names(ta$samples)) || 
  (! "sample_id" %in% names(ta$abundances)) ||
  (is.null(ta$samples$sample_id)) ||
  (is.null(ta$abundances$sample_id))) {
    stop("You cannot delete the sample_id column")
  }
}

retain_abundances <- function(ta) {
  if (! "abundance" %in% names(ta$abundances) ||
  (is.null(ta$abundances$abundance))) {
    stop("You cannot delete the abundance column")
  }
}

any_taxa_left <- function(ta,
  error_message = "No taxa left after filtering") {
    if (length(ta$taxa$taxon_id) == 0) {
    stop(error_message)
  }
}

any_samples_left <- function(ta,
  error_message = "No samples left after filtering") {
    if (length(ta$samples$sample_id) == 0) {
    stop(error_message)
  }
}