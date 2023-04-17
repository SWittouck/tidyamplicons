#' Rarefy the samples to a given number of reads
#'
#' This function performs rarefying. Make sure that all samples contain at least
#' the minimum number of reads; otherwise, an error might be thrown.
#'
#' @export
rarefy <- function(ta, n, replace = F) {

  ta$abundances <-
    ta$abundances %>%
    group_by(sample_id) %>%
    mutate(
      abundance =
        sample(x = 1:sum(abundance), size = !! n, replace = !! replace) %>%
        cut(breaks = c(0, cumsum(abundance)), labels = taxon_id) %>%
        table() %>%
        as.integer()
    ) %>%
    ungroup()

  ta %>%
    purrr::modify_at("abundances", filter, abundance > 0) %>%
    process_abundance_selection()

}

#' Change sample IDs to a given expression
#'
#' @param ta A tidyamplicons object
#' @param sample_id_new An expression that evaluates to a unique sample
#'   identifier
#'
#' @export
change_id_samples <- function(ta, sample_id_new) {

  sample_id_new <- rlang::enexpr(sample_id_new)

  ta <- mutate_samples(ta, sample_id_new = as.character(!! sample_id_new))

  if (any(duplicated(ta$samples$sample_id_new))) {
    stop("the new sample ids are not unique")
  }

  ta$abundances <-
    ta$abundances %>%
    left_join(
      ta$samples %>% select(sample_id, sample_id_new), by = "sample_id"
    ) %>%
    select(- sample_id) %>%
    rename(sample_id = sample_id_new)

  ta$samples <-
    ta$samples %>%
    select(- sample_id) %>%
    rename(sample_id = sample_id_new)

  ta

}

#' Change taxon IDs to a given expression
#'
#' @param ta A tidyamplicons object
#' @param taxon_id_new An expression that evaluates to a unique taxon
#'   identifier
#'
#' @export
change_id_taxa <- function(ta, taxon_id_new) {

  taxon_id_new <- rlang::enexpr(taxon_id_new)

  ta <- mutate_taxa(ta, taxon_id_new = as.character(!! taxon_id_new))

  if (any(duplicated(ta$taxa$taxon_id_new))) {
    stop("the new taxon ids are not unique")
  }

  ta$abundances <-
    ta$abundances %>%
    left_join(ta$taxa %>% select(taxon_id, taxon_id_new), by = "taxon_id") %>%
    select(- taxon_id) %>%
    rename(taxon_id = taxon_id_new)

  ta$taxa <-
    ta$taxa %>%
    select(- taxon_id) %>%
    rename(taxon_id = taxon_id_new)

  ta

}

#' Aggregate samples with identical values for all metadata
#' @param ta a tidyamplicons object
#' @export
aggregate_samples <- function(ta) {

  # sample table with only old and new sample names
  metadata <- setdiff(names(ta$samples), "sample_id")
  names <- ta$samples %>%
    select(-sample_id) %>%
    distinct() %>%
    mutate(sample_id_new = paste0("m", 1:n())) %>%
    right_join(ta$samples, by=metadata, multiple="all") %>%
    select(sample_id, sample_id_new)

  # adapt sample table with new names
  ta$samples <- ta$samples %>%
    left_join(names, by = "sample_id") %>%
    select(- sample_id) %>%
    rename(sample_id = sample_id_new) %>%
    distinct()

  # merge samples in abundance table and adapt with new names
  ta$abundances <- ta$abundances %>%
    left_join(names, by = "sample_id") %>%
    select(- sample_id) %>%
    group_by(sample_id_new, taxon_id) %>%
    summarize(abundance = sum(abundance)) %>%
    ungroup() %>%
    rename(sample_id = sample_id_new)

  # return ta object
  ta

}

#' Aggregate taxa on a given taxonomic rank
#'
#' There are two ways to call this function:
#'
#' * If the rank you are interested in is in the standard list, just supply it
#' as an argument.
#' * If not, delete all taxon variables except taxon_id and the ranks you are
#' still interested in prior to calling this function.
#'
#' @export
aggregate_taxa <- function(ta, rank = NULL) {

  if (! is.null(rank)) {

    rank_names <-
      rank_names(ta) %>%
      intersect(names(ta$taxa))

    if (length(rank_names) == 0) {
      stop("at least one of the taxonomic rank names should be present ",
           "in the taxon table")
    }

    if (! rank %in% rank_names) {
      stop("the rank you supplied should be one of the rank names")
    }

    rank_index <- which(rank_names == rank)
    rank_names_to_keep <- rank_names[1:rank_index]
    ta <- select_taxa(ta, taxon_id, !! rank_names_to_keep)

  }

  # this avoids some problems
  ta$taxa[is.na(ta$taxa)] <- "unknown"

  ta$taxa <-
    ta$taxa %>%
    chop(taxon_id) %>%
    mutate(taxon_id_new = paste("t", 1:n(), sep = ""))

  id_conversion <-
    ta$taxa %>%
    unnest(taxon_id) %>%
    select(taxon_id, taxon_id_new)

  ta$taxa <-
    ta$taxa %>%
    select(- taxon_id) %>%
    rename(taxon_id = taxon_id_new)

  ta$abundances <-
    ta$abundances %>%
    left_join(id_conversion, by = "taxon_id") %>%
    select(- taxon_id) %>%
    group_by(taxon_id_new, sample_id) %>%
    {
      if ("rel_abundance" %in% names(ta$abundances)) {
        summarize(
          ., abundance = sum(abundance), rel_abundance = sum(rel_abundance)
        )
      } else {
        summarize(., abundance = sum(abundance))
      }
    } %>%
    ungroup() %>%
    rename(taxon_id = taxon_id_new)

  # cleanup
  ta$taxa[ta$taxa == "unknown"] <- NA

  ta

}

#' Trim all sequences
#'
#' This function assumes that the sequence variable in the taxon table is called
#' "sequence".
#'
#' @export
trim_asvs <- function(ta, start, end) {

  ta$taxa <- ta$taxa %>%
    mutate(sequence = str_sub(sequence, start = !! start, end = !! end))
  if ("sequence" %in% names(ta$abundances)){
    ta$abundances <- ta$abundances %>%
      mutate(sequence = str_sub(sequence, start = !! start, end = !! end))
  }
  ta <- merge_redundant_taxa(ta)

  ta

}

#' Retain or remove a set of sample variables
#'
#' @export
select_samples <- function(ta, ...) {

  ta$samples <- ta$samples %>%
    select(...)

  if (! "sample_id" %in% names(ta$samples)) {
    stop("you cannot delete the sample_id column")
  }

  ta

}

#' Retain or remove a set of taxon variables
#'
#' @export
select_taxa <- function(ta, ...) {

  ta$taxa <- ta$taxa %>%
    select(...)

  retain_taxon_id(ta)

  ta

}

#' Retain or remove a set of abundance variables
#'
#' @export
select_abundances <- function(ta, ...) {

  ta$abundances <- ta$abundances %>%
    select(...)

  no_del <- c("sample_id", "taxon_id", "abundance")
  retain_sample_id(ta)
  retain_taxon_id(ta)
  retain_abundances(ta)

  ta

}

#' Create extra variables in the sample table
#'
#' @export
mutate_samples <- function(ta, ...) {

  ta$samples <- ta$samples %>%
    mutate(...)
  retain_sample_id(ta)

  ta

}

#' Create extra variables in the taxon table
#'
#' @export
mutate_taxa <- function(ta, ...) {

  ta$taxa <- ta$taxa %>%
    mutate(...)
  retain_taxon_id(ta)

  ta

}

#' Create extra variables in the abundances table
#'
#' @export
mutate_abundances <- function(ta, ...) {

  ta$abundances <- ta$abundances %>%
    mutate(...)
  retain_sample_id(ta)
  retain_taxon_id(ta)
  retain_abundances(ta)

  ta

}

#' Filter the samples
#'
#' @export
filter_samples <- function(ta, ...) {

  ta$samples <- ta$samples %>%
    filter(...)

  ta <- ta %>%
    process_sample_selection()

  ta

}

#' Filter the taxa
#'
#' @export
filter_taxa <- function(ta, ...) {

  ta$taxa <- ta$taxa %>%
    filter(...)

  ta <- ta %>%
    process_taxon_selection()

  ta

}

#' Filter the abundances
#'
#' @export
filter_abundances <- function(ta, ...) {

  ta$abundances <- ta$abundances %>%
    filter(...)

  ta <- ta %>%
    process_abundance_selection()

  ta

}
