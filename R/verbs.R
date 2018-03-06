select_samples <- function(ta, ...) {

  ta$samples <- ta$samples %>%
    select(...)

  if (! "sample" %in% names(ta$samples)) {
    stop("you cannot delete the sample column")
  }

  ta

}

select_taxa <- function(ta, ...) {

  ta$taxa <- ta$taxa %>%
    select(...)

  if (! "taxon" %in% names(ta$taxa)) {
    stop("you cannot delete the taxon column")
  }

  ta

}

select_abundances <- function(ta, ...) {

  ta$abundances <- ta$abundances %>%
    select(...)

  if (! all(c("sample", "taxon", "abundance") %in% names(ta$abundances))) {
    stop("you cannot delete the sample, taxon or abundance columns")
  }

  ta

}

mutate_samples <- function(ta, ...) {

  # to do: error if sample_id is mutated

  ta$samples <- ta$samples %>%
    mutate(...)

  ta

}

mutate_taxa <- function(ta, ...) {

  # to do: error if taxon_id is mutated

  ta$taxa <- ta$taxa %>%
    mutate(...)

  ta

}

mutate_abundances <- function(ta, ...) {

  # to do: error if sample_id or taxon_id is mutated

  ta$abundances <- ta$abundances %>%
    mutate(...)

  ta

}

filter_samples <- function(ta, ...) {

  ta$samples <- ta$samples %>%
    filter(...)

  ta <- ta %>%
    process_sample_selection()

  ta

}

filter_taxa <- function(ta, ...) {

  ta$taxa <- ta$taxa %>%
    filter(...)

  ta <- ta %>%
    process_taxon_selection()

  ta

}

filter_abundances <- function(ta, ...) {

  ta$abundances <- ta$abundances %>%
    filter(...)

  ta <- ta %>%
    process_abundance_selection()

  ta

}
