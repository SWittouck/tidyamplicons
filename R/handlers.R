
# sample_id_new should be an expression that evaluates to a unique sample
# identifier
change_id_samples <- function(ta, sample_id_new) {

  sample_id_new <- rlang::enexpr(sample_id_new)

  ta <- mutate_samples(ta, sample_id_new = as.character(!! sample_id_new))

  if (any(duplicated(ta$samples$sample_id_new))) {
    stop("the new sample ids are not unique")
  }

  ta$abundances <-
    ta$abundances %>%
    left_join(ta$samples %>% select(sample_id, sample_id_new)) %>%
    select(- sample_id) %>%
    rename(sample_id = sample_id_new)

  ta$samples <-
    ta$samples %>%
    select(- sample_id) %>%
    rename(sample_id = sample_id_new)

  ta

}

# taxon_id_new should be an expression that evaluates to a unique taxon
# identifier
change_id_taxa <- function(ta, taxon_id_new) {

  taxon_id_new <- rlang::enexpr(taxon_id_new)

  ta <- mutate_taxa(ta, taxon_id_new = as.character(!! taxon_id_new))

  if (any(duplicated(ta$taxataxon_id_new))) {
    stop("the new taxon ids are not unique")
  }

  ta$abundances <-
    ta$abundances %>%
    left_join(ta$taxa %>% select(taxon_id, taxon_id_new)) %>%
    select(- taxon_id) %>%
    rename(taxon_id = taxon_id_new)

  ta$taxa <-
    ta$taxa %>%
    select(- taxon_id) %>%
    rename(taxon_id = taxon_id_new)

  ta

}

# Preprocessing: delete all sample variables that are different within
# groups of samples that need to be merged. Keep the sample variable!
aggregate_samples <- function(ta) {

  # sample table with only old and new sample names
  names <- ta$samples %>%
    select(- sample_id) %>%
    distinct() %>%
    mutate(sample_id_new = paste("m", 1:n(), sep = "")) %>%
    right_join(ta$samples) %>%
    select(sample_id, sample_id_new)

  # adapt sample table with new names
  ta$samples <- ta$samples %>%
    left_join(names) %>%
    select(- sample_id) %>%
    rename(sample_id = sample_id_new) %>%
    distinct()

  # merge samples in abundance table and adapt with new names
  ta$abundances <- ta$abundances %>%
    left_join(names) %>%
    select(- sample_id) %>%
    group_by(sample_id_new, taxon_id) %>%
    summarize(abundance = sum(abundance)) %>%
    ungroup() %>%
    rename(sample_id = sample_id_new)

  # return ta object
  ta

}

# Two options to call this function:
# - if the rank you are interested in is in the standard list, just supply it as
# an argument
# - if not, delete all taxon variables except taxon_id and the ranks you are
# still interested in prior to calling this function
aggregate_taxa <- function(ta, rank = NULL) {

  if (! is.null(rank)) {

    ranks <- c("kingdom", "phylum", "class", "order", "family", "genus")

    if (! rank %in% ranks) {
      stop("the rank you supplied is not one of the standard ranks; ...
           it should be lower case")
    }

    rank_index <- which(ranks == rank)
    ranks_to_keep <- ranks[1:rank_index]
    ta <- select_taxa(ta, taxon_id, !! ranks_to_keep)

  }

  # this avoids some problems
  ta$taxa[is.na(ta$taxa)] <- "unknown"

  # taxon table with only old and new taxon names
  names <- ta$taxa %>%
    select(- taxon_id) %>%
    distinct() %>%
    mutate(taxon_id_new = paste("t", 1:n(), sep = "")) %>%
    right_join(ta$taxa) %>%
    select(taxon_id, taxon_id_new)

  # adapt taxon table with new names
  ta$taxa <- ta$taxa %>%
    left_join(names) %>%
    select(- taxon_id) %>%
    rename(taxon_id = taxon_id_new) %>%
    distinct()

  # merge taxa in abundance table and adapt with new names
  ta$abundances <- ta$abundances %>%
    left_join(names) %>%
    select(- taxon_id) %>%
    group_by(taxon_id_new, sample_id) %>%
    summarize(abundance = sum(abundance)) %>%
    ungroup() %>%
    rename(taxon_id = taxon_id_new)

  # this avoids some problems (part 2)
  ta$taxa[ta$taxa == "unknown"] <- NA

  # return ta object
  ta

}

trim_asvs <- function(ta, start, end) {

  ta$taxa <- ta$taxa %>%
    mutate(taxon_id = str_sub(taxon_id, start = !! start, end = !! end))
  ta$abundances <- ta$abundances %>%
    mutate(taxon_id = str_sub(taxon_id, start = !! start, end = !! end))
  ta <- merge_redundant_taxa(ta)

  ta

}

select_samples <- function(ta, ...) {

  ta$samples <- ta$samples %>%
    select(...)

  if (! "sample_id" %in% names(ta$samples)) {
    stop("you cannot delete the sample_id column")
  }

  ta

}

select_taxa <- function(ta, ...) {

  ta$taxa <- ta$taxa %>%
    select(...)

  if (! "taxon_id" %in% names(ta$taxa)) {
    stop("you cannot delete the taxon_id column")
  }

  ta

}

select_abundances <- function(ta, ...) {

  ta$abundances <- ta$abundances %>%
    select(...)

  if (! all(c("sample_id", "taxon_id", "abundance") %in% names(ta$abundances))) {
    stop("you cannot delete the sample_id, taxon_id or abundance columns")
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
