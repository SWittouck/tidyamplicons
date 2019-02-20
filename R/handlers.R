
change_ids_samples <- function(ta, sample_new) {

  if (any(duplicated(ta$samples[[sample_new]]))) {
    stop("the new sample ids are not unique")
  }

  ta$samples <-
    ta$samples %>%
    rename(sample_new = !! sample_new) %>%
    mutate_at("sample_new", as.character)

  ta$abundances <-
    ta$abundances %>%
    left_join(ta$samples %>% select(sample, sample_new)) %>%
    select(- sample) %>%
    rename(sample = sample_new)

  ta$samples <-
    ta$samples %>%
    select(- sample) %>%
    rename(sample = sample_new)

  ta

}

change_ids_taxa <- function(ta, taxon_new) {

  if (any(duplicated(ta$taxa[[taxon_new]]))) {
    stop("the new taxon ids are not unique")
  }

  ta$taxa <-
    ta$taxa %>%
    rename(taxon_new = !! taxon_new) %>%
    mutate_at("taxon_new", as.character)

  ta$abundances <-
    ta$abundances %>%
    left_join(ta$taxa %>% select(taxon, taxon_new)) %>%
    select(- taxon) %>%
    rename(taxon = taxon_new)

  ta$taxa <-
    ta$taxa %>%
    select(- taxon) %>%
    rename(taxon = taxon_new)

  ta

}

# Preprocessing: delete all sample variables that are different within
# groups of samples that need to be merged. Keep the sample variable!
aggregate_samples <- function(ta) {

  # sample table with only old and new sample names
  names <- ta$samples %>%
    select(- sample) %>%
    distinct() %>%
    mutate(sample_new = paste("m", 1:n(), sep = "")) %>%
    right_join(ta$samples) %>%
    select(sample, sample_new)

  # adapt sample table with new names
  ta$samples <- ta$samples %>%
    left_join(names) %>%
    select(- sample) %>%
    rename(sample = sample_new) %>%
    distinct()

  # merge samples in abundance table and adapt with new names
  ta$abundances <- ta$abundances %>%
    left_join(names) %>%
    select(- sample) %>%
    group_by(sample_new, taxon) %>%
    summarize(abundance = sum(abundance)) %>%
    ungroup() %>%
    rename(sample = sample_new)

  # return ta object
  ta

}

# Preprocessing: delete all taxonomic levels you do not want and all other junk,
# but keep the taxon variable!
aggregate_taxa <- function(ta) {

  # this avoids some problems
  ta$taxa[is.na(ta$taxa)] <- "unknown"

  # taxon table with only old and new taxon names
  names <- ta$taxa %>%
    select(- taxon) %>%
    distinct() %>%
    mutate(taxon_new = paste("t", 1:n(), sep = "")) %>%
    right_join(ta$taxa) %>%
    select(taxon, taxon_new)

  # adapt taxon table with new names
  ta$taxa <- ta$taxa %>%
    left_join(names) %>%
    select(- taxon) %>%
    rename(taxon = taxon_new) %>%
    distinct()

  # merge taxa in abundance table and adapt with new names
  ta$abundances <- ta$abundances %>%
    left_join(names) %>%
    select(- taxon) %>%
    group_by(taxon_new, sample) %>%
    summarize(abundance = sum(abundance)) %>%
    ungroup() %>%
    rename(taxon = taxon_new)

  # this avoids some problems (part 2)
  ta$taxa[ta$taxa == "unknown"] <- NA

  # return ta object
  ta

}

trim_asvs <- function(ta, start, end) {

  ta$taxa <- ta$taxa %>%
    mutate(taxon = str_sub(taxon, start = !! start, end = !! end))
  ta$abundances <- ta$abundances %>%
    mutate(taxon = str_sub(taxon, start = !! start, end = !! end))
  ta <- merge_redundant_taxa(ta)

  ta

}

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
