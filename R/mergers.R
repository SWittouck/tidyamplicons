
# Preprocessing: delete all sample variables that are different within
# groups of samples that need to be merged. Keep the sample variable!
merge_samples <- function(ta) {

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
merge_taxa <- function(ta) {

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

# requires that both as objects have a "run" variable in their samples table
merge_tidyamplicons <- function(ta1, ta2) {

  # make sure that sample names are unique
  ta1$samples <- ta1$samples %>%
    mutate(sample_new = paste(run, sample, sep = "_"))
  ta2$samples <- ta2$samples %>%
    mutate(sample_new = paste(run, sample, sep = "_"))
  ta1 <- process_new_sample_name(ta1)
  ta2 <- process_new_sample_name(ta2)

  # merge sample tables
  samples <- bind_rows(ta1$samples, ta2$samples)

  # merge taxa tables
  taxa <- bind_rows(ta1$taxa, ta2$taxa) %>%
    select(taxon, kingdom, phylum, class, order, family, genus, species) %>%
    group_by(taxon) %>%
    mutate(i = 1:n()) %>%
    ungroup() %>%
    filter(i == 1) %>%
    select(- i)

  # merge abundances tables
  abundances <- bind_rows(ta1$abundances, ta2$abundances)

  # make new ta object
  ta <- make_tidyamplicons(samples, taxa, abundances)

  # give new sample names in new ta object
  ta$samples <- ta$samples %>%
    mutate(sample_new = paste("m", 1:n(), sep = ""))
  ta <- process_new_sample_name(ta)

  # return ta object
  ta

}
