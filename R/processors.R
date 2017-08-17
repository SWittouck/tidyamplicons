
# Execute after samples are selected in samples
process_sample_selection <- function(ta) {

  # filter abundance table
  selected_samples <- ta$samples$sample
  ta$abundances <- ta$abundances %>%
    filter(sample %in% selected_samples)

  # filter taxon table
  selected_taxa <- ta$abundances$taxon %>% unique()
  ta$taxa <- ta$taxa %>%
    filter(taxon %in% selected_taxa)

  # return ta object
  ta

}

# Execute after taxa are selected in taxa
process_taxon_selection <- function(ta) {

  # filter abundance table
  selected_taxa <- ta$taxa$taxon
  ta$abundances <- ta$abundances %>%
    filter(taxon %in% selected_taxa)

  # return ta object
  ta

}

# Execute after selection on abundances
process_abundance_selection <- function(ta) {

  # filter taxon table
  selected_taxa <- ta$abundances$taxon %>% unique()
  ta$taxa <- ta$taxa %>%
    filter(taxon %in% selected_taxa)

  # return ta object
  ta

}

# required that you made a variable "sample_new" in samples
process_new_sample_name <- function(ta) {

  names <- ta$samples %>%
    select(sample, sample_new)

  ta$abundances <- ta$abundances %>%
    left_join(names) %>%
    mutate(sample = sample_new) %>%
    select(- sample_new)

  ta$samples <- ta$samples %>%
    mutate(sample = sample_new) %>%
    select(- sample_new)

  # return ta object
  ta

}
