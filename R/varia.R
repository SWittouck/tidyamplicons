
rank_names <- function(ta) {

  if (! is.null(ta$rank_names)) {
    ta$rank_names
  } else {
    c("kingdom", "phylum", "class", "order", "family", "genus")
  }

}

set_rank_names <- function(ta, rank_names) {

  ta$rank_names <- rank_names

  ta

}

# Execute after samples are selected in samples
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

# Execute after taxa are selected in taxa
process_taxon_selection <- function(ta) {

  # filter abundance table
  selected_taxa <- ta$taxa$taxon_id
  ta$abundances <- ta$abundances %>%
    filter(taxon_id %in% selected_taxa)

  # return ta object
  ta

}

# Execute after selection on abundances
process_abundance_selection <- function(ta) {

  # filter taxon table
  selected_taxa <- ta$abundances$taxon_id %>% unique()
  ta$taxa <- ta$taxa %>%
    filter(taxon_id %in% selected_taxa)

  # return ta object
  ta

}

# required that you made a variable "sample_new" in samples
process_new_sample_name <- function(ta) {

  names <- ta$samples %>%
    select(sample_id, sample_id_new)

  ta$abundances <- ta$abundances %>%
    left_join(names, by = "sample_id") %>%
    mutate(sample_id = sample_id_new) %>%
    select(- sample_id_new)

  ta$samples <- ta$samples %>%
    mutate(sample_id = sample_id_new) %>%
    select(- sample_id_new)

  # return ta object
  ta

}

update_lib_sizes <- function(ta, step) {

  # current lib_size in tidy table
  lib_sizes_new <- ta$abundances %>%
    group_by(sample_id) %>%
    summarize(lib_size = sum(abundance)) %>%
    mutate(step = step)

  # make lib_sizes table if it doesn't exist
  if (is.null(ta$lib_sizes)) {
    ta$lib_sizes <- lib_sizes_new %>%
      mutate(step = factor(step))
    # update lib_sizes table if it already existed
  } else {
    levels <- levels(ta$lib_sizes$step)
    levels <- c(levels, step)
    ta$lib_sizes <- ta$lib_sizes %>%
      mutate(step = as.character(step)) %>%
      bind_rows(lib_sizes_new) %>%
      mutate(step = factor(step, levels = !! levels))
  }

  # return ta object
  ta

}

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
