
report_numbers <- function(ta) {

  sprintf("samples: %i", nrow(ta$samples)) %>% message()
  sprintf("taxa: %i", nrow(ta$taxa)) %>% message()
  sprintf("reads: %i", sum(ta$abundances$abundance)) %>% message()

}

get_numbers <- function(ta) {

  c(
    n_samples = nrow(ta$samples),
    n_taxa = nrow(ta$taxa),
    n_reads = sum(ta$abundances$abundance)
  )

}

# for internal use - do not export
get_rel_abundance_matrix <- function(ta) {

  # add relative abundances if not present
  if (is.null(ta$abundances$rel_abundance)) {
    ta <- add_rel_abundance(ta)
  }

  as_abundances_matrix(ta$abundances, value = rel_abundance)

}

# get sample distances as a tidy table
get_betas <- function(ta, unique = T) {

  # make "dist" object with beta values
  rel_abundance_matrix <- get_rel_abundance_matrix(ta)
  betas_dist <- vegdist(rel_abundance_matrix, method = "bray")

  # save number of betas in betas_dist in shortcut variable
  n <- attr(betas_dist, "Size")

  # make tibble with beta values if we want only unique sample pairs
  if (unique) {

    betas <- expand.grid(i = 1:n, j = 1:n) %>%
      filter(i < j) %>%
      mutate(sample_1 = labels(betas_dist)[i]) %>%
      mutate(sample_2 = labels(betas_dist)[j]) %>%
      mutate(beta = betas_dist[n * (i - 1) - i * (i - 1) / 2 + j - i]) %>%
      select(- i, - j)

    # make tibble with beta values if we want all sample pairs (redundant!)
  } else {

    betas <- as.matrix(betas_dist) %>%
      as_tibble() %>%
      mutate(sample_1 = attr(betas_dist, "Labels")) %>%
      gather(key = "sample_2", value = "beta", - sample_1)

  }

  # add sample info to betas table
  betas <- ta$samples %>%
    `names<-`(names(.) %>% str_c("_2")) %>%
    right_join(betas)
  betas <- ta$samples %>%
    `names<-`(names(.) %>% str_c("_1")) %>%
    right_join(betas)

  # return betas table
  betas

}

get_abundances_extended <- function(ta) {

  # make and return large table
  ta$abundances %>%
    left_join(ta$samples) %>%
    left_join(ta$taxa)

}
