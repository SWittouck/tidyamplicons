
add_lib_size <- function(ta, step = "current") {

  # remove lib_size if already present
  ta$samples$lib_size <- NULL

  if (step == "current") {

    # make table with sample and library size
    lib_sizes <- ta$abundances %>%
      group_by(sample) %>%
      summarize(lib_size = sum(abundance)) %>%
      select(sample, lib_size)

  } else {

    # make table with sample and library size
    step_of_interest <- step
    lib_sizes <- ta$lib_sizes %>%
      filter(step == step_of_interest) %>%
      select(sample, lib_size)

  }

  # add library size to sample table
  ta$samples <- left_join(ta$samples, lib_sizes) %>%
    mutate(lib_size = ifelse(is.na(lib_size), 0, lib_size))

  # return ta object
  ta

}

# function will also add relative abundances if not present
add_diversity_measures <- function(ta) {

  # add relative abundances if not present
  remove_rel_abundance <- FALSE
  if (is.null(ta$abundances$rel_abundance)) {
    remove_rel_abundance <- TRUE
    ta <- add_rel_abundance(ta)
  }

  # make table with sample, divObserved and divInvSimpson
  diversities <- ta$abundances %>%
    group_by(sample) %>%
    summarize(
      div_observed = n(),
      div_inv_simpson = 1 / sum(rel_abundance ^ 2)
    ) %>%
    ungroup()

  # add diversity measure to sample table
  ta$samples = left_join(ta$samples, diversities)

  # remove relative abundances
  if (remove_rel_abundance) ta$abundances$rel_abundance <- NULL

  # return ta object
  ta

}

add_sample_clustered <- function(ta) {

  # make relative abundance matrix
  rel_abundance_matrix <- get_rel_abundance_matrix(ta)

  # make Bray-Curtis distance matrix
  dist_matrix <- vegdist(rel_abundance_matrix, method = "bray")

  # perform hierarchical clustering
  clust <- hclust(dist_matrix, method = "average")

  # make table with samples in order of clustering
  samples_clustered <- tibble(
    sample = clust$labels[clust$order],
    sample_clustered = factor(sample, levels = sample)
  )

  # add sample_clustered to samples table
  ta$samples <- ta$samples %>%
    left_join(samples_clustered)

  # return ta object
  ta

}

add_pcoa <- function(ta) {

  # make relative abundance matrix
  rel_abundance_matrix <- get_rel_abundance_matrix(ta)

  # make Bray-Curtis distance matrix
  dist_matrix = vegdist(rel_abundance_matrix, method = "bray")

  # perform PCoA
  pcoa <- cmdscale(dist_matrix, k = 2, eig = T, list = T)
  pcoa_variances <- pcoa$eig/sum(pcoa$eig)
  pcoa_dimensions <- pcoa$points %>%
    as.data.frame() %>%
    rownames_to_column(var = "sample") %>%
    rename(pcoa1 = V1, pcoa2 = V2)

  # add PCoA dimensions to sample table
  ta$samples <- ta$samples %>%
    left_join(pcoa_dimensions)

  # add PCoA variances to ta object
  ta$pcoa_variances <- pcoa_variances

  # return ta object
  ta

}
