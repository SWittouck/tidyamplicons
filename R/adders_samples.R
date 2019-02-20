
add_sample_tibble <- function(ta, sample_tibble) {

  modify_at(ta, "samples", left_join, sample_tibble)

}

add_lib_size <- function(ta, step = "current") {

  # remove lib_size if already present
  ta$samples$lib_size <- NULL

  if (step == "current") {

    # make table with sample and library size
    lib_sizes <- ta$abundances %>%
      group_by(sample_id) %>%
      summarize(lib_size = sum(abundance)) %>%
      select(sample_id, lib_size)

  } else {

    # make table with sample and library size
    step_of_interest <- step
    lib_sizes <- ta$lib_sizes %>%
      filter(step == step_of_interest) %>%
      select(sample_id, lib_size)

  }

  # add library size to sample table
  ta$samples <- left_join(ta$samples, lib_sizes) %>%
    mutate(lib_size = ifelse(is.na(lib_size), 0, lib_size))

  # return ta object
  ta

}

# function will also add relative abundances if not present
add_diversity_measures <- function(ta) {

  # if rel abundances not present: add and remove again on exit
  if (! "rel_abundance" %in% names(ta$abundances)) {
    ta <- add_rel_abundance(ta)
    on.exit(ta$abundances$rel_abundance <- NULL)
  }

  # make table with sample, divObserved and divInvSimpson
  diversities <- ta$abundances %>%
    group_by(sample_id) %>%
    summarize(
      div_observed = n(),
      div_inv_simpson = 1 / sum(rel_abundance ^ 2)
    ) %>%
    ungroup()

  # add diversity measure to sample table
  ta$samples = left_join(ta$samples, diversities)

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
    sample_id = clust$labels[clust$order],
    sample_clustered = factor(sample_id, levels = sample_id)
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
    as_tibble() %>%
    mutate(sample_id = !! rownames(pcoa$points)) %>%
    rename(pcoa1 = V1, pcoa2 = V2)

  # add PCoA dimensions to sample table
  ta$samples <- ta$samples %>%
    left_join(pcoa_dimensions)

  # add PCoA variances to ta object
  ta$pcoa_variances <- pcoa_variances

  # return ta object
  ta

}

# Credits to Wenke Smets for the idea of spiking samples prior to 16S sequencing
# (Smets et al., 2016) and the initial implementation of this function
add_spike_ratio <- function(ta, spike_taxon) {

  # if lib_size not present: add and remove again on exit
  if (! "lib_size" %in% names(ta$samples)) {
    ta <- add_lib_size(ta)
    on.exit(ta$samples$lib_size <- NULL)
  }

  # make sample table with spike abundances
  spike_abundances <- ta$abundances %>%
    filter(taxon_id == spike_taxon) %>%
    select(sample_id, spike_abundance = abundance)

  # calculate spike ratio (non-spike abundance to spike abundance)
  ta$samples <- ta$samples %>%
    left_join(spike_abundances) %>%
    mutate(spike_ratio = ( lib_size - spike_abundance ) / spike_abundance)

  # remove spike_abundance
  ta$samples$spike_abundance <- NULL

  # return ta object
  ta

}

# Adds a variable "cluster" to the samples table
# To do: merge with add_sample_clustered somehow
add_cluster <- function(ta, n_clusters) {

  # make relative abundance matrix
  rel_abundance_matrix <- get_rel_abundance_matrix(ta)

  # make Bray-Curtis distance matrix
  dist_matrix <- vegdist(rel_abundance_matrix, method = "bray")

  # perform hierarchical clustering
  clust <- hclust(dist_matrix, method = "average")

  samples_clusters <-
    tibble(
      sample_id = clust$labels,
      cluster = cutree(clust, k = n_clusters)
    ) %>%
    mutate(cluster = str_c("cluster", cluster, sep = " "))

  ta$samples <-
    left_join(ta$samples, samples_clusters)

  ta

}
