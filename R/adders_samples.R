#' Add sample table to the tidytacos object
#'
#' \code{add_sample_tibble} adds a sample tibble to the tidytacos object.
#'
#' This function adds a sample tibble containing metadata for each sample to the
#' tidytacos object. It is used after initiating a tidytacos object
#' using a numerical abundance matrix and the function
#' \code{\link{create_tidytacos}}. Also see \code{\link{add_taxon_tibble}}
#' to update the taxon data of the tidytacos object.
#'
#' @param ta tidytacos object.
#' @param sample_tibble A tibble containing sample data for each sample. samples
#'   should be rows, while sample data should be columns. At least one column
#'   name needs to be shared with the sample tibble of ta. The default shared
#'   column name is 'sample'.
#'
#' @examples
#' # Initiate abundance matrix
#' x <- matrix(
#'  c(1500, 1300, 280, 356),
#'  ncol = 2
#' )
#' rownames(x) <- c("taxon1", "taxon2")
#' colnames(x) <- c("sample1", "sample2")
#'
#' # Convert to tidytacos object
#' data <- create_tidytacos(x,
#'                      taxa_are_columns = FALSE
#'                      )
#'
#' # Initiate sample tibble
#' sample <- c("sample1", "sample2")
#' environment <- c("food fermentation", "human stool")
#' sample_tibble <- tibble::tibble(sample, environment)
#'
#' # Add sample tibble to tidytacos object
#' data <- data %>%
#' add_sample_tibble(sample_tibble)
#'
#' @export
add_sample_tibble <- function(ta, sample_tibble) {

  purrr::modify_at(ta, "samples", left_join, sample_tibble)

}

#' Add total abundance per sample
#'
#' \code{add_lib_size} adds the total abundance per sample to the samples tibble
#' of a tidytacos object.
#'
#' This function adds the total abundance per sample to the samples tibble of a
#' tidytacos object under the variable name total_abundance.
#'
#' @param ta tidytacos object.
#'
#' @examples
#' # Initiate abundance matrix
#' x <- matrix(
#'  c(1500, 1300, 280, 356),
#'  ncol = 2
#' )
#' rownames(x) <- c("taxon1", "taxon2")
#' colnames(x) <- c("sample1", "sample2")
#'
#' # Convert to tidytacos object
#' data <- create_tidytacos(x,
#'                      taxa_are_columns = FALSE
#'                      )
#'
#' # Add total abundance
#' data <- data %>%
#'  add_lib_size()
#'
#' @export
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
  ta$samples <-
    ta$samples %>%
    left_join(lib_sizes, by = "sample_id") %>%
    mutate(lib_size = ifelse(is.na(lib_size), 0, lib_size))

  # return ta object
  ta

}


#' Add alpha diversity measures
#'
#' \code{add_diversity_measures} adds two alpha diversity measures to the
#' samples tibble of a tidytacos object.
#'
#' This function adds two alpha diversity measures (observed and inverse
#' Simpson) to the samples tibble of a tidytacos object under the variable
#' names observed and inverse_simpson, respectively. This function will also
#' add relative abundances if not present using \code{\link{add_rel_abundance}}.
#'
#' @param ta tidytacos object.
#'
#' @examples
#' # Initiate abundance matrix
#' x <- matrix(
#'  c(1500, 1300, 280, 356),
#'  ncol = 2
#' )
#' rownames(x) <- c("taxon1", "taxon2")
#' colnames(x) <- c("sample1", "sample2")
#'
#' # Convert to tidytacos object
#' data <- create_tidytacos(x,
#'                      taxa_are_columns = FALSE
#'                      )
#'
#' # Add total abundance
#' data <- data %>%
#'  add_diversity_measures()
#' @export
add_alphas <- function(ta) {

  # if rel abundances not present: add temporarily
  rel_abundance_tmp <- ! "rel_abundance" %in% names(ta$abundances)
  if (rel_abundance_tmp) ta <- add_rel_abundance(ta)

  # make table with sample, divObserved and divInvSimpson
  diversities <- ta$abundances %>%
    filter(abundance > 0) %>%
    group_by(sample_id) %>%
    summarize(
      observed = n(),
      inverse_simpson = 1 / sum(rel_abundance ^ 2)
    ) %>%
    ungroup()

  # add diversity measure to sample table
  ta$samples = left_join(ta$samples, diversities, by = "sample_id")

  # cleanup
  if (rel_abundance_tmp) ta$abundances$rel_abundance <- NULL

  # return ta object
  ta

}

#' Add some alpha diversity measures
#'
#' DEPRECATED, use \code{\link{add_alphas}}
#'
#' @export
add_diversity_measures <- function(ta) {

  ta %>%
    add_alphas %>%
    purrr::modify_at(
      "samples", rename, div_observed = observed,
      div_inv_simpson = inverse_simpson
    )

}

#' Add clustered sample order
#'
#' \code{add_sample_clustered} adds a new variable defining a sample order based
#' on similarity after clustering to the samples tibble of a tidytacos
#' object.
#'
#' This function calculates the Bray-Curtis distance between samples followed by
#' hierarchical average linkage clustering of samples. It will then add a new
#' factor variable "samples_clustered" to the samples tibble of a tidytacos
#' object. This function is extremely useful if one wants to plot similar
#' samples together.
#'
#' @param ta tidytacos object.
#'
#' @examples
#' # Initiate abundance matrix
#' x <- matrix(
#'  c(1500, 1300, 280, 356),
#'  ncol = 2
#' )
#' rownames(x) <- c("taxon1", "taxon2")
#' colnames(x) <- c("sample1", "sample2")
#'
#' # Convert to tidytacos object
#' data <- create_tidytacos(x,
#'                      taxa_are_columns = FALSE
#'                      )
#'
#' # Add total abundance
#' data <- data %>%
#'  add_sample_clustered()
#'
#' @export
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
    left_join(samples_clustered, by = "sample_id")

  # return ta object
  ta

}

#' Add first two dimensions of PCOA
#'
#' \code{add_pcoa} adds the first two dimensions of a principal components
#' analysis on a Bray-Curtis dissimilarity matrix to two new variables of the
#' samples tibble of a tidytacos object.
#'
#' This function calculates the Bray-Curtis distance between samples followed by
#' a principal components analysis. It will then add the two first dimensions to
#' the samples tibble of a tidytacos object named "pcoa1" and "pcoa2". This
#' function will also add relative abundances if not present using
#' \code{\link{add_rel_abundance}}.
#'
#' @param ta tidytacos object.
#'
#' @examples
#' # Initiate abundance matrix
#' x <- matrix(
#'  c(1500, 1300, 280, 356, 456, 678),
#'  ncol = 3
#' )
#' rownames(x) <- c("taxon1", "taxon2")
#' colnames(x) <- c("sample1", "sample2", "sample3")
#'
#' # Convert to tidytacos object
#' data <- create_tidytacos(x,
#'                      taxa_are_columns = FALSE
#'                      )
#'
#' # Add total abundance
#' data <- data %>%
#'  add_pcoa()
#'
#' @export
add_pcoa <- function(ta) {

  # make relative abundance matrix
  rel_abundance_matrix <- get_rel_abundance_matrix(ta)

  # make Bray-Curtis distance matrix
  dist_matrix = vegdist(rel_abundance_matrix, method = "bray")

  # perform PCoA
  pcoa <- cmdscale(dist_matrix, k = 2, eig = T, list = T)
  pcoa_variances <- pcoa$eig / sum(pcoa$eig)
  pcoa_dimensions <- pcoa$points %>%
    `colnames<-`(c("pcoa1", "pcoa2")) %>%
    as_tibble() %>%
    mutate(sample_id = !! rownames(pcoa$points))

  # add PCoA dimensions to sample table
  ta$samples <- ta$samples %>%
    left_join(pcoa_dimensions, by = "sample_id")

  # add PCoA variances to ta object
  ta$pcoa_variances <- pcoa_variances

  # return ta object
  ta

}

#' Add spike ratio
#'
#' \code{add_spike_ratio} adds a new variable showing the ratio total abundance
#' to spike abundance to the samples tibble of a tidytacos object.
#'
#' This function calculates the spike ratio defined as the total sample
#' abundance to the spike abundance and adds this as a new variable
#' "spike_ratio" to the samples tibble of a tidytacos object. This function
#' is useful if a DNA spike was added prior to sequencing and is based on the
#' method described by
#' \href{https://doi.org/10.1016/j.soilbio.2016.02.003}{Smets et al., 2016}.
#'
#' Credits to Wenke Smets for the idea of spiking samples prior to 16S
#' sequencing and the initial implementation of this function.
#'
#' @param ta A tidytacos object.
#' @param spike_taxon The taxon_id of the spike.
#'
#' @examples
#' # Initiate abundance matrix
#' x <- matrix(
#'  c(1500, 1300, 280, 356),
#'  ncol = 2
#' )
#' rownames(x) <- c("taxon1", "taxon2")
#' colnames(x) <- c("sample1", "sample2")
#'
#' # Convert to tidytacos object
#' data <- create_tidytacos(x,
#'                      taxa_are_columns = FALSE
#'                      )
#'
#' # Add total abundance
#' data <- data %>%
#'  add_spike_ratio(spike_taxon = "t1")
#
#' @export
add_spike_ratio <- function(ta, spike_taxon) {

  # if lib_size not present: add temporarily
  lib_size_tmp <- ! "lib_size" %in% names(ta$samples)
  if (lib_size_tmp) ta <- add_lib_size(ta)

  # make sample table with spike abundances
  spike_abundances <- ta$abundances %>%
    filter(taxon_id == spike_taxon) %>%
    select(sample_id, spike_abundance = abundance)

  # calculate spike ratio (non-spike abundance to spike abundance)
  ta$samples <- ta$samples %>%
    left_join(spike_abundances, by = "sample_id") %>%
    mutate(spike_ratio = ( lib_size - spike_abundance ) / spike_abundance)

  # remove spike_abundance
  ta$samples$spike_abundance <- NULL

  # cleanup
  if (lib_size_tmp) ta$samples$lib_size <- NULL

  # return ta object
  ta

}

#' Add cluster number
#'
#' \code{add_cluster} adds a new variable to the samples tibble of a
#' tidytacos object defining to what cluster a sample belongs.
#'
#' This function calculates the Bray-Curtis distance between samples followed by
#' hierarchical average linkage clustering of samples. The user provides a
#' number of desired clusters which will be used to assign the samples to. A new
#' variable named "cluster" will be added to the samples tibble of a
#' tidytacos object defining to what cluster a sample belongs.
#'
#' @param ta tidytacos object.
#' @param n_clusters Numerical. Number of desired clusters.
#'
#' @examples
#' # Initiate abundance matrix
#' x <- matrix(
#'  c(1500, 1300, 280, 356),
#'  ncol = 2
#' )
#' rownames(x) <- c("taxon1", "taxon2")
#' colnames(x) <- c("sample1", "sample2")
#'
#' # Convert to tidytacos object
#' data <- create_tidytacos(x,
#'                      taxa_are_columns = FALSE
#'                      )
#'
#' # Add total abundance
#' data <- data %>%
#'  add_cluster(n_clusters = 2)
#'
# Adds a variable "cluster" to the samples table
# To do: merge with add_sample_clustered somehow
#
#' @export
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
    left_join(ta$samples, samples_clusters, by = "sample_id")

  ta

}
