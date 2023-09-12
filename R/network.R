sparcc <- function(ta, rarefact = 0.05, taxon_name = taxon, sample_name = sample) {
  force_optional_dependency(
    "SpiecEasi",
    "\nInstall using: install_github('zdk123/SpiecEasi')"
  )


  sample_name <- rlang::enquo(sample_name)
  taxon_name <- rlang::enquo(taxon_name)

  cutoff <- nrow(ta$samples) * rarefact
  if (!"occurrence" %in% names(ta$taxa)) {
    ta_occ <- ta %>% add_occurrences()
  } else {
    ta_occ <- ta
  }
  ta_occ <- ta_occ %>%
    filter_taxa(occurrence >= cutoff)
  counts <- ta_occ %>%
    counts_matrix(sample_name = !!sample_name, taxon_name = !!taxon_name)
  sparcc.out <- SpiecEasi::sparcc(counts)
  sparcc.out$names <- colnames(counts)
  sparcc.out
}

sparcc_to_network <- function(sparcc.out, treshold = 0.1) {
  force_optional_dependency("Matrix")
  se.sparcc.graph <- abs(sparcc.out$Cor) >= treshold
  Matrix::diag(se.sparcc.graph) <- 0
  sparcc.graph <- Matrix::Matrix(se.sparcc.graph, sparse = T)
  sparcc.elist <- as.data.frame(as.matrix(summary(sparcc.graph * sparcc.out$Cor)))
  network <- as.data.frame(as.matrix(Matrix::Matrix(se.sparcc.graph) * sparcc.out$Cor))
  rownames(network) <- sparcc.out$names
  colnames(network) <- sparcc.out$names
  as.matrix(network)
}

markov_cluster <- function(network, min_n = 3, visualize = F) {
  force_optional_dependency("MCL")
  network[network < 0] <- 0
  res <- MCL::mcl(network, addLoops = T, ESM = T)

  if (visualize) {
    force_optional_dependency("igraph")
    gd <- igraph::graph.adjacency(res$Equilibrium.state.matrix)
    plot(gd)
  }

  clusters <- tibble(taxon = colnames(network), cluster = res$Cluster)
  k <- clusters %>%
    group_by(cluster) %>%
    summarize(taxa_per_clust = n()) %>%
    filter(cluster > 0, taxa_per_clust >= min_n) %>%
    pull(cluster)

  clusters %>% filter(cluster %in% k) %>% mutate(cluster = str_c("c", cluster))
}

add_sparcc_network_clusters <- function(
    ta,
    rarefact = 0.05, network_tresh = 0.1, min_n_cluster = 3,
    taxon_name = taxon, sample_name = sample) {
  sample_name <- rlang::enquo(sample_name)
  taxon_name <- rlang::enquo(taxon_name)

  clusters <- ta %>%
    sparcc(
      rarefact = rarefact,
      sample_name = !!sample_name, taxon_name = !!taxon_name
    ) %>%
    sparcc_to_network(treshold = network_tresh) %>%
    markov_cluster(min_n = min_n_cluster)

  ta$taxa <- ta$taxa %>% left_join(clusters, by=rlang::quo_name(taxon_name))
  ta
}

pca_taxa <- function(ta, cluster_name){

    cm <- ta %>% filter_taxa(cluster==cluster_name) %>% 
        counts_matrix()

    res <- cm %>% 
        stats::prcomp(scale=T)
    colname = paste0("scaled_pca_", cluster_name)
    
    ta$samples <- ta$samples %>% left_join(
        tibble(sample=rownames(res$x), !!colname:=res$x[,1]),
        by="sample"
    )
    ta
}

add_eigentaxa <- function(ta) {

    ta_tmp <- ta
    if (! "cluster" %in% names(ta$taxa)) {
    # keep cluster ids in tax table
    ta <- ta %>% add_sparcc_network_clusters()
    ta_tmp <- ta %>% clr_transform_counts(overwrite=T)
    }

    clusters <- ta_tmp$taxa %>% drop_na() %>% pull(cluster)

    for (clust in unique(clusters)){
    tryCatch({
        ta_tmp <- pca_taxa(ta_tmp, clust)}, 
        error=function(cond){warning(cond)}
    )
    }

    ta$samples <- ta_tmp$samples
    ta
}

