
library(tidyverse)
library(stringr)
library(vegan)

# Author: Stijn Wittouck
# Last modified: 12/07/2017

# These functions are an alternative to, or at least complementary to the
# phyloseq package. They are based on the tidy data pardigm and allow for
# easier, more intuitive data access, processing and visualization. Especially
# if you like the tidyverse way of working with data.

make_tidyamplicons <- function(samples, taxa, abundances) {
  
  # remove abundances of zero
  abundances <- abundances %>%
    filter(abundance > 0)
  
  # make tidyamplicons (ta) object
  ta <- list(
    samples = samples, 
    taxa = taxa, 
    abundances = abundances
  )
  class(ta) <- "tidyamplicons"
  
  # make sure that all tables contain the same unique
  # taxa and samples, and return ta object
  ta %>%
    process_abundance_selection() %>%
    process_sample_selection() %>%
    process_taxon_selection()
  
}

# converts a phyloseq object to a tidyamplicons object
# the plyloseq object should contain absolute abundances
tidy_phyloseq <- function(ps) {
  
  # convert sample data
  samples <- sample_data(ps)@.Data %>%
    `names<-`(sample_data(ps)@names) %>%
    do.call(what = tibble) %>%
    mutate(sample = sample_data(ps)@row.names)
  
  # convert taxon table
  taxa <- tax_table(ps)@.Data %>%
    as_tibble() %>%
    mutate(taxon = tax_table(ps) %>% row.names())
  
  # make sure that taxa are rows in taxon table
  if (! taxa_are_rows(ps)) otu_table(ps) <- t(otu_table(ps))
  
  # convert taxon table
  abundances <- otu_table(ps)@.Data %>%
    as_tibble() %>%
    mutate(taxon = otu_table(ps) %>% row.names()) %>%
    gather(key = "sample", value = abundance, - taxon)
  
  # make and return tidyamplicons object 
  make_tidyamplicons(
    samples = samples,
    taxa = taxa,
    abundances = abundances
  ) 
  
}

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

report_numbers <- function(ta) {
  
  sprintf("samples: %i", nrow(ta$samples)) %>% message()
  sprintf("taxa: %i", nrow(ta$taxa)) %>% message()
  sprintf("reads: %i", sum(ta$abundances$abundance)) %>% message()
  
}

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

# OUTDATED
add_zeros <- function(ta) {
  
  # make abundance table with zeros
  abundances_with_zeros <- ta$abundances %>%
    select(taxon, sample, abundance) %>%
    spread(key = taxon, value = abundance, fill = 0) %>%
    gather(key = "taxon", value = "abundance", -sample)
  
  # re-add relative abundances if they were present
  if (is.null(ta$abundances$rel_abundance)) {
    ta$abundances <- abundances_with_zeros
  } else {
    ta$abundances <- ta$abundances %>%
      right_join(abundances_with_zeros) %>%
      mutate(rel_abundance = ifelse(is.na(rel_abundance), 0, rel_abundance))    
  }
  
  # return ta object
  ta
  
}

add_rel_abundance <- function(ta) {
  
  # add relative abundance to abundance table
  ta$abundances <- ta$abundances %>%
    group_by(sample) %>%
    mutate(rel_abundance = abundance / sum(abundance)) %>%
    ungroup()
  
  # return ta object
  ta
  
}

add_max_rel_abundance <- function(ta) {
  
  # add relative abundances if not present
  remove_rel_abundance <- FALSE
  if (is.null(ta$abundances$rel_abundance)) {
    remove_rel_abundance <- TRUE
    ta <- add_rel_abundance(ta)
  }
  
  # make table with taxon and maximum relative abundance
  max_rel_abundances <- ta$abundances %>%
    group_by(taxon) %>%
    summarize(max_rel_abundance = max(rel_abundance))
  
  # add max relative abundance to taxon table
  ta$taxa <- left_join(ta$taxa, max_rel_abundances)
  
  # remove relative abundances
  if (remove_rel_abundance) ta$abundances$rel_abundance <- NULL
  
  # return ta object
  ta
  
}

add_total_rel_abundance <- function(ta) {
  
  # add relative abundances if not present
  remove_rel_abundance <- FALSE
  if (is.null(ta$abundances$rel_abundance)) {
    remove_rel_abundance <- TRUE
    ta <- add_rel_abundance(ta)
  }
  
  # make table with taxon and total relative abundance
  total_rel_abundances <- ta$abundances %>%
    group_by(taxon) %>%
    summarize(total_rel_abundance = sum(rel_abundance)/nrow(ta$samples))
  
  # add total relative abundance to taxon table
  ta$taxa <- left_join(ta$taxa, total_rel_abundances)
  
  # remove relative abundances
  if (remove_rel_abundance) ta$abundances$rel_abundance <- NULL
  
  # return ta object
  ta
  
}

# percentage of samples in which a taxon is present
# credits to Wenke Smets for the idea and initial implementation
add_rel_occurrence <- function(ta) {
  
  # make table with taxon and relative occurrence
  rel_occurrences <- ta$abundances %>%
    group_by(taxon) %>%
    summarize(occurrence = sum(abundance > 0)) %>%
    mutate(rel_occurrence = occurrence/nrow(ta$samples)) %>%
    select(- occurrence)
  
  # add relative occurrence to taxon table
  ta$taxa <- left_join(ta$taxa, rel_occurrences)
  
  # return ta object
  ta
  
}

add_taxon_name <- function(ta, method = "max_rel_abundance") {

  # throw error if method unknown
  if (! method %in% c("max_rel_abundance", "total_rel_abundance")) {
    stop("method unknown")
  }
  
  # make quosure of method
  var <- quo(get(method))

  # add max relative abundance to taxon table
  remove_max_rel_abundance <- FALSE
  if (is.null(ta$taxa$max_rel_abundance)) {
    remove_max_rel_abundance <- TRUE
    ta <- add_max_rel_abundance(ta)
  }
  
  # add total relative abundance to taxon table
  remove_total_rel_abundance <- FALSE
  if (is.null(ta$taxa$total_rel_abundance)) {
    remove_total_rel_abundance <- TRUE
    ta <- add_total_rel_abundance(ta)
  }

  # make version of taxon table with taxonomy levels in the right order
  tax_levels <- c("kingdom", "phylum", "class", "order", "family", "genus")
  taxa <- ta$taxa[, tax_levels]
  
  # make temporary taxon name: most specific level of taxonomy available
  taxon_name_temp <- apply(taxa, 1, FUN = function(row) {
    if(is.na(row["kingdom"])) return("unclassified")
    row %>% na.omit() %>% last() 
  })
  
  # add temporary taxon name to taxon table and add numbers if 
  # temporary name is not unique
  ta$taxa <- ta$taxa %>%
    mutate(taxon_name_temp = taxon_name_temp) %>%
    group_by(taxon_name_temp) %>%
    arrange(desc(!!var)) %>%
    mutate(n_taxa = n()) %>%
    mutate(taxon_number = ifelse(n_taxa > 1, as.character(1:n()), "")) %>%
    mutate(taxon_name = paste(taxon_name_temp, taxon_number, sep = " ")) %>%
    ungroup() %>%
    select(- taxon_name_temp, - n_taxa, - taxon_number)
  
  # remove max_rel_abundance
  if(remove_max_rel_abundance) ta$taxa$max_rel_abundance <- NULL
  
  # remove total_rel_abundance
  if(remove_total_rel_abundance) ta$taxa$total_rel_abundance <- NULL
  
  # return ta object
  ta
  
}

add_taxon_name_color <- function(ta, method = "max_rel_abundance", n = 12, samples = NULL, taxa = NULL) {
  
  # throw error if method unknown
  if (! method %in% c("max_rel_abundance", "total_rel_abundance")) {
    stop("method unknown")
  }
  
  # make quosure of method
  var <- quo(get(method))
  
  # add taxon_name if not present (needs to happen before subset)
  remove_taxon_name <- FALSE
  if (is.null(ta$taxa$taxon_name)) {
    remove_taxon_name <- TRUE
    ta <- add_taxon_name(ta)
  }
  
  # make subset of ta object (with selection of samples and/or taxa
  # if requested)
  ta_subset <- ta
  if (! is.null(samples)) {
    ta_subset$samples <- filter(ta_subset$samples, sample %in% samples)
    ta_subset <- process_sample_selection(ta_subset)
  }
  if (! is.null(taxa)) {
    ta_subset$taxa <- filter(ta_subset$taxa, taxon %in% taxa)
    ta_subset <- process_taxon_selection(ta_subset)
  }
  
  # add max relative abundance to taxon table if not present
  if (is.null(ta_subset$taxa$max_rel_abundance)) {
    ta_subset <- add_max_rel_abundance(ta_subset)
  }
  
  # add total relative abundance to taxon table
  if (is.null(ta_subset$taxa$total_rel_abundance)) {
    ta_subset <- add_total_rel_abundance(ta_subset)
  }
  
  # extract taxon names to visualize, in order
  levels <- ta_subset$taxa %>%
    arrange(desc(!!var)) %>%
    pull(taxon_name) %>%
    `[`(1:(n-1))
  levels <- levels[order(levels)]
  levels <- c( "residual", levels)
  
  # add taxon_name_color factor to taxa table
  ta$taxa <- ta$taxa %>%
    mutate(taxon_name_color = ifelse(taxon_name %in% levels, taxon_name, "residual")) %>%
    mutate(taxon_name_color = factor(taxon_name_color, levels = levels))
  
  # remove taxon_name if necessary
  if(remove_taxon_name) ta$taxa$taxon_name <- NULL
  
  # return ta object
  ta
  
}

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

# for internal use - do not export
get_rel_abundance_matrix <- function(ta) {
  
  # add relative abundances if not present
  if (is.null(ta$abundances$rel_abundance)) {
    ta <- add_rel_abundance(ta)
  }
  
  # convert abundance table to wide format
  abundances_wide <- ta$abundances %>%
    select(sample, taxon, rel_abundance) %>%
    spread(key = taxon, value = rel_abundance, fill = 0)
  
  # convert wide abundance table to relative abundance matrix
  # and return
  rel_abundance_matrix <- abundances_wide %>%
    select(- sample) %>%
    as.matrix() %>%
    `row.names<-`(abundances_wide$sample)
  
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

get_abundances_extended <- function(ta) {
  
  # make and return large table
  ta$abundances %>%
    left_join(ta$samples) %>%
    left_join(ta$taxa)
  
}

update_lib_sizes <- function(ta, step) {

  # current lib_size in tidy table
  lib_sizes_new <- ta$abundances %>%
    group_by(sample) %>%
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

palette_paired <- c(
  "#e8e8e8", # light grey
  "#a6cee3", # light blue
  "#1f78b4", # dark blue
  "#b2df8a", # light green
  "#33a02c", # dark green
  "#fb9a99", # light red
  "#e31a1c", # dark red
  "#fdbf6f", # light orange
  "#ff7f00", # dark orange
  "#cab2d6", # light purple
  "#6a3d9a", # dark purple
  "#ffff99", # light brown
  "#b15928"  # dark brown
)

get_bar_plot <- function(ta, x = sample_clustered) {
  
  # convert promise to formula
  x <- substitute(x)
  
  # add sample_clustered if not present
  if (is.null(ta$samples$sample_clustered)) {
    ta <- add_sample_clustered(ta)
  }
  
  # add taxon_name_color if not present
  if (is.null(ta$taxa$taxon_name_color)) {
    ta <- add_taxon_name_color(ta)
  }
  
  # make plot and return
  get_abundances_extended(ta) %>%
    ggplot(aes_(x = x, y = ~abundance, fill = ~taxon_name_color)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_brewer(palette = "Paired", name = "Taxon") +
    xlab("sample") + ylab("relative abundance") +
    theme(
      axis.text.x = element_text(angle = 90),
      axis.ticks.x = element_blank(),
      panel.background = element_rect(fill = 'white', colour = 'white')
    )
  
}

# needs lib_sizes table
get_history_plot <- function(ta, col = NULL) {
  
  # convert promise to formula
  col <- substitute(col)
  
  # remove lib_size if present
  ta$samples$lib_size <- NULL
  
  # make plot and return
  ta$lib_sizes %>%
    left_join(ta$samples) %>%
    ggplot(aes_(x = ~step, y = ~lib_size, group = ~sample, col = col)) +
    geom_line(size = 0.5) + 
    scale_y_log10() +
    scale_color_brewer(palette = "Paired") +
    theme(axis.text.x = element_text(angle = 90, vjust = 1),
          panel.background = element_rect(fill = "white", colour = "black"))
  
}

