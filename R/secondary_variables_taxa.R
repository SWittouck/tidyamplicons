
add_max_rel_abundance <- function(ta) {

  # if rel_abundance not present: add and remove on exit
  if (! "rel_abundance" %in% names(ta$abundances)) {
    ta <- add_rel_abundance(ta)
    on.exit(ta$abundances$rel_abundance <- NULL)
  }

  # make table with taxon and maximum relative abundance
  max_rel_abundances <- ta$abundances %>%
    group_by(taxon) %>%
    summarize(max_rel_abundance = max(rel_abundance))

  # add max relative abundance to taxon table
  ta$taxa <- left_join(ta$taxa, max_rel_abundances)

  # return ta object
  ta

}

add_total_rel_abundance <- function(ta) {

  # if rel_abundance not present: add and remove on exit
  if (! "rel_abundance" %in% names(ta$abundances)) {
    ta <- add_rel_abundance(ta)
    on.exit(ta$abundances$rel_abundance <- NULL)
  }

  # make table with taxon and total relative abundance
  total_rel_abundances <- ta$abundances %>%
    group_by(taxon) %>%
    summarize(total_rel_abundance = sum(rel_abundance)/nrow(ta$samples))

  # add total relative abundance to taxon table
  ta$taxa <- left_join(ta$taxa, total_rel_abundances)

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

add_taxon_name <- function(ta, method = "max_rel_abundance", include_species = F) {

  # throw error if method unknown
  if (! method %in% c("max_rel_abundance", "total_rel_abundance")) {
    stop("method unknown")
  }

  # make quosure of method
  var <- quo(get(method))

  # if max_rel_abundance not present: add and remove on exit
  if (! "max_rel_abundance" %in% names(ta$taxa)) {
    ta <- add_max_rel_abundance(ta)
    on.exit(ta$taxa$max_rel_abundance <- NULL, add = T)
  }

  # if total_rel_abundance not present: add and remove on exit
  if (! "total_rel_abundance" %in% names(ta$taxa)) {
    ta <- add_total_rel_abundance(ta)
    on.exit(ta$taxa$total_rel_abundance <- NULL, add = T)
  }

  # make version of taxon table with taxonomy levels in the right order
  tax_levels <- c("kingdom", "phylum", "class", "order", "family", "genus")
  if (include_species) tax_levels <- c(tax_levels, "species")
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

  # if taxon_name not present: add and remove on exit
  if (! "taxon_name" %in% names(ta$taxa)) {
    ta <- add_taxon_name(ta)
    on.exit(ta$taxa$taxon_name <- NULL, add = T)
  }

  ta_subset <- ta

  # take subset of samples if requested
  if (! is.null(samples)) {
    ta_subset$samples <- filter(ta_subset$samples, sample %in% samples)
    ta_subset <- process_sample_selection(ta_subset)
  }

  # if max_rel_abundance not present: add and remove on exit
  if (! "max_rel_abundance" %in% names(ta_subset$taxa)) {
    ta_subset <- add_max_rel_abundance(ta_subset)
  }

  # if total_rel_abundance not present: add and remove on exit
  if (! "total_rel_abundance" %in% names(ta_subset$taxa)) {
    ta_subset <- add_total_rel_abundance(ta_subset)
  }

  # take subset of taxa if requested
  if (! is.null(taxa)) {
    ta_subset$taxa <- filter(ta_subset$taxa, taxon %in% taxa)
    ta_subset <- process_taxon_selection(ta_subset)
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

  # return ta object
  ta

}

add_jervis_bardy <- function(ta, min_pres = 3) {

  # if rel_abundance not present: add and remove on exit
  if (! "rel_abundance" %in% names(ta$abundances)) {
    ta <- add_rel_abundance(ta)
    on.exit(ta$abundances$rel_abundance <- NULL)
  }

  # if lib_size not present: add and remove on exit
  if (! "lib_size" %in% names(ta$samples)) {
    ta <- add_lib_size(ta)
    on.exit(ta$samples$lib_size <- NULL, add = T)
  }

  # perform jervis bardy calculation
  taxa_jb <- ta$abundances %>%
    left_join(ta$samples %>% select(sample, lib_size)) %>%
    group_by(taxon) %>%
    filter(n() >= !! min_pres) %>%
    do(jb = cor.test(x = .$rel_abundance, y = .$lib_size, alternative = "less", method = "spearman")) %>%
    mutate(jb_cor = jb$estimate, jb_p = jb$p.value) %>%
    select(- jb)

  # add jb_p and jb_cor to taxa table
  ta$taxa <- left_join(ta$taxa, taxa_jb)

  # return ta object
  ta

}

# Adds taxon presence and absence counts in sample conditions to the taxa table,
# as well as a fisher exact test for differential presence.
# Condition is a variable that should be present in the samples table.
add_presence_counts <- function(ta, condition) {

  condition <- enquo(condition)

  counts_tidy <- taxon_counts_in_conditions(ta, !! condition)

  taxa_counts <- counts_tidy %>%
    mutate(presence_in_condition = str_c(presence, !! condition, sep = "_in_")) %>%
    select(taxon, presence_in_condition, n) %>%
    spread(value = n, key = presence_in_condition)

  taxa_fisher <- counts_tidy %>%
    group_by(taxon) %>%
    arrange(!! condition, presence) %>%
    do(fisher = c(.$n) %>%
         matrix(ncol = 2, byrow = T) %>%
         fisher.test()
    ) %>%
    mutate(fisher_p = fisher$p.value) %>%
    select(- fisher)

  ta %>%
    modify_at("taxa", left_join, taxa_counts) %>%
    modify_at("taxa", left_join, taxa_fisher)

}
