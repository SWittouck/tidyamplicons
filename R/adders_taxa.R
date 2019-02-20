#' Add taxa table to the tidyamplicons object
#'
#' \code{add_taxon_tibble} adds a taxon tibble to the tidyamplicons object.
#'
#' This function adds a taxon tibble containing taxon data (e.g. taxon ranks
#' such as genus, family, ...) for each taxon to the tidyamplicons object. It is
#' used after initiating a tidyamplicons object using a numerical abundance
#' matrix and the function \code{\link{create_tidyamplicons}}. Also see
#' \code{\link{add_sample_tibble}} to update the sample data of the
#' tidyamplicons object.
#'
#' @param ta Tidyamplicons object.
#' @param taxon_tibble A tibble containing taxon data for each taxon. Taxa
#'   should be rows, while taxon data should be columns. At least one column
#'   name needs to be shared with the taxon tibble of ta. The default shared
#'   column name is 'taxon'.
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
#' # Convert to tidyamplicons object
#' data <- create_tidyamplicons(x,
#'                      taxa_are_columns = FALSE,
#'                      taxon_names_are_sequences = FALSE)
#'
#' # Initiate taxon tibble
#' taxon <- c("taxon1", "taxon2")
#' genus <- c("Salmonella", "Lactobacillus")
#' taxon_tibble <- tibble(taxa, genus)
#'
#' # Add taxon tibble to tidyamplicons object
#' data <- data %>%
#' add_taxon_tibble(taxon_tibble)
#'


add_taxon_tibble <- function(ta, taxon_tibble) {

  modify_at(ta, "taxa", left_join, taxon_tibble)

}

add_max_rel_abundance <- function(ta) {

  # if rel_abundance not present: add and remove on exit
  if (! "rel_abundance" %in% names(ta$abundances)) {
    ta <- add_rel_abundance(ta)
    on.exit(ta$abundances$rel_abundance <- NULL)
  }

  # make table with taxon and maximum relative abundance
  max_rel_abundances <- ta$abundances %>%
    group_by(taxon_id) %>%
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
    group_by(taxon_id) %>%
    summarize(total_rel_abundance = sum(rel_abundance)/nrow(ta$samples))

  # add total relative abundance to taxon table
  ta$taxa <- left_join(ta$taxa, total_rel_abundances)

  # return ta object
  ta

}

# DEPRICATED: use add_occurrences()
# percentage of samples in which a taxon is present
# credits to Wenke Smets for the idea and initial implementation
add_rel_occurrence <- function(ta) {

  # make table with taxon and relative occurrence
  rel_occurrences <- ta$abundances %>%
    group_by(taxon_id) %>%
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
  tax_levels <- tax_levels[tax_levels %in% names(ta$taxa)]
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
    ta_subset$samples <- filter(ta_subset$samples, sample_id %in% samples)
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
    ta_subset$taxa <- filter(ta_subset$taxa, taxon_id %in% taxa)
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

# Function to estimate spearman correlation between relative abundance and sample dna concentration,
# for each taxon.
# Inputs:
#   - ta: tidyamplicons object
#   - dna_conc: variable in the samples table that contains dna concetrations (unquoted)
#   - sample condition: optional extra condition that samples must pass before calculations
#   - min_pres: minimum number of samples a taxon has to be present in for its correlation to be calculated
add_jervis_bardy <- function(ta, dna_conc, sample_condition = T, min_pres = 3) {

  dna_conc <- enquo(dna_conc)
  sample_condition <- enquo(sample_condition)

  # if rel_abundance not present: add and remove on exit
  if (! "rel_abundance" %in% names(ta$abundances)) {
    ta <- add_rel_abundance(ta)
    on.exit(ta$abundances$rel_abundance <- NULL)
  }

  # if sample condition is given, use only samples that fulfill it
  if (is.null(sample_condition)) {
    ta_jb <- ta
  } else {
    ta_jb <- ta %>%
      filter_samples(!! sample_condition)
  }

  # perform jervis bardy calculation
  taxa_jb <- ta_jb$abundances %>%
    left_join(ta_jb$samples %>% select(sample_id, dna_conc = !! dna_conc)) %>%
    group_by(taxon_id) %>%
    filter(n() >= !! min_pres) %>%
    do(jb = cor.test(x = .$rel_abundance, y = .$dna_conc, alternative = "less", method = "spearman")) %>%
    mutate(jb_cor = jb$estimate, jb_p = jb$p.value) %>%
    select(- jb)

  # add jb_p and jb_cor to taxa table
  ta$taxa <- left_join(ta$taxa, taxa_jb)

  # return ta object
  ta

}

# DEPRICATED: use add_occurrences()
# Adds taxon presence and absence counts in sample conditions to the taxa table,
# as well as a fisher exact test for differential presence.
# Condition is a variable that should be present in the samples table.
add_presence_counts <- function(ta, condition) {

  condition <- enquo(condition)

  counts_tidy <- taxon_counts_in_conditions(ta, !! condition)

  taxa_counts <- counts_tidy %>%
    mutate(presence_in_condition = str_c(presence, !! condition, sep = "_in_")) %>%
    select(taxon_id, presence_in_condition, n) %>%
    spread(value = n, key = presence_in_condition)

  taxa_fisher <- counts_tidy %>%
    group_by(taxon_id) %>%
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

# Adds taxon occurrences (overall or per condition) to the taxa table.
# Condition should be a categorical variable present in the samples table.
# Supply condition as a string.
add_occurrences <- function(ta, condition = NULL, relative = F, fischer_test = F) {

  if (is.null(condition)) {

    taxa_occurrences <-
      occurrences(ta, condition = condition)

  } else if (fischer_test) {

    occurrences <-
      occurrences(ta, condition = condition, pres_abs = T)

    condition_sym <- sym(condition)

    taxa_fischer <-
      occurrences %>%
      group_by(taxon_id) %>%
      arrange(!! condition_sym, presence) %>%
      do(
        fisher = c(.$n) %>%
          matrix(ncol = 2, byrow = T) %>%
          fisher.test()
      ) %>%
      mutate(fisher_p = fisher$p.value) %>%
      select(- fisher)

    taxa_occurrences <-
      occurrences %>%
      filter(presence == "present") %>%
      select(taxon_id, !! condition_sym, occurrence = n) %>%
      mutate_at(condition, ~ str_c("occurrence_in", ., sep = "_")) %>%
      spread(value = "occurrence", key = condition) %>%
      left_join(taxa_fischer)

  } else {

    taxa_occurrences <-
      occurrences(ta, condition = condition) %>%
      mutate_at(condition, ~ str_c("occurrence_in", ., sep = "_")) %>%
      spread(value = "occurrence", key = condition)

  }

  if (relative & is.null(condition)) {

    taxa_occurrences %>%
      mutate(occurrence = occurrence / nrow(ta$samples))

  }

  if (relative & ! is.null(condition)) {

    condition_sym <- sym(condition)

    conditions <-
      ta$samples %>%
      count(!! condition_sym)

    for (con_ix in 1:nrow(conditions)) {

      con <- conditions[[condition]][con_ix]
      n_samples <- conditions[["n"]][con_ix]
      taxa_occurrences[[str_c("occurrence_in_", con)]] <-
        taxa_occurrences[[str_c("occurrence_in_", con)]] / n_samples

    }

  }

  ta %>%
    modify_at("taxa", left_join, taxa_occurrences)

}

# Adds taxon average relative abundances (overall or per condition) to the taxa
# table.
# Condition should be a categorical variable present in the samples table.
# Supply condition as a string.
add_mean_rel_abundances <- function(ta, condition = NULL, t_test = F) {

  mean_rel_abundances <- mean_rel_abundances(ta, condition = condition)

  if (is.null(condition)) {

    taxa_mean_rel_abundances <- mean_rel_abundances

  } else if(t_test) {

    condition_sym <- ensym(condition)

    # if rel_abundance not present: add and remove on exit
    if (! "rel_abundance" %in% names(ta$abundances)) {
      ta <- add_rel_abundance(ta)
      on.exit(ta$abundances$rel_abundance <- NULL)
    }

    taxa_t_tests <-
      abundances(ta) %>%
      left_join(ta$samples) %>%
      complete(nesting(sample_id, !! condition_sym), taxon_id, fill = list(rel_abundance = 0)) %>%
      group_by(taxon_id) %>%
      do(t_test = t.test(data = ., rel_abundance ~ !! condition_sym)) %>%
      mutate(t_test_p = t_test$p.value, t_test_t = t_test$statistic) %>%
      select(- t_test)

    taxa_mean_rel_abundances <-
      mean_rel_abundances %>%
      mutate_at(condition, ~ str_c("mean_rel_abundance_in", ., sep = "_")) %>%
      spread(value = mean_rel_abundance, key = condition) %>%
      left_join(taxa_t_tests)

  } else {

    taxa_mean_rel_abundances <-
      mean_rel_abundances %>%
      mutate_at(condition, ~ str_c("mean_rel_abundance_in", ., sep = "_")) %>%
      spread(value = mean_rel_abundance, key = condition)

  }

  ta %>%
    modify_at("taxa", left_join, taxa_mean_rel_abundances)

}