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

# DEPRECATED: use add_occurrences()
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

add_taxon_name <- function(ta, method = "total_rel_abundance", include_species = F) {

  if (method == "total_rel_abundance") {

    # if total_rel_abundance not present: add and remove on exit
    if (! "total_rel_abundance" %in% names(ta$taxa)) {
      ta <- add_total_rel_abundance(ta)
      on.exit(ta$taxa$total_rel_abundance <- NULL, add = T)
    }

    ta <- mutate_taxa(ta, arrange_by_me = total_rel_abundance)

  } else if (method == "max_rel_abundance") {

    # if total_rel_abundance not present: add and remove on exit
    if (! "max_rel_abundance" %in% names(ta$taxa)) {
      ta <- add_max_rel_abundance(ta)
      on.exit(ta$taxa$max_rel_abundance <- NULL, add = T)
    }

    ta <- mutate_taxa(arrange_by_me = max_rel_abundance)

  } else {

    # throw error if method unknown
    if (! method %in% c("total_rel_abundance", "max_rel_abundance")) {
      stop("method unknown")
    }

  }

  on.exit(ta$taxa$arrange_by_me <- NULL, add = T)

  rank_names <-
    rank_names(ta) %>%
    purrr::when(include_species ~ c(., "species"), ~ .) %>%
    intersect(names(ta$taxa))

  if (length(rank_names) == 0) {

    ta <- mutate_taxa(ta, best_classification = "unclassified")

  } else {

    ta$taxa <-
      ta$taxa %>%
      mutate(
        best_classification =
          pmap_chr(
            ta$taxa[, rank_names],
            function(...) {
              classification = as.character(list(...))
              if (all(is.na(classification))) return("unclassified")
              classification %>% na.omit() %>% last()
            }
          )
      )

  }

  ta$taxa <-
    ta$taxa %>%
    group_by(best_classification) %>%
    arrange(desc(arrange_by_me)) %>%
    mutate(n_taxa = n()) %>%
    mutate(taxon_number = if_else(n_taxa > 1, as.character(1:n()), "")) %>%
    mutate(taxon_name = str_c(best_classification, taxon_number, sep = " ")) %>%
    mutate_at("taxon_name", str_trim) %>%
    ungroup() %>%
    select(- best_classification, - n_taxa, - taxon_number)

  # return ta object
  ta

}

add_taxon_name_color <- function(ta, method = "total_rel_abundance", n = 12, samples = NULL, taxa = NULL) {

  # if taxon_name not present: add and remove on exit
  if (! "taxon_name" %in% names(ta$taxa)) {
    ta <- add_taxon_name(ta)
    on.exit(ta$taxa$taxon_name <- NULL, add = T)
  }

  if (method == "total_rel_abundance") {

    # if total_rel_abundance not present: add and remove on exit
    if (! "total_rel_abundance" %in% names(ta$taxa)) {
      ta <- add_total_rel_abundance(ta)
      on.exit(ta$taxa$total_rel_abundance <- NULL, add = T)
    }

    ta <- mutate_taxa(ta, arrange_by_me = total_rel_abundance)

  } else if (method == "max_rel_abundance") {

    # if total_rel_abundance not present: add and remove on exit
    if (! "max_rel_abundance" %in% names(ta$taxa)) {
      ta <- add_max_rel_abundance(ta)
      on.exit(ta$taxa$max_rel_abundance <- NULL, add = T)
    }

    ta <- mutate_taxa(arrange_by_me = max_rel_abundance)

  } else {

    # throw error if method unknown
    if (! method %in% c("total_rel_abundance", "max_rel_abundance")) {
      stop("method unknown")
    }

  }

  ta_subset <- ta

  # take subset of samples if requested
  if (! is.null(samples)) {
    ta_subset <- filter_samples(ta_subset, sample_id %in% !! samples)
  }

  # take subset of taxa if requested
  if (! is.null(taxa)) {
    ta_subset <- filter_taxa(ta_subset, taxon_id %in% !! taxa)
  }

  # extract taxon names to visualize, in order
  levels <-
    ta_subset$taxa %>%
    arrange(desc(arrange_by_me)) %>%
    pull(taxon_name) %>%
    `[`(1:(n-1)) %>%
    sort() %>%
    purrr::prepend("residual")

  # add taxon_name_color factor to taxa table
  ta$taxa <-
    ta$taxa %>%
    mutate(taxon_name_color = if_else(taxon_name %in% levels, taxon_name, "residual")) %>%
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

# DEPRECATED: use add_occurrences()
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

    taxa_occurrences <-
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

#' Add average relative abundances
#'
#' This function adds mean relative abundance values for each taxon to the taxa
#' table, overall or per sample group.
#'
#' If `condition` is specified, the mean relative abundances will be calculated
#' separately for each group defined by the condition variable. This variable
#' should be present in the sample table.
#'
#' If `condition` is specified, differential abundance testing can be performed
#' by setting the `test` argument. Options are NULL (default), "wilcox" or
#' "t-test".
#'
#' @param ta A tidyamplicons object
#' @param condition A condition variable (character)
#' @param test Differential abundance test to perform
#'
#' @return A tidyamplicons object
add_mean_rel_abundances <- function(ta, condition = NULL, test = NULL) {

  mean_rel_abundances <- mean_rel_abundances(ta, condition = condition)

  if (is.null(condition)) {

    taxa_mean_rel_abundances <- mean_rel_abundances

  } else if (! is.null(test)) {

    condition_sym <- ensym(condition)

    # if rel_abundance not present: add and remove on exit
    if (! "rel_abundance" %in% names(ta$abundances)) {
      ta <- add_rel_abundance(ta)
      on.exit(ta$abundances$rel_abundance <- NULL)
    }

    taxa_test <-
      abundances(ta) %>%
      left_join(ta$samples) %>%
      complete(
        nesting(sample_id, !! condition_sym), taxon_id,
        fill = list(rel_abundance = 0)
      ) %>%
      group_by(taxon_id)

    if (test == "wilcox") {
      taxa_test <-
        taxa_test %>%
        do(result = wilcox.test(
          data = ., rel_abundance ~ !! condition_sym
        )) %>%
        mutate(wilcox_p = result$p.value, wilcox_stat = result$statistic) %>%
        select(- result)
    } else if (test == "t-test") {
      taxa_test <-
        taxa_test  %>%
        do(result = t.test(
          data = ., rel_abundance ~ !! condition_sym
        )) %>%
        mutate(t_test_p = result$p.value, t_test_stat = result$statistic) %>%
        select(- result)
    }

    taxa_mean_rel_abundances <-
      mean_rel_abundances %>%
      mutate_at(condition, ~ str_c("mean_rel_abundance_in", ., sep = "_")) %>%
      spread(value = mean_rel_abundance, key = condition) %>%
      left_join(taxa_test, by = "taxon_id")

  } else {

    taxa_mean_rel_abundances <-
      mean_rel_abundances %>%
      mutate_at(condition, ~ str_c("mean_rel_abundance_in", ., sep = "_")) %>%
      spread(value = mean_rel_abundance, key = condition)

  }

  ta %>%
    modify_at("taxa", left_join, taxa_mean_rel_abundances)

}
