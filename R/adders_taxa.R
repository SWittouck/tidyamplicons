#' (Re)classify amplicon sequences
#'
#' This function requires the DADA2 package to be installed.
#'
#' This function will (re)classify either all or a subset of the taxa, given
#' that a variable is present in the taxon table that contains (representative)
#' sequences of the taxa.
#'
#' Ranks can be supplied as a named integer vector, where the names represent
#' taxonomic ranks and the integers represent positions of these ranks in the
#' taxonomy strings present in the reference database. Ranks can also be
#' supplied as just a character vector with the rank names; in that case, it is
#' assumed that the database taxonomy string follows the default order {domain,
#' phylum, class, order, family, genus, species}. If no ranks are supplied, taxa
#' will be (re)classified at all default ranks.
#'
#' @param ta A tidyamplicons object.
#' @param refdb The path to a DADA2-compatible reference database.
#' @param taxa An expression that specifies which taxa to (re)classify.
#' @param ranks A vector that specifies which ranks to (re)classify.
#' @param sequence_var The (quoted) name of a variable within the taxa table
#'   that contains (representative) sequences of the taxa.
#' @param multithread A boolean indicating whether to use multiple threads.
#' @param min_boot The minimum bootstrap value for taxonomy assignment.
#' @param n_ranks The number of ranks present in the reference database.
#'
#' @return An updated tidyamplicons object.
#'
#' @export
classify_taxa <- function(
  ta, refdb, taxa = rep(T, times = length(taxon_id)), ranks = "default",
  sequence_var = "sequence", multithread = T, min_boot = 50, n_ranks = 7
) {

  # throw error if sequence_var doesn't exist
  if (! sequence_var %in% names(ta$taxa)) {
    stop(paste0("variable '", sequence_var, "' not found in taxon table"))
  }

  # convert taxa argument to an expression
  taxa <- rlang::enexpr(taxa)

  # determine which taxa to (re)classify
  to_reclassify <- eval(taxa, ta$taxa)
  to_reclassify <- to_reclassify & ! is.na(to_reclassify)

  # throw error if no taxa to (re)classify
  if (sum(to_reclassify) == 0) stop("zero taxa obey the given criteria")

  # lookup default ranks and/or rank numbers if not given
  if (! is.numeric(ranks)) {
    default_ranks <- 1:7
    names(default_ranks) <-
      c("domain", "phylum", "class", "order", "family", "genus", "species")
    if (ranks[1] == "default") ranks <- names(default_ranks)
    ranks <- default_ranks[ranks]
    ranks <- ranks[! is.na(ranks)]
  }

  # throw error if a rank number exceeds the number of ranks in the db
  if (max(ranks) > n_ranks) stop("not enough ranks in the database")

  # perform (re)classification
  # tryRC is important; see <https://github.com/benjjneb/dada2/issues/1441>
  taxa_reclassified <-
    ta$taxa[to_reclassify, ][[sequence_var]] %>%
    dada2::assignTaxonomy(
      refFasta = refdb, multithread = multithread, minBoot = min_boot,
      tryRC = T, taxLevels = letters[1:n_ranks]
    )
  ta$taxa[to_reclassify, names(ranks)] <- taxa_reclassified[ , ranks]

  ta

}

#' Add taxon metadata to the tidyamplicons object
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
#' @export
add_taxon_tibble <- function(ta, taxon_tibble) {

  modify_at(ta, "taxa", left_join, taxon_tibble)

}

#' Add the maximum relative abundance of taxa to the taxon table
#'
#' @export
add_max_rel_abundance <- function(ta) {

  # if rel_abundance not present: add temporarily
  rel_abundance_tmp <- ! "rel_abundance" %in% names(ta$abundances)
  if (rel_abundance_tmp) ta <- add_rel_abundance(ta)

  # make table with taxon and maximum relative abundance
  max_rel_abundances <- ta$abundances %>%
    group_by(taxon_id) %>%
    summarize(max_rel_abundance = max(rel_abundance))

  # add max relative abundance to taxon table
  ta$taxa <- left_join(ta$taxa, max_rel_abundances, by = "taxon_id")

  # cleanup
  if (rel_abundance_tmp) ta$abundances$rel_abundance <- NULL

  # return ta object
  ta

}

#' Add the total relative abundance of taxa to the taxon table
#'
#' @export
add_total_rel_abundance <- function(ta) {

  # if rel_abundance not present: add temporarily
  rel_abundance_tmp <- ! "rel_abundance" %in% names(ta$abundances)
  if (rel_abundance_tmp) ta <- add_rel_abundance(ta)

  # make table with taxon and total relative abundance
  total_rel_abundances <- ta$abundances %>%
    group_by(taxon_id) %>%
    summarize(total_rel_abundance = sum(rel_abundance)/nrow(ta$samples))

  # add total relative abundance to taxon table
  ta$taxa <- left_join(ta$taxa, total_rel_abundances, by = "taxon_id")

  # cleanup
  if (rel_abundance_tmp) ta$abundances$rel_abundance <- NULL

  # return ta object
  ta

}

#' Add the relative occurrence of taxa to the taxon table
#'
#' DEPRECATED, use \code{\link{add_occurrences}}
#'
#' Credits to Wenke Smets for the idea and initial implementation.
#'
#' @export
add_rel_occurrence <- function(ta) {

  # make table with taxon and relative occurrence
  rel_occurrences <- ta$abundances %>%
    group_by(taxon_id) %>%
    summarize(occurrence = sum(abundance > 0)) %>%
    mutate(rel_occurrence = occurrence/nrow(ta$samples)) %>%
    select(- occurrence)

  # add relative occurrence to taxon table
  ta$taxa <- left_join(ta$taxa, rel_occurrences, by = "taxon_id")

  # return ta object
  ta

}

#' Create sensible names for the taxa and add to taxon table
#'
#' @export
add_taxon_name <- function(
  ta, method = "total_rel_abundance", include_species = F
  ) {

  if (method == "total_rel_abundance") {

    # if total_rel_abundance not present: add temporarily
    total_rel_ab_tmp <- ! "total_rel_abundance" %in% names(ta$taxa)
    if (total_rel_ab_tmp) ta <- add_total_rel_abundance(ta)

    ta <- mutate_taxa(ta, arrange_by_me = total_rel_abundance)

  } else if (method == "max_rel_abundance") {

    # if max_rel_abundance not present: add temporarily
    max_rel_ab_tmp <- ! "max_rel_abundance" %in% names(ta$taxa)
    if (max_rel_ab_tmp) ta <- add_max_rel_abundance(ta)

    ta <- mutate_taxa(ta, arrange_by_me = max_rel_abundance)

  } else {

    # throw error if method unknown
    if (! method %in% c("total_rel_abundance", "max_rel_abundance")) {
      stop("method unknown")
    }

  }

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

  # cleanup
  if (exists("total_rel_ab_tmp")) ta$taxa$total_rel_abundance <- NULL
  if (exists("max_rel_ab_tmp")) ta$taxa$max_rel_abundance <- NULL
  ta$taxa$arrange_by_me <- NULL

  # return ta object
  ta

}

#' Create taxon names suitable for visualization with color
#'
#' @export
add_taxon_name_color <- function(
  ta, method = "total_rel_abundance", n = 12, samples = NULL, taxa = NULL
  ) {

  # if taxon_name not present: add temporarily
  taxon_name_tmp <- ! "taxon_name" %in% names(ta$taxa)
  if (taxon_name_tmp) ta <- add_taxon_name(ta)

  if (method == "total_rel_abundance") {

    # if total_rel_abundance not present: add temporarily
    total_rel_ab_tmp <- ! "total_rel_abundance" %in% names(ta$taxa)
    if (total_rel_ab_tmp) ta <- add_total_rel_abundance(ta)

    ta <- mutate_taxa(ta, arrange_by_me = total_rel_abundance)

  } else if (method == "max_rel_abundance") {

    # if max_rel_abundance not present: add temporarily
    max_rel_ab_tmp <- ! "max_rel_abundance" %in% names(ta$taxa)
    if (max_rel_ab_tmp) ta <- add_max_rel_abundance(ta)

    ta <- mutate_taxa(ta, arrange_by_me = max_rel_abundance)

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

  # cleanup
  if (taxon_name_tmp) ta$taxa$taxon_name <- NULL
  if (exists("total_rel_ab_tmp")) ta$taxa$total_rel_abundance <- NULL
  if (exists("max_rel_ab_tmp")) ta$taxa$max_rel_abundance <- NULL
  ta$taxa$arrange_by_me <- NULL

  # return ta object
  ta

}

#' Apply the taxon QC method of Jervis-Bardy
#'
#' Function to estimate spearman correlation between relative abundance and
#' sample dna concentration, for each taxon.
#'
#' See:
#' J. Jervis-Bardy et al., “Deriving accurate microbiota profiles from
#' human samples with low bacterial content through post-sequencing processing
#' of Illumina MiSeq data,” Microbiome, vol. 3, no. 1, Art. no. 1, 2015, doi:
#' 10.1186/s40168-015-0083-8.
#'
#' @param ta A tidyamplicons object.
#' @param dna_conc A variable in the samples table that contains dna
#'   concetrations (unquoted).
#' @param sample_condition An optional extra condition that samples must pass
#'   before calculations.
#' @param min_pres The minimum number of samples a taxon has to be present in
#'   for its correlation to be calculated.
#'
#' @export
add_jervis_bardy <- function(ta, dna_conc, sample_condition = T, min_pres = 3) {

  dna_conc <- enquo(dna_conc)
  sample_condition <- enquo(sample_condition)

  # if rel_abundance not present: add temporarily
  rel_abundance_tmp <- ! "rel_abundance" %in% names(ta$abundances)
  if (rel_abundance_tmp) ta <- add_rel_abundance(ta)

  # if sample condition is given, use only samples that fulfill it
  if (is.null(sample_condition)) {
    ta_jb <- ta
  } else {
    ta_jb <- ta %>%
      filter_samples(!! sample_condition)
  }

  # perform jervis bardy calculation
  taxa_jb <- ta_jb$abundances %>%
    left_join(
      ta_jb$samples %>% select(sample_id, dna_conc = !! dna_conc),
      by = "sample_id"
    ) %>%
    group_by(taxon_id) %>%
    filter(n() >= !! min_pres) %>%
    do(
      jb = cor.test(
        x = .$rel_abundance, y = .$dna_conc, alternative = "less",
        method = "spearman"
      )
    ) %>%
    mutate(jb_cor = jb$estimate, jb_p = jb$p.value) %>%
    select(- jb)

  # add jb_p and jb_cor to taxa table
  ta$taxa <- left_join(ta$taxa, taxa_jb, by = "taxon_id")

  # cleanup
  if (rel_abundance_tmp) ta$abundances$rel_abundance <- NULL

  # return ta object
  ta

}

#' Add absolute occurrences of taxa to the taxon table
#'
#' Adds taxon presence and absence counts in sample conditions to the taxa
#' table, as well as a fisher exact test for differential presence.
#'
#' Condition is a variable that should be present in the samples table.
#'
#' DEPRECATED, use \code{\link{add_occurrences}}
#'
#' @export
add_presence_counts <- function(ta, condition) {

  condition <- enquo(condition)

  counts_tidy <- taxon_counts_in_conditions(ta, !! condition)

  taxa_counts <- counts_tidy %>%
    mutate(
      presence_in_condition = str_c(presence, !! condition, sep = "_in_")
    ) %>%
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
    modify_at("taxa", left_join, taxa_counts, by = "taxon_id") %>%
    modify_at("taxa", left_join, taxa_fisher, by = "taxon_id")

}

#' Add taxon occurrences to the taxon table
#'
#' Adds taxon occurrences (overall or per condition) to the taxa table.
#'
#' Condition should be a categorical variable present in the samples table.
#' Supply condition as a string.
#'
#' @export
add_occurrences <- function(
  ta, condition = NULL, relative = F, fischer_test = F
  ) {

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
      left_join(taxa_fischer, by = "taxon_id")

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
    modify_at("taxa", left_join, taxa_occurrences, by = "taxon_id")

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
#'
#' @export
add_mean_rel_abundances <- function(ta, condition = NULL, test = NULL) {

  mean_rel_abundances <- mean_rel_abundances(ta, condition = condition)

  if (is.null(condition)) {

    taxa_mean_rel_abundances <- mean_rel_abundances

  } else if (! is.null(test)) {

    condition_sym <- ensym(condition)

    # if rel_abundance not present: add temporarily
    rel_abundance_tmp <- ! "rel_abundance" %in% names(ta$abundances)
    if (rel_abundance_tmp) ta <- add_rel_abundance(ta)

    rel_abundances_complete <-
      abundances(ta) %>%
      left_join(ta$samples, by = "sample_id") %>%
      select(sample_id, !! condition_sym, taxon_id, rel_abundance) %>%
      complete(
        nesting(sample_id, !! condition_sym), taxon_id,
        fill = list(rel_abundance = 0)
      )

    if (test == "wilcox") {
      taxa_test <-
        rel_abundances_complete %>%
        group_by(taxon_id) %>%
        do(result = wilcox.test(
          data = ., rel_abundance ~ !! condition_sym
        )) %>%
        mutate(wilcox_p = result$p.value, wilcox_stat = result$statistic) %>%
        select(- result)
    } else if (test == "t-test") {
      taxa_test <-
        rel_abundances_complete %>%
        group_by(taxon_id) %>%
        do(result = t.test(
          data = ., rel_abundance ~ !! condition_sym
        )) %>%
        mutate(t_test_p = result$p.value, t_test_stat = result$statistic) %>%
        select(- result)
    } else {
      stop("please supply a valid test")
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

  # cleanup
  if (exists("rel_abundance_tmp")) ta$abundances$rel_abundance <- NULL

  ta %>%
    modify_at("taxa", left_join, taxa_mean_rel_abundances, by = "taxon_id")

}
