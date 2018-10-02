# Performs adonis function of the vegan package and returns output.
# Predictors should be a character vector.
# Samples where one of the predictors is NA are removed.
perform_adonis <- function(ta, predictors, permutations = 999) {

  abundances_matrix <- ta %>%
    modify_at("samples", drop_na, one_of(predictors)) %>%
    process_sample_selection() %>%
    add_rel_abundance() %>%
    abundances() %>%
    as_abundances_matrix(value = "rel_abundance")

  formula_RHS <- paste0(predictors, collapse = " + ")

  metadata <- tibble(sample = rownames(abundances_matrix)) %>%
    left_join(ta$samples)

  adonis(
    as.formula(paste("abundances_matrix", formula_RHS, sep = " ~ ")),
    metadata,
    permutations = permutations
  )

}
