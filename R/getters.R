#' Report numbers
#'
#' \code{report_numbers} returns the number of samples, taxa and reads in a
#' tidyamplicons object.
#'
#' This function prints the number of samples, taxa and reads in a tidyamplicons
#' object. To retrieve the numbers stored in named numeric vector, use
#' \code{\link{get_numbers}} instead.
#'
#' @param ta Tidyamplicons object.
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
#'                      taxa_are_columns = FALSE
#'                      )
#' # Report numbers
#' data %>%
#'  report_numbers()

report_numbers <- function(ta) {

  sprintf("samples: %i", nrow(ta$samples)) %>% message()
  sprintf("taxa: %i", nrow(ta$taxa)) %>% message()
  sprintf("reads: %i", sum(ta$abundances$abundance)) %>% message()

}

#' Get numbers
#'
#' \code{get_numbers} returns the number of samples, taxa and reads in a
#' tidyamplicons object.
#'
#' This function returns the number of samples, taxa and reads in a tidyamplicons
#' object, stored in a named numeric vector. To print the numbers, use
#' \code{\link{report_numbers}} instead.
#'
#' @param ta Tidyamplicons object.
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
#'                      taxa_are_columns = FALSE
#'                      )
#' # Report numbers
#' numbers <- data %>%
#'  get_numbers()
get_numbers <- function(ta) {

  c(
    n_samples = nrow(ta$samples),
    n_taxa = nrow(ta$taxa),
    n_reads = sum(ta$abundances$abundance)
  )

}

# for internal use - do not export
get_rel_abundance_matrix <- function(ta) {

  # add relative abundances if not present
  if (! "rel_abundance" %in% names(ta$abundances)) {
    ta <- add_rel_abundance(ta)
  }

  as_abundances_matrix(ta$abundances, value = rel_abundance)

}

#' Get betas
#'
#' \code{betas} returns a tidy tibble with the beta diversity for each
#' combination of samples.
#'
#' This function calculates the beta diversity using the
#' \code{\link[vegan]{vegdist}} function of Vegan. It will report one diversity
#' estimate for each combination of samples.
#'
#'
#' @param ta Tidyamplicons object.
#' @param unique A logical scalar. Avoid redundancy by removing all self sample
#'   comparisons and keep only one of two pairwise comparisons? Default is TRUE.
#' @param method The dissimilarity index. See \code{\link[vegan]{vegdist}} for
#'   all options. Default is "bray".
#' @param binary A logical scalar. Perform presence/absence standardization
#'   before analysis. See \code{\link[vegan]{vegdist}}. Default is FALSE.
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
#'                      taxa_are_columns = FALSE
#'                      )
#' # Report numbers
#' numbers <- data %>%
#'  get_betas()
betas <- function(ta, unique = T, method = "bray", binary = F) {

  # make "dist" object with beta values
  rel_abundance_matrix <- get_rel_abundance_matrix(ta)
  betas_dist <- vegdist(rel_abundance_matrix, method = method, binary = binary)

  # save number of betas in betas_dist in shortcut variable
  n <- attr(betas_dist, "Size")

  # make tibble with beta values if we want only unique sample pairs
  if (unique) {

    betas <- expand.grid(i = 1:n, j = 1:n) %>%
      filter(i < j) %>%
      mutate(sample_id_1 = labels(betas_dist)[i]) %>%
      mutate(sample_id_2 = labels(betas_dist)[j]) %>%
      mutate(beta = betas_dist[n * (i - 1) - i * (i - 1) / 2 + j - i]) %>%
      select(- i, - j)

    # make tibble with beta values if we want all sample pairs (redundant!)
  } else {

    betas <- as.matrix(betas_dist) %>%
      as_tibble() %>%
      mutate(sample_id_1 = attr(betas_dist, "Labels")) %>%
      gather(key = "sample_id_2", value = "beta", - sample_id_1)

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

# synonym
#' Get betas
#'
#' \code{get_betas} returns a tidy tibble with the beta diversity for each
#' combination of samples.
#'
#' This function calculates the beta diversity using the
#' \code{\link[vegan]{vegdist}} function of Vegan. It will report one diversity
#' estimate for each combination of samples.
#'
#'
#' @param ta Tidyamplicons object.
#' @param unique A logical scalar. Avoid redundancy by removing all self sample
#'   comparisons and keep only one of two pairwise comparisons? Default is TRUE.
#' @param method The dissimilarity index. See \code{\link[vegan]{vegdist}} for
#'   all options. Default is "bray".
#' @param binary A logical scalar. Perform presence/absence standardization
#'   before analysis. See \code{\link[vegan]{vegdist}}. Default is FALSE.
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
#'                      taxa_are_columns = FALSE
#'                      )
#' # Report numbers
#' numbers <- data %>%
#'  get_betas()
get_betas <- betas

# DEPRICATED: use occurrences()
# Returns a tidy table of taxon presence and absence counts in sample conditions.
# Condition is a variable that should be present in the samples table.
taxon_counts_in_conditions <- function(ta, condition) {

  condition <- enquo(condition)
  condition_name <- quo_name(condition)

  ta$abundances %>%
    filter(abundance > 0) %>%
    left_join(ta$samples) %>%
    select(taxon_id, sample_id, condition = !! condition) %>%
    mutate(presence = "present") %>%
    complete(nesting(condition, sample_id), taxon_id, fill = list(presence = "absent")) %>%
    count(taxon_id, condition, presence) %>%
    complete(taxon_id, condition, presence, fill = list(n = 0)) %>%
    rename(!! condition_name := condition)

}

# Returns a tidy table of occurrences: taxon presence counts in samples, overall
# or per condition.
# Condition should be a categorical variable present in the samples table.
# Supply condition as a string.
occurrences <- function(ta, condition = NULL, pres_abs = F) {

  abundances_extended <-
    ta$abundances %>%
    filter(abundance > 0) %>%
    left_join(ta$samples)

  if (is.null(condition)) {

    abundances_extended %>%
      count(taxon_id) %>%
      rename(occurrence = n)

  } else if (pres_abs) {

    condition <- sym(condition)

    abundances_extended %>%
      select(taxon_id, sample_id, !! condition) %>%
      mutate(presence = "present") %>%
      complete(nesting(!! condition, sample_id), taxon_id, fill = list(presence = "absent")) %>%
      count(taxon_id, !! condition, presence) %>%
      complete(taxon_id, !! condition, presence, fill = list(n = 0))

  } else {

    condition <- sym(condition)

    abundances_extended %>%
      count(taxon_id, !! condition) %>%
      rename(occurrence = n) %>%
      complete(taxon_id, !! condition, fill = list(occurrence = 0))

  }

}

# Returns tidy table with average relatively abundances of taxa, overall or per
# condition.
# Condition should be a categorical variable present in the samples table.
# Supply condition as a string.
mean_rel_abundances <- function(ta, condition = NULL) {

  # if rel_abundance not present: add and remove on exit
  if (! "rel_abundance" %in% names(ta$abundances)) {
    ta <- add_rel_abundance(ta)
  }

  if (is.null(condition)) {

    ta$abundances %>%
      select(sample_id, taxon_id, rel_abundance) %>%
      complete(sample_id, taxon_id, fill = list(rel_abundance = 0)) %>%
      group_by(taxon_id) %>%
      summarize(mean_rel_abundance = mean(rel_abundance)) %>%
      ungroup()

  } else {

    condition <- sym(condition)

    ta$abundances %>%
      left_join(ta$samples) %>%
      select(!! condition, sample_id, taxon_id, rel_abundance) %>%
      complete(nesting(!! condition, sample_id), taxon_id, fill = list(rel_abundance = 0)) %>%
      group_by(taxon_id, !! condition) %>%
      summarize(mean_rel_abundance = mean(rel_abundance)) %>%
      ungroup()

  }

}

# Returns joined tibble with all abundances, samples and taxa variables
everything <- function(ta) {

  # make and return large table
  ta$abundances %>%
    left_join(ta$samples) %>%
    left_join(ta$taxa)

}

# Returns samples tibble
samples <- function(ta) ta$samples

# Returns taxa tibble
taxa <- function(ta) ta$taxa

# Returns abundances
abundances <- function(ta) ta$abundances

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

  metadata <- tibble(sample_id = rownames(abundances_matrix)) %>%
    left_join(ta$samples)

  adonis(
    as.formula(paste("abundances_matrix", formula_RHS, sep = " ~ ")),
    metadata,
    permutations = permutations
  )

}
