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
#'
#' @export
report_numbers <- function(ta) {

  sprintf("samples: %i", nrow(ta$samples)) %>% message()
  sprintf("taxa: %i", nrow(ta$taxa)) %>% message()
  sprintf("reads: %i", sum(ta$abundances$abundance)) %>% message()

}

#' Return some descriptive numbers
#'
#' \code{get_numbers} returns the number of samples, taxa and reads in a
#' tidyamplicons object.
#'
#' This function returns the number of samples, taxa and reads in a
#' tidyamplicons object, stored in a named numeric vector. To print the numbers,
#' use \code{\link{report_numbers}} instead.
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
#' data <- create_tidyamplicons(x, taxa_are_columns = FALSE)
#'
#' # Report numbers
#' numbers <- data %>% get_numbers()
#'
#' @export
numbers <- function(ta) {

  c(
    n_samples = nrow(ta$samples),
    n_taxa = nrow(ta$taxa),
    n_reads = sum(ta$abundances$abundance)
  )

}

#' Return some descriptive numbers
#'
#' DEPRECATED, use \code{\link{numbers}}
#'
#' @export
get_numbers <- numbers

#' Return a relative abundance matrix
#'
#' DEPRECATED, use \code{\link{abundances_matrix}}
#'
#' @export
get_rel_abundance_matrix <- function(ta) {

  # add relative abundances if not present
  if (! "rel_abundance" %in% names(ta$abundances)) {
    ta <- add_rel_abundance(ta)
  }

  as_abundances_matrix(ta$abundances, value = rel_abundance)

}

#' Get beta diversity table
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
#' data <-
#'   create_tidyamplicons(x, taxa_are_columns = FALSE)
#'
#' # Report numbers
#' numbers <- data %>% betas()
#'
#' @export
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
    right_join(betas, by = "sample_id_2")
  betas <- ta$samples %>%
    `names<-`(names(.) %>% str_c("_1")) %>%
    right_join(betas, by = "sample_id_1")

  # return betas table
  betas

}

#' Get beta diversity table
#'
#' DEPRECATED, use \code{\link{betas}}
#'
#' @export
get_betas <- betas

#' Get occurrences of taxa in general or in conditions
#'
#' DEPRECATED, use \code{\link{occurrences}}
#'
#' @export
taxon_counts_in_conditions <- function(ta, condition) {

  condition <- enquo(condition)
  condition_name <- quo_name(condition)

  ta$abundances %>%
    filter(abundance > 0) %>%
    left_join(ta$samples, by = "sample_id") %>%
    select(taxon_id, sample_id, condition = !! condition) %>%
    mutate(presence = "present") %>%
    complete(nesting(condition, sample_id), taxon_id, fill = list(presence = "absent")) %>%
    count(taxon_id, condition, presence) %>%
    complete(taxon_id, condition, presence, fill = list(n = 0)) %>%
    rename(!! condition_name := condition)

}

#' Get occurrences of taxa in general or per condition
#'
#' Returns a tidy table of occurrences: taxon presence counts in samples,
#' overall or per condition.
#'
#' Condition should be a categorical variable present in the samples table.
#' Supply condition as a string.
#'
#' @export
occurrences <- function(ta, condition = NULL, pres_abs = F) {

  abundances_extended <-
    ta$abundances %>%
    filter(abundance > 0) %>%
    left_join(ta$samples, by = "sample_id")

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

#' Get mean relative abundances of taxa in general or per condition
#'
#' Returns tidy table with average relatively abundances of taxa, overall or per
#' condition.
#'
#' Condition should be a categorical variable present in the samples table.
#' Supply condition as a string.
#'
#' @export
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
      left_join(ta$samples, by = "sample_id") %>%
      select(!! condition, sample_id, taxon_id, rel_abundance) %>%
      complete(nesting(!! condition, sample_id), taxon_id, fill = list(rel_abundance = 0)) %>%
      group_by(taxon_id, !! condition) %>%
      summarize(mean_rel_abundance = mean(rel_abundance)) %>%
      ungroup()

  }

}

#' Get all data in one single table
#'
#' @export
everything <- function(ta) {

  # make and return large table
  ta$abundances %>%
    left_join(ta$samples, by = "sample_id") %>%
    left_join(ta$taxa, by = "taxon_id")

}

#' Extract the sample table
#'
#' @export
samples <- function(ta) ta$samples

#' Extract the taxon table
#'
#' @export
taxa <- function(ta) ta$taxa

#' Extract the abundance table
#'
#' @export
abundances <- function(ta) ta$abundances

#' Perform an adonis test
#'
#' This function executes the \link[vegan]{adonis} function of the vegan package
#' and returns the result.
#'
#' Samples where one or more predictors are NA are removed.
#'
#' @param ta A tidyamplicons object.
#' @param predictors A character vector with predictors to include in the model.
#' @param permutations The number of permutations (more permutations takes
#'   longer but gives a more accurate p-value).
#'
#' @return An object of class "adonis" (see \link[vegan]{adonis}).
#'
#' @export
perform_adonis <- function(ta, predictors, permutations = 999) {

  abundances_matrix <- ta %>%
    purrr::modify_at("samples", drop_na, one_of(predictors)) %>%
    process_sample_selection() %>%
    add_rel_abundance() %>%
    abundances() %>%
    as_abundances_matrix(value = "rel_abundance")

  formula_RHS <- paste0(predictors, collapse = " + ")

  metadata <- tibble(sample_id = rownames(abundances_matrix)) %>%
    left_join(ta$samples, by = "sample_id")

  adonis(
    as.formula(paste("abundances_matrix", formula_RHS, sep = " ~ ")),
    metadata,
    permutations = permutations
  )

}

#' Return an abundances matrix
#'
#' This function returns the abundances table in the form of a matrix where the
#' rows are samples and the column are taxa.
#'
#' @param ta A tidyamplicons object.
#' @param value The name of a variable in the abundances table that contains the
#'   abundances (unquoted). Could be relative abundances, if present.
#' @param sample_name The name of the variable in the sample table to use as row
#'   names (unquoted).
#' @param taxon_name The name of the variable in the taxon table to use as
#'   column names (unquoted).
#' @return A matrix with abundance values.
#'
#' @export
abundances_matrix <-
  function(ta, value = abundance, sample_name = sample, taxon_name = taxon) {

  if (
    ! "tidyamplicons" %in% class(ta)
  ) stop("first argument should be a tidyamplicons object")

  value <- rlang::enquo(value)
  sample_name <- rlang::enquo(sample_name)
  taxon_name <- rlang::enquo(taxon_name)

  ta %>%
    change_id_samples(sample_id_new = {{sample_name}}) %>%
    change_id_taxa(taxon_id_new = {{taxon_name}}) %>%
    abundances() %>%
    as_abundances_matrix(value = {{value}})

}
