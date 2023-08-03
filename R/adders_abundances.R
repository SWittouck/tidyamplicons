#' Add relative abundance to abundance table
#'
#' \code{add_rel_abundance} adds relative abundance to the abundance table of a
#' tidytacos object.
#'
#' This function adds the relative abundance per sample to the abundance table
#' of a tidytacos object under the variable name "rel_abundance".
#'
#' @param ta tidytacos object.
#'
#' @examples
#' # Initiate count matrix
#' x <- matrix(
#'  c(1500, 1300, 280, 356),
#'  ncol = 2
#' )
#' rownames(x) <- c("taxon1", "taxon2")
#' colnames(x) <- c("sample1", "sample2")
#'
#' # Convert to tidytacos object
#' data <- create_tidytacos(x,
#'                      taxa_are_columns = FALSE)
#'
#' # Add relative abundance
#' data <- data %>%
#'  add_rel_abundance()
#'
#' @export
add_rel_abundance <- function(ta) {

  # add relative abundance to abundance table
  ta$counts <- ta$counts %>%
    group_by(sample_id) %>%
    mutate(rel_abundance = readcount / sum(readcount)) %>%
    ungroup()

  # return ta object
  ta

}
