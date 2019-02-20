#' Add relative abundance to abundance table
#'
#' \code{add_rel_abundance} adds relative abundance to the abundance table of a
#' tidyamplicons object.
#'
#' This function adds the relative abundance per sample to the abundance table
#' of a tidyamplicons object under the variable name "rel_abundance".
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
#'                      taxa_are_columns = FALSE)
#'
#' # Add relative abundance
#' data <- data %>%
#'  add_rel_abundance()
#' 
add_rel_abundance <- function(ta) {

  # add relative abundance to abundance table
  ta$abundances <- ta$abundances %>%
    group_by(sample_id) %>%
    mutate(rel_abundance = abundance / sum(abundance)) %>%
    ungroup()

  # return ta object
  ta

}
