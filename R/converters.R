#' Initiate tidyamplicons object
#'
#' \code{tidyamplicons} returns a tidyamplicons object given a numeric matrix.
#'
#' This function initiates a tidyamplicons object based on a numeric matrix. It
#' will automatically create a dummy taxa and sample table which will need to be
#' updated using the functions \code{\link{add_taxon_tibble}} and
#' \code{\link{add_sample_tibble}}.
#'
#' @param abundance_matrix Numerical matrix containing the abundance data.
#' @param taxa_are_columns A logical scalar. Are the taxa defined in columns?
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
#'
#' \dontrun{
#' tidyamplicons("a")
#' }
#'
#' @export
create_tidyamplicons <- function(abundance_matrix, taxa_are_columns = TRUE) {

  if (
    ! is.matrix(abundance_matrix) |
    ! is.numeric(abundance_matrix)
  ) stop("first argument should be a numeric matrix")

  if (! taxa_are_columns) abundance_matrix = t(abundance_matrix)

  abundance_matrix <-
    abundance_matrix[, colSums(abundance_matrix) != 0]

  ta <- list()
  class(ta) <- "tidyamplicons"

  n_samples <- nrow(abundance_matrix)
  n_taxa <- ncol(abundance_matrix)

  sample_ids <- str_c("s", 1:n_samples)
  taxon_ids <- str_c("t", 1:n_taxa)

  ta$abundances <-
    abundance_matrix %>%
    as.vector() %>%
    tibble(abundance = .) %>%
    mutate(sample_id = rep(!! sample_ids, times = !! n_taxa)) %>%
    mutate(taxon_id = rep(!! taxon_ids, each = !! n_samples)) %>%
    filter(abundance > 0)

  ta$samples <-
    abundance_matrix %>%
    rownames() %>%
    tibble(sample = .) %>%
    mutate(sample_id = !! sample_ids)

  ta$taxa <-
    abundance_matrix %>%
    colnames() %>%
    tibble(taxon = .) %>%
    mutate(taxon_id = !! taxon_ids)

  ta

}

#' Write community data in tidyamplicons format
#' @importFrom readr write_csv
#' @param ta a tidyamplicons object
#' @param dout the directory to store the three tidyamplicons tables in
#' @export
write_tidyamplicons <- function(ta, dout) {
  dir.create(dout)
  write_csv(ta$samples, paste0(dout, "/samples.csv"))
  write_csv(ta$taxa, paste0(dout, "/taxa.csv"))
  write_csv(ta$abundances, paste0(dout, "/abundances.csv"))
}

#' Read community data written by tidyamplicons
#' @importFrom readr read_csv
#' @param din directory containing the a sample, taxa and abundances table in csv format
#' @param samples the name of the samples table, defaults to samples.csv
#' @param taxa the name of the taxa table, defaults to taxa.csv
#' @param abundances the name of the abundances table, defaults to abundances.csv
#' @export
read_tidyamplicons <- function(din, samples = "samples.csv", taxa = "taxa.csv",
                               abundances = "abundances.csv") {
  samples <- read_csv(paste0(din, "/", samples), col_types = readr::cols())
  taxa <- read_csv(paste0(din, "/", taxa), col_types = readr::cols())
  abundances <- read_csv(paste0(din, "/", abundances), col_types = readr::cols())
  make_tidyamplicons(
    samples, taxa, abundances, sample_name = sample_id, taxon_name = taxon_id
  )
}

#' Update old tidyamplicons object to new one.
#'
#' \code{update_tidyamplicons} updates an old tidyamplicons object to a new one.
#'
#' This function will update a tidyamplicons object created prior to version
#' 0.1.0 to a tidyamplicons object compatible with version 0.1.0.
#'
#' @param ta Old tidyamplicons object.
#'
#' @export
update_tidyamplicons <- function(ta) {

  ta %>%
    update_id_names() %>%
    mutate_taxa(sequence = taxon_id) %>%
    reset_ids()

}

#' Reset the taxon and sample IDs
#' @param ta a tidyamplicons object
#' @export
reset_ids <- function(ta, keep_prev = F) {

  if (keep_prev) {

    ta <-
      ta %>%
      mutate_samples(sample_id_prev = sample_id) %>%
      mutate_taxa(taxon_id_prev = taxon_id)

  }

  ta %>%
    change_id_samples(sample_id_new = str_c("s", seq_len(n()))) %>%
    change_id_taxa(taxon_id_new = str_c("t", seq_len(n())))

}

#' Rename the "sample" and "taxon" columns to "sample_id" and "taxon_id"
#' @param ta a tidyamplicons object
#' @export
update_id_names <- function(ta) {

  ta %>%
    purrr::modify_at("samples", rename, sample_id = sample) %>%
    purrr::modify_at("taxa", rename, taxon_id = taxon) %>%
    purrr::modify_at("abundances", rename, sample_id = sample, taxon_id = taxon)

}

#' Convert tidyamplicons object to phyloseq object
#'
#' \code{as_phyloseq} returns a phyloseq object given a tidyamplicons object.
#'
#' This function will convert a tidyamplicons object into a phyloseq object for
#' alternative processing using the phyloseq package. To convert from a phyloseq
#' object to a tidyamplicons object use \code{\link{as_tidyamplicons}}.
#'
#' @param ta Tidyamplicons object.
#' @param sample  The sample names required for a phyloseq object. Default is
#'   "sample" column in sample tibble of the tidyamplicons object.
#' @param taxon The taxon names required for a phyloseq object. Default is
#'   "taxon" column in taxon tibble of the tidyamplicons object.
#'
#' @export
as_phyloseq <- function(ta, sample = sample, taxon = taxon) {

  if ("phyloseq" %in% class(ta)) return(ta)

  sample <- rlang::enexpr(sample)
  taxon <- rlang::enexpr(taxon)

  ta <- change_id_samples(ta, sample_id_new = !! sample)
  ta <- change_id_taxa(ta, taxon_id_new = !! taxon)

  otu_table <-
    ta$abundances %>%
    spread(key = taxon_id, value = abundance, fill = 0) %>%
    `attr<-`("class", "data.frame") %>%
    `rownames<-`(.$sample_id) %>%
    select(- sample_id) %>%
    as.matrix() %>%
    phyloseq::otu_table(taxa_are_rows = F)

  if (ncol(ta$samples) == 1) {
    ta <- mutate_samples(ta, dummy = as.character(1:nrow(ta$samples)))
  }

  if (ncol(ta$taxa) == 1) {
    ta <- mutate_taxa(ta, dummy = as.character(1:nrow(ta$taxa)))
  }

  sample_data <-
    ta$samples %>%
    `attr<-`("class", "data.frame") %>%
    `rownames<-`(.$sample_id) %>%
    select(- sample_id) %>%
    phyloseq::sample_data()

  tax_table <-
    ta$taxa %>%
    `attr<-`("class", "data.frame") %>%
    `rownames<-`(.$taxon_id) %>%
    select(- taxon_id) %>%
    as.matrix() %>%
    phyloseq::tax_table()

  phyloseq::phyloseq(otu_table, sample_data, tax_table)

}

#' Convert phyloseq object to tidyamplicons object
#'
#' \code{as_tidyamplicons} returns a tidyamplicons object given a phyloseq
#' object.
#'
#' This function will convert a phyloseq object into a tidyamplicons object. To
#' convert from a tidyamplicons object to a phyloseq object use
#' \code{\link{as_phyloseq}}.
#'
#' @param ps Phyloseq object.
#'
#' @export
as_tidyamplicons <- function(ps) {

  if ("tidyamplicons" %in% class(ps)) return(ps)

  # convert sample data to tibble
  samples <-
    phyloseq::sample_data(ps)@.Data %>%
    `names<-`(phyloseq::sample_data(ps)@names) %>%
    do.call(what = tibble) %>%
    mutate(sample = phyloseq::sample_data(ps)@row.names)

  # convert taxon table to tibble
  taxa <-
    phyloseq::tax_table(ps)@.Data %>%
    as_tibble() %>%
    mutate(taxon = phyloseq::tax_table(ps) %>% row.names()) %>%
    `names<-`(names(.) %>% str_to_lower())

  # make sure that taxa are columns in abundances table
  if (phyloseq::taxa_are_rows(ps)) {
    phyloseq::otu_table(ps) <- phyloseq::t(phyloseq::otu_table(ps))
  }

  phyloseq::otu_table(ps)@.Data %>%
    create_tidyamplicons() %>%
    add_sample_tibble(samples) %>%
    add_taxon_tibble(taxa)

}

#' Convert matrix with abundances to tidy data frame
#'
#' \code{as_abundances} returns a tidy data frame given a numerical abundance
#' matrix.
#'
#' This function will convert a numerical abundance matrix into a tidy data
#' frame. To convert a tidy data frame into a numerical abundance matrix
#' use \code{\link{as_abundances_matrix}}.
#'
#' @param abundances_matrix The ambundance matrix that will be converted.
#' @param taxa_are_columns A logical scalar. Are the taxa defined in columns?
#'   Default is TRUE.
#' @param value Name of resulting colum containing the abundance data. Default
#'   is "abundance".
#'
#' @export
as_abundances <- function(abundances_matrix, taxa_are_columns = TRUE,
                          value = "abundance") {

  if (
    ! is.matrix(abundances_matrix) |
    ! is.numeric(abundances_matrix)
  ) stop("first argument should be an abundances matrix")

  if (! taxa_are_columns) abundances_matrix = t(abundances_matrix)

  abundances_matrix %>%
    as_tibble() %>%
    mutate(sample_id = row.names(abundances_matrix)) %>%
    gather(key = "taxon_id", value = !! value, - sample_id) %>%
    filter(!! value > 0)

}

#' Convert abundances tidy data frame to matrix.
#'
#' \code{as_abundances_matrix} returns a numerical matrix given a tidy
#' abundances data frame.
#'
#' This function will convert a abundances tidy data frame into a numerlical
#' abundance matrix. To convert a numerical abundance matrix into a abundances
#' tidy data frame use \code{\link{as_abundances_matrix}}.
#'
#' @param abundances The abundance tidy data frame that will be converted.
#' @param value Name of colum containing the abundance data. Default is
#'   "abundance".
#'
#' @export
as_abundances_matrix <- function(abundances, value = abundance) {

  if (
    ! is.data.frame(abundances) |
    is.null(abundances$taxon_id) |
    is.null(abundances$sample_id)
  ) stop("first argument should be an abundances table (data frame)")

  value <- enquo(value)

  abundances_wide <- abundances %>%
    select(sample_id, taxon_id, !! value) %>%
    spread(key = taxon_id, value = !! value, fill = 0)

  abundances_wide %>%
    select(- sample_id) %>%
    as.matrix() %>%
    `row.names<-`(abundances_wide$sample_id)

}

#' Merge two tidyamplicons objects.
#'
#' \code{merge_tidyamplicons} merges two tidyamplicons objects and returns one
#' single tidyamplicons object.
#'
#' This function will merge two tidyamplicons objects into one. It is useful if
#' one wants to merge data obtained from different sequencing runs. Therefore,
#' this function requirers that both tidyamplicons objects contain a "run"
#' variable in their samples table, indicating their origin.
#'
#' @param ta1 The first tidyamplicons object.
#' @param ta2 The second tidyamplicons object.
#'
#' @export
merge_tidyamplicons <- function(ta1, ta2, taxon_identifier = sequence) {

  taxon_identifier <- rlang::ensym(taxon_identifier)

  ti <- rlang::as_string(taxon_identifier)
  if (! (ti %in% names(ta1$taxa) & ti %in% names(ta2$taxa))) {
    stop("the taxon identifier was not found in one or both of the ta objects")
  }

  # make sure that sample names are unique
  ta1 <- change_id_samples(ta1, paste("ta1", sample_id, sep = "_"))
  ta2 <- change_id_samples(ta2, paste("ta2", sample_id, sep = "_"))

  # merge sample tables
  samples <- bind_rows(ta1$samples, ta2$samples)

  # change taxon ids to something meaningful across ta objects
  ta1 <- change_id_taxa(ta1, taxon_id_new = !! taxon_identifier)
  ta2 <- change_id_taxa(ta2, taxon_id_new = !! taxon_identifier)

  # merge taxa tables
  taxa <-
    bind_rows(ta1$taxa, ta2$taxa) %>%
    group_by(taxon_id) %>%
    summarize_all(function(x) {
      x <- unique(x)
      x <- x[! is.na(x)]
      if (length(x) == 1) return(x)
      as.character(NA)
    })

  # merge abundances tables
  abundances <- bind_rows(ta1$abundances, ta2$abundances)

  # make new ta object
  ta <- list(samples = samples, taxa = taxa, abundances = abundances)
  class(ta) <- "tidyamplicons"

  # give new sample names in new ta object
  ta <- reset_ids(ta)

  # return ta object
  ta

}

#' Create a tidyamplicons object from three tidy tables
#'
#' @export
make_tidyamplicons <- function(samples, taxa, abundances,
                               sample_name = sample, taxon_name = taxon) {

  sample_name <- rlang::enexpr(sample_name)
  taxon_name <- rlang::enexpr(taxon_name)

  list(samples = samples, taxa = taxa, abundances = abundances) %>%
    purrr::modify_at("samples", mutate, sample_id = !! sample_name) %>%
    purrr::modify_at("taxa", mutate, taxon_id = !! taxon_name) %>%
    purrr::modify_at("abundances", rename, sample_id = !! sample_name) %>%
    purrr::modify_at("abundances", rename, taxon_id = !! taxon_name) %>%
    purrr::modify_at("abundances", filter, abundance > 0) %>%
    purrr::modify_at("abundances", filter, sample_id %in% .$samples$sample_id) %>%
    purrr::modify_at("abundances", filter, taxon_id %in% .$taxa$taxon_id) %>%
    `class<-`("tidyamplicons") %>%
    reset_ids()

}

#' Convert phyloseq object to tidyamplicons object
#'
#' DEPRECATED, use \code{\link{as_tidyamplicons}}
#'
#' @export
tidy_phyloseq <- function(ps) {

  # convert sample data
  samples <-
    phyloseq::sample_data(ps)@.Data %>%
    `names<-`(phyloseq::sample_data(ps)@names) %>%
    do.call(what = tibble) %>%
    mutate(sample_id = phyloseq::sample_data(ps)@row.names)

  # convert taxon table
  taxa <- phyloseq::tax_table(ps)@.Data %>%
    as_tibble() %>%
    mutate(taxon_id = phyloseq::tax_table(ps) %>% row.names()) %>%
    set_names(names(.) %>% str_to_lower())

  # make sure that taxa are columns in taxon table
  if (phyloseq::taxa_are_rows(ps)) phyloseq::otu_table(ps) <- phyloseq::t(phyloseq::otu_table(ps))

  # convert taxon table
  abundances <- phyloseq::otu_table(ps)@.Data %>%
    as_abundances(taxa_are_columns = T)

  # make and return tidyamplicons object
  make_tidyamplicons(
    samples = samples,
    taxa = taxa,
    abundances = abundances
  )

}
