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
#' @param taxon_names_are_sequences A logical scalar. Are the taxon names the
#'   full length sequences they present?
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
#'
#' \dontrun{
#' tidyamplicons("a")
#' }

create_tidyamplicons <- function(abundance_matrix, taxa_are_columns = TRUE, taxon_names_are_sequences = TRUE) {

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

  sample_names <- str_c("s", 1:n_samples)
  taxon_names <- str_c("t", 1:n_taxa)

  ta$abundances <-
    abundance_matrix %>%
    as.vector() %>%
    tibble(abundance = .) %>%
    mutate(sample = rep(!! sample_names, times = !! n_taxa)) %>%
    mutate(taxon = rep(!! taxon_names, each = !! n_samples)) %>%
    filter(abundance > 0)

  ta$samples <-
    abundance_matrix %>%
    rownames() %>%
    tibble(sample_name = .) %>%
    mutate(sample = !! sample_names)

  ta$taxa <-
    abundance_matrix %>%
    colnames() %>%
    tibble(sequence = .) %>%
    mutate(taxon = !! taxon_names)

  if (! taxon_names_are_sequences) {
    ta$taxa <- rename(ta$taxa, taxon_name = sequence)
  }

  ta

}

reset_ids <- function(ta) {

  ta %>%
    mutate_samples(sample_prev = sample) %>%
    mutate_taxa(taxon_prev = taxon) %>%
    mutate_samples(sample_new = str_c("s", 1:n())) %>%
    mutate_taxa(taxon_new = str_c("t", 1:n())) %>%
    change_ids_samples(sample_new = "sample_new") %>%
    change_ids_taxa(taxon_new = "taxon_new")

}

as_phyloseq <- function(ta, sample_id = "sample_name", taxon_id = "sequence") {

  if ("phyloseq" %in% class(ta)) return(ta)

  if (! is.null(sample_id)) {
    ta <- change_ids_samples(ta, sample_new = sample_id)
  }

  if (! is.null(taxon_id)) {
    ta <- change_ids_taxa(ta, taxon_new = taxon_id)
  }

  otu_table <-
    ta$abundances %>%
    spread(key = taxon, value = abundance, fill = 0) %>%
    `attr<-`("class", "data.frame") %>%
    `rownames<-`(.$sample) %>%
    select(- sample) %>%
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
    `rownames<-`(.$sample) %>%
    select(- sample) %>%
    phyloseq::sample_data()

  tax_table <-
    ta$taxa %>%
    `attr<-`("class", "data.frame") %>%
    `rownames<-`(.$taxon) %>%
    select(- taxon) %>%
    as.matrix() %>%
    phyloseq::tax_table()

  phyloseq::phyloseq(otu_table, sample_data, tax_table)

}

# converts a phyloseq object to a tidyamplicons object
# the plyloseq object should contain absolute abundances
as_tidyamplicons <- function(ps) {

  if ("tidyamplicons" %in% class(ps)) return(ps)

  # convert sample data to tibble
  samples <-
    phyloseq::sample_data(ps)@.Data %>%
    `names<-`(phyloseq::sample_data(ps)@names) %>%
    do.call(what = tibble) %>%
    mutate(sample_name = phyloseq::sample_data(ps)@row.names)

  # convert taxon table to tibble
  taxa <-
    phyloseq::tax_table(ps)@.Data %>%
    as_tibble() %>%
    mutate(sequence = phyloseq::tax_table(ps) %>% row.names()) %>%
    set_names(names(.) %>% str_to_lower())

  # make sure that taxa are columns in abundances table
  if (phyloseq::taxa_are_rows(ps)) {
    phyloseq::otu_table(ps) <- phyloseq::t(phyloseq::otu_table(ps))
  }

  phyloseq::otu_table(ps)@.Data %>%
    tidyamplicons() %>%
    add_sample_tibble(samples) %>%
    add_taxon_tibble(taxa)

}

# convert matrix with abundances to tidy data frame
as_abundances <- function(abundances_matrix, taxa_are_columns = TRUE, value = "abundance") {

  if (
    ! is.matrix(abundances_matrix) |
    ! is.numeric(abundances_matrix)
  ) stop("first argument should be an abundances matrix")

  if (! taxa_are_columns) abundances_matrix = t(abundances_matrix)

  abundances_matrix %>%
    as_tibble() %>%
    mutate(sample = row.names(abundances_matrix)) %>%
    gather(key = "taxon", value = !! value, - sample) %>%
    filter(!! value > 0)

}

# convert abundances tidy data frame to matrix
as_abundances_matrix <- function(abundances, value = abundance) {

  if (
    ! is.data.frame(abundances) |
    is.null(abundances$taxon) |
    is.null(abundances$sample)
  ) stop("first argument should be an abundances table (data frame)")

  value <- enquo(value)

  abundances_wide <- abundances %>%
    select(sample, taxon, !! value) %>%
    spread(key = taxon, value = !! value, fill = 0)

  abundances_wide %>%
    select(- sample) %>%
    as.matrix() %>%
    `row.names<-`(abundances_wide$sample)

}

# requires that both as objects have a "run" variable in their samples table
merge_tidyamplicons <- function(ta1, ta2) {

  # make sure that sample names are unique
  ta1$samples <- ta1$samples %>%
    mutate(sample_new = paste(run, sample, sep = "_"))
  ta2$samples <- ta2$samples %>%
    mutate(sample_new = paste(run, sample, sep = "_"))
  ta1 <- process_new_sample_name(ta1)
  ta2 <- process_new_sample_name(ta2)

  # merge sample tables
  samples <- bind_rows(ta1$samples, ta2$samples)

  # merge taxa tables
  taxa <- bind_rows(ta1$taxa, ta2$taxa) %>%
    select(taxon, kingdom, phylum, class, order, family, genus, species) %>%
    group_by(taxon) %>%
    summarize_all(function(x) {
      x <- unique(x)
      x <- x[! is.na(x)]
      if (length(x) == 1) return(x)
      as.character(NA)
    })

  # merge abundances tables
  abundances <- bind_rows(ta1$abundances, ta2$abundances)

  # make new ta object
  ta <- make_tidyamplicons(samples, taxa, abundances)

  # give new sample names in new ta object
  ta$samples <- ta$samples %>%
    mutate(sample_new = paste("m", 1:n(), sep = ""))
  ta <- process_new_sample_name(ta)

  # return ta object
  ta

}

# DEPRECATED
make_tidyamplicons <- function(samples, taxa, abundances) {

  # remove abundances of zero
  abundances <- abundances %>%
    filter(abundance > 0)

  # make tidyamplicons (ta) object
  ta <- list(
    samples = samples,
    taxa = taxa,
    abundances = abundances
  )
  class(ta) <- "tidyamplicons"

  # make sure that all tables contain the same unique
  # taxa and samples, and return ta object
  ta %>%
    process_abundance_selection() %>%
    process_sample_selection() %>%
    process_taxon_selection()

}

# DEPRECATED
# converts a phyloseq object to a tidyamplicons object
# the plyloseq object should contain absolute abundances
tidy_phyloseq <- function(ps) {

  # convert sample data
  samples <-
    phyloseq::sample_data(ps)@.Data %>%
    `names<-`(phyloseq::sample_data(ps)@names) %>%
    do.call(what = tibble) %>%
    mutate(sample = phyloseq::sample_data(ps)@row.names)

  # convert taxon table
  taxa <- phyloseq::tax_table(ps)@.Data %>%
    as_tibble() %>%
    mutate(taxon = phyloseq::tax_table(ps) %>% row.names()) %>%
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
