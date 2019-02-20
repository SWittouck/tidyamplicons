
tidyamplicons <- function(abundance_matrix, taxa_are_columns = TRUE, taxon_names_are_sequences = TRUE) {

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

add_sample_tibble <- function(ta, sample_tibble) {

  modify_at(ta, "samples", left_join, sample_tibble)

}

add_taxon_tibble <- function(ta, taxon_tibble) {

  modify_at(ta, "taxa", left_join, taxon_tibble)

}

change_sample_ids <- function(ta, sample_new) {

  if (any(duplicated(ta$samples[[sample_new]]))) {
    stop("the new sample ids are not unique")
  }

  ta$samples <-
    ta$samples %>%
    rename(sample_new = !! sample_new) %>%
    mutate_at("sample_new", as.character)

  ta$abundances <-
    ta$abundances %>%
    left_join(ta$samples %>% select(sample, sample_new)) %>%
    select(- sample) %>%
    rename(sample = sample_new)

  ta$samples <-
    ta$samples %>%
    select(- sample) %>%
    rename(sample = sample_new)

  ta

}

change_taxon_ids <- function(ta, taxon_new) {

  if (any(duplicated(ta$taxa[[taxon_new]]))) {
    stop("the new taxon ids are not unique")
  }

  ta$taxa <-
    ta$taxa %>%
    rename(taxon_new = !! taxon_new) %>%
    mutate_at("taxon_new", as.character)

  ta$abundances <-
    ta$abundances %>%
    left_join(ta$taxa %>% select(taxon, taxon_new)) %>%
    select(- taxon) %>%
    rename(taxon = taxon_new)

  ta$taxa <-
    ta$taxa %>%
    select(- taxon) %>%
    rename(taxon = taxon_new)

  ta

}

as_phyloseq <- function(ta, sample_id = "sample_name", taxon_id = "sequence") {

  if ("phyloseq" %in% class(ta)) return(ta)

  if (! is.null(sample_id)) {
    ta <- change_sample_ids(ta, sample_new = sample_id)
  }

  if (! is.null(taxon_id)) {
    ta <- change_taxon_ids(ta, taxon_new = taxon_id)
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

