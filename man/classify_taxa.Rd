% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adders_taxa.R
\name{classify_taxa}
\alias{classify_taxa}
\title{(Re)classify amplicon sequences}
\usage{
classify_taxa(
  ta,
  refdb,
  taxa = rep(T, times = length(taxon_id)),
  ranks = "default",
  sequence_var = "sequence",
  multithread = T,
  min_boot = 50,
  n_ranks = 7
)
}
\arguments{
\item{ta}{A tidyamplicons object.}

\item{refdb}{The path to a DADA2-compatible reference database.}

\item{taxa}{An expression that specifies which taxa to (re)classify.}

\item{ranks}{A vector that specifies which ranks to (re)classify.}

\item{sequence_var}{The (quoted) name of a variable within the taxa table
that contains (representative) sequences of the taxa.}

\item{multithread}{A boolean indicating whether to use multiple threads.}

\item{min_boot}{The minimum bootstrap value for taxonomy assignment.}

\item{n_ranks}{The number of ranks present in the reference database.}
}
\value{
An updated tidyamplicons object.
}
\description{
This function requires the DADA2 package to be installed.
}
\details{
This function will (re)classify either all or a subset of the taxa, given
that a variable is present in the taxon table that contains (representative)
sequences of the taxa.

Ranks can be supplied as a named integer vector, where the names represent
taxonomic ranks and the integers represent positions of these ranks in the
taxonomy strings present in the reference database. Ranks can also be
supplied as just a character vector with the rank names; in that case, it is
assumed that the database taxonomy string follows the default order {domain,
phylum, class, order, family, genus, species}. If no ranks are supplied, taxa
will be (re)classified at all default ranks.
}
