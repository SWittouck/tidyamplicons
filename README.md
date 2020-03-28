# tidyamplicons

## Overview

Tidyamplicons is an R package for the analysis of amplicon abundance data: counts of amplicon sequences (either clustered in OTUs or [exact variants](https://www.ncbi.nlm.nih.gov/pubmed/28731476)) in samples. The package builds on the [tidyverse](https://www.tidyverse.org/) created by [Hadley Wickham](http://hadley.nz/): the data are stored in tidy tables where each row is an observation and each column a variable. In addition, the package supplies a set of "verbs": functions that take a tidyamplicons object as first argument and also return a tidyamplicons object.

## Installation

Run the following (only install devtools if not yet installed): 

```R
# install.packages("devtools")
devtools::install_github("SWittouck/tidyamplicons", ref = "v0.2.1")
```

## Documentation

Documentation is available for most functions but not all, and it is still very preliminary. A vignette is included with a basic explanation on the package design philosophy and an illustration of some common workflows on an example dataset. To make sure the vignette is built when you install the package, install it in the following way:

```R
devtools::install_github(
  "SWittouck/tidyamplicons", ref = "v0.2.1", build_vignettes = TRUE
)
```
