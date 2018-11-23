# tidyamplicons: Functions to Manipulate and Visualize Amplicon Count Data

Tidyamplicons is an R package for the analysis of amplicon count data: abundances of amplicon sequences (either clustered in OTUs or [exact variants](https://www.ncbi.nlm.nih.gov/pubmed/28731476)) in samples. The package builds on the [tidyverse](https://www.tidyverse.org/) created by [Hadley Wickham](http://hadley.nz/): the data are stored in tidy tables where each row is an observation and each column a variable. In addition, the package supplies a set of "verbs": functions that take a tidyamplicons object as first argument and also return a tidyamplicons object. Not all functionality is currently implemented in the form of verbs, but this will soon be remediated. 

# Installation

Run the following (only install devtools if not yet installed): 

```R
# install.packages("devtools")
devtools::install_github("SWittouck/tidyamplicons")
```

# Documentation

Full documentation for this package is not yet available. A vignette is included with a basic explanation on the package design philosophy and an illustration of some common workflows on an example dataset. To make sure the vignette is built when you install the package, install it in the following way:

```R
devtools::install_github("SWittouck/tidyamplicons", build_vignettes = TRUE)
```
