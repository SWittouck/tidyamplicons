# SparCC - MCL
test_that("Can run sparcc-MCL worfklow", {

urt_clust <- urt %>% add_sparcc_network_clusters(taxon_name=sequence)
expect_true("cluster" %in% names(urt_clust$taxa))

})

# Eigentaxa
test_that("Can add eigentaxa", {
    urt_eigen <- urt %>% add_eigentaxa(taxon_name=sequence)
    expect_gt(length(urt_eigen$samples), length(urt$samples))
})
