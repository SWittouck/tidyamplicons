ta.test <- read_tidyamplicons(test_path("data/urt"))

# RAREFY
test_that("Raise error on rarification higher 
than lowest abundance when replace=FALSE.", {
expect_error(ta.test %>% rarefy(100))
})

test_that("Max abundance equals n", {
    ab <- ta.test %>% rarefy(10) %>% 
    abundances %>% dplyr::pull(abundance)
    expect_equal(max(ab), 10)
})

test_that("Max abundance equals n with replace=TRUE", {
    ab <- ta.test %>% rarefy(100, replace=TRUE) %>% 
    abundances %>% dplyr::pull(abundance)
    expect_equal(max(ab), 100)
})

# CHANGE_ID_SAMPLES
test_that("Change sample ID based on unique column",{
    uq_part <- ta.test %>% filter_samples(location == "NF", method=="A")
    # remove empty samples, as they will cause a mismatch
    uq_part <- uq_part %>% remove_empty_samples()
    uq_part_id <- uq_part %>% change_id_samples(participant)
    expect_true(identical(
        sort(unique(uq_part$samples$participant)), 
        sort(unique(uq_part_id$abundances$sample_id))
    ))
})

test_that("Raise error on non-unique sample_id",{
    expect_error(ta.test %>% change_id_samples(condition))
})

# CHANGE_ID_TAXA

test_that("Change sample ID based on unique column",{
    
    ta.test_id <- ta.test %>% change_id_taxa(sequence)
    expect_true(identical(
        sort(unique(ta.test_id$abundances$taxon_id)), 
        sort(unique(ta.test$taxa$sequence))
    ))
})

test_that("Raise error on non-unique taxon_id",{
    expect_error(ta.test %>% change_id_taxa(genus))
})

# AGG SAMPLES
nsamples <- 217
nsamples_agg <- 139

test_that("No aggregation of sample table when metadata is unique for every row", {
    sample_id_before <- ta.test$samples$sample_id
    ta.agg <- ta.test %>% aggregate_samples()
    sample_id_after <- ta.agg$samples$sample_id
    expect_length(sample_id_after, nsamples)
    expect_length(sample_id_before, nsamples)
})

test_that("Aggregation of sample table succeeds", {
    sample_id_before <- ta.test$samples$sample_id
    # Make one column not unique to test agg
    ta.agg <- ta.test %>% 
        mutate_samples(location = "NF") %>% 
        aggregate_samples()
    sample_id_after <- ta.agg$samples$sample_id
    expect_length(sample_id_before, nsamples)
    expect_length(sample_id_after, nsamples_agg)
    # test abundances table for irregularities
    expect_lte(
        length(unique(ta.agg$abundances$sample_id)), 
        nsamples_agg
    )
})

# AGG TAXA
ntaxa <- 1957
ntaxa_agg <- 1066
ntaxa_agg_phylum <- 32

test_that("No aggregation of taxa table when metadata is unique for every row", {
    taxon_id_before <- ta.test$taxa$taxon_id
    ta.agg <- ta.test %>% aggregate_taxa()
    taxon_id_after <- ta.agg$taxa$taxon_id
    expect_length(taxon_id_after, ntaxa)
    expect_length(taxon_id_before, ntaxa)
})

test_that("Aggregation of taxa table succeeds when no rank is specified", {
    taxon_id_before <- ta.test$taxa$taxon_id
    # Make one column not unique to test agg
    ta.agg <- ta.test %>% 
        mutate_taxa(sequence = 0) %>% 
        aggregate_taxa()
    taxon_id_after <- ta.agg$taxa$taxon_id
    expect_length(taxon_id_before, ntaxa)
    expect_length(taxon_id_after, ntaxa_agg)
    # test abundances table for irregularities
    expect_lte(
        length(unique(ta.agg$abundances$taxon_id)), 
        ntaxa_agg
    )
})

test_that("Aggregation of taxa table succeeds when rank is specified", {
    taxon_id_before <- ta.test$taxa$taxon_id
    # Make one column not unique to test agg
    ta.agg <- ta.test %>% aggregate_taxa(rank="phylum")
    taxon_id_after <- ta.agg$taxa$taxon_id
    expect_length(taxon_id_before, ntaxa)
    expect_length(taxon_id_after, ntaxa_agg_phylum)
    # test abundances table for irregularities
    expect_lte(
        length(unique(ta.agg$abundances$taxon_id)), 
        ntaxa_agg_phylum
    )
})

# TRIM ASVS

test_that("Trim ASVs works", {
    ta.trim <- ta.test %>% trim_asvs(1, 5)
    longest_seq <- max(str_length(ta.trim$taxa$sequence))
    expect_lte(longest_seq, 5)
})


# SELECT SAMPLES

test_that("Removing sample_id throws error", {
expect_error(ta.test %>% select_samples(-sample_id))
})

test_that("Removing a column from samples is succesful", {
    expect_true("method" %in% names(ta.test$samples))
    ta.test <- ta.test %>% select_samples(-method)
    expect_true(!"method" %in% names(ta.test$samples))
})

# SELECT TAXA

test_that("Removing taxon_id throws error", {
expect_error(ta.test %>% select_taxa(-taxon_id))
})

test_that("Removing a column from taxa is succesful", {
    expect_true("phylum" %in% names(ta.test$taxa))
    ta.test <- ta.test %>% select_taxa(-phylum)
    expect_true(!"phylum" %in% names(ta.test$taxa))
})

# SELECT ABUNDANCES

# Add removable column
ta.test$abundances <-ta.test$abundances %>% mutate(remove_me = 0)

test_that("Removing sample_id throws error", {
expect_error(ta.test %>% select_abundances(-sample_id))
})

test_that("Removing taxon_id throws error", {
expect_error(ta.test %>% select_abundances(-taxon_id))
})

test_that("Removing abundance column throws error", {
expect_error(ta.test %>% select_abundances(-abundance))
})

test_that("Removing a column from abundances is succesful", {
    expect_true("remove_me" %in% names(ta.test$abundances))
    ta.test <- ta.test %>% select_abundances(-remove_me)
    expect_true(!"remove_me" %in% names(ta.test$abundances))
})
