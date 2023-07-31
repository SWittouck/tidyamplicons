ta_test <- read_tidytacos(test_path("data/urt"))

# RAREFY
test_that("Raise error on rarification higher 
than lowest abundance when replace=FALSE.", {
expect_error(ta_test %>% rarefy(100))
})

test_that("Max abundance equals n", {
    ab <- ta_test %>% rarefy(10) %>% 
    abundances %>% dplyr::pull(abundance)
    expect_equal(max(ab), 10)
})

test_that("Max abundance equals n with replace=TRUE", {
    ab <- ta_test %>% rarefy(100, replace=TRUE) %>% 
    abundances %>% dplyr::pull(abundance)
    expect_equal(max(ab), 100)
})

# CHANGE_ID_SAMPLES
test_that("Change sample ID based on unique column",{
    uq_part <- ta_test %>% filter_samples(location == "NF", method=="A")
    # remove empty samples, as they will cause a mismatch
    uq_part <- uq_part %>% remove_empty_samples()
    uq_part_id <- uq_part %>% change_id_samples(participant)
    expect_true(identical(
        sort(unique(uq_part$samples$participant)), 
        sort(unique(uq_part_id$abundances$sample_id))
    ))
})

test_that("Raise error on non-unique sample_id",{
    expect_error(ta_test %>% change_id_samples(condition))
})

# CHANGE_ID_TAXA

test_that("Change sample ID based on unique column",{
    
    ta_test_id <- ta_test %>% change_id_taxa(sequence)
    expect_true(identical(
        sort(unique(ta_test_id$abundances$taxon_id)), 
        sort(unique(ta_test$taxa$sequence))
    ))
})

test_that("Raise error on non-unique taxon_id",{
    expect_error(ta_test %>% change_id_taxa(genus))
})

# AGG SAMPLES
nsamples <- 217
nsamples_agg <- 139

test_that("No aggregation of sample table when 
    metadata is unique for every row", {
    sample_id_before <- ta_test$samples$sample_id
    ta_agg <- ta_test %>% aggregate_samples()
    sample_id_after <- ta_agg$samples$sample_id
    expect_length(sample_id_after, nsamples)
    expect_length(sample_id_before, nsamples)
})

test_that("Aggregation of sample table succeeds", {
    sample_id_before <- ta_test$samples$sample_id
    # Make one column not unique to test agg
    ta_agg <- ta_test %>% 
        mutate_samples(location = "NF") %>% 
        aggregate_samples()
    sample_id_after <- ta_agg$samples$sample_id
    expect_length(sample_id_before, nsamples)
    expect_length(sample_id_after, nsamples_agg)
    # test abundances table for irregularities
    expect_lte(
        length(unique(ta_agg$abundances$sample_id)), 
        nsamples_agg
    )
})

# AGG TAXA
ntaxa <- 1957
ntaxa_agg <- 1066
ntaxa_agg_phylum <- 32

test_that("No aggregation of taxa table when metadata 
    is unique for every row", {
    taxon_id_before <- ta_test$taxa$taxon_id
    ta_agg <- ta_test %>% aggregate_taxa()
    taxon_id_after <- ta_agg$taxa$taxon_id
    expect_length(taxon_id_after, ntaxa)
    expect_length(taxon_id_before, ntaxa)
})

test_that("Aggregation of taxa table succeeds when no rank is specified", {
    taxon_id_before <- ta_test$taxa$taxon_id
    # Make one column not unique to test agg
    ta_agg <- ta_test %>% 
        mutate_taxa(sequence = 0) %>% 
        aggregate_taxa()
    taxon_id_after <- ta_agg$taxa$taxon_id
    expect_length(taxon_id_before, ntaxa)
    expect_length(taxon_id_after, ntaxa_agg)
    # test abundances table for irregularities
    expect_lte(
        length(unique(ta_agg$abundances$taxon_id)), 
        ntaxa_agg
    )
})

test_that("Aggregation of taxa table succeeds when rank is specified", {
    taxon_id_before <- ta_test$taxa$taxon_id
    # Make one column not unique to test agg
    ta_agg <- ta_test %>% aggregate_taxa(rank="phylum")
    taxon_id_after <- ta_agg$taxa$taxon_id
    expect_length(taxon_id_before, ntaxa)
    expect_length(taxon_id_after, ntaxa_agg_phylum)
    # test abundances table for irregularities
    expect_lte(
        length(unique(ta_agg$abundances$taxon_id)), 
        ntaxa_agg_phylum
    )
})

# TRIM ASVS

test_that("Trim ASVs works", {
    ta_trim <- ta_test %>% trim_asvs(1, 5)
    longest_seq <- max(str_length(ta_trim$taxa$sequence))
    expect_lte(longest_seq, 5)
})


# SELECT SAMPLES

test_that("Removing sample_id throws error", {
expect_error(ta_test %>% select_samples(-sample_id))
})

test_that("Removing a column from samples is succesful", {
    expect_true("method" %in% names(ta_test$samples))
    ta_test <- ta_test %>% select_samples(-method)
    expect_true(!"method" %in% names(ta_test$samples))
})

# SELECT TAXA

test_that("Removing taxon_id throws error", {
expect_error(ta_test %>% select_taxa(-taxon_id))
})

test_that("Removing a column from taxa is succesful", {
    expect_true("phylum" %in% names(ta_test$taxa))
    ta_test <- ta_test %>% select_taxa(-phylum)
    expect_true(!"phylum" %in% names(ta_test$taxa))
})

# SELECT ABUNDANCES

# Add removable column
ta_test$abundances <- ta_test$abundances %>% mutate(remove_me = 0)

test_that("Removing sample_id throws error", {
expect_error(ta_test %>% select_abundances(-sample_id))
})

test_that("Removing taxon_id throws error", {
expect_error(ta_test %>% select_abundances(-taxon_id))
})

test_that("Removing abundance column throws error", {
expect_error(ta_test %>% select_abundances(-abundance))
})

test_that("Removing a column from abundances is succesful", {
    expect_true("remove_me" %in% names(ta_test$abundances))
    ta_test <- ta_test %>% select_abundances(-remove_me)
    expect_true(!"remove_me" %in% names(ta_test$abundances))
})

# MUTATE TAXA
test_that("Mutating taxa", {
    ta_gs <- ta_test %>% mutate_taxa(new_col = paste(genus, species))
    expect_true("new_col" %in% names(ta_gs$taxa))
})

test_that("Mutating taxa, where taxon_id is removed throws error", {
    expect_error(ta_test %>% mutate_taxa(taxon_id = NULL))
})

# MUTATE SAMPLES
test_that("Mutating samples", {
    ta_con <- ta_test %>% mutate_samples(new_col = paste(condition, location))
    expect_true("new_col" %in% names(ta_con$samples))
})

test_that("Mutating samples, where sample_id is removed throws error", {
    expect_error(ta_test %>% mutate_samples(sample_id = NULL))
})

# MUTATE SAMPLES
test_that("Mutating abundances", {
    ta_st <- ta_test %>% mutate_abundances(new_col = paste(sample_id, taxon_id))
    expect_true("new_col" %in% names(ta_st$abundances))
})

test_that("Mutating abundances, where sample_id is removed throws error", {
    expect_error(ta_test %>% mutate_abundances(sample_id = NULL))
})

test_that("Mutating abundances, where taxon_id is removed throws error", {
    expect_error(ta_test %>% mutate_abundances(taxon_id = NULL))
})

test_that("Mutating abundances, where abundances is removed throws error", {
    expect_error(ta_test %>% mutate_abundances(abundance = NULL))
})

test_that("Filtering out all samples raises error", {
    expect_error(ta_test %>% filter_samples(condition == "Non-existent"))
})

test_that("Filtering returns the expected amount of samples", {
    ta_f <- ta_test %>% filter_samples(method == "S", location == "NF")
    expect_equal(length(ta_f$samples$sample_id), 96)
})

test_that("Filtering out all taxa raises error", {
    expect_error(ta_test %>% filter_taxa(species == "Impossiblum"))
})

test_that("Filtering returns the expected amount of taxa", {
    ta_f <- ta_test %>% filter_taxa(species == "coli")
    expect_equal(length(ta_f$taxa$taxon_id), 1)
})

test_that("Filtering out all abundances raises error", {
    expect_error(ta_test %>% filter_abundances(abundance == "Impossible"))
})

test_that("Filtering returns the expected amount of samples", {
    ta_f <- ta_test %>% filter_abundances(abundance >= 10000)
    expect_equal(length(ta_f$abundances$sample_id), 70)
    expect_equal(length(ta_f$taxa$taxon_id), 18)
})