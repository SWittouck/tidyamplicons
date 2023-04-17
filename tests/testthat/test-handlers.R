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

# AGG TAXA

# TRIM ASVS

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