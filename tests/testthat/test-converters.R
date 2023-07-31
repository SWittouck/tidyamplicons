ta_test <- read_tidyamplicons(test_path("data/urt"))

test_that("Can read tidyamplicons object.", {
  ta <- read_tidyamplicons(test_path("data/urt"))
  expect_equal(attr(ta, "class"), "tidyamplicons")
})

test_that("Can update old tidyamplicons object to new format.", {
  load(test_path("data/urt.rda"))
  expect_false("sample_id" %in% colnames(urt$samples))
  ta <- update_tidyamplicons(urt)
  expect_true("sample_id" %in% colnames(ta$samples))
})

test_that("Can save a tidyamplicons object.", {
  expect_no_warning(write_tidyamplicons(ta_test, "test"))
  on.exit(unlink("test", recursive=TRUE), add=TRUE, after=FALSE)
})

test_that("Can convert to phyloseq object.", {
  ta_phylo <- as_phyloseq(ta_test, sample=sample_id, taxon=taxon_id)
  expect_true(class(ta_phylo) == "phyloseq")
})

test_that("Can convert abundances to abundances matrix", {
  ab_mat <- ta_test %>% abundances() %>% as_abundances_matrix()
  expected_width <- 1957
  expected_height <- 214
  width <- dim(ab_mat)[2]
  height <- dim(ab_mat)[1]
  expect_equal(width, expected_width)
  expect_equal(height, expected_height)
})

test_that("Can merge two tidyamplicons", {
  ta_merged <- merge_tidyamplicons(ta_test, ta_test)
  
  final_sample_id <- "s434"
  expect_equal(
    ta_merged$samples$sample_id[length(ta_merged$samples$sample_id)], 
    final_sample_id
  )
})