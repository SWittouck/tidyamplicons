
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
  ta <- read_tidyamplicons(test_path("data/urt"))
  expect_no_warning(write_tidyamplicons(ta, "test"))
  on.exit(unlink("test", recursive=TRUE), add=TRUE, after=FALSE)
})