ta.test <- read_tidytacos(test_path("data/urt"))

test_that("Gives default rank names when not provided", {
  expect_equal(ta.test %>% rank_names(), 
  c("kingdom", "phylum", "class", "order", "family", "genus")
  )
})

test_that("Gives custom rank names when provided", {
  ta.test <- ta.test %>% set_rank_names(
    c("domain","phylum","order","family", "genus", "species")
  ) 
  expect_equal(ta.test %>% rank_names(), 
  c("domain","phylum","order","family", "genus", "species")
  )
})