ta_test <- read_tidyamplicons(test_path("data/urt"))

test_that("Correct numbers reported", {
    expect_snapshot(ta_test %>% report_numbers())
})

test_that("Correct numbers reported", {
    expect_snapshot(ta_test %>% numbers())
})

test_that("Relative abundance matrix can be created", {
    relabun_mat <- ta_test %>% get_rel_abundance_matrix()
    expect_lte(max(relabun_mat), 1)
    expect_equal(dim(relabun_mat), c(214, 1957))
})

test_that("Beta diversity table can be generated", {
    beta_table <- ta_test %>% get_betas()
    expect_lte(max(beta_table$beta), 1)
    expect_equal(dim(beta_table), c(22791, 17))
})