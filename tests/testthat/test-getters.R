ta_test <- read_tidytacos(test_path("data/urt"))
len_taxa <- length(ta_test$taxa$taxon_id)
len_samples <- 214 #some empty samples in the df

test_that("Correct numbers reported", {
    expect_snapshot(ta_test %>% report_numbers())
})

test_that("Correct numbers reported", {
    expect_snapshot(ta_test %>% numbers())
})

test_that("Relative abundance matrix can be created", {
    relabun_mat <- ta_test %>% get_rel_abundance_matrix()
    expect_lte(max(relabun_mat), 1)
    expect_equal(dim(relabun_mat), c(len_samples, len_taxa))
})

test_that("Beta diversity table can be generated", {
    beta_table <- ta_test %>% betas()
    expect_lte(max(beta_table$beta), 1)
    expect_equal(dim(beta_table), c(22791, 17))
})

test_that("Occurences are calculated",{
    occ <- ta_test %>% occurrences("location")
    expect_equal(dim(occ), c(2*len_taxa, 3))
    expect_true(all(c("location","occurrence") %in% names(occ)))
    expect_equal(sum(occ$occurrence), 7693)
})

test_that("Presence/absense is calculated", {
    pres_ab <- ta_test %>% occurrences("location", pres_abs=TRUE)
    expect_equal(dim(pres_ab), c(4*len_taxa, 4))
    expect_true(all(c("location", "presence", "n") %in% names(pres_ab)))
    expect_equal(sum(pres_ab$n), len_samples * len_taxa)
})

test_that("Mean relative abundance is computed", {
    mra <- ta_test %>% mean_rel_abundances()
    expect_equal(sum(mra$mean_rel_abundance), 1)
    expect_equal(length(mra$taxon_id), len_taxa)
})

test_that("Mean relative abundance is computed for a condition", {
    mra <- ta_test %>% mean_rel_abundances("location")
    cond_opt <- unique(ta_test$samples$location)
    expect_equal(sum(mra$mean_rel_abundance), length(cond_opt))
    expect_equal(length(mra$taxon_id), len_taxa*length(cond_opt))
})

test_that("Extract sample tibble", {
    expect_equal(ta_test %>% samples(), ta_test$samples)
})

test_that("Extract taxon tibble", {
    expect_equal(ta_test %>% taxa(), ta_test$taxa)
})

test_that("Extract abundances tibble", {
    expect_equal(ta_test %>% abundances(), ta_test$abundances)
})

test_that("Extract all three tables in one large tibble",{
    wide_ta <- ta_test %>% everything()
    expect_equal(dim(wide_ta), c(7693, 18))
})

test_that("Perform adonis shows stable output", {
    set.seed(42)
    adn <- ta_test %>% perform_adonis(c("location", "method"))
    expect_snapshot(adn)
})