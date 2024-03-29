x <- matrix(
  c(1500, 1300, 280, 356),
  ncol = 2
)
rownames(x) <- c("taxon1", "taxon2")
colnames(x) <- c("sample1", "sample2")
test_data <- create_tidyamplicons(x, 
    taxa_are_columns=FALSE)

test_that("Can add sample tibble to ta object",{
    sample <- c("sample1","sample2")
    environment <- c("food fermentation", "human stool")
    smp_tibble <- tibble::tibble(sample, environment)
    suppressMessages(
        test_data <- test_data %>% add_sample_tibble(smp_tibble)
    )
    expect_true("environment" %in% names(test_data$samples))
})

test_that("Can add lib sizes", {
    ta_lib <- test_data %>% add_lib_size()
    expect_equal(ta_lib$samples$lib_size, c(2800, 636))
})

test_that("Can add alpha diversity metrics", {
    ta_alpha <- test_data %>% add_alphas()
    expect_equal(ta_alpha$samples$inverse_simpson, 
                    c(1.989, 1.971), tolerance=1e-3)
})

