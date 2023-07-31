# Tests could be expanded with proto testing of diff layers,
# but for now a heuristic test
# to see if the same plot is generate will suffice

ta_test <- read_tidytacos(test_path("data/urt"))

test_that("Barplot returns identical plot", {
    bp  <- ta_test %>% bar_plot()
    vdiffr::expect_doppelganger("Default barplot", bp)
})

test_that("Barplot returns identical plot when different arguments are used", {
    bp  <- ta_test %>% bar_plot(n=5, x=participant)
    vdiffr::expect_doppelganger("Custom barplot", bp)
})