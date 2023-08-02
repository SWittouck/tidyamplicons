# Tests could be expanded with proto testing of diff layers,
# but for now a heuristic test
# to see if the same plot is generate will suffice

ta_test <- read_tidytacos(test_path("data/urt"))

test_that("Barplot returns identical plot", {
    bp  <- ta_test %>% bar_plot()
    vdiffr::expect_doppelganger("Default barplot", bp)
})

test_that("Barplot raises warning when aggregating samples", {
    expect_warning(bp <- ta_test %>% bar_plot(n=5, x=participant))
    vdiffr::expect_doppelganger("Custom barplot", bp)
})

test_that("Barplot raises error when providing non-existant label", {
    expect_error(ta_test %>% bar_plot(x=imagined))
})

test_that("Can create venndiagram", {
    skip_if_not_installed("ggVenDiagram")
    venn <- ta_test %>% tacoplot_venn(location)
    vdiffr::expect_doppelganger("Venndiagram", venn)
})
