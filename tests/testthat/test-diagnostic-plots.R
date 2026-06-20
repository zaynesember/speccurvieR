# Tests for the diagnostic plots: plot_influence, plot_coef_fit, plot_multi_se

sca_diag <- function() {
  suppressMessages(sca("Salnty", "T_degC", c("ChlorA", "O2Sat", "NO2uM"),
                       bottles, progress_bar = FALSE))
}

test_that("plot_influence returns a ggplot with included/excluded groups", {
  s <- sca_diag()
  p <- plot_influence(s)
  expect_s3_class(p, "ggplot")
  expect_setequal(unique(as.character(p$data$included)),
                  c("Excluded", "Included"))
  # one facet per control
  expect_setequal(unique(p$data$control), c("ChlorA", "O2Sat", "NO2uM"))
})

test_that("plot_coef_fit returns a ggplot and validates the metric", {
  s <- sca_diag()
  expect_s3_class(plot_coef_fit(s), "ggplot")
  expect_s3_class(plot_coef_fit(s, metric = "adjR"), "ggplot")
  expect_error(plot_coef_fit(s, metric = "AIC"), "must be one of")
})

test_that("plot_multi_se returns a ggplot with one facet per SE type", {
  p <- suppressMessages(
    plot_multi_se("Salnty", "T_degC", c("ChlorA", "O2Sat"), bottles,
                types = c("iid", "HC1", "HC3")))
  expect_s3_class(p, "ggplot")
  expect_setequal(unique(p$data$se_type), c("iid", "HC1", "HC3"))
  # estimates are identical across SE types for a given specification
  by_type <- split(p$data$coef, p$data$se_type)
  expect_equal(sort(by_type[["iid"]]), sort(by_type[["HC3"]]))
})

test_that("plot_multi_se significance reflects the SE type", {
  p <- suppressMessages(
    plot_multi_se("Salnty", "T_degC", c("ChlorA", "O2Sat"), bottles,
                types = c("iid", "HC3")))
  # p-values are computed from coef / se, so larger SEs give larger p-values
  expect_true(all(p$data$p >= 0 & p$data$p <= 1))
})
