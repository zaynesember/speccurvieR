# Tests for plot_se()

test_that("plot_se returns a ggplot for a basic se_compare output", {
  s <- suppressMessages(se_compare("Salnty ~ T_degC + STheta", bottles,
                                   types = c("iid", "HC0", "HC3")))
  p <- plot_se(s)
  expect_s3_class(p, "ggplot")
})

test_that("plot_se handles fixed-effects + clustered output", {
  s <- suppressMessages(se_compare("Salnty ~ T_degC + STheta | Sta_ID", bottles,
                                   types = "all", cluster = "Depth_ID"))
  p <- plot_se(s)
  expect_s3_class(p, "ggplot")
})

test_that("plot_se excludes the intercept by default and can include it", {
  s <- suppressMessages(se_compare("Salnty ~ T_degC + STheta", bottles,
                                   types = c("iid", "HC0")))
  without <- plot_se(s)
  with_int <- plot_se(s, intercept = TRUE)
  terms_without <- unique(without$data$term)
  terms_with <- unique(with_int$data$term)
  expect_false("(Intercept)" %in% terms_without)
  expect_true("(Intercept)" %in% terms_with)
})

test_that("plot_se significance flag reflects the confidence level", {
  s <- suppressMessages(se_compare("Salnty ~ T_degC + STheta", bottles,
                                   types = c("iid", "HC0", "HC3")))
  p <- plot_se(s, level = 0.95)
  # the built data frame carries lower/upper interval bounds and a sig flag
  expect_true(all(c("lower", "upper", "sig") %in% names(p$data)))
  expect_equal(p$data$sig, (p$data$lower > 0) | (p$data$upper < 0))
})

test_that("plot_se returns NULL with a message when there are no SE columns", {
  # a data frame with only an estimate column has nothing to compare
  df <- data.frame(estimate = c(1, 2), row.names = c("a", "b"))
  expect_message(res <- plot_se(df), "No standard error columns")
  expect_null(res)
})
