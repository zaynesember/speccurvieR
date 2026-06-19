# Tests for plotSE()

test_that("plotSE returns a ggplot for a basic se_compare output", {
  s <- suppressMessages(se_compare("Salnty ~ T_degC + STheta", bottles,
                                   types = c("iid", "HC0", "HC3")))
  p <- plotSE(s)
  expect_s3_class(p, "ggplot")
})

test_that("plotSE handles fixed-effects + clustered output", {
  s <- suppressMessages(se_compare("Salnty ~ T_degC + STheta | Sta_ID", bottles,
                                   types = "all", cluster = "Depth_ID"))
  p <- plotSE(s)
  expect_s3_class(p, "ggplot")
})

test_that("plotSE excludes the intercept by default and can include it", {
  s <- suppressMessages(se_compare("Salnty ~ T_degC + STheta", bottles,
                                   types = c("iid", "HC0")))
  without <- plotSE(s)
  with_int <- plotSE(s, intercept = TRUE)
  terms_without <- unique(without$data$term)
  terms_with <- unique(with_int$data$term)
  expect_false("(Intercept)" %in% terms_without)
  expect_true("(Intercept)" %in% terms_with)
})

test_that("plotSE significance flag reflects the confidence level", {
  s <- suppressMessages(se_compare("Salnty ~ T_degC + STheta", bottles,
                                   types = c("iid", "HC0", "HC3")))
  p <- plotSE(s, level = 0.95)
  # the built data frame carries lower/upper interval bounds and a sig flag
  expect_true(all(c("lower", "upper", "sig") %in% names(p$data)))
  expect_equal(p$data$sig, (p$data$lower > 0) | (p$data$upper < 0))
})

test_that("plotSE returns NULL with a message when there are no SE columns", {
  # a data frame with only an estimate column has nothing to compare
  df <- data.frame(estimate = c(1, 2), row.names = c("a", "b"))
  expect_message(res <- plotSE(df), "No standard error columns")
  expect_null(res)
})
