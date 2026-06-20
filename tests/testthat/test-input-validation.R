# Tests for input validation and family/link handling

test_that("sca() gives a clear error for columns missing from data", {
  expect_error(
    suppressMessages(sca("Salnty", "TYPO", "STheta", bottles, progress_bar = FALSE)),
    "not found in data"
  )
  expect_error(
    suppressMessages(sca("Salnty", "T_degC", c("STheta", "NOPE"), bottles,
                         progress_bar = FALSE)),
    "not found in data"
  )
  expect_error(
    suppressMessages(sca("Salnty", "T_degC", "STheta", bottles,
                         fixed_effects = "NOPE", progress_bar = FALSE)),
    "Fixed-effects"
  )
  expect_error(
    suppressMessages(sca("Salnty", "T_degC", "STheta", bottles,
                         weights = "NOPE", progress_bar = FALSE)),
    "Weights"
  )
})

test_that("sca() validation does not interfere with return_formulae", {
  # return_formulae does not touch data, so it should not require valid columns.
  expect_silent(
    f <- sca("y", "x", c("a", "b"), data = data.frame(z = 1),
             return_formulae = TRUE)
  )
  expect_length(f, 3)
})

test_that("sca() treats family = 'gaussian' as OLS", {
  g <- suppressMessages(sca("Salnty", "T_degC", "STheta", bottles,
                            family = "gaussian", progress_bar = FALSE))
  l <- suppressMessages(sca("Salnty", "T_degC", "STheta", bottles,
                            family = "linear", progress_bar = FALSE))
  expect_equal(g, l)
  expect_true("RMSE" %in% names(g))
})

test_that("sca() defaults the link for a glm family and errors on a bad family", {
  d <- bottles
  d$bin <- as.integer(d$Salnty > median(d$Salnty, na.rm = TRUE))
  # binomial with no link supplied should default to the canonical logit link.
  g <- suppressMessages(sca("bin", "T_degC", "STheta", d, family = "binomial",
                            progress_bar = FALSE))
  expect_true("AIC" %in% names(g))
  expect_error(
    suppressMessages(sca("bin", "T_degC", "STheta", d, family = "binomail",
                         progress_bar = FALSE)),
    "not a recognised model family"
  )
})

test_that("se_compare() validates weights and formula columns", {
  expect_error(
    suppressMessages(se_compare("Salnty ~ T_degC", bottles, weights = "NOPE")),
    "Weights"
  )
  expect_error(
    suppressMessages(se_compare("Salnty ~ TYPO", bottles, types = "HC0")),
    "not found in data"
  )
})
