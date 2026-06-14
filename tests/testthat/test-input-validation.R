# Tests for input validation and family/link handling

test_that("sca() gives a clear error for columns missing from data", {
  expect_error(
    suppressMessages(sca("Salnty", "TYPO", "STheta", bottles, progressBar = FALSE)),
    "not found in data"
  )
  expect_error(
    suppressMessages(sca("Salnty", "T_degC", c("STheta", "NOPE"), bottles,
                         progressBar = FALSE)),
    "not found in data"
  )
  expect_error(
    suppressMessages(sca("Salnty", "T_degC", "STheta", bottles,
                         fixedEffects = "NOPE", progressBar = FALSE)),
    "Fixed-effects"
  )
  expect_error(
    suppressMessages(sca("Salnty", "T_degC", "STheta", bottles,
                         weights = "NOPE", progressBar = FALSE)),
    "Weights"
  )
})

test_that("sca() validation does not interfere with returnFormulae", {
  # returnFormulae does not touch data, so it should not require valid columns.
  expect_silent(
    f <- sca("y", "x", c("a", "b"), data = data.frame(z = 1),
             returnFormulae = TRUE)
  )
  expect_length(f, 3)
})

test_that("sca() treats family = 'gaussian' as OLS", {
  g <- suppressMessages(sca("Salnty", "T_degC", "STheta", bottles,
                            family = "gaussian", progressBar = FALSE))
  l <- suppressMessages(sca("Salnty", "T_degC", "STheta", bottles,
                            family = "linear", progressBar = FALSE))
  expect_equal(g, l)
  expect_true("RMSE" %in% names(g))
})

test_that("sca() defaults the link for a glm family and errors on a bad family", {
  d <- bottles
  d$bin <- as.integer(d$Salnty > median(d$Salnty, na.rm = TRUE))
  # binomial with no link supplied should default to the canonical logit link.
  g <- suppressMessages(sca("bin", "T_degC", "STheta", d, family = "binomial",
                            progressBar = FALSE))
  expect_true("AIC" %in% names(g))
  expect_error(
    suppressMessages(sca("bin", "T_degC", "STheta", d, family = "binomail",
                         progressBar = FALSE)),
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
