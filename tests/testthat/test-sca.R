# Tests for sca(), the workhorse of the package

test_that("sca() returns one row per specification with the expected columns", {
  s <- suppressMessages(sca(y = "Salnty", x = "T_degC",
                            controls = c("O2Sat", "STheta"),
                            data = bottles, progressBar = FALSE,
                            parallel = FALSE))
  expect_s3_class(s, "data.frame")
  # 2 controls -> 3 specifications
  expect_equal(nrow(s), 3)
  expect_true(all(c("coef", "se", "statistic", "p", "RMSE", "adjR",
                    "sig.level", "index", "terms", "control_coefs") %in%
                    names(s)))
  # index should be a 1..n sequence
  expect_equal(sort(s$index), seq_len(nrow(s)))
  # dummy indicator columns for each control are added
  expect_true(all(c("O2Sat", "STheta") %in% names(s)))
})

test_that("sca(returnFormulae = TRUE) returns formulae instead of estimates", {
  f <- sca(y = "Salnty", x = "T_degC", controls = c("O2Sat", "STheta"),
           data = bottles, returnFormulae = TRUE)
  expect_length(f, 3)
})

test_that("sca() supports fixed effects via feols", {
  fe <- suppressMessages(sca(y = "Salnty", x = "T_degC", controls = "STheta",
                             data = bottles, fixedEffects = "Sta_ID",
                             progressBar = FALSE))
  expect_s3_class(fe, "data.frame")
  expect_equal(nrow(fe), 1)
  expect_true(all(c("RMSE", "adjR") %in% names(fe)))
})

test_that("sca() supports glm families and reports AIC/deviance", {
  d <- bottles
  d$bin <- as.integer(d$Salnty > median(d$Salnty, na.rm = TRUE))
  g <- suppressMessages(sca(y = "bin", x = "T_degC", controls = "STheta",
                            data = d, family = "binomial", link = "logit",
                            progressBar = FALSE))
  expect_s3_class(g, "data.frame")
  expect_true(all(c("AIC", "deviance") %in% names(g)))
})

test_that("sca() warns for fixed effects + non-linear family and then ignores them", {
  d <- bottles
  d$bin <- as.integer(d$Salnty > median(d$Salnty, na.rm = TRUE))
  expect_warning(
    res <- suppressMessages(sca(y = "bin", x = "T_degC", controls = "STheta",
                                data = d, family = "binomial", link = "logit",
                                fixedEffects = "Sta_ID", progressBar = FALSE)),
    "Fixed effects unsupported"
  )
  # The warning promises to ignore the fixed effects, so estimation should fall
  # back to the glm path and return a valid result rather than erroring.
  expect_s3_class(res, "data.frame")
  expect_true("AIC" %in% names(res))
})
