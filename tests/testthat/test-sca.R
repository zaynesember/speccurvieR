# Tests for sca(), the workhorse of the package

test_that("sca() returns one row per specification with the expected columns", {
  s <- suppressMessages(sca(y = "Salnty", x = "T_degC",
                            controls = c("O2Sat", "STheta"),
                            data = bottles, progress_bar = FALSE,
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

test_that("sca() control-indicator columns use exact term membership", {
  # Regression test: indicators were built with substring matching on a
  # deparsed terms list, so a control name that is a substring of another term
  # (e.g. "O2" inside "O2Sat") set a false 1. They must be exact membership.
  d <- bottles
  d$O2 <- d$O2Sat * 0.5
  s <- suppressMessages(sca(y = "Salnty", x = "T_degC",
                            controls = c("O2", "O2Sat"), data = d,
                            progress_bar = FALSE, parallel = FALSE))
  for(i in seq_len(nrow(s))){
    terms_i <- s$terms[[i]]
    expect_equal(s[["O2"]][i], as.integer("O2" %in% terms_i))
    expect_equal(s[["O2Sat"]][i], as.integer("O2Sat" %in% terms_i))
  }
  # Every indicator cell is exactly 0 or 1 (no NA from coercion).
  ind <- s[, c("O2", "O2Sat")]
  expect_true(all(unlist(ind) %in% c(0L, 1L)))
})

test_that("sca() reports n_obs, which can vary across specifications by default", {
  s <- suppressMessages(sca(y = "Salnty", x = "T_degC",
                            controls = c("ChlorA", "O2Sat", "NO2uM"),
                            data = bottles, progress_bar = FALSE))
  expect_true("n_obs" %in% names(s))
  expect_true(all(s$n_obs > 0))
  # bottles has differential missingness, so specs are fit on different samples.
  expect_gt(length(unique(s$n_obs)), 1)
  # n_obs is metadata, not a control indicator.
  expect_false("n_obs" %in% sca_control_cols(s))
})

test_that("sca(common_sample = TRUE) fits every specification on one sample", {
  s <- suppressMessages(sca(y = "Salnty", x = "T_degC",
                            controls = c("ChlorA", "O2Sat", "NO2uM"),
                            data = bottles, common_sample = TRUE,
                            progress_bar = FALSE))
  expect_equal(length(unique(s$n_obs)), 1)
  cc <- sum(stats::complete.cases(
    bottles[, c("Salnty", "T_degC", "ChlorA", "O2Sat", "NO2uM")]))
  expect_equal(unique(s$n_obs), cc)
})

test_that("plot_samplesizes() returns a ggplot and needs n_obs", {
  s <- suppressMessages(sca(y = "Salnty", x = "T_degC",
                            controls = c("ChlorA", "O2Sat"), data = bottles,
                            progress_bar = FALSE))
  expect_s3_class(plot_samplesizes(s), "ggplot")
  s$n_obs <- NULL
  expect_error(plot_samplesizes(s), "n_obs")
})

test_that("sca(return_formulae = TRUE) returns formulae instead of estimates", {
  f <- sca(y = "Salnty", x = "T_degC", controls = c("O2Sat", "STheta"),
           data = bottles, return_formulae = TRUE)
  expect_length(f, 3)
})

test_that("sca() supports fixed effects via feols", {
  fe <- suppressMessages(sca(y = "Salnty", x = "T_degC", controls = "STheta",
                             data = bottles, fixed_effects = "Sta_ID",
                             progress_bar = FALSE))
  expect_s3_class(fe, "data.frame")
  expect_equal(nrow(fe), 1)
  expect_true(all(c("RMSE", "adjR") %in% names(fe)))
})

test_that("sca() supports glm families and reports AIC/deviance", {
  d <- bottles
  d$bin <- as.integer(d$Salnty > median(d$Salnty, na.rm = TRUE))
  g <- suppressMessages(sca(y = "bin", x = "T_degC", controls = "STheta",
                            data = d, family = "binomial", link = "logit",
                            progress_bar = FALSE))
  expect_s3_class(g, "data.frame")
  expect_true(all(c("AIC", "deviance") %in% names(g)))
})

test_that("sca() errors informatively when the focal variable is not a single coefficient", {
  # A factor focal variable expands to multiple model terms (e.g. "regionwarm"),
  # so there is no single focal coefficient to extract; this previously crashed
  # with a cryptic "subscript out of bounds".
  d <- bottles
  d$region <- factor(ifelse(d$T_degC > stats::median(d$T_degC, na.rm = TRUE),
                            "warm", "cold"))
  expect_error(
    suppressMessages(sca("Salnty", "region", "STheta", d, progress_bar = FALSE)),
    "single (model|focal) coefficient"
  )
})

test_that("sca() warns for fixed effects + non-linear family and then ignores them", {
  d <- bottles
  d$bin <- as.integer(d$Salnty > median(d$Salnty, na.rm = TRUE))
  expect_warning(
    res <- suppressMessages(sca(y = "bin", x = "T_degC", controls = "STheta",
                                data = d, family = "binomial", link = "logit",
                                fixed_effects = "Sta_ID", progress_bar = FALSE)),
    "Fixed effects unsupported"
  )
  # The warning promises to ignore the fixed effects, so estimation should fall
  # back to the glm path and return a valid result rather than erroring.
  expect_s3_class(res, "data.frame")
  expect_true("AIC" %in% names(res))
})

test_that("sca() keeps the focal variable when a control name contains it", {
  # Regression: paste_factory() used str_detect (a regex substring test), so a
  # control whose name contains the focal name -- or a focal name with regex
  # metacharacters -- dropped the focal term and crashed sca() with "subscript
  # out of bounds". Now matched by exact equality.
  set.seed(1); n <- 120
  d <- data.frame(y = stats::rnorm(n), Temp = stats::rnorm(n),
                  TempX = stats::rnorm(n), Z = stats::rnorm(n))
  s <- suppressMessages(sca("y", "Temp", c("Z", "TempX"), d,
                            progress_bar = FALSE))
  expect_s3_class(s, "sca")
  expect_true(all(vapply(s$terms, function(t) "Temp" %in% t, logical(1))))
  expect_false(anyNA(s$coef))
  # Helper-level checks of the exact-match behaviour.
  expect_equal(paste_factory(c("Z", "TempX"), "Temp"), "Temp + Z + TempX")
  expect_equal(paste_factory(c("Temp", "Z"), "Temp"), "Temp + Z")
  # A focal variable that is an interaction partner is dropped as a duplicate.
  expect_equal(duplicate_remover(c("x*b", "b"), "x"), "x*b")
  # An interaction not involving the focal is left untouched.
  expect_setequal(duplicate_remover(c("a*b", "a", "b"), "x"), c("a*b", "a", "b"))
})
