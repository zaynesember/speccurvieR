# Tests for sca_variance() and plot_variance()

ctrls <- c("ChlorA", "O2Sat", "NO2uM")

curve <- function(controls = ctrls, data = bottles, ...) {
  suppressMessages(sca(y = "Salnty", x = "T_degC", controls = controls,
                       data = data, progress_bar = FALSE, parallel = FALSE, ...))
}

test_that("sca_variance() returns a tidy choice/variance/percent table summing to 100", {
  s <- curve()
  v <- sca_variance(s)
  expect_s3_class(v, "data.frame")
  expect_named(v, c("choice", "variance", "percent"))
  expect_true("Residual" %in% v$choice)
  expect_true(all(ctrls %in% v$choice))
  expect_equal(sum(v$percent), 100, tolerance = 1e-8)
  # Sorted by descending share (Residual is small here, lands last).
  non_resid <- v$percent[v$choice != "Residual"]
  expect_false(is.unsorted(rev(non_resid)))
})

test_that("LMG shares sum to the full-model R-squared", {
  s <- curve()
  v <- sca_variance(s)
  r2 <- summary(stats::lm(coef ~ ChlorA + O2Sat + NO2uM, data = s))$r.squared
  terms_pct <- sum(v$percent[v$choice != "Residual"])
  expect_equal(terms_pct, 100 * r2, tolerance = 1e-6)
  expect_equal(v$percent[v$choice == "Residual"], 100 * (1 - r2),
               tolerance = 1e-6)
})

test_that("residual = FALSE drops the Residual row and rescales to 100", {
  s <- curve()
  v <- sca_variance(s, residual = FALSE)
  expect_false("Residual" %in% v$choice)
  expect_equal(sum(v$percent), 100, tolerance = 1e-8)
})

test_that("'shapley' is an alias for 'lmg'", {
  s <- curve()
  expect_equal(sca_variance(s, method = "shapley"),
               sca_variance(s, method = "lmg"))
})

test_that("method = 'type2' sums to 100 and differs from lmg", {
  s <- curve()
  v2 <- sca_variance(s, method = "type2")
  expect_equal(sum(v2$percent), 100, tolerance = 1e-8)
  expect_false(isTRUE(all.equal(v2$percent, sca_variance(s)$percent)))
})

test_that("the variance column partitions var(coef) for both methods", {
  s <- curve()
  for(m in c("lmg", "type2")){
    v <- sca_variance(s, method = m)
    expect_equal(sum(v$variance), stats::var(s$coef), tolerance = 1e-8)
    # variance is exactly percent/100 of the total variance.
    expect_equal(v$variance, v$percent / 100 * stats::var(s$coef),
                 tolerance = 1e-8)
  }
})

test_that("sca_variance() ignores stray non-binary columns", {
  s <- curve()
  s$junk <- stats::rnorm(nrow(s))
  expect_warning(v <- sca_variance(s), "non-binary")
  expect_false("junk" %in% v$choice)
})

test_that("sca_variance() does not reject a small but real spread at large location", {
  # The near-constant guard must key off the centered spread, not the
  # coefficient's magnitude.
  s <- curve()
  s$coef <- s$coef - mean(s$coef) + 1e8
  expect_silent(v <- sca_variance(s))
  expect_equal(sum(v$percent), 100, tolerance = 1e-8)
})

test_that("sca_variance() handles an interaction control without rank deficiency", {
  # The interaction is a single choice (its term column), not three collinear
  # term columns -- so this stays well-posed even though distinct terms exceed
  # the spec count.
  s <- curve(controls = c("ChlorA*O2Sat", "NO2uM"))
  v <- sca_variance(s)
  expect_true("ChlorA:O2Sat" %in% v$choice)
  expect_equal(sum(v$percent), 100, tolerance = 1e-8)
  expect_true(all(is.finite(v$percent)))
})

test_that("sca_variance() errors with a single control (too few specifications)", {
  s <- curve(controls = "ChlorA")
  expect_error(sca_variance(s), "two specifications")
})

test_that("sca_variance() errors when the estimate is (near-)constant", {
  s <- curve()
  s$coef <- 0.5
  expect_error(sca_variance(s), "near-.constant|no variance")
})

test_that("sca_variance() drops specifications with a missing estimate", {
  s <- curve()
  s_na <- s
  s_na$coef[1] <- NA
  expect_warning(v <- sca_variance(s_na), "missing")
  # Result matches physically removing that row.
  v_ref <- sca_variance(s[-1, , drop = FALSE])
  expect_equal(v$percent, v_ref$percent, tolerance = 1e-8)
})

test_that("sca_variance() works on a glm curve", {
  d <- bottles
  d$bin <- as.integer(d$Salnty > stats::median(d$Salnty, na.rm = TRUE))
  s <- suppressWarnings(suppressMessages(
    sca(y = "bin", x = "T_degC", controls = c("ChlorA", "O2Sat"), data = d,
        family = "binomial", progress_bar = FALSE, parallel = FALSE)))
  v <- sca_variance(s)
  expect_true(all(is.finite(v$percent)))
  expect_equal(sum(v$percent), 100, tolerance = 1e-8)
})

test_that("sca_variance() validates its inputs", {
  s <- curve()
  expect_error(sca_variance(s, estimate = "nope"), "not found")
  s$terms <- NULL
  s$bad <- letters[seq_len(nrow(s))]
  expect_error(sca_variance(s, estimate = "bad"), "numeric")
  expect_error(sca_variance(42), "data frame")
})

test_that("plot_variance() returns a single ggplot", {
  s <- curve()
  p <- plot_variance(s)
  expect_s3_class(p, "ggplot")
  expect_false(inherits(p, "patchwork"))
  expect_s3_class(plot_variance(s, method = "type2", residual = FALSE), "ggplot")
})
