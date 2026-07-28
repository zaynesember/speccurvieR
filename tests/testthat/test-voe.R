# Tests for the vibration-of-effects summary and volcano plot

make_curve <- function(){
  sca(y = "Salnty", x = "T_degC",
      controls = c("O2Sat", "STheta", "ChlorA"),
      data = bottles, progress_bar = FALSE)
}

test_that("sca_voe() reproduces the Patel et al. definitions", {
  s <- make_curve()
  v <- sca_voe(s)

  expect_s3_class(v, "data.frame")
  expect_equal(nrow(v), 1)
  expect_equal(v$n_specs, nrow(s))

  # Percentile spread of the estimate, difference-type for OLS
  q <- unname(stats::quantile(s$coef, c(0.01, 0.99)))
  expect_equal(v$estimate_lo, q[1])
  expect_equal(v$estimate_hi, q[2])
  expect_equal(v$relative_effect, q[2] - q[1])
  expect_identical(v$relative_effect_type, "difference")

  # RP = difference of the -log10(p) percentiles
  qp <- unname(stats::quantile(-log10(s$p), c(0.01, 0.99)))
  expect_equal(v$RP, qp[2] - qp[1])

  # Janus flag matches the percentile signs
  expect_identical(v$janus, q[1] < 0 && q[2] > 0)

  expect_equal(v$median_estimate, stats::median(s$coef))
  expect_equal(v$median_p, stats::median(s$p))
})

test_that("sca_voe() reports a hazard-ratio ratio for cox curves", {
  s <- sca(y = c("time", "status"), x = "age",
           controls = c("sex", "ph.ecog"),
           data = survival::lung, family = "cox", progress_bar = FALSE)
  v <- sca_voe(s)
  q <- unname(stats::quantile(s$coef, c(0.01, 0.99)))
  expect_equal(v$relative_effect, exp(q[2]) / exp(q[1]))
  expect_identical(v$relative_effect_type, "ratio (RHR)")
})

test_that("sca_voe() detects a Janus effect when the spread straddles zero", {
  s <- make_curve()
  # Force a straddling curve by centering the estimates around zero
  s$coef <- s$coef - stats::median(s$coef)
  s$coef[1] <- -abs(s$coef[nrow(s)]) - 1   # guarantee both signs at extremes
  v <- sca_voe(s)
  expect_true(v$janus)
})

test_that("sca_voe() validates its inputs", {
  s <- make_curve()
  expect_error(sca_voe(data.frame(a = 1)), "does not look like sca")
  expect_error(sca_voe(s, probs = c(0.5)), "two increasing")
  expect_error(sca_voe(s, probs = c(0.9, 0.1)), "two increasing")
  expect_error(sca_voe(s, probs = c(0, 0.99)), "two increasing")
})

test_that("plot_voe() returns a ggplot for linear and cox curves", {
  s <- make_curve()
  p <- plot_voe(s)
  expect_s3_class(p, "ggplot")

  s_cox <- sca(y = c("time", "status"), x = "age",
               controls = c("sex", "ph.ecog"),
               data = survival::lung, family = "cox", progress_bar = FALSE)
  p_cox <- plot_voe(s_cox)
  expect_s3_class(p_cox, "ggplot")

  expect_error(plot_voe(data.frame(a = 1)), "does not look like sca")
})
