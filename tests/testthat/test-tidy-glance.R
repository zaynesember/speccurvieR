# Tests for tidy()/glance()/as.data.frame() methods and the sca() class

curve <- function() {
  suppressMessages(sca(y = "Salnty", x = "T_degC",
                       controls = c("ChlorA", "O2Sat", "NO2uM"),
                       data = bottles, progress_bar = FALSE, parallel = FALSE))
}
a_test <- function() {
  suppressMessages(sca_test(y = "Salnty", x = "T_degC",
                            controls = c("ChlorA", "O2Sat"), data = bottles,
                            n_permutations = 50, seed = 1, progress_bar = FALSE))
}

test_that("sca() output is classed but still a data frame (back-compat)", {
  s <- curve()
  expect_s3_class(s, "sca")
  expect_s3_class(s, "data.frame")     # inherits, exact = FALSE
  expect_true(is.data.frame(s))
  # The class is a pure dispatch tag: downstream consumers still work.
  expect_s3_class(plot_curve(s, plot_vars = FALSE), "ggplot")
  expect_s3_class(plot_samplesizes(s), "ggplot")
  expect_s3_class(sca_variance(s), "data.frame")
  expect_setequal(sca_control_cols(s), c("ChlorA", "O2Sat", "NO2uM"))
  # Survives a dplyr round-trip.
  s2 <- dplyr::mutate(dplyr::filter(s, p < 1), keep = TRUE)
  expect_true(is.data.frame(s2))
})

test_that("the focal variable resolves for a single-control curve", {
  s <- suppressMessages(sca(y = "Salnty", x = "T_degC", controls = "ChlorA",
                            data = bottles, progress_bar = FALSE,
                            parallel = FALSE))
  expect_equal(nrow(s), 1)
  t <- tidy(s)
  expect_equal(unique(t$term), "T_degC")
  expect_false("T_degC" %in% t$controls)   # focal must not leak into controls
  expect_equal(glance(s)$focal, "T_degC")
})

test_that("glance.sca reports the specific glm family", {
  d <- bottles
  d$bin <- as.integer(d$Salnty > stats::median(d$Salnty, na.rm = TRUE))
  s <- suppressWarnings(suppressMessages(
    sca(y = "bin", x = "T_degC", controls = c("ChlorA", "O2Sat"), data = d,
        family = "binomial", progress_bar = FALSE, parallel = FALSE)))
  expect_equal(glance(s)$family, "binomial")
  expect_equal(unique(tidy(s)$fit_stat_name), "deviance")
})

test_that("tidy.sca returns broom-style per-specification rows", {
  s <- curve()
  t <- tidy(s)
  expect_named(t, c("term", "estimate", "std.error", "statistic", "p.value",
                    "spec_id", "controls", "n_obs", "fit_stat", "fit_stat_name"))
  expect_equal(nrow(t), nrow(s))
  expect_equal(t$estimate, s$coef[order(s$index)])
  expect_equal(unique(t$term), "T_degC")
  expect_equal(unique(t$fit_stat_name), "RMSE")
})

test_that("glance.sca summarises the curve and surfaces the N range", {
  s <- curve()
  g <- glance(s)
  expect_equal(nrow(g), 1)
  expect_equal(g$n_specs, nrow(s))
  expect_equal(g$n_controls, 3)
  expect_equal(g$median_estimate, stats::median(s$coef))
  # The default sample varies across specifications -- locks the footgun.
  expect_gt(g$n_obs_max, g$n_obs_min)
  expect_false(g$common_sample)
})

test_that("tidy/glance for sca_test match the observed statistics", {
  skip_on_cran()
  r <- a_test()
  t <- tidy(r)
  expect_equal(t$term, r$params$test_stats)
  expect_equal(t$estimate, as.numeric(r$observed[r$params$test_stats]))
  expect_equal(t$p.value, as.numeric(r$p_values[r$params$test_stats]))
  expect_equal(unique(t$null_type), "shuffle_x")

  g <- glance(r)
  expect_equal(nrow(g), 1)
  expect_equal(g$n_used, r$params$n_used)
  expect_equal(g$p_resolution, 1 / (r$params$n_used + 1))
  expect_equal(g$null_type, "shuffle_x")
})

test_that("as.data.frame.sca_test returns a plain data frame with what=", {
  skip_on_cran()
  r <- a_test()
  d <- as.data.frame(r)
  expect_identical(class(d), "data.frame")
  expect_named(d, c("statistic", "observed", "p_value"))
  expect_equal(d$observed, as.numeric(r$observed[r$params$test_stats]))
  expect_equal(as.data.frame(r, what = "null"), r$null_distribution)
  expect_equal(nrow(as.data.frame(r, what = "params")), 1)
})

test_that("the generics broom dispatches on reach our methods", {
  skip_on_cran()
  # broom::tidy/glance ARE generics::tidy/glance, so dispatching the generic
  # exercises the same path broom (and modelsummary) would use.
  r <- a_test()
  expect_identical(generics::tidy(r), tidy(r))
  expect_identical(generics::glance(r), glance(r))
})
