# Tests for sca_test() null mechanisms (shuffle_x / freedman_lane /
# residual_bootstrap)

ct <- c("ChlorA", "O2Sat", "NO2uM")

# A seeded collinear true-null data set: x is correlated with confounder z1 but
# has no effect on y given the controls.
null_dgp <- function(n = 300, seed = 7, beta_x = 0) {
  set.seed(seed)
  z1 <- stats::rnorm(n)
  z2 <- stats::rnorm(n)
  x <- 0.8 * z1 + sqrt(0.36) * stats::rnorm(n)
  y <- beta_x * x + z1 - 0.5 * z2 + stats::rnorm(n)
  data.frame(y = y, x = x, z1 = z1, z2 = z2)
}

test_that("null_type defaults to shuffle_x and equals the explicit value", {
  d <- bottles
  a <- sca_test("Salnty", "T_degC", ct, d, n_permutations = 60, seed = 1,
                progress_bar = FALSE)
  b <- sca_test("Salnty", "T_degC", ct, d, n_permutations = 60, seed = 1,
                null_type = "shuffle_x", progress_bar = FALSE)
  expect_equal(a$params$null_type, "shuffle_x")
  # Back-compat lock: passing the default explicitly changes nothing.
  expect_equal(a$p_values, b$p_values)
  expect_equal(a$null_distribution, b$null_distribution)
  # Including the per-spec null curves.
  ak <- sca_test("Salnty", "T_degC", ct, d, n_permutations = 40, seed = 1,
                 keep_curves = TRUE, progress_bar = FALSE)
  bk <- sca_test("Salnty", "T_degC", ct, d, n_permutations = 40, seed = 1,
                 keep_curves = TRUE, null_type = "shuffle_x",
                 progress_bar = FALSE)
  expect_equal(ak$null_curves, bk$null_curves)
})

test_that("sca_data is ignored (with a warning) for the design-preserving nulls", {
  d <- bottles
  curve <- suppressMessages(sca("Salnty", "T_degC", ct, d, progress_bar = FALSE))
  # The observed curve must come from the forced common sample, not the supplied
  # (possibly per-spec-sample) sca_data, so it matches the null's design.
  expect_warning(
    r <- sca_test("Salnty", "T_degC", ct, d, sca_data = curve,
                  null_type = "freedman_lane", n_permutations = 40, seed = 1,
                  progress_bar = FALSE),
    "ignored")
  r2 <- suppressMessages(
    sca_test("Salnty", "T_degC", ct, d, null_type = "freedman_lane",
             n_permutations = 40, seed = 1, progress_bar = FALSE))
  expect_equal(r$observed, r2$observed)
})

test_that("sca_test() validates null_type and its scope", {
  d <- bottles
  expect_error(sca_test("Salnty", "T_degC", ct, d, null_type = "bogus",
                        progress_bar = FALSE))
  # glm families are out of scope for the design-preserving nulls.
  d$bin <- as.integer(d$Salnty > stats::median(d$Salnty, na.rm = TRUE))
  expect_error(
    sca_test("bin", "T_degC", "ChlorA", d, family = "binomial",
             null_type = "freedman_lane", n_permutations = 10,
             progress_bar = FALSE),
    "linear models only")
  # A reserved-column collision errors before any work.
  d2 <- bottles
  d2$.sca_test_response <- 1
  expect_error(
    sca_test("Salnty", "T_degC", "ChlorA", d2, null_type = "residual_bootstrap",
             n_permutations = 10, progress_bar = FALSE),
    "already contains")
})

test_that("freedman_lane and residual_bootstrap force a common sample (message)", {
  d <- bottles
  expect_message(
    r <- sca_test("Salnty", "T_degC", ct, d, n_permutations = 40, seed = 1,
                  null_type = "freedman_lane", progress_bar = FALSE),
    "requires a common sample")
  expect_true(r$params$common_sample)
  expect_equal(r$params$reduced_model, "control_superset")
})

test_that("each null_type is reproducible for a fixed seed", {
  d <- bottles
  for(nt in c("shuffle_x", "freedman_lane", "residual_bootstrap")){
    a <- suppressMessages(sca_test("Salnty", "T_degC", ct, d,
                                   n_permutations = 50, seed = 3, null_type = nt,
                                   progress_bar = FALSE))
    b <- suppressMessages(sca_test("Salnty", "T_degC", ct, d,
                                   n_permutations = 50, seed = 3, null_type = nt,
                                   progress_bar = FALSE))
    expect_equal(a$p_values, b$p_values)
    expect_equal(a$null_distribution, b$null_distribution)
  }
})

test_that("serial and parallel agree for freedman_lane and residual_bootstrap", {
  # Locks the closure-capture contract: perm_fun reaches the PSOCK workers with
  # the master-side null objects (data_cs, reduced) intact, so the design-
  # preserving nulls reproduce the serial path exactly. Covers a fixed-effects
  # run so the kept-row / block-group machinery is exercised on the workers too.
  skip_on_cran()
  cases <- list(
    list(controls = ct, fe = NULL),
    list(controls = "STheta", fe = "Sta_ID"))
  for(case in cases){
    for(nt in c("freedman_lane", "residual_bootstrap")){
      s <- suppressMessages(sca_test("Salnty", "T_degC", case$controls, bottles,
                                     fixed_effects = case$fe, n_permutations = 30,
                                     seed = 5, null_type = nt, keep_curves = TRUE,
                                     progress_bar = FALSE))
      p <- suppressMessages(sca_test("Salnty", "T_degC", case$controls, bottles,
                                     fixed_effects = case$fe, n_permutations = 30,
                                     seed = 5, null_type = nt, keep_curves = TRUE,
                                     parallel = TRUE, workers = 2,
                                     progress_bar = FALSE))
      expect_equal(s$p_values, p$p_values)
      expect_equal(s$null_distribution, p$null_distribution)
      expect_equal(s$null_curves$null_coef, p$null_curves$null_coef)
    }
  }
})

test_that("the new nulls produce a well-formed object and support keep_curves", {
  d <- bottles
  for(nt in c("freedman_lane", "residual_bootstrap")){
    r <- suppressMessages(sca_test("Salnty", "T_degC", ct, d,
                                   n_permutations = 40, seed = 1, null_type = nt,
                                   keep_curves = TRUE, progress_bar = FALSE))
    expect_s3_class(r, "sca_test")
    expect_equal(dim(r$null_curves$null_coef),
                 c(r$params$n_specs, r$params$n_used))
    expect_true(all(r$p_values >= 1 / (r$params$n_used + 1) - 1e-9 &
                      r$p_values <= 1))
  }
})

test_that("residual_bootstrap imposes the null exactly on the full spec", {
  d <- null_dgp()
  ynull <- speccurvieR:::sca_test_ynull(d, "y", "x", c("z1", "z2"), NULL, NULL)
  d$yn <- ynull
  expect_lt(abs(stats::coef(stats::lm(yn ~ x + z1 + z2, d))[["x"]]), 1e-8)
})

test_that("freedman_lane reduced fit reconstructs y on the kept rows", {
  d <- null_dgp()
  # Non-FE: every row kept.
  r <- speccurvieR:::sca_test_reduced_fit(d, "y", c("z1", "z2"), NULL, NULL, NULL)
  expect_equal(r$n_kept, nrow(d))
  expect_lt(max(abs(r$fitted + r$resid - d$y[r$kept])), 1e-8)
  # FE with a singleton level: that row is dropped, the rest reconstruct exactly.
  d$fe <- rep(seq_len(60), 5)
  d$fe[1] <- 999
  rfe <- speccurvieR:::sca_test_reduced_fit(d, "y", c("z1", "z2"), "fe", NULL,
                                            "fe")
  expect_equal(rfe$n_kept, nrow(d) - 1)
  expect_lt(max(abs(rfe$fitted + rfe$resid - d$y[rfe$kept])), 1e-8)
})

test_that("freedman_lane has the right conditional-null centering", {
  d <- null_dgp(n = 400, seed = 7)
  r <- suppressMessages(sca_test("y", "x", c("z1", "z2"), d,
                                 n_permutations = 200, seed = 11,
                                 null_type = "freedman_lane", keep_curves = TRUE,
                                 progress_bar = FALSE))
  nc <- r$null_curves
  full_mean <- mean(nc$null_coef[nc$spec$spec == "z1 + z2", ], na.rm = TRUE)
  z2_mean <- mean(nc$null_coef[nc$spec$spec == "z2", ], na.rm = TRUE)
  # No partial effect given the full controls -> centred at 0.
  expect_lt(abs(full_mean), 0.05)
  # Omitting the confounder z1 -> the confounded association is retained.
  expect_gt(z2_mean, 0.5)
})

test_that("freedman_lane is better calibrated than shuffle_x under collinearity", {
  skip_on_cran()
  n_sets <- 20L
  fl_p <- numeric(n_sets)
  sh_p <- numeric(n_sets)
  for(i in seq_len(n_sets)){
    d <- null_dgp(n = 150, seed = 100 + i)  # true null
    fl_p[i] <- suppressMessages(sca_test("y", "x", c("z1", "z2"), d,
                 n_permutations = 99, seed = 1, null_type = "freedman_lane",
                 progress_bar = FALSE))$p_values[["median"]]
    sh_p[i] <- sca_test("y", "x", c("z1", "z2"), d, n_permutations = 99,
                 seed = 1, null_type = "shuffle_x",
                 progress_bar = FALSE)$p_values[["median"]]
  }
  # FL median p-values should sit near uniform (mean ~ 0.5); shuffle_x is
  # anti-conservative under collinearity, so its p-values run smaller.
  expect_gt(mean(fl_p), 0.3)
  expect_lt(mean(fl_p), 0.7)
  expect_lt(mean(sh_p), mean(fl_p))
})
