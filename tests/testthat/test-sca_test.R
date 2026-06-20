# Tests for sca_test() (joint-inference test) and plot_sca_test()

# A small, fast configuration reused across tests.
ctrls <- c("ChlorA", "O2Sat")

test_that("sca_test() returns a well-formed sca_test object", {
  r <- sca_test("Salnty", "T_degC", ctrls, bottles,
                n_permutations = 50, seed = 1, progress_bar = FALSE)
  expect_s3_class(r, "sca_test")
  expect_named(r, c("observed", "null_distribution", "p_values", "params"))
  stats <- c("median", "share_significant", "stouffer")
  expect_named(r$observed, stats)
  expect_named(r$p_values, stats)
  expect_s3_class(r$null_distribution, "data.frame")
  expect_equal(nrow(r$null_distribution), r$params$n_used)
  expect_true(all(c(stats, "n_valid") %in% names(r$null_distribution)))
  # Every p-value is a valid permutation p-value in [1/(n+1), 1].
  expect_true(all(r$p_values >= 1 / (r$params$n_used + 1) - 1e-9))
  expect_true(all(r$p_values <= 1))
})

test_that("sca_test() is reproducible for a fixed seed (serial)", {
  a <- sca_test("Salnty", "T_degC", ctrls, bottles,
                n_permutations = 50, seed = 99, progress_bar = FALSE)
  b <- sca_test("Salnty", "T_degC", ctrls, bottles,
                n_permutations = 50, seed = 99, progress_bar = FALSE)
  expect_equal(a$null_distribution, b$null_distribution)
  expect_equal(a$p_values, b$p_values)
})

test_that("sca_test() serial and parallel paths agree for a fixed seed", {
  skip_on_cran()
  s <- sca_test("Salnty", "T_degC", ctrls, bottles,
                n_permutations = 30, seed = 5, progress_bar = FALSE)
  p <- sca_test("Salnty", "T_degC", ctrls, bottles,
                n_permutations = 30, seed = 5, parallel = TRUE, workers = 2,
                progress_bar = FALSE)
  # The permutations are generated in the master, so both paths are identical.
  expect_equal(s$p_values, p$p_values)
  expect_equal(s$null_distribution, p$null_distribution)
})

test_that("sca_test() p-values respect the 1/(n+1) floor and are never zero", {
  r <- suppressWarnings(sca_test("Salnty", "T_degC", "ChlorA", bottles,
                                 n_permutations = 10, seed = 2,
                                 progress_bar = FALSE))
  floor <- 1 / (r$params$n_used + 1)
  expect_true(all(r$p_values >= floor - 1e-9))
  expect_true(all(r$p_values > 0))
})

test_that("sca_test() warns when n_permutations is too small to resolve alpha", {
  expect_warning(
    sca_test("Salnty", "T_degC", "ChlorA", bottles,
             n_permutations = 5, seed = 1, progress_bar = FALSE),
    "small"
  )
})

test_that("sca_test() errors immediately on a non-single-coefficient focal x", {
  d <- bottles
  d$region <- factor(ifelse(d$T_degC > stats::median(d$T_degC, na.rm = TRUE),
                            "warm", "cold"))
  expect_error(
    sca_test("Salnty", "region", "ChlorA", d, n_permutations = 20,
             progress_bar = FALSE),
    "single (model|focal) coefficient"
  )
})

test_that("sca_test() validates test_stats and direction", {
  expect_error(
    sca_test("Salnty", "T_degC", ctrls, bottles, test_stats = "bogus",
             progress_bar = FALSE),
    "Invalid `test_stats`"
  )
  expect_error(
    sca_test("Salnty", "T_degC", ctrls, bottles, direction = "up",
             progress_bar = FALSE),
    "should be one of"
  )
})

test_that("sca_test() supports all three direction modes", {
  for(d in c("two.sided", "positive", "negative")){
    r <- sca_test("Salnty", "T_degC", ctrls, bottles, n_permutations = 50,
                  seed = 4, direction = d, progress_bar = FALSE)
    expect_equal(r$params$direction, d)
    expect_true(all(r$p_values >= 1 / (r$params$n_used + 1) - 1e-9 &
                      r$p_values <= 1))
  }
})

test_that("sca_test() formula interface matches the argument interface", {
  a <- sca_test("Salnty", "T_degC", ctrls, bottles, n_permutations = 40,
                seed = 8, progress_bar = FALSE)
  b <- sca_test(Salnty ~ T_degC + ChlorA + O2Sat, data = bottles,
                n_permutations = 40, seed = 8, progress_bar = FALSE)
  expect_equal(a$observed, b$observed)
  expect_equal(a$p_values, b$p_values)
})

test_that("sca_test() blocks the permutation within fixed effects", {
  r <- sca_test(Salnty ~ T_degC + ChlorA | Sta_ID, data = bottles,
                n_permutations = 20, seed = 3, progress_bar = FALSE)
  expect_true(r$params$blocked)
  expect_equal(r$params$fixed_effects, "Sta_ID")
})

test_that("sca_test_perm_index keeps x within its fixed-effect block", {
  grp <- rep(1:5, each = 6)
  n <- length(grp)
  block_groups <- unname(split(seq_len(n), grp))
  set.seed(1)
  idx <- sca_test_perm_index(n, block_groups)
  # The permutation only reorders rows within a group, so group membership is
  # preserved: grp[idx] == grp.
  expect_equal(grp[idx], grp)
  # And it is a genuine permutation.
  expect_setequal(idx, seq_len(n))
  # A global permutation (no blocks) need not preserve grp.
  set.seed(1)
  gidx <- sca_test_perm_index(n, NULL)
  expect_setequal(gidx, seq_len(n))
})

test_that("sca_test_statistics drops NA specifications from the denominators", {
  curve <- data.frame(coef = c(1, 2, NA, -1),
                      p = c(0.01, 0.20, NA, 0.30))
  s <- sca_test_statistics(curve, c("median", "share_significant"),
                           direction = "two.sided", alpha = 0.05)
  # Median over non-NA coefs: median(c(1, 2, -1)) = 1.
  expect_equal(unname(s[["median"]]), 1)
  # Share significant: 1 of 3 valid (non-NA p) specs has p < 0.05.
  expect_equal(unname(s[["share_significant"]]), 1 / 3)
  expect_equal(attr(s, "n_valid"), 3)
})

test_that("sca_test() accepts a precomputed observed curve via sca_data", {
  curve <- suppressMessages(sca("Salnty", "T_degC", ctrls, bottles,
                                progress_bar = FALSE))
  r1 <- sca_test("Salnty", "T_degC", ctrls, bottles, sca_data = curve,
                 n_permutations = 30, seed = 6, progress_bar = FALSE)
  r2 <- sca_test("Salnty", "T_degC", ctrls, bottles,
                 n_permutations = 30, seed = 6, progress_bar = FALSE)
  # Observed statistics come from the supplied curve; identical either way.
  expect_equal(r1$observed, r2$observed)
  # Inconsistent controls are rejected.
  expect_error(
    sca_test("Salnty", "T_degC", "ChlorA", bottles, sca_data = curve,
             n_permutations = 20, progress_bar = FALSE),
    "inconsistent"
  )
})

test_that("sca_test() accounts for every permutation (used + failed)", {
  r <- sca_test("Salnty", "T_degC", ctrls, bottles, n_permutations = 40,
                seed = 7, progress_bar = FALSE)
  expect_equal(r$params$n_used + r$params$n_failed, r$params$n_permutations)
  expect_gt(r$params$n_used, 0)
})

test_that("sca_test() runs with weights and with a glm family", {
  d <- bottles
  d$w <- 1
  rw <- sca_test("Salnty", "T_degC", "ChlorA", d, weights = "w",
                 n_permutations = 20, seed = 1, progress_bar = FALSE)
  expect_s3_class(rw, "sca_test")

  d$bin <- as.integer(d$Salnty > stats::median(d$Salnty, na.rm = TRUE))
  rg <- suppressWarnings(sca_test("bin", "T_degC", "ChlorA", d,
                                  family = "binomial", n_permutations = 20,
                                  seed = 1, progress_bar = FALSE))
  expect_s3_class(rg, "sca_test")
})

test_that("plot_sca_test() returns a single ggplot", {
  r <- sca_test("Salnty", "T_degC", ctrls, bottles, n_permutations = 40,
                seed = 1, progress_bar = FALSE)
  g <- plot_sca_test(r)
  expect_s3_class(g, "ggplot")
  expect_false(inherits(g, "patchwork"))
  expect_s3_class(plot_sca_test(r, type = "density"), "ggplot")
  expect_error(plot_sca_test(mtcars), "sca_test")
})

test_that("print.sca_test returns its input invisibly", {
  r <- sca_test("Salnty", "T_degC", ctrls, bottles, n_permutations = 20,
                seed = 1, progress_bar = FALSE)
  expect_output(print(r), "joint-inference test")
  expect_invisible(print(r))
})
