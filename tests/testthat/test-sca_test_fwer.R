# Tests for per-specification FWER (min-P) inference: sca_minp(), the
# auto-attach path in sca_test(), and the plot/reporting integration.

ctrls <- c("ChlorA", "O2Sat", "NO2uM")

# A seeded collinear true-null DGP (same shape as test-sca_test_nulls.R): x is
# correlated with the confounder z1 but has no effect on y given the controls,
# so the under-controlled "z2" specification has an OFF-CENTRE null.
null_dgp <- function(n = 300, seed = 7, beta_x = 0) {
  set.seed(seed)
  z1 <- stats::rnorm(n)
  z2 <- stats::rnorm(n)
  x <- 0.8 * z1 + sqrt(0.36) * stats::rnorm(n)
  y <- beta_x * x + z1 - 0.5 * z2 + stats::rnorm(n)
  data.frame(y = y, x = x, z1 = z1, z2 = z2)
}

# Convenience: a Freedman-Lane sca_test with retained curves + auto FWER.
fl_fit <- function(d, controls = c("z1", "z2"), n = 199, seed = 1, beta = 0,
                   ...) {
  suppressMessages(sca_test("y", "x", controls, d, null_type = "freedman_lane",
                            n_permutations = n, keep_curves = TRUE, seed = seed,
                            progress_bar = FALSE, ...))
}

test_that("sca_test() auto-attaches single-step FWER under a valid null", {
  r <- fl_fit(null_dgp(300, 7, 0.4))
  fwer <- r$null_curves$fwer
  expect_false(is.null(fwer))
  expect_equal(fwer$summary$method, "single_step")
  expect_equal(fwer$summary$fwer_type, "weak")
  expect_equal(r$params$fwer_method, "single_step")
  # One row per original spec; the per-spec columns are present.
  expect_equal(nrow(fwer$specs), r$params$n_specs)
  expect_true(all(c("spec", "observed", "p_raw", "p_adj", "significant_adj")
                  %in% names(fwer$specs)))
  # The auto-attached result equals an explicit single-step sca_minp() call.
  expect_equal(fwer, sca_minp(r, method = "single_step")$null_curves$fwer)
})

test_that("FWER adjusted p-values respect the core invariants", {
  r <- fl_fit(null_dgp(300, 7, 0.4))
  sp <- r$null_curves$fwer$specs
  floor <- 1 / (r$params$n_used + 1)
  # The assertion that actually guards the subtle includes-self / grid bug:
  # correction can only make a p-value larger, never smaller. (This fails under
  # an exclude-self or no-+1 grid, or a min(1,2*min) two-sided p.)
  expect_true(all(sp$p_adj >= sp$p_raw - 1e-12, na.rm = TRUE))
  expect_true(all(sp$p_adj >= floor - 1e-12, na.rm = TRUE))
  expect_true(all(sp$p_adj <= 1 + 1e-12, na.rm = TRUE))
  expect_true(all(sp$p_raw >= floor - 1e-12, na.rm = TRUE))
})

test_that("the per-spec null p-values are marginally uniform, even off-centre", {
  skip_on_cran()
  # Large B so the empirical rates are sharp. The 'z2' spec (controls omit the
  # confounder z1) has a null centred far from zero; its own-null p must still
  # be uniform -- the property a zero-referenced marginal p would violate.
  r <- fl_fit(null_dgp(400, 7, 0), n = 2000, seed = 11)
  nc <- r$null_curves
  for(key in c("z2", "z1 + z2")){
    col <- nc$null_coef[which(nc$spec$spec == key)[1], ]
    col <- col[!is.na(col)]
    ps <- speccurvieR:::sca_fwer_null_pstar(col, "two.sided")
    ps <- ps[!is.na(ps)]
    expect_gt(mean(ps <= 0.05), 0.035)
    expect_lt(mean(ps <= 0.05), 0.065)
    expect_gt(mean(ps <= 0.10), 0.080)
    expect_lt(mean(ps <= 0.10), 0.120)
    expect_lt(abs(mean(ps) - 0.5), 0.03)
  }
})

test_that("the off-centre under-controlled spec is NOT falsely flagged", {
  # Under a true null the confounded z2 spec sits at the centre of its own
  # (off-centre) null band, so it must not be FWER-significant; a marginal-p
  # procedure would (wrongly) flag it.
  r <- fl_fit(null_dgp(400, 7, 0), n = 499, seed = 3)
  sp <- r$null_curves$fwer$specs
  expect_false(isTRUE(sp$significant_adj[sp$spec == "z2"][1]))
})

test_that("single-step controls weak FWER under the complete null", {
  skip_on_cran()
  n_sets <- 200L
  rej05 <- rej10 <- logical(n_sets)
  for(i in seq_len(n_sets)){
    d <- null_dgp(n = 150, seed = 1000 + i, beta_x = 0)   # true null
    pa <- fl_fit(d, n = 199, seed = 1)$null_curves$fwer$specs$p_adj
    rej05[i] <- any(pa < 0.05, na.rm = TRUE)
    rej10[i] <- any(pa < 0.10, na.rm = TRUE)
  }
  se05 <- sqrt(0.05 * 0.95 / n_sets)
  se10 <- sqrt(0.10 * 0.90 / n_sets)
  # At most alpha (allowing Monte-Carlo slack); and not grossly conservative.
  expect_lt(mean(rej05), 0.05 + 3 * se05)
  expect_lt(mean(rej10), 0.10 + 3 * se10)
  expect_gt(mean(rej10), 0.01)
})

test_that("the procedure has power against a real effect", {
  skip_on_cran()
  n_sets <- 80L
  hit <- logical(n_sets)
  for(i in seq_len(n_sets)){
    d <- null_dgp(n = 300, seed = 2000 + i, beta_x = 0.4)
    sp <- fl_fit(d, n = 199, seed = 1)$null_curves$fwer$specs
    hit[i] <- isTRUE(sp$significant_adj[sp$spec == "z1 + z2"][1])
  }
  expect_gt(mean(hit), 0.5)
})

test_that("step_down is monotone, floored, and at least as powerful", {
  r <- fl_fit(null_dgp(300, 7, 0.4))
  ss <- r$null_curves$fwer$specs
  sd <- sca_minp(r, method = "step_down")$null_curves$fwer
  expect_equal(sd$summary$method, "step_down")
  expect_equal(sd$summary$fwer_type, "strong")
  d <- sd$specs
  floor <- 1 / (r$params$n_used + 1)
  o <- order(d$p_raw)
  expect_true(all(diff(d$p_adj[o]) >= -1e-12))            # monotone in p_raw
  expect_true(all(d$p_adj >= d$p_raw - 1e-12, na.rm = TRUE))
  expect_true(all(d$p_adj >= floor - 1e-12, na.rm = TRUE))
  # Free step-down is uniformly at least as powerful as single-step.
  expect_true(all(d$p_adj <= ss$p_adj + 1e-12, na.rm = TRUE))
})

test_that("sca_minp() recomputes live with a new alpha", {
  r <- fl_fit(null_dgp(300, 7, 0.4))
  loose <- sca_minp(r, alpha = 0.2)$null_curves$fwer
  tight <- sca_minp(r, alpha = 0.001)$null_curves$fwer
  # p_adj does not depend on alpha; only the significance flag and tally do.
  expect_equal(loose$specs$p_adj, tight$specs$p_adj)
  expect_equal(loose$specs$significant_adj, loose$specs$p_adj < 0.2)
  expect_gte(loose$summary$n_significant, tight$summary$n_significant)
  expect_equal(loose$summary$alpha, 0.2)
})

test_that("sca_minp() gates on the null type and keep_curves", {
  d <- null_dgp(300, 7, 0.4)
  # shuffle_x: per-spec FWER is invalid -> hard error naming the valid nulls.
  rs <- sca_test("y", "x", c("z1", "z2"), d, n_permutations = 40, seed = 1,
                 keep_curves = TRUE, progress_bar = FALSE)
  expect_null(rs$null_curves$fwer)               # auto path is a silent no-op
  expect_error(sca_minp(rs), "freedman_lane")
  # keep_curves = FALSE: nothing to compute from.
  r0 <- suppressMessages(sca_test("y", "x", c("z1", "z2"), d,
                                  null_type = "freedman_lane",
                                  n_permutations = 40, seed = 1,
                                  progress_bar = FALSE))
  expect_null(r0$null_curves)
  expect_error(sca_minp(r0), "keep_curves")
  expect_error(sca_minp(mtcars), "sca_test")
})

test_that("shuffle_x keep_curves object is unchanged (back-compat)", {
  d <- null_dgp(300, 7, 0)
  a <- sca_test("y", "x", c("z1", "z2"), d, n_permutations = 40, seed = 1,
                keep_curves = TRUE, progress_bar = FALSE)
  # No fwer slot, and only the original two elements.
  expect_named(a$null_curves, c("spec", "null_coef"))
  expect_true(is.na(a$params$fwer_method))
})

test_that("duplicate-key specifications collapse to one tested model", {
  # An interaction control plus its components fit the identical model, so
  # several control subsets share a terms key. Inference must count each
  # distinct model once and assign rows sharing a key identical values.
  r <- suppressMessages(sca_test("Salnty", "T_degC",
                                 c("ChlorA*O2Sat", "ChlorA", "O2Sat"), bottles,
                                 null_type = "freedman_lane", n_permutations = 60,
                                 keep_curves = TRUE, seed = 1,
                                 progress_bar = FALSE))
  fwer <- r$null_curves$fwer
  sp <- fwer$specs
  expect_lt(fwer$summary$n_specs_tested, nrow(sp))      # fewer distinct models
  expect_equal(fwer$summary$n_specs_tested, length(unique(sp$spec)))
  dup <- unique(sp$spec[duplicated(sp$spec)])[1]
  rows <- sp[sp$spec == dup, ]
  expect_equal(rows$p_raw, rep(rows$p_raw[1], nrow(rows)))
  expect_equal(rows$p_adj, rep(rows$p_adj[1], nrow(rows)))
})

test_that("the engine drops all-NA specs and excludes them from the tally", {
  # Drive the engine directly with a crafted matrix: key 'b' is non-estimable
  # in every permutation, and key 'a' appears twice (duplicate model).
  null_coef <- rbind(a1 = c(1, 2, 3, 4, 5),
                     b  = c(NA, NA, NA, NA, NA),
                     a2 = c(1, 2, 3, 4, 5))
  spec <- data.frame(spec = c("a", "b", "a"), observed = c(6, 0.5, 6),
                     stringsAsFactors = FALSE)
  out <- speccurvieR:::sca_fwer_compute(null_coef, spec, direction = "two.sided",
                                        method = "single_step", alpha = 0.05,
                                        null_type = "freedman_lane")
  # Only the distinct, estimable key 'a' is tested.
  expect_equal(out$summary$n_specs_tested, 1L)
  expect_true(is.na(out$specs$p_raw[2]))               # all-NA 'b'
  expect_true(is.na(out$specs$p_adj[2]))
  expect_false(isTRUE(out$specs$significant_adj[2]))
  # Rows sharing key 'a' are identical and finite.
  expect_equal(out$specs$p_adj[1], out$specs$p_adj[3])
  expect_false(is.na(out$specs$p_adj[1]))
})

test_that("ragged per-spec usable-draw counts keep p_adj >= p_raw", {
  # A specification estimable in only SOME permutations has a shorter, coarser
  # own-null grid than its peers. Without the engine's clamp, the pooled
  # minimum-p null could hand it an anti-conservative p_adj < p_raw. Spec 'a' is
  # usable in 7 of 20 permutations, 'b' in all 20, and 'a' is observed beyond
  # all its usable draws -- the configuration that exposed the bug.
  set.seed(42)
  B <- 20L
  a_col <- c(stats::rnorm(7), rep(NA_real_, B - 7))
  b_col <- stats::rnorm(B)
  null_coef <- rbind(a = a_col, b = b_col)
  spec <- data.frame(spec = c("a", "b"),
                     observed = c(max(a_col, na.rm = TRUE) + 5, 0),
                     stringsAsFactors = FALSE)
  for(dir in c("two.sided", "positive", "negative")){
    for(meth in c("single_step", "step_down")){
      out <- speccurvieR:::sca_fwer_compute(null_coef, spec, direction = dir,
                                            method = meth, alpha = 0.05,
                                            null_type = "freedman_lane")
      sp <- out$specs
      lbl <- paste(dir, meth)
      expect_true(all(sp$p_adj >= sp$p_raw - 1e-12, na.rm = TRUE), info = lbl)
      expect_true(all(sp$p_adj <= 1 + 1e-12, na.rm = TRUE), info = lbl)
    }
  }
})

test_that("step_down does not crash with a single permutation (B == 1)", {
  # B == 1 made the old t(vapply()) construction collapse orientation and
  # overrun the step-down indexing.
  null_coef <- matrix(c(1, 2), nrow = 2, ncol = 1,
                      dimnames = list(c("a", "b"), NULL))
  spec <- data.frame(spec = c("a", "b"), observed = c(3, 0.5),
                     stringsAsFactors = FALSE)
  out <- speccurvieR:::sca_fwer_compute(null_coef, spec, direction = "two.sided",
                                        method = "step_down", alpha = 0.05,
                                        null_type = "freedman_lane")
  expect_equal(nrow(out$specs), 2L)
  expect_true(all(out$specs$p_adj >= out$specs$p_raw - 1e-12, na.rm = TRUE))
  expect_true(all(out$specs$p_adj <= 1 + 1e-12, na.rm = TRUE))
})

test_that("auto FWER is identical across serial and parallel", {
  skip_on_cran()
  d <- null_dgp(300, 7, 0.4)
  s <- fl_fit(d, n = 60, seed = 5)
  p <- suppressMessages(sca_test("y", "x", c("z1", "z2"), d,
                                 null_type = "freedman_lane", n_permutations = 60,
                                 keep_curves = TRUE, seed = 5, parallel = TRUE,
                                 workers = 2, progress_bar = FALSE))
  expect_equal(s$null_curves$fwer, p$null_curves$fwer)
})

test_that("residual_bootstrap also auto-attaches FWER", {
  r <- suppressMessages(sca_test("y", "x", c("z1", "z2"), null_dgp(300, 7, 0.4),
                                 null_type = "residual_bootstrap",
                                 n_permutations = 99, keep_curves = TRUE,
                                 seed = 1, progress_bar = FALSE))
  expect_false(is.null(r$null_curves$fwer))
  expect_equal(r$null_curves$fwer$summary$null_type, "residual_bootstrap")
})

test_that("reporting surfaces expose the FWER results", {
  r <- fl_fit(null_dgp(300, 7, 0.4))
  # print
  expect_output(print(r), "After correcting for searching")
  # as.data.frame(what = "specs")
  specs <- as.data.frame(r, what = "specs")
  expect_true(all(c("p_raw", "p_adj", "significant_adj") %in% names(specs)))
  # glance carries the FWER columns
  g <- glance(r)
  expect_true(all(c("n_significant_adj", "min_p_adj", "fwer_method") %in%
                    names(g)))
  expect_equal(g$fwer_method, "single_step")
  # sca_report / sca_table mention the correction
  expect_match(sca_report(r), "family-wise error")
  tab <- sca_table(r)
  expect_true(any(grepl("Significant after correction", tab$label)))
  # A result without FWER omits the columns gracefully (stable schema: NA).
  rs <- sca_test("y", "x", c("z1", "z2"), null_dgp(300, 7, 0), seed = 1,
                 n_permutations = 40, progress_bar = FALSE)
  expect_true(is.na(glance(rs)$min_p_adj))
  expect_error(as.data.frame(rs, what = "specs"), "no per-specification")
})

test_that("plot_sca_test_specs() draws the 3-tier highlight under FWER", {
  r <- fl_fit(null_dgp(300, 7, 0.4))
  g <- plot_sca_test_specs(r)
  expect_s3_class(g, "ggplot")
  expect_false(inherits(g, "patchwork"))
  # The 2-tier fallback still works when no FWER is attached.
  rs <- sca_test("y", "x", c("z1", "z2"), null_dgp(300, 7, 0), seed = 1,
                 n_permutations = 40, keep_curves = TRUE, progress_bar = FALSE)
  expect_s3_class(plot_sca_test_specs(rs), "ggplot")
})
