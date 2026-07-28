# Tests for Cox proportional-hazards support (family = "cox")

lung_df <- survival::lung

test_that("sca(family = 'cox') matches a hand-fit coxph() exactly", {
  s <- sca(y = c("time", "status"), x = "age",
           controls = c("sex", "ph.ecog"),
           data = lung_df, family = "cox", progress_bar = FALSE)

  # 2 controls -> 3 specifications, cox-specific columns present
  expect_s3_class(s, "sca")
  expect_equal(nrow(s), 3)
  expect_true(all(c("coef", "HR", "se", "statistic", "p", "AIC",
                    "concordance", "n_obs", "n_events", "sig.level",
                    "terms", "control_coefs") %in% names(s)))
  expect_identical(attr(s, "family"), "cox")
  expect_identical(attr(s, "x"), "age")

  # The full specification must reproduce summary(coxph(...)) exactly
  m <- survival::coxph(survival::Surv(time, status) ~ age + sex + ph.ecog,
                       data = lung_df)
  sm <- summary(m)
  full <- s[vapply(s$terms,
                   function(t) all(c("sex", "ph.ecog") %in% t), logical(1)), ]
  expect_equal(full$coef, unname(sm$coefficients["age", "coef"]))
  expect_equal(full$HR, unname(sm$coefficients["age", "exp(coef)"]))
  expect_equal(full$se, unname(sm$coefficients["age", "se(coef)"]))
  expect_equal(full$statistic, unname(sm$coefficients["age", "z"]))
  expect_equal(full$p, unname(sm$coefficients["age", "Pr(>|z|)"]))
  expect_equal(full$AIC, stats::AIC(m))
  expect_equal(full$concordance, unname(sm$concordance[[1]]))
  expect_equal(full$n_obs, unname(sm$n))
  expect_equal(full$n_events, unname(sm$nevent))

  # HR is exp(coef) throughout
  expect_equal(s$HR, exp(s$coef))
})

test_that("a Surv() formula response implies family = 'cox'", {
  s_vec <- sca(y = c("time", "status"), x = "age",
               controls = c("sex", "ph.ecog"),
               data = lung_df, family = "cox", progress_bar = FALSE)
  s_fml <- sca(Surv(time, status) ~ age + sex + ph.ecog,
               data = lung_df, progress_bar = FALSE)
  expect_identical(attr(s_fml, "family"), "cox")
  expect_equal(s_fml[, c("coef", "HR", "se", "p")],
               s_vec[, c("coef", "HR", "se", "p")])
})

test_that("cox + weights uses coxph's robust standard errors", {
  d <- lung_df
  set.seed(1)
  d$w <- runif(nrow(d), 0.5, 2)
  s <- sca(y = c("time", "status"), x = "age",
           controls = c("sex", "ph.ecog"),
           data = d, family = "cox", weights = "w", progress_bar = FALSE)
  m <- survival::coxph(survival::Surv(time, status) ~ age + sex + ph.ecog,
                       data = d, weights = w)
  sm <- summary(m)
  full <- s[vapply(s$terms,
                   function(t) all(c("sex", "ph.ecog") %in% t), logical(1)), ]
  expect_equal(full$se, unname(sm$coefficients["age", "robust se"]))
  expect_equal(full$p, unname(sm$coefficients["age", "Pr(>|z|)"]))
})

test_that("cox respects common_sample", {
  s <- sca(y = c("time", "status"), x = "age",
           controls = c("sex", "ph.ecog"),
           data = lung_df, family = "cox", common_sample = TRUE,
           progress_bar = FALSE)
  expect_length(unique(s$n_obs), 1)
})

test_that("cox return_formulae builds Surv() responses", {
  f <- sca(y = c("time", "status"), x = "age", controls = c("sex", "ph.ecog"),
           data = lung_df, family = "cox", return_formulae = TRUE)
  lhs <- vapply(f, function(fm) deparse(fm[[2]]), character(1))
  expect_true(all(lhs == "Surv(time, status)"))
})

test_that("cox outcome validation errors are informative", {
  # length-2 y without family = "cox"
  expect_error(
    sca(y = c("time", "status"), x = "age", controls = "sex",
        data = lung_df, progress_bar = FALSE),
    "single variable name")
  # cox with a scalar y
  expect_error(
    sca(y = "time", x = "age", controls = "sex",
        data = lung_df, family = "cox", progress_bar = FALSE),
    "time and event")
  # Surv() response with a conflicting family
  expect_error(
    sca(Surv(time, status) ~ age + sex, data = lung_df,
        family = "binomial", progress_bar = FALSE),
    "requires family")
  # malformed Surv() response
  expect_error(
    sca(Surv(time) ~ age + sex, data = lung_df, progress_bar = FALSE),
    "exactly two columns")
})

test_that("cox drops fixed effects with a warning", {
  expect_warning(
    s <- sca(y = c("time", "status"), x = "age", controls = "sex",
             fixed_effects = "inst", data = lung_df, family = "cox",
             progress_bar = FALSE),
    "Ignoring fixed effects")
  expect_equal(nrow(s), 1)
})

test_that("se_compare() rejects family = 'cox' cleanly", {
  expect_error(
    se_compare(formula = "Surv(time, status) ~ age", data = lung_df,
               family = "cox"),
    "does not support family = \"cox\"")
})

test_that("reporting and plotting layers handle cox results", {
  s <- sca(y = c("time", "status"), x = "age",
           controls = c("sex", "ph.ecog"),
           data = lung_df, family = "cox", progress_bar = FALSE)

  td <- tidy(s)
  expect_identical(unique(td$fit_stat_name), "concordance")
  expect_equal(td$fit_stat[order(td$spec_id)],
               s$concordance[order(s$index)])

  gl <- glance(s)
  expect_identical(gl$family, "cox")

  expect_s3_class(plot_curve(s), "ggplot")

  v <- sca_variance(s)
  expect_true(all(c("sex", "ph.ecog") %in% v$choice))

  # HR / concordance / n_events are metadata, not control indicators
  expect_identical(sort(sca_control_cols(s)), sort(c("sex", "ph.ecog")))
})

test_that("sca_test() runs a shuffle_x joint test on cox curves", {
  tt <- suppressWarnings(
    sca_test(y = c("time", "status"), x = "age",
             controls = c("sex", "ph.ecog"),
             data = lung_df, family = "cox", n_permutations = 9,
             progress_bar = FALSE, seed = 42))
  expect_s3_class(tt, "sca_test")
  expect_true(all(is.finite(unlist(tt$p_values))))

  # design-preserving nulls remain linear-only
  expect_error(
    sca_test(y = c("time", "status"), x = "age", controls = "sex",
             data = lung_df, family = "cox", null_type = "freedman_lane",
             n_permutations = 9, progress_bar = FALSE),
    "linear")
})
