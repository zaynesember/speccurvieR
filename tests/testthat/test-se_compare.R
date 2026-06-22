# Tests for se_compare() and se_boot()

test_that("se_compare() returns one column per requested HC type", {
  r <- suppressMessages(se_compare("Salnty ~ T_degC + STheta", bottles,
                                   types = c("HC0", "HC1", "HC3")))
  expect_s3_class(r, "data.frame")
  expect_equal(rownames(r), c("(Intercept)", "T_degC", "STheta"))
  expect_true(all(c("estimate", "HC0", "HC1", "HC3") %in% colnames(r)))
})

test_that("se_compare() warns on an invalid type without a coercion warning", {
  # Capture every warning that fires.
  warns <- character(0)
  withCallingHandlers(
    suppressMessages(se_compare("Salnty ~ T_degC", bottles,
                                types = c("HC0", "BOGUS"))),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("not a valid type", warns)))
  # Regression test for the misplaced "!= 0" inside length(): the old code
  # produced a spurious "NAs introduced by coercion" warning.
  expect_false(any(grepl("NAs introduced by coercion", warns)))
})

test_that("se_compare() expands a vector of boot_sample_size into multiple columns", {
  # Regression test for length(boot_sample_size == 1): a length-1 boot_samples
  # combined with a length-2 boot_sample_size must yield one column per size.
  set.seed(1)
  r <- suppressMessages(se_compare("Salnty ~ T_degC + STheta", bottles,
                                   types = "bootstrapped",
                                   boot_samples = 4,
                                   boot_sample_size = c(250, 300)))
  expect_true(all(c("bootstrap_k4n250", "bootstrap_k4n300") %in% colnames(r)))
})

test_that("se_compare() handles FE + cluster with specific (non-'all') types", {
  # Regression test: previously errored with "object 'types_CL' not found"
  # because the FE clustering branch referenced a variable defined only in the
  # non-FE branch. FE clustered SEs do not depend on `types`.
  warns <- character(0)
  r <- withCallingHandlers(
    suppressMessages(se_compare("Salnty ~ T_degC + STheta | Sta_ID", bottles,
                                types = c("CL_FE"), cluster = "Depth_ID",
                                fixed_effects_only = TRUE)),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_s3_class(r, "data.frame")
  expect_true(all(c("estimate_FE", "CL_FE", "CL_Depth_ID_FE") %in% colnames(r)))
  expect_length(warns, 0)
})

test_that("se_compare() does not flag CL_FE as invalid when the formula has FE", {
  # With the default fixed_effects_only = FALSE both an FE and a non-FE model are
  # fit. "CL_FE" is a fixed-effects-only type, so the non-FE branch should not
  # warn that it is invalid.
  warns <- character(0)
  withCallingHandlers(
    suppressMessages(se_compare("Salnty ~ T_degC + STheta | Sta_ID", bottles,
                                types = c("CL_FE"), cluster = "Depth_ID")),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_false(any(grepl("CL_FE", warns)))
})

test_that("se_compare() still flags CL_FE as invalid for a non-FE formula", {
  # Guard against over-suppression: without fixed effects, CL_FE is genuinely
  # not a valid type and should still warn.
  expect_warning(
    suppressMessages(se_compare("Salnty ~ T_degC", bottles, types = c("CL_FE"))),
    "not a valid type"
  )
})

test_that("se_boot() returns a named standard error per coefficient", {
  set.seed(1)
  b <- suppressMessages(se_boot(data = bottles,
                                formula = "Salnty ~ T_degC + STheta",
                                n_x = 2, n_samples = 4, sample_size = 250))
  expect_length(b, 3)
  expect_named(b, c("(Intercept)", "T_degC", "STheta"))
  expect_true(all(b >= 0))
})

# glm-family support ----------------------------------------------------------

# A binary outcome for logistic-regression tests.
bin_bottles <- within(bottles, {
  bin <- as.integer(Salnty > stats::median(Salnty, na.rm = TRUE))
})

test_that("se_compare() estimates a glm and matches glm()/sandwich exactly", {
  r <- suppressWarnings(suppressMessages(
    se_compare("bin ~ T_degC + STheta", bin_bottles, family = "binomial",
               types = c("iid", "HC0", "HC3"))))
  expect_s3_class(r, "data.frame")
  expect_true(all(c("estimate", "iid", "HC0", "HC3") %in% colnames(r)))

  # The coefficients and iid/HC SEs should reproduce a direct glm + sandwich fit.
  m <- suppressWarnings(glm(bin ~ T_degC + STheta, data = bin_bottles,
                            family = binomial()))
  expect_equal(r$estimate, unname(coef(m)))
  expect_equal(r$iid, unname(summary(m)$coefficients[, 2]))
  expect_equal(r$HC0,
               unname(lmtest::coeftest(m, vcov. = sandwich::vcovHC,
                                       type = "HC0")[, 2]))
})

test_that("se_compare() supports clustered SEs for a glm", {
  r <- suppressWarnings(suppressMessages(
    se_compare("bin ~ T_degC + STheta", bin_bottles, family = "binomial",
               types = "HC1", cluster = "Sta_ID")))
  expect_true("HC1_Sta_ID" %in% colnames(r))
})

test_that("se_compare() bootstraps SEs for a glm", {
  set.seed(3)
  r <- suppressWarnings(suppressMessages(
    se_compare("bin ~ T_degC + STheta", bin_bottles, family = "binomial",
               types = "bootstrapped", boot_samples = 4,
               boot_sample_size = 300)))
  expect_true("bootstrap_k4n300" %in% colnames(r))
})

test_that("se_compare() warns and drops fixed effects for a glm family", {
  expect_warning(
    suppressMessages(se_compare("bin ~ T_degC | Sta_ID", bin_bottles,
                                family = "binomial", types = c("iid", "HC0"))),
    "Fixed effects unsupported"
  )
  # After dropping the FE the result is a plain non-FE glm: no FE columns.
  r <- suppressWarnings(suppressMessages(
    se_compare("bin ~ T_degC | Sta_ID", bin_bottles, family = "binomial",
               types = c("iid", "HC0"))))
  expect_false(any(grepl("_FE", colnames(r))))
  expect_true(all(c("estimate", "iid", "HC0") %in% colnames(r)))
})

test_that("se_compare() treats family = 'gaussian' as linear OLS", {
  a <- suppressMessages(se_compare("Salnty ~ T_degC + STheta", bottles,
                                   family = "gaussian", types = "iid"))
  b <- suppressMessages(se_compare("Salnty ~ T_degC + STheta", bottles,
                                   types = "iid"))
  expect_equal(a, b)
})

test_that("se_compare() errors on an unrecognised family", {
  expect_error(
    se_compare("bin ~ T_degC", bin_bottles, family = "notafamily"),
    "not a recognised model family"
  )
})

test_that("se_boot() fits glm models when given a family object", {
  set.seed(4)
  b <- suppressWarnings(suppressMessages(
    se_boot(data = bin_bottles, formula = "bin ~ T_degC + STheta",
            n_x = 2, n_samples = 4, sample_size = 300,
            fam_obj = binomial())))
  expect_length(b, 3)
  expect_named(b, c("(Intercept)", "T_degC", "STheta"))
})

# Statistical-correctness regression tests ------------------------------------

# Complete cases for the bootstrap tests, so that resampling and the analytic
# reference operate on the same observations.
cc_bottles <- bottles[stats::complete.cases(
  bottles[, c("Salnty", "T_degC", "STheta", "Sta_ID")]), ]

test_that("se_compare() CL_FE is clustered by the first fixed effect, not iid", {
  # Regression test: modern fixest defaults feols() to IID standard errors, so
  # the "CL_FE" column must cluster by the first fixed effect explicitly to
  # match its name and documentation.
  r <- suppressMessages(se_compare("Salnty ~ T_degC + STheta | Sta_ID", bottles,
                                   types = "CL_FE", fixed_effects_only = TRUE))
  m <- fixest::feols(Salnty ~ T_degC + STheta | Sta_ID, data = bottles)
  clustered <- summary(m, cluster = ~Sta_ID)$coeftable[, 2]
  iid <- summary(m, vcov = "iid")$coeftable[, 2]

  clfe <- r$CL_FE[!is.na(r$CL_FE)]  # drop the NA intercept row
  expect_equal(unname(clfe), unname(clustered))
  # And it is genuinely cluster-robust, i.e. distinct from the iid default.
  expect_false(isTRUE(all.equal(unname(clfe), unname(iid))))
})

test_that("se_boot() does not collapse to zero when sample_size == n", {
  # Regression test for the old without-replacement sampling, which made every
  # resample identical (hence SE ~ 0) when sample_size equalled the row count.
  set.seed(1)
  n <- nrow(cc_bottles)
  b <- suppressMessages(se_boot(cc_bottles, "Salnty ~ T_degC + STheta",
                                n_x = 2, n_samples = 30, sample_size = n))
  expect_true(all(b > 1e-6))
})

test_that("se_boot() full-n bootstrap recovers the heteroskedasticity-robust SE", {
  # The nonparametric pairs bootstrap estimates the robust (not iid) SE. With
  # replacement at sample_size = n this should track HC1 closely.
  set.seed(20)
  n <- nrow(cc_bottles)
  b <- suppressMessages(se_boot(cc_bottles, "Salnty ~ T_degC + STheta",
                                n_x = 2, n_samples = 800, sample_size = n))
  m <- lm(Salnty ~ T_degC + STheta, cc_bottles)
  hc1 <- lmtest::coeftest(m, vcov. = sandwich::vcovHC, type = "HC1")[, 2]
  expect_equal(unname(b), unname(hc1), tolerance = 0.2)
})

test_that("se_boot() rescales an m-out-of-n bootstrap to the full-sample SE", {
  # With sqrt(m/n) rescaling, a quarter-sample bootstrap targets the same SE as
  # a full-sample bootstrap (previously it overstated it by ~sqrt(n/m) = 2x).
  n <- nrow(cc_bottles)
  set.seed(30)
  full <- suppressMessages(se_boot(cc_bottles, "Salnty ~ T_degC + STheta",
                                   n_x = 2, n_samples = 800, sample_size = n))
  set.seed(31)
  quarter <- suppressMessages(se_boot(cc_bottles, "Salnty ~ T_degC + STheta",
                                      n_x = 2, n_samples = 800,
                                      sample_size = round(n / 4)))
  expect_equal(unname(quarter), unname(full), tolerance = 0.25)
})

# --- Two-way / multiway clustered standard errors -----------------------------

test_that("se_compare() two-way clustering matches sandwich::vcovCL (non-FE)", {
  r <- suppressWarnings(se_compare("Salnty ~ T_degC + ChlorA", bottles,
                                   types = "HC1",
                                   cluster = list(c("Sta_ID", "Depth_ID"))))
  expect_true("HC1_Depth_ID_BY_Sta_ID" %in% colnames(r))
  m <- lm(Salnty ~ T_degC + ChlorA, bottles)
  ref <- lmtest::coeftest(m, vcov. = sandwich::vcovCL, type = "HC1",
                          cluster = ~ Sta_ID + Depth_ID)[, 2]
  expect_equal(unname(r[["HC1_Depth_ID_BY_Sta_ID"]]), unname(ref))
})

test_that("se_compare() two-way clustering matches fixest (FE)", {
  r <- suppressMessages(suppressWarnings(
    se_compare("Salnty ~ T_degC + ChlorA | Sta_ID", bottles,
               types = "CL_FE", cluster = list(c("Sta_ID", "Depth_ID")))))
  expect_true("CL_Depth_ID_BY_Sta_ID_FE" %in% colnames(r))
  m <- suppressMessages(feols(Salnty ~ T_degC + ChlorA | Sta_ID, bottles))
  ref <- summary(m, cluster = ~ Sta_ID + Depth_ID)$coeftable[, 2]
  got <- r[["CL_Depth_ID_BY_Sta_ID_FE"]]
  got <- got[!is.na(got)]
  expect_equal(unname(got), unname(ref))
})

test_that("se_compare() two-way clustering works for a glm", {
  d <- bottles
  d$saline <- as.integer(d$Salnty > stats::median(d$Salnty, na.rm = TRUE))
  r <- suppressWarnings(se_compare("saline ~ T_degC + ChlorA", d,
                                   family = "binomial", types = "HC0",
                                   cluster = list(c("Sta_ID", "Depth_ID"))))
  m <- glm(saline ~ T_degC + ChlorA, d, family = stats::binomial())
  ref <- lmtest::coeftest(m, vcov. = sandwich::vcovCL, type = "HC0",
                          cluster = ~ Sta_ID + Depth_ID)[, 2]
  expect_equal(unname(r[["HC0_Depth_ID_BY_Sta_ID"]]), unname(ref))
})

test_that("se_compare() mixes one-way and multiway clustering in one call", {
  r <- suppressWarnings(se_compare("Salnty ~ T_degC + ChlorA", bottles,
                  types = "HC1",
                  cluster = list("Sta_ID", "Depth_ID", c("Sta_ID", "Depth_ID"))))
  expect_true(all(c("HC1_Sta_ID", "HC1_Depth_ID", "HC1_Depth_ID_BY_Sta_ID")
                  %in% colnames(r)))
})

test_that("se_compare() clustering is back-compatible (vector = one-way each)", {
  # A character vector still produces one separate one-way clustering per
  # element, with the historical column names and values.
  vec <- suppressWarnings(se_compare("Salnty ~ T_degC + ChlorA", bottles,
                                     types = "HC1",
                                     cluster = c("Sta_ID", "Depth_ID")))
  expect_equal(colnames(vec),
               c("estimate", "HC1", "HC1_Sta_ID", "HC1_Depth_ID"))
  m <- lm(Salnty ~ T_degC + ChlorA, bottles)
  expect_equal(unname(vec[["HC1_Sta_ID"]]),
               unname(lmtest::coeftest(m, vcov. = sandwich::vcovCL, type = "HC1",
                                       cluster = ~ Sta_ID)[, 2]))
  # A length-1 list element is identical to the bare-vector one-way (dedup).
  lst <- suppressWarnings(se_compare("Salnty ~ T_degC + ChlorA", bottles,
                                     types = "HC1", cluster = list("Sta_ID")))
  bare <- suppressWarnings(se_compare("Salnty ~ T_degC + ChlorA", bottles,
                                      types = "HC1", cluster = "Sta_ID"))
  expect_equal(colnames(lst), colnames(bare))
  expect_equal(lst, bare)
})

test_that("se_compare() collapses duplicate clustering specifications", {
  # c("a","b") and c("b","a") are the same joint clustering -> one column.
  r <- suppressWarnings(se_compare("Salnty ~ T_degC", bottles, types = "HC1",
                  cluster = list(c("Sta_ID", "Depth_ID"),
                                 c("Depth_ID", "Sta_ID"))))
  expect_equal(sum(grepl("_BY_", colnames(r))), 1L)
})

test_that("se_compare() multiway clustering is row-aligned under NA dropping", {
  set.seed(7)
  d <- data.frame(y = stats::rnorm(60), x = stats::rnorm(60),
                  g1 = rep(1:6, 10), g2 = rep(1:10, each = 6))
  d$y[c(3, 17, 40)] <- NA           # NAs in a model variable
  r <- suppressWarnings(se_compare("y ~ x", d, types = "HC1",
                                   cluster = list(c("g1", "g2"))))
  m <- lm(y ~ x, d)                 # drops the NA rows
  ref <- lmtest::coeftest(m, vcov. = sandwich::vcovCL, type = "HC1",
                          cluster = ~ g1 + g2)[, 2]
  expect_equal(unname(r[["HC1_g1_BY_g2"]]), unname(ref))
})

test_that("se_compare() warns and drops an unknown dimension in a list element", {
  expect_warning(
    r <- se_compare("Salnty ~ T_degC", bottles, types = "HC1",
                    cluster = list(c("Sta_ID", "not_a_col"))),
    "not a valid clustering variable")
  # The unknown dim is dropped, leaving a one-way clustering by Sta_ID.
  expect_true("HC1_Sta_ID" %in% colnames(r))
  expect_false(any(grepl("not_a_col", colnames(r))))
})

test_that("se_compare() errors on a non-character cluster specification", {
  expect_error(se_compare("Salnty ~ T_degC", bottles, cluster = list(1:2)),
               "character")
})

test_that("plot_se() handles a multiway-clustered se_compare() result", {
  r <- suppressWarnings(se_compare("Salnty ~ T_degC + ChlorA", bottles,
                                   types = c("HC1"),
                                   cluster = list("Sta_ID",
                                                  c("Sta_ID", "Depth_ID"))))
  expect_s3_class(plot_se(r), "ggplot")
})

test_that("normalize_cluster_spec() cleans specs", {
  cols <- c("a", "b", "c")
  # Character vector -> list of one-way specs.
  expect_equal(speccurvieR:::normalize_cluster_spec(c("a", "b"), cols),
               list("a", "b"))
  # Dims de-duplicated and sorted within a spec; empty specs dropped.
  expect_equal(speccurvieR:::normalize_cluster_spec(list(c("b", "a", "a")), cols),
               list(c("a", "b")))
  expect_null(speccurvieR:::normalize_cluster_spec(list(character(0)), cols))
  expect_null(speccurvieR:::normalize_cluster_spec(NULL, cols))
  # Unknown dims dropped (with a warning); duplicate specs collapsed.
  expect_warning(
    out <- speccurvieR:::normalize_cluster_spec(list(c("a", "z"), c("a")), cols),
    "not a valid clustering variable")
  expect_equal(out, list("a"))
})

# --- formula-object interface + clustered-HC guard ----------------------------

test_that("se_compare() accepts a formula object as well as a string", {
  # Regression: a formula object used to crash at `has_fe <- grepl("|", formula)`
  # ("the condition has length > 1") because as.character() on a formula is a
  # length-3 vector. The formula and string forms must agree.
  expect_equal(
    suppressMessages(se_compare(Salnty ~ T_degC + ChlorA, bottles,
                                types = c("iid", "HC1", "HC3"))),
    suppressMessages(se_compare("Salnty ~ T_degC + ChlorA", bottles,
                                types = c("iid", "HC1", "HC3"))))
  # Fixed-effects formula object (the `| fe` must survive deparsing).
  expect_equal(
    suppressMessages(suppressWarnings(
      se_compare(Salnty ~ T_degC + ChlorA | Sta_ID, bottles, types = "CL_FE"))),
    suppressMessages(suppressWarnings(
      se_compare("Salnty ~ T_degC + ChlorA | Sta_ID", bottles,
                 types = "CL_FE"))))
  # Clustered formula object.
  expect_equal(
    suppressWarnings(se_compare(Salnty ~ T_degC, bottles, types = "HC1",
                                cluster = list(c("Sta_ID", "Depth_ID")))),
    suppressWarnings(se_compare("Salnty ~ T_degC", bottles, types = "HC1",
                                cluster = list(c("Sta_ID", "Depth_ID")))))
})

test_that("se_compare() skips a clustered SE type that fails, keeping the rest", {
  # A clustered HC type (e.g. HC3) can fail with a singular-matrix LAPACK error
  # on some designs. The failing column is skipped with a clear warning rather
  # than crashing the whole call. Simulate the failure by mocking coeftest().
  skip_if_not_installed("lmtest")
  real_ct <- lmtest::coeftest
  testthat::local_mocked_bindings(
    coeftest = function(x, vcov., type, ...){
      if(!missing(type) && identical(type, "HC3")){
        stop("system is computationally singular")
      }
      real_ct(x, vcov. = vcov., type = type, ...)
    },
    .package = "speccurvieR")
  expect_warning(
    r <- se_compare("Salnty ~ T_degC", bottles, types = c("HC1", "HC3"),
                    cluster = "Sta_ID", clustered_only = TRUE),
    "HC3 standard errors clustered by Sta_ID could not be computed")
  expect_true("HC1_Sta_ID" %in% colnames(r))   # the stable type survives
  expect_false("HC3_Sta_ID" %in% colnames(r))  # the failing type is dropped
})
