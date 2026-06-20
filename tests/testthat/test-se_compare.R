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
