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

test_that("se_compare() expands a vector of bootSampleSize into multiple columns", {
  # Regression test for length(bootSampleSize == 1): a length-1 bootSamples
  # combined with a length-2 bootSampleSize must yield one column per size.
  set.seed(1)
  r <- suppressMessages(se_compare("Salnty ~ T_degC + STheta", bottles,
                                   types = "bootstrapped",
                                   bootSamples = 4,
                                   bootSampleSize = c(250, 300)))
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
                                fixedEffectsOnly = TRUE)),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_s3_class(r, "data.frame")
  expect_true(all(c("estimate_FE", "CL_FE", "CL_Depth_ID_FE") %in% colnames(r)))
  expect_length(warns, 0)
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
