# Tests for the deprecated camelCase aliases

test_that("deprecated function aliases warn and forward to the new names", {
  s_new <- suppressMessages(sca("Salnty", "T_degC", c("ChlorA", "O2Sat"),
                                bottles, progress_bar = FALSE))

  # Old function name warns ...
  expect_warning(p <- plotCurve(s_new), "deprecated")
  # ... and returns the same thing as the new name
  expect_equal(class(suppressWarnings(plotCurve(s_new))),
               class(plot_curve(s_new)))

  expect_warning(plotRMSE(s_new), "deprecated")
  expect_warning(unAsIs(I(1:3)), "deprecated")
})

test_that("deprecated camelCase arguments still work and warn", {
  # old function name + old argument name (the alias warns once and remaps)
  s_new <- suppressMessages(
    sca("Salnty", "T_degC", c("ChlorA", "O2Sat"), bottles,
        progress_bar = FALSE))
  expect_warning(s_old <- plotCurve(s_new, plotVars = FALSE), "deprecated")
  expect_s3_class(s_old, "ggplot")

  # deprecated argument on a function whose name did not change (sca)
  expect_warning(
    f <- sca("Salnty", "T_degC", "O2Sat", bottles, returnFormulae = TRUE),
    "deprecated"
  )
  expect_length(f, 1)
})

test_that("old and new APIs produce identical results", {
  new <- suppressMessages(
    sca("Salnty", "T_degC", c("O2Sat", "STheta"), bottles,
        progress_bar = FALSE))
  old <- suppressWarnings(suppressMessages(
    sca("Salnty", "T_degC", c("O2Sat", "STheta"), bottles,
        progressBar = FALSE)))
  expect_equal(new, old)
})
