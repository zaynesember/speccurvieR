# Tests for sca_report()

curve <- function(controls = c("ChlorA", "O2Sat")) {
  suppressMessages(sca(y = "Salnty", x = "T_degC", controls = controls,
                       data = bottles, progress_bar = FALSE, parallel = FALSE))
}
a_test <- function(test_stats = c("median", "share_significant", "stouffer"),
                   direction = "two.sided") {
  suppressMessages(sca_test(y = "Salnty", x = "T_degC",
                            controls = c("ChlorA", "O2Sat"), data = bottles,
                            n_permutations = 50, seed = 1,
                            test_stats = test_stats, direction = direction,
                            progress_bar = FALSE))
}

test_that("sca_report on a test is a one-paragraph string with a verdict", {
  skip_on_cran()
  txt <- sca_report(a_test())
  expect_length(txt, 1L)
  expect_type(txt, "character")
  expect_match(txt, "T_degC")
  expect_match(txt, "Stouffer")
  expect_match(txt, "rejected|did not reject")
})

test_that("sca_report withholds a verdict when the trio is incomplete", {
  skip_on_cran()
  txt <- sca_report(a_test(test_stats = "median"))
  expect_false(grepl("rejected the null", txt))
  expect_match(txt, "median estimate")
})

test_that("sca_report uses directional wording when a direction is set", {
  skip_on_cran()
  txt <- sca_report(a_test(direction = "positive"))
  expect_match(txt, "predicted positive direction")
})

test_that("sca_report on a curve is descriptive", {
  txt <- sca_report(curve())
  expect_length(txt, 1L)
  expect_match(txt, "descriptive")
  expect_match(txt, "across 3 specifications")
})

test_that("sca_report.default errors helpfully", {
  expect_error(sca_report(1:3), "from sca\\(\\) or sca_test\\(\\)")
})
