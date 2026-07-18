# Tests for the formula interface to sca()

test_that("formula interface matches the vector interface (linear)", {
  v <- suppressMessages(sca("Salnty", "T_degC", c("ChlorA", "O2Sat"),
                            bottles, progress_bar = FALSE))
  f <- suppressMessages(sca(Salnty ~ T_degC + ChlorA + O2Sat,
                            data = bottles, progress_bar = FALSE))
  expect_equal(f, v)
})

test_that("data can be passed positionally with a formula", {
  named <- suppressMessages(sca(Salnty ~ T_degC + ChlorA + O2Sat,
                                data = bottles, progress_bar = FALSE))
  positional <- suppressMessages(sca(Salnty ~ T_degC + ChlorA + O2Sat,
                                     bottles, progress_bar = FALSE))
  expect_equal(positional, named)
})

test_that("formula interface handles fixed effects after `|`", {
  v <- suppressMessages(sca("Salnty", "T_degC", "ChlorA", bottles,
                            fixed_effects = "Sta_ID", progress_bar = FALSE))
  f <- suppressMessages(sca(Salnty ~ T_degC + ChlorA | Sta_ID,
                            data = bottles, progress_bar = FALSE))
  expect_equal(f, v)
})

test_that("formula interface keeps interaction terms as single control units", {
  v <- suppressMessages(sca("Salnty", "T_degC",
                            c("ChlorA", "O2Sat", "ChlorA*O2Sat"),
                            bottles, progress_bar = FALSE))
  f <- suppressMessages(sca(Salnty ~ T_degC + ChlorA + O2Sat + ChlorA*O2Sat,
                            data = bottles, progress_bar = FALSE))
  expect_equal(f, v)
})

test_that("formula interface works with return_formulae", {
  v <- sca("Salnty", "T_degC", c("ChlorA", "O2Sat"), bottles,
           return_formulae = TRUE)
  f <- sca(Salnty ~ T_degC + ChlorA + O2Sat, data = bottles,
           return_formulae = TRUE)
  expect_equal(f, v)
})

test_that("formula interface errors helpfully on bad input", {
  # one-sided formula
  expect_error(suppressMessages(sca(~ T_degC + ChlorA, data = bottles)),
               "two-sided")
  # focal variable but no controls
  expect_error(suppressMessages(sca(Salnty ~ T_degC, data = bottles)),
               "at least one control")
  # a typo'd column still triggers the column-existence check
  expect_error(suppressMessages(sca(Salnty ~ T_degC + NOPE, data = bottles)),
               "not found in data")
})

test_that("supplying x/controls alongside a formula warns", {
  expect_warning(
    suppressMessages(sca(Salnty ~ T_degC + ChlorA, x = "ignored",
                         controls = "ignored", data = bottles,
                         progress_bar = FALSE)),
    "ignored when `y` is a formula"
  )
})
