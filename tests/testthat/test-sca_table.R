# Tests for sca_table()

curve <- function() {
  suppressMessages(sca(y = "Salnty", x = "T_degC",
                       controls = c("ChlorA", "O2Sat"), data = bottles,
                       progress_bar = FALSE, parallel = FALSE))
}
a_test <- function() {
  suppressMessages(sca_test(y = "Salnty", x = "T_degC",
                            controls = c("ChlorA", "O2Sat"), data = bottles,
                            n_permutations = 50, seed = 1, progress_bar = FALSE))
}

test_that("sca_table() on a test returns an sca_table results block", {
  skip_on_cran()
  tbl <- sca_table(a_test())
  expect_s3_class(tbl, "sca_table")
  expect_s3_class(tbl, "data.frame")
  expect_named(tbl, c("label", "value"))
  expect_true(all(c("Specifications", "Median estimate", "Stouffer Z",
                    "Null hypothesis") %in% tbl$label))
  expect_match(attr(tbl, "note"), "resolution floor")
})

test_that("sca_table() on a curve is descriptive (no inference rows)", {
  tbl <- sca_table(curve())
  expect_true(all(c("Estimate range", "Sign agreement") %in% tbl$label))
  expect_false("Stouffer Z" %in% tbl$label)
  expect_match(attr(tbl, "note"), "Descriptive only")
})

test_that("print.sca_table renders aligned text", {
  expect_output(print(sca_table(curve())), "Median estimate")
})

test_that("sca_table() renders markdown and latex via knitr", {
  skip_if_not_installed("knitr")
  expect_type(sca_table(curve(), format = "markdown"), "character")
  expect_type(sca_table(curve(), format = "latex"), "character")
})

test_that("sca_table.default errors helpfully", {
  expect_error(sca_table(1:3), "from sca\\(\\) or sca_test\\(\\)")
})

test_that("Suggests-gated renderers work when installed", {
  s <- curve()
  skip_if_not_installed("gt")
  expect_s3_class(sca_table(s, format = "gt"), "gt_tbl")
})

test_that("flextable rendering works when installed", {
  s <- curve()
  skip_if_not_installed("flextable")
  expect_s3_class(sca_table(s, format = "flextable"), "flextable")
})

test_that("modelsummary consumes the tidy/glance methods", {
  skip_on_cran()
  skip_if_not_installed("modelsummary")
  r <- a_test()
  ms <- modelsummary::modelsummary(r, output = "data.frame")
  expect_s3_class(ms, "data.frame")
})
