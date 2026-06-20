# Tests for the helper functions in R/helpers.R

test_that("formula_builder returns the full powerset of controls", {
  # 2 controls -> 2^2 - 1 = 3 specifications
  f2 <- formula_builder("y", "x", c("a", "b"))
  expect_length(f2, 3)
  # 3 controls -> 2^3 - 1 = 7 specifications
  f3 <- formula_builder("y", "x", c("a", "b", "c"))
  expect_length(f3, 7)
  # Each element should be a formula
  expect_true(all(vapply(f3, inherits, logical(1), "formula")))
})

test_that("formula_builder appends fixed effects with a pipe", {
  f <- formula_builder("y", "x", c("a"), fixed_effects = "fe")
  expect_true(any(grepl("|", as.character(f), fixed = TRUE)))
})

test_that("paste_factory builds the RHS string", {
  expect_equal(paste_factory(c("a", "b"), "x"), "x + a + b")
  # When the independent variable already appears (e.g. an interaction),
  # it should not be duplicated.
  expect_equal(paste_factory(c("x*a"), "x"), "x*a")
})

test_that("duplicate_remover drops standalone copies of interacted controls", {
  # "b" is already present inside the "x*b" interaction, so it is removed.
  expect_equal(duplicate_remover(c("a", "b", "x*b"), "x"), c("a", "x*b"))
  # With no interactions the input is returned unchanged.
  expect_equal(duplicate_remover(c("a", "b"), "x"), c("a", "b"))
})

test_that("control_extractor returns control coefs without intercept or x", {
  m <- summary(lm(Salnty ~ STheta + T_degC, bottles))
  ce <- control_extractor(m, "STheta")
  expect_named(ce, c("coef", "term"))
  expect_false("(Intercept)" %in% ce$term)
  expect_false("STheta" %in% ce$term)
  expect_true("T_degC" %in% ce$term)
})

test_that("un_as_is strips the AsIs class", {
  x <- I(1:4)
  expect_true(inherits(x, "AsIs"))
  expect_false(inherits(un_as_is(x), "AsIs"))
})

test_that("scp returns a data frame and label vector", {
  s <- suppressMessages(sca(y = "Salnty", x = "T_degC",
                            controls = c("O2Sat", "STheta"),
                            data = bottles, progress_bar = FALSE))
  out <- scp(s)
  expect_length(out, 2)
  expect_true(is.data.frame(out[[1]]))
  expect_true("controlID" %in% names(out[[1]]))
})
