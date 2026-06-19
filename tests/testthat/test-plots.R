# Tests for the shared theme and the visualization upgrades

sca_fixture <- function() {
  suppressMessages(sca("Salnty", "T_degC", c("ChlorA", "O2Sat"), bottles,
                       progressBar = FALSE))
}

test_that("theme_sca returns a ggplot theme", {
  expect_s3_class(theme_sca(), "theme")
})

test_that("plotCurve returns patchwork with the panel, ggplot without", {
  s <- sca_fixture()
  expect_s3_class(plotCurve(s), "patchwork")
  expect_s3_class(plotCurve(s, plotVars = FALSE), "ggplot")
})

test_that("plotCurve accepts medianLine and pointSize", {
  s <- sca_fixture()
  p <- plotCurve(s, plotVars = FALSE, medianLine = TRUE, pointSize = 3)
  expect_s3_class(p, "ggplot")
  # medianLine adds a horizontal reference line
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_true(sum(geoms == "GeomHline") >= 2)  # zero line + median line
})

test_that("model-fit plots combine panels with patchwork", {
  s <- sca_fixture()
  expect_s3_class(plotRMSE(s), "patchwork")
  expect_s3_class(plotRMSE(s, plotVars = FALSE), "ggplot")
})

test_that("plotControlDistributions draws a zero line only when requested", {
  s <- sca_fixture()
  has_vline <- function(p) {
    "GeomVline" %in% vapply(p$layers, function(l) class(l$geom)[1], character(1))
  }
  expect_true(has_vline(plotControlDistributions(s, zeroLine = TRUE)))
  expect_false(has_vline(plotControlDistributions(s, zeroLine = FALSE)))
})
