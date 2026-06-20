# Tests for the shared theme and the visualization upgrades

sca_fixture <- function() {
  suppressMessages(sca("Salnty", "T_degC", c("ChlorA", "O2Sat"), bottles,
                       progress_bar = FALSE))
}

test_that("theme_sca returns a ggplot theme", {
  expect_s3_class(theme_sca(), "theme")
})

test_that("plot_curve returns patchwork with the panel, ggplot without", {
  s <- sca_fixture()
  expect_s3_class(plot_curve(s), "patchwork")
  expect_s3_class(plot_curve(s, plot_vars = FALSE), "ggplot")
})

test_that("plot_curve accepts median_line and point_size", {
  s <- sca_fixture()
  p <- plot_curve(s, plot_vars = FALSE, median_line = TRUE, point_size = 3)
  expect_s3_class(p, "ggplot")
  # median_line adds a horizontal reference line
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_true(sum(geoms == "GeomHline") >= 2)  # zero line + median line
})

test_that("model-fit plots combine panels with patchwork", {
  s <- sca_fixture()
  expect_s3_class(plot_rmse(s), "patchwork")
  expect_s3_class(plot_rmse(s, plot_vars = FALSE), "ggplot")
})

test_that("plot_control_distributions draws a zero line only when requested", {
  s <- sca_fixture()
  has_vline <- function(p) {
    "GeomVline" %in% vapply(p$layers, function(l) class(l$geom)[1], character(1))
  }
  expect_true(has_vline(plot_control_distributions(s, zero_line = TRUE)))
  expect_false(has_vline(plot_control_distributions(s, zero_line = FALSE)))
})
