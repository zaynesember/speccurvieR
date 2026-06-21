# speccurvieR (development version)

* `sca_test()` gains a `null_type` argument with two design-preserving nulls
  alongside the default `"shuffle_x"`. `"freedman_lane"` (Freedman & Lane 1983)
  permutes the residuals of a reduced model `y ~ controls` (the focal variable
  omitted) and refits the curve, preserving the focal variable's correlation
  with the controls -- so it stays calibrated under collinearity, where
  shuffling the focal variable is anti-conservative. `"residual_bootstrap"`
  imposes the null on the response (`y - betahat * x`) and resamples rows, per
  Simonsohn, Simmons, and Nelson's (2020) observational scheme. Both are for
  linear models (including fixed effects) and fit every specification on a
  common sample; a new `common_sample` argument is also exposed.

* `sca()` now reports `n_obs`, the number of observations each specification was
  fit on, and gains a `common_sample` argument. With default listwise deletion,
  specifications with different control sets can be fit on different samples
  (which conflates control effects with sample changes); `n_obs` makes this
  visible and `common_sample = TRUE` fits every specification on the rows that
  are complete across all model variables. New `plot_samplesizes()` plots the
  per-specification sample sizes.

* New `sca_variance()` and `plot_variance()`: decompose the variance of the
  focal coefficient across the specification curve into the share attributable
  to each control choice (plus a residual). Uses an LMG / Shapley
  decomposition of R-squared by default (order-invariant shares that sum to the
  model R-squared), with a Type II partial sums-of-squares option; adds no new
  dependency. Analogous to `specr::icc_specs()` but suited to speccurvieR's
  binary control-inclusion choices.
* Bug fix: `sca()`'s control-indicator columns (used by `plot_vars()` and the
  `plot_curve()` control panel) are now built by exact term membership instead
  of substring matching on the model terms. Previously a control whose name was
  a substring of another term (for example `"O2"` inside `"O2Sat"`) was marked
  present in specifications that did not contain it.

* New `sca_test()`: a permutation-based joint-inference test for the whole
  specification curve, following Simonsohn, Simmons, and Nelson (2020). It
  tests the sharp null that the focal variable has no effect in any
  specification by repeatedly permuting that variable (blocked within fixed
  effects when present) and re-estimating the entire curve, then comparing the
  observed curve to the resulting null distribution via three statistics: the
  median estimate, the share of specifications significant in the predicted
  direction, and a Stouffer combination of the per-specification p-values.
  Supports the same model interface as `sca()` (formula or arguments, glm
  families, fixed effects, weights), one- and two-sided tests, parallelisation,
  and reproducible seeding. `plot_sca_test()` visualises each statistic's null
  distribution against the observed value, and `plot_sca_test_specs()` (with
  `sca_test(keep_curves = TRUE)`) draws the specification curve against each
  specification's own permutation null band.
* The package now ships a vignette, "Specification Curve Analysis with
  speccurvieR", touring the full workflow from `sca()` through the
  joint-inference test.
* The package now uses snake_case throughout for both function names and
  arguments (e.g. `plotCurve()` -> `plot_curve()`, `fixedEffects` ->
  `fixed_effects`). The previous camelCase names are kept as deprecated aliases
  that warn and forward to the new names, so existing code keeps working; they
  will be removed in a future release.
* `sca()` gains a formula interface:
  `sca(y ~ x + control1 + control2 | fixed_effect, data)`. The first
  right-hand-side term is the focal independent variable, the remaining terms
  are controls, and anything after `|` is treated as fixed effects. The original
  `y` / `x` / `controls` / `fixed_effects` argument interface is unchanged.
* Bug fix: bootstrapped standard errors are now estimated with a proper
  bootstrap. Resamples are drawn *with* replacement and the standard deviation
  of the resampled coefficients is rescaled by `sqrt(sample_size / nrow(data))`
  (an m-out-of-n bootstrap). Previously resamples were drawn without
  replacement and were not rescaled, which overstated the standard error by
  roughly `sqrt(nrow(data) / sample_size)` and collapsed it to ~0 when
  `boot_sample_size` equalled the number of rows.
* Bug fix: the `"CL_FE"` column of `se_compare()` is now clustered by the first
  fixed effect, as its name and documentation describe. Modern `fixest`
  (>= 0.10) defaults `feols()` to IID standard errors, so the previous code
  (which read off the `feols()` default) was returning IID rather than
  cluster-robust standard errors under that label.
* Bug fix: `sca()` now raises an informative error when the focal variable `x`
  does not correspond to a single model coefficient (for example a factor or
  interaction term), instead of failing with a cryptic "subscript out of
  bounds" error or silently returning `NA` for every specification.
* `se_compare()` gains `family` and `link` arguments mirroring `sca()`, so the
  full range of `glm()` model families (e.g. logistic, Poisson) can be compared
  across standard error types. The iid, heteroskedasticity-consistent,
  clustered, and bootstrapped standard error machinery all carry over to glm
  fits. Fixed effects remain OLS-only and are dropped with a warning when a
  non-linear family is supplied, matching `sca()`'s behaviour.
* New `plot_se()`: visualises the output of `se_compare()`, plotting each
  coefficient's estimate with a confidence interval for every standard error
  type so the sensitivity of inference to the choice of standard error is easy
  to see.
* Three new diagnostic plots: `plot_influence()` shows how including vs
  excluding each control shifts the independent variable's coefficient;
  `plot_coef_fit()` plots the coefficient against model fit (are the
  best-fitting specifications outliers?); and `plot_multi_se()` draws the
  specification curve faceted by standard error type, so you can see which
  specifications stay significant under each choice of standard error.
* Visualization overhaul: a shared, colour-blind-safe palette and an exported
  `theme_sca()` are now applied across all plots; `plot_curve()` and the
  model-fit plots combine their panels with `patchwork` (so the combined plot
  is now customisable and panels align precisely); `plot_curve()` gains
  `median_line` and `point_size` arguments; and `plot_control_distributions()`
  gains a `zero_line` argument and a cohesive fill.

# speccurvieR 0.5.0

* `sca()` and `se_compare()` now validate that the supplied column names exist
  in `data`, failing with an informative message instead of a cryptic
  deep-stack error.
* `sca()` accepts `family = "gaussian"` as an alias for OLS and defaults to a
  glm family's canonical link when `link` is not supplied (e.g.
  `family = "binomial"` no longer requires explicitly specifying the logit
  link). An unrecognised family now raises a clear error.
* `se_compare()` no longer errors on a fixed-effects model when clustering with
  specific (non-`"all"`) standard error types, and no longer warns that the
  fixed-effects-only `"CL_FE"` type is invalid when the formula has fixed
  effects.
* Fixed latent bugs in `se_compare()` standard-error type validation and
  bootstrap sample-size expansion, a stray debug print, and a crash when
  combining fixed effects with a non-linear family in `sca()`.
* Added a `testthat` test suite and refactored the internals of `sca()` and
  `se_compare()`; no change to their output for valid input.

# speccurvieR 0.4.2

Added support for weights in se_compare().

# speccurvieR 0.4.1

Added support for weights in sca().

# speccurvieR 0.4.0

* Added functionality to compare different types of standard error estimates.

# speccurvieR 0.3.0

* Minor documentation tweaks.

# speccurvieR 0.2.0

* Fixed effects estimation switched from lfe::felm() to fixest::feols().
* Minor improvements to model estimate extraction.

# speccurvieR 0.1.0

* Initial CRAN submission.
