# speccurvieR (development version)

* New `sca_voe()` and `plot_voe()`: the vibration-of-effects summary and
  volcano plot of Patel, Burford, and Ioannidis (2015) computed over an
  `sca()` curve -- the spread of the estimate between the 1st and 99th
  percentile specifications (reported as a relative hazard ratio for Cox
  models), the relative p-value `RP`, and a Janus-effect flag for when the
  direction of the effect depends on the controls chosen.

* `sca()` supports survival outcomes: `family = "cox"` estimates every
  specification with a Cox proportional-hazards model via `survival::coxph()`.
  The outcome is given as `y = c("time", "status")` or as a `Surv(time,
  status)` response in the formula interface (which implies `family = "cox"`).
  The returned curve is on the log-hazard-ratio scale with the exponentiated
  hazard ratio alongside as `HR`, model fit reported via `AIC` and
  `concordance`, and an `n_events` column next to `n_obs` (in survival data,
  power tracks the number of events rather than rows). Weighted specifications
  report `coxph()`'s robust standard errors. `sca_test()` supports Cox curves
  under the default `shuffle_x` null; the design-preserving nulls and
  `se_compare()` remain linear/glm-only and reject `family = "cox"` with a
  clear message.

# speccurvieR 1.0.0

First stable release. Since 0.5.0 the package has grown from a curve-plotting
tool into a full specification-curve workflow -- estimate the curve, compare
standard errors, test it, decompose it, and report it -- and the public API is
now considered stable. The deprecated camelCase function and argument names
introduced this cycle continue to work (with a warning).

## Inference

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
* Per-specification family-wise-error-rate-adjusted p-values. When `sca_test()`
  is run with `keep_curves = TRUE` and a confound-preserving null
  (`null_type = "freedman_lane"` or `"residual_bootstrap"`), it now attaches a
  multiple-comparison-corrected p-value for *every* specification, answering
  *which* specifications are more extreme than chance once you account for having
  searched all of them. It uses the min-P / max-statistic permutation method of
  Westfall and Young (1993), reusing the permutation null already computed (no
  extra permutations), and compares each specification to its *own* null -- which
  keeps it calibrated under the confound-preserving nulls, where an
  under-controlled specification's null is not centred on zero. The new exported
  `sca_minp()` recomputes the adjustment from an existing result with a different
  method (`"single_step"`, the default, controls weak FWER; `"step_down"` is the
  more powerful Westfall-Young free step-down) or significance threshold.
  `plot_sca_test_specs()` highlights the specifications that survive correction,
  and `print()`, `sca_report()`, `sca_table()`, `glance()`, and
  `as.data.frame(what = "specs")` surface the results.

## Reporting and export

* New reporting and export tools. `tidy()` and `glance()` methods (the generics
  that `broom` and `modelsummary` dispatch on, so no `broom` dependency is
  added) are provided for `sca()` curves and `sca_test()` results, along with
  `as.data.frame()` for `sca_test`. `sca_table()` renders a publication results
  block -- a dependency-free data frame by default, or Markdown / LaTeX / `gt` /
  `kableExtra` / `flextable` -- and `sca_report()` writes a one-paragraph,
  manuscript-ready summary. `sca()` output now carries an `"sca"` class (it
  remains a data frame in every other respect) so these methods can dispatch.

## New analyses and plots

* New `sca_variance()` and `plot_variance()`: decompose the variance of the
  focal coefficient across the specification curve into the share attributable
  to each control choice (plus a residual). Uses an LMG / Shapley
  decomposition of R-squared by default (order-invariant shares that sum to the
  model R-squared), with a Type II partial sums-of-squares option; adds no new
  dependency. Analogous to `specr::icc_specs()` but suited to speccurvieR's
  binary control-inclusion choices.
* Three new diagnostic plots: `plot_influence()` shows how including vs
  excluding each control shifts the independent variable's coefficient;
  `plot_coef_fit()` plots the coefficient against model fit (are the
  best-fitting specifications outliers?); and `plot_multi_se()` draws the
  specification curve faceted by standard error type, so you can see which
  specifications stay significant under each choice of standard error.
* New `plot_se()`: visualises the output of `se_compare()`, plotting each
  coefficient's estimate with a confidence interval for every standard error
  type so the sensitivity of inference to the choice of standard error is easy
  to see.
* Visualization overhaul: a shared significance palette (a deep-purple to
  teal/amber to muted-grey ramp, so the strength of evidence reads from the
  colour) and an exported `theme_sca()` are now applied across all plots.
  Specification-curve and standard-error points are filled by significance with
  a thin white outline so they read as distinct markers against their own error
  bars; the zero reference line is a softened red and gridlines are lighter.
  `plot_curve()` and the model-fit plots combine their panels with `patchwork`
  (so the combined plot is now customisable and panels align precisely);
  `plot_curve()` gains `median_line` and `point_size` arguments; and
  `plot_control_distributions()` gains a `zero_line` argument and a cohesive
  fill. `theme_sca()` gains a `base_family` argument (defaulting to the
  `speccurvieR.base_family` option, else the device font) for setting the plot
  font without giving up portability.

## Modelling and interface

* `se_compare()` now supports two-way and multiway clustered standard errors.
  In addition to a character vector (each element a separate one-way
  clustering, unchanged), `cluster` accepts a list whose elements are each a
  character vector of dimensions to cluster on *jointly*, so
  `cluster = list("a", "b", c("a", "b"))` returns one-way SEs by `a`, by `b`,
  and two-way clustered by both. Multiway estimates use the native
  Cameron-Gelbach-Miller computation in `sandwich::vcovCL` (non-fixed-effects,
  including glm) and `fixest` (fixed effects); multiway columns are labelled
  with the dimensions joined by `_BY_` (e.g. `"HC1_a_BY_b"`).
* `se_compare()` gains `family` and `link` arguments mirroring `sca()`, so the
  full range of `glm()` model families (e.g. logistic, Poisson) can be compared
  across standard error types. The iid, heteroskedasticity-consistent,
  clustered, and bootstrapped standard error machinery all carry over to glm
  fits. Fixed effects remain OLS-only and are dropped with a warning when a
  non-linear family is supplied, matching `sca()`'s behaviour.
* `sca()` gains a formula interface:
  `sca(y ~ x + control1 + control2 | fixed_effect, data)`. The first
  right-hand-side term is the focal independent variable, the remaining terms
  are controls, and anything after `|` is treated as fixed effects. The original
  `y` / `x` / `controls` / `fixed_effects` argument interface is unchanged.
* `sca()` now reports `n_obs`, the number of observations each specification was
  fit on, and gains a `common_sample` argument. With default listwise deletion,
  specifications with different control sets can be fit on different samples
  (which conflates control effects with sample changes); `n_obs` makes this
  visible and `common_sample = TRUE` fits every specification on the rows that
  are complete across all model variables. New `plot_samplesizes()` plots the
  per-specification sample sizes.
* The package now uses snake_case throughout for both function names and
  arguments (e.g. `plotCurve()` -> `plot_curve()`, `fixedEffects` ->
  `fixed_effects`). The previous camelCase names are kept as deprecated aliases
  that warn and forward to the new names, so existing code keeps working; they
  are scheduled for removal in version 2.0.0.
* The package now ships a vignette, "Specification Curve Analysis with
  speccurvieR", touring the full workflow from `sca()` through the
  joint-inference test.

## Bug fixes

* Bug fix: `plot_sca_test_specs()` no longer shows blank legend keys for tiers
  with no specifications (e.g. when every specification survives correction);
  unused tiers are dropped from the legend.
* Bug fix: `sca()` no longer drops the focal variable (and crash with
  "subscript out of bounds") when a control variable's name contains the focal
  variable's name as a substring (e.g. focal `"Temp"` with a control `"TempX"`)
  or the focal name contains regular-expression metacharacters. The focal term
  is now matched by exact equality rather than a substring/regex test.
* Bug fix: `se_compare()` accepts a formula *object* (`y ~ x | fe`), not only a
  formula string (`"y ~ x | fe"`). Previously a formula object crashed with
  "the condition has length > 1", because `as.character()` on a formula returns
  a length-3 vector that broke the internal fixed-effects check. `se_compare()`
  also now skips an individual clustered standard-error type that fails to
  compute (e.g. a numerically singular clustered `HC3`) with a clear warning,
  instead of letting the whole call error on a raw LAPACK message.
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
* Bug fix: `sca()`'s control-indicator columns (used by `plot_vars()` and the
  `plot_curve()` control panel) are now built by exact term membership instead
  of substring matching on the model terms. Previously a control whose name was
  a substring of another term (for example `"O2"` inside `"O2Sat"`) was marked
  present in specifications that did not contain it.

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
