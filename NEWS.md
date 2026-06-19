# speccurvieR (development version)

* `sca()` gains a formula interface:
  `sca(y ~ x + control1 + control2 | fixedEffect, data)`. The first
  right-hand-side term is the focal independent variable, the remaining terms
  are controls, and anything after `|` is treated as fixed effects. The original
  `y` / `x` / `controls` / `fixedEffects` argument interface is unchanged.
* New `plotSE()`: visualises the output of `se_compare()`, plotting each
  coefficient's estimate with a confidence interval for every standard error
  type so the sensitivity of inference to the choice of standard error is easy
  to see.
* Visualization overhaul: a shared, colour-blind-safe palette and an exported
  `theme_sca()` are now applied across all plots; `plotCurve()` and the
  model-fit plots combine their panels with `patchwork` (so the combined plot
  is now customisable and panels align precisely); `plotCurve()` gains
  `medianLine` and `pointSize` arguments; and `plotControlDistributions()` gains
  a `zeroLine` argument and a cohesive fill.

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
