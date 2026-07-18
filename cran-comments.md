## Release summary

This is a major release (1.0.0). Since 0.5.0 it adds a permutation-based
joint-inference test (`sca_test()`) with design-preserving nulls and a
per-specification family-wise-error-rate correction (`sca_minp()`), a variance
decomposition (`sca_variance()`), a reporting layer (`tidy()` / `glance()`
methods, `sca_table()`, `sca_report()`), generalized-linear-model support and
two-way / multiway clustered standard errors in `se_compare()`, a formula
interface for `sca()`, and several bug fixes. The full list of changes is in
NEWS.md. The public API is now considered stable; the camelCase function and
argument names introduced in earlier versions remain available as deprecated
aliases. The package has no reverse dependencies.

## Test environments

* Local: macOS, R-release.
* win-builder: R-release and R-devel.

## R CMD check results

0 errors | 0 warnings | 0 notes
