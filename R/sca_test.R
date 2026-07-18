# Joint-inference test for specification curve analysis----------------------

# Internal: permutation index vector for the focal variable. With no block a
# global permutation of 1:n; with a block (a list of per-group row-index
# vectors) the rows are permuted only within their group, so exchangeability is
# respected within fixed-effect blocks. Groups of size one are left untouched.
sca_test_perm_index <- function(n, block_groups){
  if(is.null(block_groups)) return(sample.int(n))
  idx <- seq_len(n)
  for(g in block_groups){
    if(length(g) > 1L) idx[g] <- g[sample.int(length(g))]
  }
  idx
}

# Internal: a canonical key identifying each specification by the set of model
# terms it includes (the focal x and the intercept are dropped, since they are
# common to every specification). Because sca() reorders its rows by coefficient
# magnitude, this terms-based key -- not row position -- is what aligns the same
# specification across the observed curve and the permuted curves. `terms_list`
# is sca()'s `terms` list-column.
#
# The key is unique per *distinct fitted model*, not necessarily per sca() row:
# a redundant control set (e.g. an interaction control together with its
# components, controls = c("a*b", "a", "b")) yields several control subsets that
# expand to the same model terms and hence the same key. That is exactly why
# matching by key (taking the first match) is safe in the keep_curves block --
# the colliding rows are the identical model, so they share the same coefficient.
sca_spec_key <- function(terms_list, x){
  unname(vapply(terms_list,
                function(t) paste(sort(setdiff(t, c("(Intercept)", x))),
                                  collapse = " + "),
                character(1)))
}

# Internal: compute the joint-inference test statistics for a single
# specification curve (the data frame returned by sca()). Returns a named
# numeric vector keyed by the requested statistics, with an "n_valid" attribute
# giving the number of specifications with a non-missing focal coefficient. NA
# specifications (a focal coefficient or p-value that failed to estimate) are
# dropped from the relevant denominators, identically for the observed and the
# permuted curves.
sca_test_statistics <- function(curve, test_stats, direction, alpha){
  coef <- curve$coef
  p <- curve$p
  ok_p <- !is.na(p)
  n_valid <- sum(!is.na(coef))

  out <- stats::setNames(rep(NA_real_, length(test_stats)), test_stats)

  if("median" %in% test_stats){
    out[["median"]] <- stats::median(coef, na.rm = TRUE)
  }

  if("share_significant" %in% test_stats){
    num <- if(direction == "positive"){
             sum(p < alpha & coef > 0, na.rm = TRUE)
           } else if(direction == "negative"){
             sum(p < alpha & coef < 0, na.rm = TRUE)
           } else {
             sum(p < alpha, na.rm = TRUE)
           }
    denom <- sum(ok_p)
    out[["share_significant"]] <- if(denom > 0) num / denom else NA_real_
  }

  if("stouffer" %in% test_stats){
    # Per-spec signed z from the two-sided p-value, oriented to the reference
    # direction; combined via Stouffer's method. The null p-value for this
    # statistic is taken from the permutation distribution (see sca_test()),
    # not pnorm(), because specifications are dependent. The upper tail is
    # capped so an astronomically significant specification (p so small that
    # 1 - p/2 rounds to 1) yields a finite z (qnorm(1 - 1e-16) ~ 8.2) rather
    # than Inf.
    q <- pmin(1 - p[ok_p] / 2, 1 - 1e-16)
    z <- stats::qnorm(q) * sign(coef[ok_p])
    if(direction == "negative") z <- -z
    n_z <- sum(ok_p)
    out[["stouffer"]] <- if(n_z > 0) sum(z) / sqrt(n_z) else NA_real_
  }

  if("share_sign" %in% test_stats){
    num <- if(direction == "positive"){
             sum(coef > 0, na.rm = TRUE)
           } else if(direction == "negative"){
             sum(coef < 0, na.rm = TRUE)
           } else {
             max(sum(coef > 0, na.rm = TRUE), sum(coef < 0, na.rm = TRUE))
           }
    out[["share_sign"]] <- if(n_valid > 0) num / n_valid else NA_real_
  }

  attr(out, "n_valid") <- n_valid
  out
}

# Internal: finite-sample permutation p-values, one per statistic, using the
# (1 + #{null as extreme as observed}) / (n_used + 1) convention.
sca_test_pvalues <- function(observed, null_df, test_stats, direction){
  pv <- stats::setNames(rep(NA_real_, length(test_stats)), test_stats)
  for(s in test_stats){
    o <- observed[[s]]
    nulls <- null_df[[s]]
    nulls <- nulls[!is.na(nulls)]
    # Effective null count for this statistic: permutations whose value for
    # this statistic is non-missing. Keeping the numerator and denominator over
    # the same valid-null set means a curve that is NA for one statistic does
    # not silently inflate that statistic's p-value.
    n_eff <- length(nulls)
    count <- if(s == "median"){
               if(direction == "positive") sum(nulls >= o)
               else if(direction == "negative") sum(nulls <= o)
               else sum(abs(nulls) >= abs(o))
             } else if(s == "stouffer"){
               if(direction == "two.sided") sum(abs(nulls) >= abs(o))
               else sum(nulls >= o)
             } else {
               # share_significant, share_sign: one-sided greater
               sum(nulls >= o)
             }
    pv[[s]] <- if(n_eff > 0) (1 + count) / (n_eff + 1) else NA_real_
  }
  pv
}

# Internal: the common-sample subset of `data` -- the rows complete across every
# model variable (response, focal, controls, fixed effects, weights). This is
# the same union sca() uses for common_sample = TRUE, and the fixed design the
# Freedman-Lane and residual-bootstrap nulls require.
sca_test_common_sample <- function(data, y, x, controls, fixed_effects, weights){
  vars <- unique(trimws(unlist(strsplit(c(y, x, controls), "[*:]"))))
  union_vars <- unique(c(vars, fixed_effects, weights))
  data[stats::complete.cases(data[, union_vars, drop = FALSE]), , drop = FALSE]
}

# Internal: fit the Freedman-Lane reduced (nuisance) model `y ~ controls (+ FE)`
# with the focal x omitted, on the common sample, once. Returns the pieces needed
# to rebuild y* per permutation: the fitted values and raw residuals on the rows
# the model kept, those row indices into `data_cs`, the per-row weights and their
# square roots (for the weighted-residual permutation), the FE block groups over
# the kept rows, and the kept-row count.
sca_test_reduced_fit <- function(data_cs, y, controls, fixed_effects, weights,
                                 block){
  rhs <- paste(controls, collapse = " + ")
  w <- if(is.null(weights)) NULL else data_cs[[weights]]
  if(is.null(fixed_effects)){
    fml <- stats::as.formula(paste(y, "~", rhs))
    fit <- if(is.null(w)) stats::lm(fml, data = data_cs)
           else stats::lm(fml, data = data_cs, weights = w)
    kept <- seq_len(nrow(data_cs))
  } else {
    fml <- stats::as.formula(paste(y, "~", rhs, "|",
                                   paste(fixed_effects, collapse = " + ")))
    fit <- if(is.null(w)) feols(fml, data = data_cs)
           else feols(fml, data = data_cs, weights = w)
    kept <- seq_len(nrow(data_cs))
    if(!is.null(fit$obs_selection$obsRemoved)){
      kept <- kept[fit$obs_selection$obsRemoved]  # obsRemoved are negative indices
    }
  }
  w_kept <- if(is.null(w)) rep(1, length(kept)) else w[kept]
  block_groups <- if(!is.null(block)){
    unname(split(seq_along(kept), data_cs[[block]][kept]))
  } else NULL
  list(fitted = as.numeric(stats::fitted(fit)),
       resid = as.numeric(stats::residuals(fit)),
       kept = kept, w_kept = w_kept, sqrt_w = sqrt(w_kept),
       block_groups = block_groups, n_kept = length(kept))
}

# Internal: impose the residual-bootstrap null on the response. Fits the full
# control-superset spec `y ~ x + controls (+ FE)` once on the common sample and
# returns `y - betahat_super * x` over all common-sample rows.
sca_test_ynull <- function(data_cs, y, x, controls, fixed_effects, weights){
  rhs <- paste(c(x, controls), collapse = " + ")
  w <- if(is.null(weights)) NULL else data_cs[[weights]]
  if(is.null(fixed_effects)){
    fml <- stats::as.formula(paste(y, "~", rhs))
    fit <- if(is.null(w)) stats::lm(fml, data = data_cs)
           else stats::lm(fml, data = data_cs, weights = w)
  } else {
    fml <- stats::as.formula(paste(y, "~", rhs, "|",
                                   paste(fixed_effects, collapse = " + ")))
    fit <- if(is.null(w)) feols(fml, data = data_cs)
           else feols(fml, data = data_cs, weights = w)
  }
  beta <- stats::coef(fit)[[x]]
  data_cs[[y]] - beta * data_cs[[x]]
}

#' Joint-inference test for a specification curve
#'
#' @description
#' `sca_test()` performs the permutation-based joint-inference test of
#' Simonsohn, Simmons, and Nelson (2020). It tests the sharp null hypothesis
#' that the focal independent variable `x` has no effect on `y` in any
#' specification. The focal variable is repeatedly permuted (shuffled) and the
#' entire specification curve is re-estimated on each permuted data set,
#' building a null distribution for curve-level summary statistics. The
#' observed statistics are then compared to that null distribution to obtain
#' permutation p-values.
#'
#' The permutation is *blocked within the first fixed effect* when
#' `fixed_effects` are supplied, because exchangeability of `x` only holds
#' within fixed-effect groups; otherwise it is a global permutation.
#'
#' @param y A string with the dependent variable column name, or a two-sided
#'          formula specifying the whole model
#'          (`y ~ x + control1 + control2 | fixed_effect`), exactly as accepted
#'          by [sca()]. When a formula is supplied `data` may be passed
#'          positionally.
#' @param x A string with the focal independent variable column name.
#' @param controls A vector of strings with the control variable column names.
#' @param data A data frame containing the variables.
#' @param weights Optional string naming a weights column in `data`. Weights are
#'                held fixed to their row and never permuted; this assumes the
#'                weights are exogenous to `x`.
#' @param family A string giving the model family, as in [sca()]. Defaults to
#'               `"linear"`.
#' @param link A string giving the link function, as in [sca()]. Defaults to
#'             `NULL`.
#' @param fixed_effects A string (or vector) naming fixed-effect variables, as
#'                      in [sca()]. When supplied, `x` is permuted within levels
#'                      of the first fixed effect.
#' @param n_permutations An integer number of permutations used to build the
#'                       null distribution. Defaults to `500`. The smallest
#'                       attainable p-value is `1 / (n_permutations + 1)`.
#' @param test_stats A character vector of statistics to compute. Any of
#'                   `"median"` (median focal coefficient across
#'                   specifications), `"share_significant"` (share of
#'                   specifications statistically significant in the reference
#'                   direction), `"stouffer"` (Stouffer's combined Z over the
#'                   per-specification p-values), and the descriptive
#'                   `"share_sign"` (share of specifications with the reference
#'                   sign). Defaults to
#'                   `c("median", "share_significant", "stouffer")`.
#' @param direction One of `"two.sided"` (default), `"positive"`, or
#'                  `"negative"`, giving the a-priori predicted direction of the
#'                  effect. The package never chooses the direction from the
#'                  data: with `"two.sided"` the tests are direction-agnostic;
#'                  use `"positive"`/`"negative"` only for a genuinely a-priori
#'                  prediction.
#' @param alpha The significance threshold used by `"share_significant"`.
#'              Defaults to `0.05`.
#' @param sca_data Optional data frame previously returned by [sca()] for the
#'                 *observed* model, used purely to skip recomputing the
#'                 observed curve. The null distribution is always recomputed
#'                 from `data`, so `data` remains required. Its control
#'                 indicator columns must be consistent with `controls`. Ignored
#'                 (with a warning) for `null_type = "freedman_lane"` or
#'                 `"residual_bootstrap"`, where the observed curve is recomputed
#'                 on the forced common sample so it matches the null.
#' @param parallel A boolean indicating whether to parallelise the permutations.
#'                 Defaults to `FALSE`. The inner [sca()] call is always run
#'                 serially to avoid nested parallelism.
#' @param workers An integer number of parallel workers. Defaults to `2`.
#' @param seed Optional integer seed. When supplied, results are reproducible
#'             and identical for the serial and parallel paths (the permutations
#'             are generated up front).
#' @param progress_bar A boolean indicating whether to show a progress bar for
#'                     the permutations. Defaults to `TRUE`.
#'
#' @return An object of class `"sca_test"`: a list with elements `observed` (a
#'         named vector of observed statistics), `null_distribution` (a data
#'         frame with one row per usable permutation and one column per
#'         statistic, plus `n_valid`), `p_values` (a named vector of
#'         permutation p-values), `params` (a list of run metadata), and, when
#'         `keep_curves = TRUE`, `null_curves` (a list with `spec`, a data frame
#'         of each specification's key and observed coefficient, and
#'         `null_coef`, a specifications-by-permutations matrix of the permuted
#'         focal coefficients). When `keep_curves = TRUE` *and* the null is
#'         confound-preserving (`null_type = "freedman_lane"` or
#'         `"residual_bootstrap"`), `null_curves` also gains `fwer`, the
#'         per-specification family-wise-error-rate-adjusted p-values attached
#'         automatically by [sca_minp()] (see there for its structure).
#'
#' @param keep_curves A boolean indicating whether to retain, for every
#'                    specification, the focal coefficient from each permuted
#'                    curve (aligned across permutations by specification, not
#'                    row order). Required by [plot_sca_test_specs()]; increases
#'                    the size of the returned object. Defaults to `FALSE`.
#' @param null_type How the null is generated. One of:
#'                  \describe{
#'                    \item{`"shuffle_x"`}{(default) Simonsohn-Simmons-Nelson
#'                      sharp-null permutation: shuffle the focal variable
#'                      (blocked within the first fixed effect). Appropriate for
#'                      experimental / as-if-randomly-assigned `x`. It breaks
#'                      `x`'s correlation with the controls, so it is
#'                      miscalibrated (anti-conservative) under collinearity --
#'                      the observational case.}
#'                    \item{`"freedman_lane"`}{Freedman-Lane (1983) partial
#'                      permutation. A reduced model `y ~ controls` (the full
#'                      control superset, `x` omitted, fixed effects retained) is
#'                      fit once; its residuals are permuted (blocked within the
#'                      first fixed effect) and added back to its fitted values
#'                      to form a null response, on which the whole curve is
#'                      refit. This preserves `x`'s partial correlation with the
#'                      controls. It tests the sharp null of no partial effect of
#'                      `x` given the superset controls, so under-controlled
#'                      specifications may have non-zero null centres by design.}
#'                    \item{`"residual_bootstrap"`}{The SSN (2020) observational
#'                      scheme: impose the null on the response
#'                      (`y - betahat * x`, with `betahat` from the full
#'                      control-superset specification), then resample rows with
#'                      replacement and refit. A null-imposed case bootstrap
#'                      (nearly equivalent to Flachaire 1999); robust to
#'                      heteroskedasticity but its p-values carry extra
#'                      Monte-Carlo variability.}
#'                  }
#'                  The `"freedman_lane"` and `"residual_bootstrap"` nulls are
#'                  defined for `family = "linear"` only (including fixed
#'                  effects) and force `common_sample = TRUE`.
#' @param common_sample A boolean passed through to [sca()]: fit every
#'                      specification on the rows complete across all model
#'                      variables. Defaults to `FALSE`; forced to `TRUE` for the
#'                      `"freedman_lane"` and `"residual_bootstrap"` nulls.
#'
#' @seealso [plot_sca_test()] to visualise the null distributions of the test
#'   statistics, and [plot_sca_test_specs()] for the per-specification null-band
#'   plot (requires `keep_curves = TRUE`).
#'
#' @references
#' Simonsohn, U., Simmons, J. P., & Nelson, L. D. (2020). Specification curve
#' analysis. \emph{Nature Human Behaviour}, 4, 1208-1214.
#' \doi{10.1038/s41562-020-0912-z}
#'
#' Freedman, D., & Lane, D. (1983). A nonstochastic interpretation of reported
#' significance levels. \emph{Journal of Business & Economic Statistics}, 1(4),
#' 292-298.
#'
#' Winkler, A. M., Ridgway, G. R., Webster, M. A., Smith, S. M., & Nichols, T. E.
#' (2014). Permutation inference for the general linear model. \emph{NeuroImage},
#' 92, 381-397.
#'
#' @importFrom stats median qnorm quantile complete.cases lm fitted residuals coef
#' @export
#'
#' @examples
#' \donttest{
#' # Test whether temperature robustly predicts salinity across specifications.
#' result <- sca_test(y = "Salnty", x = "T_degC",
#'                    controls = c("ChlorA", "O2Sat"),
#'                    data = bottles, n_permutations = 100, progress_bar = FALSE)
#' result
#' }
sca_test <- function(y, x, controls, data, weights = NULL,
                     family = "linear", link = NULL, fixed_effects = NULL,
                     n_permutations = 500,
                     test_stats = c("median", "share_significant", "stouffer"),
                     direction = "two.sided", alpha = 0.05, sca_data = NULL,
                     keep_curves = FALSE,
                     null_type = c("shuffle_x", "freedman_lane",
                                   "residual_bootstrap"),
                     common_sample = FALSE,
                     parallel = FALSE, workers = 2, seed = NULL,
                     progress_bar = TRUE){

  # Formula interface, mirroring sca(): sca_test(y ~ x + c1 | fe, data).
  if(inherits(y, "formula")){
    if(missing(data) && !missing(x) && is.data.frame(x)){
      data <- x
    }
    parsed <- formula_to_args(y)
    y <- parsed$y
    x <- parsed$x
    controls <- parsed$controls
    if(!is.null(parsed$fixed_effects)) fixed_effects <- parsed$fixed_effects
  }

  # Quiet fixest's per-fit singleton/collinearity NOTEs: sca_test refits the
  # whole curve up to n_permutations times, so they would otherwise flood the
  # console (the estimates are unaffected). Restored on exit.
  old_fixest_notes <- getOption("fixest_notes")
  options(fixest_notes = FALSE)
  on.exit(options(fixest_notes = old_fixest_notes), add = TRUE)

  # Validate arguments before doing any expensive work.
  direction <- match.arg(direction, c("two.sided", "positive", "negative"))
  null_type <- match.arg(null_type)

  # sca() (and hence sca_test) supports a single fixed-effects variable; more
  # than one currently fails inside formula_builder. Fail with a clear message.
  if(length(fixed_effects) > 1){
    stop("sca_test() supports a single `fixed_effects` variable; ",
         length(fixed_effects), " were supplied.", call. = FALSE)
  }

  # The design-preserving nulls (Freedman-Lane, residual bootstrap) are defined
  # for linear models on a fixed sample. Restrict to the linear path and force a
  # common sample (they construct a single null response indexed to fixed rows).
  resp_col <- ".sca_test_response"
  if(null_type != "shuffle_x"){
    if(family != "linear"){
      stop("null_type = \"", null_type, "\" is defined for linear models only ",
           "(family = \"linear\"); use null_type = \"shuffle_x\" for glm ",
           "families.", call. = FALSE)
    }
    if(resp_col %in% names(data)){
      stop("`data` already contains a column named \"", resp_col,
           "\", which null_type = \"", null_type, "\" needs to build the null ",
           "response. Please rename it.", call. = FALSE)
    }
    if(!isTRUE(common_sample)){
      message("null_type = \"", null_type, "\" requires a common sample; ",
              "computing the observed curve and the null on the rows complete ",
              "across all model variables.")
      common_sample <- TRUE
    }
    # The null is always built on the forced common sample, so a precomputed
    # `sca_data` (which may have been fit per-spec, on different samples) would
    # put the observed statistics on a different design than the null. Recompute
    # the observed curve instead.
    if(!is.null(sca_data)){
      warning("`sca_data` is ignored for null_type = \"", null_type,
              "\"; the observed curve is recomputed on the common sample so it ",
              "matches the null.", call. = FALSE)
      sca_data <- NULL
    }
  }

  valid_stats <- c("median", "share_significant", "stouffer", "share_sign")
  if(!all(test_stats %in% valid_stats)){
    stop("Invalid `test_stats`: ",
         paste(setdiff(test_stats, valid_stats), collapse = ", "),
         ". Choose from ", paste(valid_stats, collapse = ", "), ".",
         call. = FALSE)
  }
  if(!is.numeric(alpha) || alpha <= 0 || alpha >= 1){
    stop("`alpha` must be between 0 and 1.", call. = FALSE)
  }
  if(!is.numeric(workers) || workers < 1){
    stop("`workers` must be a positive integer.", call. = FALSE)
  }
  if(!is.numeric(n_permutations) || n_permutations < 1){
    stop("`n_permutations` must be a positive integer.", call. = FALSE)
  }
  n_permutations <- as.integer(n_permutations)
  # Warn when the smallest attainable p-value, 1 / (n_permutations + 1), is not
  # strictly below alpha, since nothing can then reach significance under the
  # strict `p < alpha` rule. (This mirrors that rule exactly, unlike a
  # ceiling()-based gate which is off by one when 1/alpha is an integer.)
  if((1 / (n_permutations + 1)) >= alpha){
    warning("`n_permutations` (", n_permutations, ") is small: the smallest ",
            "attainable p-value is ", signif(1 / (n_permutations + 1), 3),
            ", which is not below alpha = ", alpha,
            ", so no statistic can reach significance.", call. = FALSE)
  }

  # Observed curve. When sca_data is supplied it is used only here; computing it
  # via sca() also triggers sca()'s single-coefficient guard on the real data,
  # so a factor/interaction focal variable errors immediately.
  if(is.null(sca_data)){
    observed_curve <- sca(y = y, x = x, controls = controls, data = data,
                          weights = weights, family = family, link = link,
                          fixed_effects = fixed_effects,
                          common_sample = common_sample,
                          parallel = FALSE, progress_bar = FALSE)
  } else {
    if(!all(c("coef", "p") %in% names(sca_data))){
      stop("`sca_data` does not look like sca() output (missing coef/p).",
           call. = FALSE)
    }
    expected <- str_replace(controls, fixed("*"), fixed(":"))
    if(!setequal(sca_control_cols(sca_data), expected)){
      stop("`sca_data` control columns are inconsistent with `controls`.",
           call. = FALSE)
    }
    observed_curve <- sca_data
  }

  observed <- sca_test_statistics(observed_curve, test_stats, direction, alpha)
  n_specs <- nrow(observed_curve)

  # The observed curve must yield a usable value for every requested statistic.
  # This is unreachable through sca() (its single-coefficient guard errors
  # first), but a malformed `sca_data` could otherwise produce silent NA
  # p-values; fail loudly instead.
  bad_stats <- test_stats[is.na(observed[test_stats])]
  if(length(bad_stats) > 0){
    stop("The observed specification curve produced no usable focal estimate ",
         "for: ", paste(bad_stats, collapse = ", "),
         ". p-values cannot be computed (check `sca_data`).", call. = FALSE)
  }

  # Block variable for the permutation: the first fixed effect, if any.
  block <- if(!is.null(fixed_effects)) fixed_effects[1] else NULL
  blocked <- !is.null(block)
  n <- nrow(data)
  block_groups <- if(blocked) unname(split(seq_len(n), data[[block]])) else NULL

  # Master-side null construction for the design-preserving nulls. These are
  # deterministic given the data (no RNG), computed once, and reused across all
  # permutations. Freedman-Lane caches the reduced-model fit; the residual
  # bootstrap injects the null-imposed response into the common-sample data.
  data_cs <- NULL
  reduced <- NULL
  if(null_type != "shuffle_x"){
    data_cs <- sca_test_common_sample(data, y, x, controls, fixed_effects,
                                      weights)
    if(null_type == "freedman_lane"){
      reduced <- tryCatch(
        sca_test_reduced_fit(data_cs, y, controls, fixed_effects, weights,
                             block),
        error = function(e)
          stop("Freedman-Lane reduced model failed to fit: ",
               conditionMessage(e), call. = FALSE))
    } else {
      ynull <- tryCatch(
        sca_test_ynull(data_cs, y, x, controls, fixed_effects, weights),
        error = function(e)
          stop("Residual-bootstrap null model failed to fit: ",
               conditionMessage(e), call. = FALSE))
      data_cs[[resp_col]] <- ynull
    }
  }

  # Pre-generate the permutations in the master under `seed`. This makes the
  # serial and parallel paths produce identical, reproducible results: all RNG
  # happens here and the workers do deterministic work. Setting the seed mutates
  # the caller's global RNG stream, so the previous state is restored on exit
  # (CRAN policy) -- safe because everything after this is deterministic.
  if(!is.null(seed)){
    if(exists(".Random.seed", envir = .GlobalEnv)){
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    } else {
      on.exit(if(exists(".Random.seed", envir = .GlobalEnv))
                rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    }
    set.seed(seed)
  }
  # The index each permutation uses depends on the null mechanism: shuffle_x
  # permutes the focal column over all rows (FE-blocked); freedman_lane permutes
  # the reduced-model residuals over the kept rows (FE-blocked); residual
  # bootstrap resamples common-sample rows with replacement.
  perm_indices <- if(null_type == "shuffle_x"){
    lapply(seq_len(n_permutations),
           function(i) sca_test_perm_index(n, block_groups))
  } else if(null_type == "freedman_lane"){
    lapply(seq_len(n_permutations),
           function(i) sca_test_perm_index(reduced$n_kept, reduced$block_groups))
  } else {
    lapply(seq_len(n_permutations),
           function(i) sample.int(nrow(data_cs), nrow(data_cs), replace = TRUE))
  }

  # One permutation: apply the permuted focal column, re-estimate the curve,
  # and return the per-specification focal coefficient and p-value (or NULL on
  # failure). The summary statistics are computed in the master so the workers
  # only ever call the exported sca(), never the package internals.
  perm_fun <- function(idx){
    if(null_type == "shuffle_x"){
      work <- data
      work[[x]] <- data[[x]][idx]
      resp <- y
    } else if(null_type == "freedman_lane"){
      # y* = reduced fitted + permuted (weighted) residuals, scattered into the
      # kept rows; rows the reduced model dropped (FE singletons) stay NA and are
      # dropped uniformly by every specification.
      estar <- (reduced$sqrt_w * reduced$resid)[idx] / reduced$sqrt_w
      ystar <- rep(NA_real_, nrow(data_cs))
      ystar[reduced$kept] <- reduced$fitted + estar
      work <- data_cs
      work[[resp_col]] <- ystar
      resp <- resp_col
    } else {
      # Resample common-sample rows (with the null-imposed response) jointly.
      work <- data_cs[idx, , drop = FALSE]
      resp <- resp_col
    }
    curve <- tryCatch(
      sca(y = resp, x = x, controls = controls, data = work,
          weights = weights, family = family, link = link,
          fixed_effects = fixed_effects, common_sample = common_sample,
          parallel = FALSE, progress_bar = FALSE),
      error = function(e) NULL)
    if(is.null(curve)) return(NULL)
    out <- list(coef = curve$coef, p = curve$p)
    # Retain the per-specification term sets so the master can align each
    # specification's coefficient across permutations (sca() reorders rows).
    if(keep_curves) out$terms <- curve$terms
    out
  }

  # Run the permutations, parallelising the outer loop if requested.
  if(parallel){
    cl <- makePSOCKcluster(rep("localhost", workers))
    on.exit(stopCluster(cl), add = TRUE)
    clusterEvalQ(cl, library(speccurvieR))
    clusterEvalQ(cl, library(fixest))
    clusterEvalQ(cl, options(fixest_notes = FALSE))
    # NOTE on how the workers actually get their data: perm_fun is a closure, so
    # pblapply()/parLapply() serialise it together with this sca_test() frame --
    # that closure is the real mechanism that carries `data`, `data_cs`,
    # `reduced`, and everything else perm_fun references. This clusterExport is
    # therefore belt-and-suspenders: it makes the worker dependencies explicit
    # and would keep feeding them if perm_fun were ever refactored into a
    # top-level function (which would break the closure capture). It lists every
    # object perm_fun references, including the freedman_lane/residual_bootstrap
    # objects `data_cs`/`reduced`/`resp_col` (NULL for shuffle_x). The
    # serial == parallel reproducibility this contract guarantees is locked in by
    # tests in test-sca_test_nulls.R.
    clusterExport(cl,
                  c("data", "x", "y", "controls", "weights", "family", "link",
                    "fixed_effects", "null_type", "common_sample", "keep_curves",
                    "resp_col", "data_cs", "reduced"),
                  envir = environment())
    null_list <- if(progress_bar){
      pblapply(perm_indices, perm_fun, cl = cl)
    } else {
      parLapply(cl, perm_indices, perm_fun)
    }
  } else {
    null_list <- if(progress_bar){
      pblapply(perm_indices, perm_fun)
    } else {
      lapply(perm_indices, perm_fun)
    }
  }

  # Drop permutations whose entire curve failed; the p-value denominator uses
  # only the usable permutations.
  ok <- !vapply(null_list, is.null, logical(1))
  null_list <- null_list[ok]
  n_used <- length(null_list)
  n_failed <- n_permutations - n_used
  if(n_used == 0){
    stop("All permutations failed to estimate; cannot build a null ",
         "distribution.", call. = FALSE)
  }

  # Compute each permuted curve's statistics in the master.
  null_rows <- lapply(null_list, function(cp){
    s <- sca_test_statistics(data.frame(coef = cp$coef, p = cp$p),
                             test_stats, direction, alpha)
    c(s, n_valid = attr(s, "n_valid"))
  })
  null_df <- as.data.frame(do.call(rbind, null_rows))

  p_values <- sca_test_pvalues(observed, null_df, test_stats, direction)

  # When requested, assemble the per-specification null coefficients: a matrix
  # with one row per specification (keyed by its term set, in the observed
  # curve's order) and one column per usable permutation. Each permuted curve's
  # coefficient is matched to its specification by key, so reordering by sca()
  # cannot misalign specifications; a specification absent from a permuted curve
  # (dropped/failed) is left NA.
  null_curves <- NULL
  if(keep_curves){
    observed_keys <- sca_spec_key(observed_curve$terms, x)
    null_coef <- matrix(NA_real_, nrow = n_specs, ncol = n_used,
                        dimnames = list(observed_keys, NULL))
    for(j in seq_len(n_used)){
      cp <- null_list[[j]]
      null_coef[, j] <- cp$coef[match(observed_keys,
                                      sca_spec_key(cp$terms, x))]
    }
    null_curves <- list(
      spec = data.frame(spec = observed_keys,
                        observed = observed_curve$coef,
                        stringsAsFactors = FALSE),
      null_coef = null_coef)

    # Per-specification FWER (min-P) inference comes for free here: it is pure
    # arithmetic on null_coef, no extra permutations. Attach the single-step
    # adjustment automatically, but ONLY under a confound-preserving null --
    # shuffle_x is severely anti-conservative for per-specification inference,
    # so it is a silent no-op there (the descriptive null band still applies and
    # the shuffle_x null_curves stay byte-identical to the pre-feature object).
    if(sca_fwer_valid_null(null_type)){
      fwer <- sca_fwer_compute(
        null_curves$null_coef, null_curves$spec, direction = direction,
        method = "single_step", alpha = alpha, null_type = null_type)
      # Only attach (and announce) when at least one specification had a usable
      # null distribution to test; otherwise leave the slot absent so the
      # reporting surfaces stay silent rather than printing "0 of 0".
      if(fwer$summary$n_specs_tested > 0){
        null_curves$fwer <- fwer
        s <- fwer$summary
        message("Attached per-specification FWER-adjusted p-values ",
                "(single-step min-P, ", s$n_used, " permutations, alpha = ",
                format(alpha), "): ", s$n_significant, " of ", s$n_specs_tested,
                " specifications significant after correction.")
      }
    }
  }

  out <- list(
    observed = observed[test_stats],
    null_distribution = null_df,
    p_values = p_values)
  # Only present when requested, so the default return is unchanged.
  if(keep_curves) out$null_curves <- null_curves
  out$params <- list(n_permutations = n_permutations, n_used = n_used,
                     n_failed = n_failed, direction = direction, alpha = alpha,
                     x = x, seed = seed, parallel = parallel, workers = workers,
                     family = family, fixed_effects = fixed_effects,
                     n_specs = n_specs, blocked = blocked,
                     test_stats = test_stats, keep_curves = keep_curves,
                     null_type = null_type, common_sample = common_sample,
                     reduced_model = if(null_type != "shuffle_x")
                       "control_superset" else NA_character_,
                     fwer_method = if(keep_curves && sca_fwer_valid_null(null_type))
                       "single_step" else NA_character_)
  structure(out, class = "sca_test")
}

# Internal: human-readable labels for the statistics.
sca_test_labels <- function(){
  c(median = "Median estimate", share_significant = "Share significant",
    stouffer = "Stouffer Z", share_sign = "Share dominant sign")
}

#' Print a specification curve joint-inference test
#'
#' @param x An object of class `"sca_test"` returned by [sca_test()].
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#'
#' @export
print.sca_test <- function(x, ...){
  p <- x$params
  labels <- sca_test_labels()
  floor <- 1 / (p$n_used + 1)
  # The smallest attainable p-value equals the resolution floor (it is not
  # below it), so print exact values; the resolution-floor line below conveys
  # which p-values sit at the floor.
  fmt_p <- function(pv) formatC(pv, format = "f", digits = 4)

  cat("Specification curve joint-inference test",
      "(Simonsohn, Simmons & Nelson 2020)\n\n")
  cat(sprintf("Focal variable:   %s\n", p$x))
  cat(sprintf("Specifications:   %d\n", p$n_specs))
  cat(sprintf("Permutations:     %d used (%d failed)   |  blocked within FE: %s\n",
              p$n_used, p$n_failed, if(p$blocked) "yes" else "no"))
  null_label <- c(shuffle_x = "shuffle x",
                  freedman_lane = "Freedman-Lane (control superset, common sample)",
                  residual_bootstrap = "residual bootstrap (common sample)")
  null_type <- if(is.null(p$null_type)) "shuffle_x" else p$null_type
  cat(sprintf("Null:             %s\n", null_label[[null_type]]))
  cat(sprintf("Direction:        %s   alpha = %s\n\n",
              p$direction, format(p$alpha)))

  cat(sprintf("  %-20s %12s %10s\n", "Statistic", "Observed", "p-value"))
  for(s in p$test_stats){
    cat(sprintf("  %-20s %12s %10s\n", labels[[s]],
                formatC(x$observed[[s]], format = "f", digits = 4),
                fmt_p(x$p_values[[s]])))
  }
  cat(sprintf("\np-values are permutation-based; resolution floor = %s.\n",
              formatC(floor, format = "f", digits = 4)))

  # Per-specification FWER summary, when it was attached (keep_curves + a
  # confound-preserving null). Plain language: no min-P / Westfall-Young jargon.
  fwer <- x$null_curves$fwer
  if(!is.null(fwer)){
    s <- fwer$summary
    if(s$n_significant > 0){
      cat(sprintf(paste0("\nAfter correcting for searching %d specifications, ",
                         "%d remain statistically significant\n",
                         "(smallest corrected p = %s).\n"),
                  s$n_specs_tested, s$n_significant,
                  formatC(s$min_p_adj, format = "f", digits = 4)))
    } else {
      cat(sprintf(paste0("\nAfter correcting for searching %d specifications, ",
                         "no individual specification remains significant\n",
                         "(smallest corrected p = %s).\n"),
                  s$n_specs_tested,
                  formatC(s$min_p_adj, format = "f", digits = 4)))
    }
  }

  # SSN's joint-inference decision rule only applies when the full canonical
  # trio was computed; omit it for any other subset of statistics.
  if(all(c("median", "share_significant", "stouffer") %in% p$test_stats)){
    cat("Interpretation (SSN): conclude a robust effect when the median test AND\n",
        "at least one of {share significant, Stouffer} are significant.\n",
        sep = "")
  }
  invisible(x)
}

#' Plot the null distribution of a joint-inference test
#'
#' @description
#' `plot_sca_test()` visualises the output of [sca_test()]. For each test
#' statistic it draws the permutation null distribution with the observed value
#' marked and the permutation p-value annotated, making it easy to see how
#' extreme the observed specification curve is relative to the sharp null.
#'
#' @param test_result An object of class `"sca_test"` returned by [sca_test()].
#' @param type A string, `"histogram"` (default) or `"density"`, selecting how
#'             the null distribution is drawn.
#' @param title A string used as the plot title. Defaults to `""`.
#'
#' @return A ggplot object with one facet per test statistic.
#'
#' @export
#'
#' @examples
#' \donttest{
#' result <- sca_test(y = "Salnty", x = "T_degC", controls = c("ChlorA", "O2Sat"),
#'                    data = bottles, n_permutations = 100, progress_bar = FALSE)
#' plot_sca_test(result)
#' }
plot_sca_test <- function(test_result, type = "histogram", title = ""){
  if(!inherits(test_result, "sca_test")){
    stop("`test_result` must be an object returned by sca_test().",
         call. = FALSE)
  }
  stats <- test_result$params$test_stats
  labels <- sca_test_labels()

  # Long form of the null distribution (base stack avoids a tidyselect dep).
  nd <- test_result$null_distribution[, stats, drop = FALSE]
  long <- utils::stack(nd)
  names(long) <- c("value", "statistic")
  long$statistic <- factor(long$statistic, levels = stats,
                           labels = labels[stats])

  obs <- data.frame(
    statistic = factor(stats, levels = stats, labels = labels[stats]),
    observed = as.numeric(test_result$observed[stats]),
    p = as.numeric(test_result$p_values[stats]))
  obs$label <- paste0("p = ", formatC(obs$p, format = "f", digits = 3))

  null_fill <- sca_sig_colors()[["p >= .1"]]
  obs_color <- sca_sig_colors()[["p < .005"]]

  ggplot(long, aes(x = value)) +
    {if(tolower(type) == "density")
       geom_density(fill = null_fill, color = "grey20", alpha = .85)
     else
       geom_histogram(fill = null_fill, color = "white", bins = 30)} +
    geom_vline(data = obs, aes(xintercept = observed),
               color = obs_color, linewidth = .8) +
    geom_text(data = obs, aes(x = Inf, y = Inf, label = label),
              hjust = 1.1, vjust = 1.5, size = 3, color = "grey20",
              inherit.aes = FALSE) +
    facet_wrap(~statistic, scales = "free") +
    labs(x = "Statistic under the null", y = "", title = title) +
    theme_sca()
}

#' Plot the specification curve against its per-specification null band
#'
#' @description
#' `plot_sca_test_specs()` draws the observed specification curve (each
#' specification's focal coefficient, ranked) together with a shaded band giving
#' that specification's own null distribution under the permutation test. This
#' is the specification-curve view of the joint-inference test of Simonsohn,
#' Simmons, and Nelson (2020): specifications whose observed estimate falls
#' outside their null band are highlighted, so it is easy to see which parts of
#' the curve are more extreme than chance.
#'
#' When the result carries per-specification family-wise-error-rate-adjusted
#' p-values (i.e. it was computed with `keep_curves = TRUE` and a
#' confound-preserving null, so [sca_minp()] has run -- automatically or
#' explicitly), the points are coloured in three tiers -- within the chance
#' band, beyond it but not significant after correction, and significant after
#' multiple-comparison correction (enlarged) -- so the specifications that
#' survive correction stand out. Otherwise the usual two-tier
#' within/outside-band colouring is used.
#'
#' It requires an [sca_test()] result computed with `keep_curves = TRUE`.
#'
#' @param test_result An object of class `"sca_test"` returned by [sca_test()]
#'                    with `keep_curves = TRUE`.
#' @param level The width of the null band, as a probability. Defaults to `0.95`
#'              (the 2.5th to 97.5th percentile of each specification's null
#'              coefficients).
#' @param title A string used as the plot title. Defaults to `""`.
#'
#' @return A ggplot object.
#'
#' @references
#' Simonsohn, U., Simmons, J. P., & Nelson, L. D. (2020). Specification curve
#' analysis. \emph{Nature Human Behaviour}, 4, 1208-1214.
#' \doi{10.1038/s41562-020-0912-z}
#'
#' @export
#'
#' @examples
#' \donttest{
#' result <- sca_test(y = "Salnty", x = "T_degC", controls = c("ChlorA", "O2Sat"),
#'                    data = bottles, n_permutations = 100, keep_curves = TRUE,
#'                    progress_bar = FALSE)
#' plot_sca_test_specs(result)
#' }
plot_sca_test_specs <- function(test_result, level = 0.95, title = ""){
  if(!inherits(test_result, "sca_test")){
    stop("`test_result` must be an object returned by sca_test().",
         call. = FALSE)
  }
  nc <- test_result$null_curves
  if(is.null(nc)){
    stop("This sca_test object has no per-specification null curves. Re-run ",
         "sca_test() with `keep_curves = TRUE`.", call. = FALSE)
  }
  if(!is.numeric(level) || level <= 0 || level >= 1){
    stop("`level` must be between 0 and 1.", call. = FALSE)
  }

  a <- (1 - level) / 2
  null_coef <- nc$null_coef
  df <- nc$spec

  # Drop specifications that produced no usable null coefficient in any
  # permutation (dropped or non-estimable on every permuted data set); they
  # have no null band to compare against.
  keep <- rowSums(!is.na(null_coef)) > 0
  if(any(!keep)){
    warning(sum(!keep), " specification(s) had no usable null coefficients ",
            "and were omitted from the plot.", call. = FALSE)
    null_coef <- null_coef[keep, , drop = FALSE]
    df <- df[keep, , drop = FALSE]
  }
  if(nrow(df) == 0){
    stop("No specifications have a usable null distribution to plot.",
         call. = FALSE)
  }

  df$lower <- apply(null_coef, 1, quantile, probs = a, na.rm = TRUE)
  df$upper <- apply(null_coef, 1, quantile, probs = 1 - a, na.rm = TRUE)

  # Redundant control sets (e.g. an interaction control plus its components)
  # produce the same fitted model and hence the same key; draw each distinct
  # specification once rather than as overlapping points.
  df <- df[!duplicated(df$spec), ]

  # Rank specifications by observed estimate, as in plot_curve().
  df <- df[order(df$observed), ]
  df$index <- seq_len(nrow(df))
  df$outside <- df$observed < df$lower | df$observed > df$upper

  band_fill <- sca_sig_colors()[["p >= .1"]]

  # When per-specification FWER p-values are attached (keep_curves + a
  # confound-preserving null), upgrade the two-tier "within / outside band"
  # colouring to three tiers so the user can see which specifications survive
  # multiple-comparison correction, not just which fall outside their raw band.
  fwer <- test_result$null_curves$fwer
  if(!is.null(fwer)){
    sig <- fwer$specs$significant_adj[match(df$spec, fwer$specs$spec)]
    sig[is.na(sig)] <- FALSE
    df$band_status <- factor(
      ifelse(sig, "fwer_significant",
             ifelse(df$outside, "outside_uncorrected", "within")),
      levels = c("within", "outside_uncorrected", "fwer_significant"))
    status_cols <- c(within = sca_sig_colors()[["p >= .1"]],
                     outside_uncorrected = sca_sig_colors()[["p < .1"]],
                     fwer_significant = sca_sig_colors()[["p < .005"]])
    status_labs <- c(within = "within chance band",
                     outside_uncorrected = "beyond chance (uncorrected)",
                     fwer_significant = "significant after correction")
    return(
      ggplot(df, aes(x = index)) +
        geom_ribbon(aes(ymin = lower, ymax = upper), fill = band_fill,
                    alpha = .45) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red",
                   linewidth = .5, alpha = .4) +
        # Significant-after-correction points get a size lift and an outline so
        # they read even in greyscale.
        geom_point(aes(y = observed, fill = band_status, size = band_status),
                   shape = 21, colour = "grey20", stroke = .3) +
        # drop = TRUE so a tier with no specifications (e.g. when every
        # specification survives correction) is omitted from the legend rather
        # than shown as a blank key.
        scale_fill_manual(values = status_cols, labels = status_labs,
                          drop = TRUE, name = "Specification") +
        scale_size_manual(values = c(within = 1.3, outside_uncorrected = 1.3,
                                     fwer_significant = 2.6),
                          guide = "none") +
        labs(x = "Specification (ranked by estimate)",
             y = "Focal coefficient", title = title) +
        theme_sca())
  }

  inside_col <- sca_sig_colors()[["p >= .1"]]
  outside_col <- sca_sig_colors()[["p < .005"]]

  ggplot(df, aes(x = index)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = band_fill,
                alpha = .45) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red",
               linewidth = .5, alpha = .4) +
    geom_point(aes(y = observed, color = outside), size = 1.1) +
    scale_color_manual(values = c("FALSE" = inside_col, "TRUE" = outside_col),
                       labels = c("FALSE" = "within null band",
                                  "TRUE" = "outside null band"),
                       drop = TRUE,
                       name = "Observed estimate") +
    labs(x = "Specification (ranked by estimate)",
         y = "Focal coefficient",
         title = title) +
    theme_sca()
}
