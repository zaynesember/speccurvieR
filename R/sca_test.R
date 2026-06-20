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
#'                 indicator columns must be consistent with `controls`.
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
#'         permutation p-values), and `params` (a list of run metadata).
#'
#' @references
#' Simonsohn, U., Simmons, J. P., & Nelson, L. D. (2020). Specification curve
#' analysis. \emph{Nature Human Behaviour}, 4, 1208-1214.
#' \doi{10.1038/s41562-020-0912-z}
#'
#' @seealso [plot_sca_test()] to visualise the result.
#'
#' @importFrom stats median qnorm
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

  # Validate arguments before doing any expensive work.
  direction <- match.arg(direction, c("two.sided", "positive", "negative"))
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
  perm_indices <- lapply(seq_len(n_permutations),
                         function(i) sca_test_perm_index(n, block_groups))

  # One permutation: apply the permuted focal column, re-estimate the curve,
  # and return the per-specification focal coefficient and p-value (or NULL on
  # failure). The summary statistics are computed in the master so the workers
  # only ever call the exported sca(), never the package internals.
  perm_fun <- function(idx){
    perm_data <- data
    perm_data[[x]] <- data[[x]][idx]
    curve <- tryCatch(
      sca(y = y, x = x, controls = controls, data = perm_data,
          weights = weights, family = family, link = link,
          fixed_effects = fixed_effects, parallel = FALSE, progress_bar = FALSE),
      error = function(e) NULL)
    if(is.null(curve)) return(NULL)
    list(coef = curve$coef, p = curve$p)
  }

  # Run the permutations, parallelising the outer loop if requested.
  if(parallel){
    cl <- makePSOCKcluster(rep("localhost", workers))
    on.exit(stopCluster(cl), add = TRUE)
    clusterEvalQ(cl, library(speccurvieR))
    clusterEvalQ(cl, library(fixest))
    clusterExport(cl,
                  c("data", "x", "y", "controls", "weights", "family", "link",
                    "fixed_effects"),
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

  structure(
    list(
      observed = observed[test_stats],
      null_distribution = null_df,
      p_values = p_values,
      params = list(n_permutations = n_permutations, n_used = n_used,
                    n_failed = n_failed, direction = direction, alpha = alpha,
                    x = x, seed = seed, parallel = parallel, workers = workers,
                    family = family, fixed_effects = fixed_effects,
                    n_specs = n_specs, blocked = blocked,
                    test_stats = test_stats)
    ),
    class = "sca_test")
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
