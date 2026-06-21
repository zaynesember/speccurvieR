# broom-style tidiers and generic re-exports-----------------------------------

# Re-export the generics so a user who has attached only speccurvieR can call
# tidy()/glance() directly. They are the same generics broom::tidy()/glance()
# (and modelsummary) dispatch on, so no broom dependency is needed.

#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance

# Internal formatters shared by the tidiers, sca_table(), and sca_report().
sca_fmt_num <- function(x, digits = 3){
  if(!is.finite(x)) return("NA")
  format(signif(x, digits), trim = TRUE, scientific = FALSE)
}
sca_fmt_pct <- function(x, digits = 1){
  if(!is.finite(x)) return("NA")
  paste0(formatC(100 * x, format = "f", digits = digits), "%")
}
sca_fmt_p <- function(p){
  if(is.na(p)) return("NA")
  if(p < 1e-4) return("<0.0001")
  formatC(p, format = "f", digits = 4)
}

#' Tidy and glance methods for speccurvieR objects
#'
#' @description
#' [generics::tidy()] and [generics::glance()] methods for the data frame
#' returned by [sca()] (class `"sca"`) and the object returned by [sca_test()].
#' These are the same generics that `broom::tidy()` / `broom::glance()` and
#' `modelsummary` dispatch on, so the results can feed those tools.
#'
#' Note that the tidied column names follow the broom convention
#' (`estimate` / `std.error` / `p.value`), which differs from the raw column
#' names of [sca()] output (`coef` / `se` / `p`). For an `sca` curve every row
#' shares the same focal `term`, so a tool that keys on `term` (such as
#' `modelsummary`) needs each row distinguished (e.g. by `spec_id`); the
#' per-statistic rows of `tidy.sca_test()` are already distinct. The `controls`
#' column lists the model terms included (so an interaction control `a*b`
#' appears as its expanded terms), which can differ from the user-facing control
#' count in `glance()$n_controls`.
#'
#' @param x An object returned by [sca()] (class `"sca"`) or by [sca_test()].
#' @param ... Ignored.
#' @param row.names,optional Passed through from the [as.data.frame()] generic;
#'   ignored.
#' @param what For `as.data.frame.sca_test()`, which component to return:
#'   `"summary"` (default; observed statistics and p-values), `"null"` (the raw
#'   null distribution), `"params"` (the run metadata), or `"specs"` (the
#'   per-specification family-wise-error-rate-adjusted p-values, available only
#'   when the result was computed with `keep_curves = TRUE` and a
#'   confound-preserving null; see [sca_minp()]).
#'
#' @return A data frame. For `tidy()`, one row per specification (`sca`) or per
#'   test statistic (`sca_test`); for `glance()`, a one-row summary.
#'
#' @name sca_tidiers
#' @seealso [sca_table()] and [sca_report()] for formatted reporting.
#'
#' @examples
#' s <- sca(y = "Salnty", x = "T_degC", controls = c("ChlorA", "O2Sat"),
#'          data = bottles, progress_bar = FALSE, parallel = FALSE)
#' tidy(s)
#' glance(s)
NULL

#' @rdname sca_tidiers
#' @export
tidy.sca <- function(x, ...){
  focal <- sca_focal_var(x)
  controls <- vapply(x$terms,
                     function(t) paste(sort(setdiff(t, c("(Intercept)", focal))),
                                       collapse = " + "),
                     character(1))
  is_glm <- !("RMSE" %in% names(x))
  out <- data.frame(
    term = focal,
    estimate = x$coef,
    std.error = x$se,
    statistic = x$statistic,
    p.value = x$p,
    spec_id = x$index,
    controls = controls,
    n_obs = x$n_obs,
    fit_stat = if(is_glm) x$deviance else x$RMSE,
    fit_stat_name = if(is_glm) "deviance" else "RMSE",
    stringsAsFactors = FALSE)
  out <- out[order(out$spec_id), ]
  rownames(out) <- NULL                       # broom-conventional integer rows
  out
}

#' @rdname sca_tidiers
#' @export
glance.sca <- function(x, ...){
  is_glm <- !("RMSE" %in% names(x))
  n_obs <- x$n_obs
  fam <- attr(x, "family", exact = TRUE)
  if(is.null(fam)) fam <- if(is_glm) "glm" else "linear"
  data.frame(
    focal = sca_focal_var(x),
    n_specs = nrow(x),
    n_controls = length(sca_control_cols(x)),
    family = fam,
    median_estimate = stats::median(x$coef, na.rm = TRUE),
    estimate_min = min(x$coef, na.rm = TRUE),
    estimate_max = max(x$coef, na.rm = TRUE),
    share_significant = mean(x$p < 0.05, na.rm = TRUE),
    share_positive = mean(x$coef > 0, na.rm = TRUE),
    share_negative = mean(x$coef < 0, na.rm = TRUE),
    n_obs_min = min(n_obs, na.rm = TRUE),
    n_obs_max = max(n_obs, na.rm = TRUE),
    common_sample = min(n_obs, na.rm = TRUE) == max(n_obs, na.rm = TRUE),
    stringsAsFactors = FALSE)
}

#' @rdname sca_tidiers
#' @export
tidy.sca_test <- function(x, ...){
  stats <- x$params$test_stats
  data.frame(
    term = stats,
    estimate = as.numeric(x$observed[stats]),
    p.value = as.numeric(x$p_values[stats]),
    direction = x$params$direction,
    null_type = if(is.null(x$params$null_type)) "shuffle_x" else x$params$null_type,
    stringsAsFactors = FALSE)
}

#' @rdname sca_tidiers
#' @export
glance.sca_test <- function(x, ...){
  p <- x$params
  # FWER columns are always emitted (NA when no per-specification adjustment was
  # attached) so the glance schema is stable across results that do and do not
  # carry it.
  fwer <- x$null_curves$fwer
  data.frame(
    n_specs = p$n_specs,
    n_permutations = p$n_permutations,
    n_used = p$n_used,
    n_failed = p$n_failed,
    direction = p$direction,
    alpha = p$alpha,
    null_type = if(is.null(p$null_type)) "shuffle_x" else p$null_type,
    blocked = p$blocked,
    family = p$family,
    focal = p$x,
    common_sample = isTRUE(p$common_sample),
    p_resolution = 1 / (p$n_used + 1),
    n_significant_adj = if(is.null(fwer)) NA_integer_
                        else as.integer(fwer$summary$n_significant),
    min_p_adj = if(is.null(fwer)) NA_real_ else fwer$summary$min_p_adj,
    fwer_method = if(is.null(fwer)) NA_character_ else fwer$summary$method,
    stringsAsFactors = FALSE)
}

#' @rdname sca_tidiers
#' @export
as.data.frame.sca_test <- function(x, row.names = NULL, optional = FALSE, ...,
                                   what = c("summary", "null", "params",
                                            "specs")){
  what <- match.arg(what)
  if(what == "null") return(as.data.frame(x$null_distribution))
  if(what == "params") return(glance.sca_test(x))
  if(what == "specs"){
    fwer <- x$null_curves$fwer
    if(is.null(fwer)){
      stop("This sca_test object has no per-specification FWER p-values. ",
           "Re-run sca_test() with keep_curves = TRUE and a confound-preserving ",
           "null (null_type = \"freedman_lane\" or \"residual_bootstrap\"), or ",
           "call sca_minp() on such a result.", call. = FALSE)
    }
    out <- fwer$specs
    rownames(out) <- NULL
    return(out)
  }
  stats <- x$params$test_stats
  data.frame(
    statistic = stats,
    observed = as.numeric(x$observed[stats]),
    p_value = as.numeric(x$p_values[stats]),
    stringsAsFactors = FALSE)
}
