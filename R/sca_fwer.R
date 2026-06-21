# Per-specification family-wise-error-rate (FWER) inference--------------------

# Internal: the per-specification empirical p-value of a single `value` against
# a specification's own permutation null vector `null_col`. NA null entries are
# dropped; B_j is the number of usable null draws for this specification.
#
# The (1 + count)/(B_j + 1) convention -- the same one sca_test_pvalues() uses
# for the curve-level statistics -- is applied to BOTH the observed estimate
# (here) and, when building the null distribution of the minimum p-value, to
# each permuted estimate scored against its OWN column including itself (see
# sca_fwer_null_pstar()). That shared grid is load-bearing for calibration:
# scoring the observed on (1 + count)/(B + 1) but the permuted points on count/B
# would make the adjusted p-values anti-conservative.
#
# Crucially the comparison is to the specification's OWN null, not to zero.
# Under the Freedman-Lane and residual-bootstrap nulls an under-controlled
# specification's null is centred on its confounded value, not on zero, so a
# zero-referenced (marginal) p-value would be systematically tiny for those
# specifications and would dominate the minimum in every permutation -- which is
# why a naive minP over the marginal regression p-values is wrong here.
#
# `direction` orients the test with a SINGLE non-negative extremeness statistic
# (no two-tail doubling): "positive"/"negative" rank the signed estimate;
# "two.sided" ranks the absolute deviation from the null's own median,
# |value - median(null)|, which is self-centring (so it handles the off-centre
# Freedman-Lane / residual-bootstrap nulls) and, being a single statistic scored
# by the same (1 + count)/(B + 1) rule for the observed and the permuted draws,
# guarantees the adjusted p-value never falls below the raw one. (A min(1,
# 2*min(tails)) two-sided p would put the observed on a finer grid than the
# self-counted permuted draws and could yield p_adj < p_raw.) Returns NA if
# `value` is NA or the column has no usable draws.
sca_fwer_perspec_p <- function(value, null_col, direction){
  nn <- null_col[!is.na(null_col)]
  B <- length(nn)
  if(B == 0L || is.na(value)) return(NA_real_)
  if(direction == "positive"){
    (1 + sum(nn >= value)) / (B + 1)
  } else if(direction == "negative"){
    (1 + sum(nn <= value)) / (B + 1)
  } else {
    m <- stats::median(nn)
    (1 + sum(abs(nn - m) >= abs(value - m))) / (B + 1)
  }
}

# Internal: the vector of own-null empirical p-values for every permutation draw
# in a single specification's null column, aligned to the column's positions
# (NA where the column entry is NA). Each draw is scored against its own column
# INCLUDING itself, identically to how sca_fwer_perspec_p() scores the observed
# estimate -- this is the rank-based vectorisation of
# sca_fwer_perspec_p(null_col[b], null_col, direction) over b. Keeping the
# observed and the permuted draws on this one shared grid is what makes the
# minimum-p null exchangeable with the observed statistic; do not "optimise" it
# to exclude-self or drop the +1 (that breaks calibration, detectably at small
# B).
sca_fwer_null_pstar <- function(null_col, direction){
  out <- rep(NA_real_, length(null_col))
  ok <- !is.na(null_col)
  nn <- null_col[ok]
  B <- length(nn)
  if(B == 0L) return(out)
  # rank(., "min")[i] = 1 + #{strictly less}, so #{>= x_i} = B - rank_min + 1;
  # rank(., "max")[i] = #{<= x_i}. Both count the focal draw itself, matching
  # sca_fwer_perspec_p()'s scoring of a permuted draw against its own column.
  if(direction == "positive"){
    n_ge <- B - rank(nn, ties.method = "min") + 1
    out[ok] <- (1 + n_ge) / (B + 1)
  } else if(direction == "negative"){
    n_le <- rank(nn, ties.method = "max")
    out[ok] <- (1 + n_le) / (B + 1)
  } else {
    # Absolute deviation from the null median, ranked descending.
    a <- abs(nn - stats::median(nn))
    n_ge <- B - rank(a, ties.method = "min") + 1
    out[ok] <- (1 + n_ge) / (B + 1)
  }
  out
}

# Internal: Westfall-Young free step-down adjusted p-values. `p_raw` is the
# vector of per-specification observed own-null p-values (distinct specs);
# `pstar` is the distinct-specs x permutations matrix of own-null p-values for
# the permuted draws (NA where a draw was non-estimable). Specifications are
# ordered by p_raw (most to least significant); the reference for the k-th
# ordered specification is the minimum permuted p-value over the not-yet-
# rejected (k-th and less significant) specifications, so the multiplicity
# penalty shrinks as we step down the curve. Monotonicity is then enforced by a
# cumulative maximum, and every value is floored at 1/(n_eff + 1). This controls
# strong FWER under subset pivotality (exact for Freedman-Lane via the common
# control-superset conditioning; approximate for the residual bootstrap and
# under partial alternatives), so it is opt-in rather than the default.
sca_fwer_stepdown <- function(p_raw, pstar, n_eff){
  m <- length(p_raw)
  floor_p <- 1 / (n_eff + 1)
  o <- order(p_raw)                       # most significant first
  Pord <- pstar[o, , drop = FALSE]
  Pord[is.na(Pord)] <- Inf                # NA draws never lower the running min
  # Suffix (reverse cumulative) minimum down the ordered specifications, per
  # permutation column: Q[k, b] = min over ordered specs k..m of Pord[., b].
  Q <- apply(Pord, 2, function(col) rev(cummin(rev(col))))
  if(is.null(dim(Q))) Q <- matrix(Q, nrow = m)   # m == 1 guard
  pr_ord <- p_raw[o]
  padj_ord <- vapply(seq_len(m), function(k){
    qk <- Q[k, ]
    qk <- qk[is.finite(qk)]
    (1 + sum(qk <= pr_ord[k])) / (n_eff + 1)
  }, numeric(1))
  padj_ord <- cummax(padj_ord)            # enforce monotone non-decreasing
  padj_ord <- pmax(padj_ord, floor_p)
  out <- numeric(m)
  out[o] <- padj_ord
  out
}

# Internal: the shared FWER engine behind both the auto-augment path in
# sca_test() and the user-facing sca_minp(). It is pure arithmetic on the
# already-stored per-specification null-coefficient matrix (no RNG, no
# refitting), so the serial and parallel sca_test() paths produce identical
# results for free.
#
# `null_coef` is the n_rows x n_used matrix from keep_curves; `spec` is the
# matching data frame (`spec` key + `observed` coefficient), one row per sca()
# row. Returns a list with `specs` (one row per ORIGINAL sca() row, so it stays
# aligned with `spec`/`null_coef`) and `summary` (the plain-language tallies).
sca_fwer_compute <- function(null_coef, spec, direction, method, alpha,
                             null_type = NA_character_){
  keys_all <- spec$spec
  obs_all <- spec$observed
  B <- ncol(null_coef)

  res <- data.frame(spec = keys_all, observed = obs_all,
                    p_raw = NA_real_, p_adj = NA_real_,
                    significant_adj = NA, stringsAsFactors = FALSE)
  empty_summary <- function(){
    list(method = method, alpha = alpha,
         fwer_type = if(method == "step_down") "strong" else "weak",
         statistic = "per_spec_min_p", n_specs_tested = 0L,
         n_significant = 0L, min_p_adj = NA_real_, null_type = null_type,
         n_used = B, direction = direction)
  }
  if(B == 0L) return(list(specs = res, summary = empty_summary()))

  # Distinct, usable specifications. Drop rows that are NA in every permutation
  # (no null band to test against), then collapse redundant duplicate-key rows
  # (interaction controls plus their components fit the identical model and
  # share a key) so each distinct fitted model contributes ONCE to the minimum
  # -- counting it several times would inflate the multiplicity penalty and the
  # reported N. Mirrors plot_sca_test_specs()'s reduction so the plot and the
  # inference agree on which specifications are tested. A specification with no
  # usable observed estimate cannot be tested either.
  usable <- rowSums(!is.na(null_coef)) > 0 & !is.na(obs_all)
  d_idx <- which(usable)
  if(length(d_idx) == 0L) return(list(specs = res, summary = empty_summary()))
  d_idx <- d_idx[!duplicated(keys_all[d_idx])]
  d_keys <- keys_all[d_idx]
  d_obs <- obs_all[d_idx]
  d_coef <- null_coef[d_idx, , drop = FALSE]
  m <- length(d_idx)

  # Observed per-specification own-null p-values.
  p_raw_d <- vapply(seq_len(m),
                    function(j) sca_fwer_perspec_p(d_obs[j], d_coef[j, ],
                                                   direction),
                    numeric(1))

  if(method == "step_down"){
    # Materialise the distinct-specs x permutations own-null p-value matrix.
    # Build it row by row (not t(vapply())) so the orientation is m x B even
    # when there is a single permutation (B == 1), where vapply()/t() would
    # otherwise collapse to a 1 x m matrix and misalign the step-down.
    pstar <- matrix(NA_real_, nrow = m, ncol = B)
    for(j in seq_len(m)){
      pstar[j, ] <- sca_fwer_null_pstar(d_coef[j, ], direction)
    }
    minP <- apply(pstar, 2, function(col){
      col <- col[!is.na(col)]
      if(length(col) == 0L) NA_real_ else min(col)
    })
    n_eff <- sum(!is.na(minP))
    p_adj_d <- sca_fwer_stepdown(p_raw_d, pstar, n_eff)
  } else {
    # Single-step: stream the per-permutation minimum so peak memory stays at
    # the already-resident null_coef plus an O(n_used) vector.
    minP <- rep(Inf, B)
    any_ok <- rep(FALSE, B)
    for(j in seq_len(m)){
      ps <- sca_fwer_null_pstar(d_coef[j, ], direction)
      ok <- !is.na(ps)
      any_ok <- any_ok | ok
      minP[ok] <- pmin(minP[ok], ps[ok])
    }
    minP[!any_ok] <- NA_real_
    minP_ok <- minP[!is.na(minP)]
    n_eff <- length(minP_ok)
    floor_p <- 1 / (n_eff + 1)
    p_adj_d <- vapply(p_raw_d, function(pr){
      if(is.na(pr)) return(NA_real_)
      max(floor_p, (1 + sum(minP_ok <= pr)) / (n_eff + 1))
    }, numeric(1))
  }

  # Multiple-comparison correction can never make a p-value smaller than its
  # uncorrected value. The shared-grid construction guarantees p_adj >= p_raw
  # when every tested specification has the same number of usable null draws;
  # this clamp also enforces it in the ragged case (a specification estimable in
  # the observed fit but non-estimable in only some permuted data sets has a
  # shorter, coarser own-null grid than its peers, against which the pooled
  # minimum-p null could otherwise yield a slightly anti-conservative p_adj).
  p_adj_d <- pmax(p_adj_d, p_raw_d)

  # Back-map the distinct results to every original row by key, so duplicate-key
  # rows receive identical values and the table stays aligned with `spec`.
  map <- match(keys_all, d_keys)
  res$p_raw <- p_raw_d[map]
  res$p_adj <- p_adj_d[map]
  res$significant_adj <- res$p_adj < alpha

  n_sig <- sum(p_adj_d < alpha, na.rm = TRUE)
  summary <- list(
    method = method, alpha = alpha,
    fwer_type = if(method == "step_down") "strong" else "weak",
    statistic = "per_spec_min_p",
    n_specs_tested = m, n_significant = n_sig,
    min_p_adj = if(m > 0L) min(p_adj_d, na.rm = TRUE) else NA_real_,
    null_type = null_type, n_used = B, direction = direction)
  list(specs = res, summary = summary)
}

# Internal: TRUE when a null_type supports valid per-specification FWER
# inference. shuffle_x breaks the focal variable's correlation with the controls
# and is severely anti-conservative for per-specification inference, so it is
# excluded; the design-preserving nulls (which also force a common sample) are
# the valid ones.
sca_fwer_valid_null <- function(null_type){
  isTRUE(null_type %in% c("freedman_lane", "residual_bootstrap"))
}

#' Family-wise-error-rate-adjusted p-values for each specification
#'
#' @description
#' `sca_minp()` adds, to an existing [sca_test()] result, a multiple-comparison-
#' corrected p-value for *every* specification in the curve. It answers the
#' question the per-specification null band only shows descriptively: *which*
#' individual specifications are more extreme than chance once you account for
#' having searched all of them?
#'
#' It reuses the permutation null [sca_test()] already computed (so it never
#' re-runs the permutations) and requires that result to have been produced with
#' `keep_curves = TRUE` and a confound-preserving null
#' (`null_type = "freedman_lane"` or `"residual_bootstrap"`). When those
#' prerequisites are met, [sca_test()] already attaches the single-step
#' adjustment automatically; call `sca_minp()` to recompute it with a different
#' method or significance threshold.
#'
#' @details
#' Reporting a specification as "significant" because its ordinary p-value is
#' below 0.05 is misleading when the curve contains many correlated
#' specifications: searching dozens of them all but guarantees some will look
#' significant by chance. `sca_minp()` corrects for this with the min-P /
#' max-statistic permutation method of Westfall and Young (1993). For each
#' permutation it records the most extreme specification anywhere in the curve,
#' building the null distribution of "the best result a search could turn up by
#' chance", and compares each observed specification to that distribution.
#'
#' Each specification is compared to its *own* permutation null rather than to
#' zero. This matters under the confound-preserving nulls: an under-controlled
#' specification's null distribution is centred on its confounded value, not on
#' zero, so a specification is flagged only when its estimate is more extreme
#' than that specification could itself produce under the null. Consequently a
#' flagged under-controlled specification indicates an association beyond the
#' superset-conditional sharp null -- it is **not** a causal effect of size
#' (estimate minus null centre), because omitting a confounder re-routes part of
#' the focal variable's coefficient.
#'
#' `method = "single_step"` (the default, and what [sca_test()] attaches
#' automatically) controls the family-wise error rate in the weak sense, under
#' the global null of no effect in any specification, with no further
#' assumptions. `method = "step_down"` applies Westfall and Young's free
#' step-down refinement, which is more powerful and controls the family-wise
#' error rate in the strong sense under subset pivotality (exact for
#' `"freedman_lane"` via the common control-superset conditioning; approximate
#' for `"residual_bootstrap"` and under partial alternatives). Adjusted
#' p-values have a resolution floor of `1 / (n_eff + 1)`, where `n_eff` is the
#' number of permutations that estimated at least one specification (equal to
#' `n_used` unless some permutations were wholly non-estimable).
#'
#' @param result An object of class `"sca_test"` returned by [sca_test()] with
#'               `keep_curves = TRUE` and a confound-preserving `null_type`.
#' @param method One of `"single_step"` (default) or `"step_down"`; see Details.
#' @param alpha The significance threshold for `significant_adj`. Defaults to
#'              `NULL`, which inherits the `alpha` used by [sca_test()].
#'
#' @return The input `result`, augmented with a `fwer` element under
#'         `null_curves`: a list with `specs` (a data frame with one row per
#'         specification giving its term-set key `spec`, `observed` coefficient,
#'         uncorrected own-null p-value `p_raw`, FWER-adjusted p-value `p_adj`,
#'         and logical `significant_adj`) and `summary` (the method, alpha,
#'         whether weak or strong FWER is controlled, the number of distinct
#'         specifications tested, the number significant after correction, and
#'         the smallest adjusted p-value).
#'
#' @seealso [sca_test()], which attaches the single-step adjustment
#'   automatically; [plot_sca_test_specs()] to see which specifications survive
#'   correction.
#'
#' @references
#' Westfall, P. H., & Young, S. S. (1993). \emph{Resampling-based multiple
#' testing: Examples and methods for p-value adjustment}. Wiley.
#'
#' @export
#'
#' @examples
#' \donttest{
#' result <- sca_test(y = "Salnty", x = "T_degC",
#'                    controls = c("ChlorA", "O2Sat", "NO2uM"), data = bottles,
#'                    null_type = "freedman_lane", n_permutations = 199,
#'                    keep_curves = TRUE, seed = 1, progress_bar = FALSE)
#' # The single-step adjustment is already attached; recompute it step-down.
#' result <- sca_minp(result, method = "step_down")
#' as.data.frame(result, what = "specs")
#' }
sca_minp <- function(result, method = c("single_step", "step_down"),
                     alpha = NULL){
  if(!inherits(result, "sca_test")){
    stop("`result` must be an object returned by sca_test().", call. = FALSE)
  }
  method <- match.arg(method)
  nc <- result$null_curves
  if(is.null(nc) || is.null(nc$null_coef)){
    stop("Per-specification FWER needs the retained null curves; re-run ",
         "sca_test(..., keep_curves = TRUE).", call. = FALSE)
  }
  null_type <- if(is.null(result$params$null_type)) "shuffle_x"
               else result$params$null_type
  if(!sca_fwer_valid_null(null_type)){
    stop("Per-specification FWER p-values require a confound-preserving null ",
         "(null_type = \"freedman_lane\" or \"residual_bootstrap\"); ",
         "null_type = \"", null_type, "\" breaks the focal variable's ",
         "correlation with the controls and is severely anti-conservative for ",
         "per-specification inference.", call. = FALSE)
  }
  if(is.null(alpha)) alpha <- result$params$alpha
  if(!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1){
    stop("`alpha` must be a single number between 0 and 1.", call. = FALSE)
  }
  fwer <- sca_fwer_compute(
    nc$null_coef, nc$spec, direction = result$params$direction,
    method = method, alpha = alpha, null_type = null_type)
  if(fwer$summary$n_specs_tested == 0L){
    stop("No specification has a usable null distribution to adjust ",
         "(every specification was non-estimable across the permutations).",
         call. = FALSE)
  }
  result$null_curves$fwer <- fwer
  result$params$fwer_method <- method
  result
}
