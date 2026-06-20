# Variance decomposition of the specification curve--------------------------

# Internal: R^2 of regressing y on the intercept plus the columns of Xm indexed
# by `cols` (an integer vector); 0 for the empty set. Uses lm.fit for speed.
sca_r2 <- function(y, Xm, cols, SST){
  if(length(cols) == 0) return(0)
  fit <- stats::lm.fit(cbind(1, Xm[, cols, drop = FALSE]), y)
  1 - sum(fit$residuals^2) / SST
}

# Internal: LMG / Shapley decomposition of R^2 over the columns of Xm. Returns a
# named vector of shares, one per column, that sum to R^2 of the full model
# (the Shapley efficiency axiom), each order-invariant and non-negative.
lmg_r2_shares <- function(y, Xm){
  m <- ncol(Xm)
  SST <- sum((y - mean(y))^2)
  bits <- 2^(0:(m - 1))

  # R^2 for every subset of columns, indexed by bitmask (mask + 1).
  r2 <- numeric(2^m)
  for(mask in seq_len(2^m) - 1){
    r2[mask + 1] <- sca_r2(y, Xm, which(bitwAnd(mask, bits) > 0), SST)
  }

  # Shapley value: average marginal R^2 gain of each column over all subsets.
  shares <- numeric(m)
  for(mask in seq_len(2^m) - 1){
    in_s <- bitwAnd(mask, bits) > 0
    free <- which(!in_s)
    if(length(free) == 0) next               # full set: no marginal gains left
    size_s <- sum(in_s)
    w <- factorial(size_s) * factorial(m - size_s - 1) / factorial(m)
    for(v in free){
      shares[v] <- shares[v] + w * (r2[bitwOr(mask, bits[v]) + 1] - r2[mask + 1])
    }
  }
  names(shares) <- colnames(Xm)
  pmax(shares, 0)
}

# Internal: Type II "drop-one" partial sums of squares for each column, plus the
# full-model residual SS. These do NOT sum to the model SS under a non-orthogonal
# design, so the caller normalises them.
type2_ss <- function(y, Xm){
  m <- ncol(Xm)
  rss_full <- sum(stats::lm.fit(cbind(1, Xm), y)$residuals^2)
  ss <- numeric(m)
  for(v in seq_len(m)){
    rss_v <- sum(stats::lm.fit(cbind(1, Xm[, -v, drop = FALSE]), y)$residuals^2)
    ss[v] <- max(rss_v - rss_full, 0)
  }
  names(ss) <- colnames(Xm)
  list(ss = ss, rss_full = rss_full)
}

#' Decompose the variance of the specification curve by analytic choice
#'
#' @description
#' `sca_variance()` quantifies how much of the variation in the focal
#' coefficient across a specification curve is attributable to each analytic
#' choice -- here, the inclusion or exclusion of each control variable. It is the
#' descriptive complement to [sca_test()]: where the joint-inference test asks
#' whether the curve as a whole is more extreme than chance, the variance
#' decomposition asks which modelling choices *drive* the spread of estimates.
#'
#' The focal coefficient (the `coef` column of [sca()] output) is regressed on
#' the 0/1 control-inclusion indicators, and the variance explained is
#' partitioned across controls. By default this uses the LMG (Shapley-value)
#' decomposition of R-squared, whose shares are order-invariant and sum exactly
#' to the model R-squared; the remaining `1 - R^2` is reported as a `"Residual"`
#' row capturing interactions among choices and unexplained variation.
#'
#' @details
#' **Method.** With `method = "lmg"` (default) each control's share is the
#' Shapley value of R-squared: its average marginal contribution to R-squared
#' over all orderings of the controls. These shares are non-negative,
#' order-invariant, and sum to the full-model R-squared, so the reported
#' percentages (controls plus `"Residual"`) sum to 100. `method = "type2"`
#' instead reports each control's "drop-one" partial sum of squares, normalised
#' to sum to 100; because partial sums of squares do not add up under a
#' non-orthogonal design, those percentages are unique contributions rather than
#' an exact partition. (The regression is purely additive -- each control,
#' including an interaction control, enters as a single 0/1 indicator with no
#' product terms -- so Type II and Type III partial sums of squares coincide;
#' the `"type2"` label is kept only for familiarity.) `"shapley"` is accepted as
#' an alias for `"lmg"`. The LMG cost grows as 2 to the power of the number of
#' varying controls (like `sca()`'s own combinatorial growth) and is capped at
#' 20 controls; for larger sets use `"type2"`, which is linear.
#'
#' **Choice unit.** The decomposition is over *controls* (one binary in/out
#' choice each), read from the indicator columns of [sca()] output. An
#' interaction control such as `"a*b"` is a single choice (its column `"a:b"`).
#' Within one [sca()] call the model family, fixed effects, and focal variable
#' are constant, so only control choices (plus `"Residual"`) appear.
#'
#' **Relationship to specr.** This is analogous to `specr::icc_specs()` /
#' `specr::plot_variance()` and returns a comparable choice-by-percent table and
#' bar chart. It deliberately differs in method: `icc_specs()` fits a
#' random-effects model and reads intraclass-correlation-based variance
#' components, which treats each analytic choice as one grouping factor with many
#' levels. speccurvieR's only enumerated choice is each control's binary status,
#' for which two-level random-effect variance components are unreliable, so an
#' LMG R-squared decomposition is used instead and no `lme4` dependency is added.
#' The output is intentionally *not* labelled "ICC".
#'
#' **Caveats.** The shares are descriptive, not inferential: they summarise an
#' already-estimated set of point estimates and carry no sampling-uncertainty
#' interpretation. A share measures how much a control's in/out status *moves*
#' the focal coefficient, not the control's predictive importance and not the
#' direction of any shift. For glm or fixed-effects curves the coefficient is on
#' the model's native (e.g. log-odds) scale, so shares are not comparable across
#' families. At least two specifications with a varying control are required; with
#' exactly two controls the three specifications are saturated by the two main
#' effects, so the residual is zero. If an interaction control is supplied
#' *together with* its component main-effect controls, the controls are not
#' cleanly separable and the shares should be interpreted with care.
#'
#' @param sca_data A data frame returned by [sca()].
#' @param estimate A string naming the numeric column to decompose. Defaults to
#'                 `"coef"` (the focal coefficient).
#' @param method One of `"lmg"` (default; LMG/Shapley decomposition of
#'               R-squared) or `"type2"` (normalised Type II partial sums of
#'               squares). `"shapley"` is an alias for `"lmg"`.
#' @param residual A boolean indicating whether to append a `"Residual"` row for
#'                 the unexplained (and choice-interaction) variance. Defaults to
#'                 `TRUE`. When `FALSE`, the control percentages are rescaled to
#'                 sum to 100.
#'
#' @return A data frame with one row per control (and a final `"Residual"` row
#'         when `residual = TRUE`), sorted by descending share, with columns:
#'         `choice` (the control name, or `"Residual"`), `variance` (the portion
#'         of the focal coefficient's variance attributed to that choice), and
#'         `percent` (its share as a percentage; the column sums to 100).
#'
#' @seealso [plot_variance()] to visualise the decomposition; [sca_test()] for
#'   the joint-inference test.
#'
#' @export
#'
#' @examples
#' s <- sca(y = "Salnty", x = "T_degC",
#'          controls = c("ChlorA", "O2Sat", "NO2uM"),
#'          data = bottles, progress_bar = FALSE, parallel = FALSE)
#' sca_variance(s)
#' sca_variance(s, method = "type2", residual = FALSE)
sca_variance <- function(sca_data, estimate = "coef",
                         method = c("lmg", "type2"), residual = TRUE){
  if(length(method) == 1 && method == "shapley") method <- "lmg"
  method <- match.arg(method)

  if(!is.data.frame(sca_data)){
    stop("`sca_data` must be a data frame returned by sca().", call. = FALSE)
  }
  if(!estimate %in% names(sca_data)){
    stop("Column `", estimate, "` not found in `sca_data`.", call. = FALSE)
  }
  y <- sca_data[[estimate]]
  if(!is.numeric(y)){
    stop("Column `", estimate, "` must be numeric.", call. = FALSE)
  }

  cols <- sca_control_cols(sca_data)
  if(length(cols) == 0){
    stop("No control-indicator columns found in `sca_data`; nothing to ",
         "decompose.", call. = FALSE)
  }
  # Only genuine 0/1 control indicators are choices; ignore any extra numeric
  # columns the user may have added to the sca() output.
  is_binary <- vapply(cols, function(cc) all(sca_data[[cc]] %in% c(0, 1, NA)),
                      logical(1))
  if(any(!is_binary)){
    warning("Ignoring non-binary column(s) not produced by sca(): ",
            paste(cols[!is_binary], collapse = ", "), call. = FALSE)
    cols <- cols[is_binary]
  }
  if(length(cols) == 0){
    stop("No binary control-indicator columns found in `sca_data`.",
         call. = FALSE)
  }

  # Drop specifications with a missing estimate.
  keep <- !is.na(y)
  if(any(!keep)){
    warning(sum(!keep), " specification(s) with a missing `", estimate,
            "` were dropped.", call. = FALSE)
    y <- y[keep]
    sca_data <- sca_data[keep, , drop = FALSE]
  }
  if(length(y) < 2){
    stop("Variance decomposition needs at least two specifications; with a ",
         "single control sca() returns one model.", call. = FALSE)
  }

  Xm <- as.matrix(sca_data[, cols, drop = FALSE])

  # Drop controls that do not vary across the retained specifications: a
  # constant indicator explains nothing and is collinear with the intercept.
  varies <- apply(Xm, 2, function(z) length(unique(z)) > 1)
  if(!any(varies)){
    stop("No control varies across the supplied specifications; there is ",
         "nothing to decompose.", call. = FALSE)
  }
  Xm <- Xm[, varies, drop = FALSE]

  n <- length(y)
  SST <- sum((y - mean(y))^2)
  # Reject only a coefficient that is constant up to floating-point round-off.
  # Use a centered round-off floor (not the uncentered sum of squares, which
  # would scale with the coefficient's magnitude and could reject a small but
  # genuine spread when the coefficient's location is large).
  if(SST <= n * (.Machine$double.eps * max(abs(y)))^2){
    stop("The `", estimate, "` values are (near-)constant across ",
         "specifications; there is no variance to decompose.", call. = FALSE)
  }

  if(method == "lmg"){
    m <- ncol(Xm)
    if(m > 20L){
      stop("LMG decomposition over ", m, " varying controls is infeasible ",
           "(2^", m, " model fits); use method = \"type2\".", call. = FALSE)
    }
    if(m > 14L){
      message("LMG decomposition over ", m, " controls may be slow (2^", m,
              " fits); method = \"type2\" is much faster.")
    }
    shares <- lmg_r2_shares(y, Xm)
    r2_full <- sum(shares)
    percent_terms <- 100 * shares
    variance_terms <- shares * SST / (n - 1)
    resid_percent <- 100 * (1 - r2_full)
    resid_variance <- (1 - r2_full) * SST / (n - 1)
  } else {
    t2 <- type2_ss(y, Xm)
    denom <- sum(t2$ss) + t2$rss_full
    percent_terms <- 100 * t2$ss / denom
    resid_percent <- 100 * t2$rss_full / denom
    # Scale `variance` to a true partition of var(coef), consistent with the
    # `percent` column and the lmg path (Type II partial SS do not add up to
    # SST under a non-orthogonal design, so the raw SS/(n-1) would not).
    variance_terms <- (percent_terms / 100) * SST / (n - 1)
    resid_variance <- (resid_percent / 100) * SST / (n - 1)
  }

  out <- data.frame(choice = names(percent_terms),
                    variance = as.numeric(variance_terms),
                    percent = as.numeric(percent_terms),
                    stringsAsFactors = FALSE)
  out <- out[order(-out$percent), ]

  if(residual){
    out <- rbind(out, data.frame(choice = "Residual",
                                 variance = resid_variance,
                                 percent = resid_percent))
  } else {
    out$percent <- 100 * out$percent / sum(out$percent)
  }
  rownames(out) <- NULL
  out
}

#' Plot the variance decomposition of a specification curve
#'
#' @description
#' `plot_variance()` visualises [sca_variance()] as a bar chart of the share of
#' the focal coefficient's variance attributable to each control choice (plus a
#' `"Residual"` bar), making it easy to see which modelling choices drive the
#' spread of estimates across the curve.
#'
#' @inheritParams sca_variance
#' @param title A string used as the plot title. Defaults to `""`.
#'
#' @return A ggplot object.
#'
#' @seealso [sca_variance()] for the underlying table.
#'
#' @export
#'
#' @examples
#' s <- sca(y = "Salnty", x = "T_degC",
#'          controls = c("ChlorA", "O2Sat", "NO2uM"),
#'          data = bottles, progress_bar = FALSE, parallel = FALSE)
#' plot_variance(s)
plot_variance <- function(sca_data, estimate = "coef",
                          method = c("lmg", "type2"), residual = TRUE,
                          title = ""){
  vd <- sca_variance(sca_data, estimate = estimate, method = method,
                     residual = residual)

  vd$is_residual <- vd$choice == "Residual"
  # Order bars by share (largest at the top after coord_flip()).
  vd$choice <- factor(vd$choice, levels = vd$choice[order(vd$percent)])

  bar_col <- sca_sig_colors()[["p < .005"]]
  resid_col <- sca_sig_colors()[["p >= .1"]]

  ggplot(vd, aes(x = choice, y = percent, fill = is_residual)) +
    geom_col() +
    geom_text(aes(label = sprintf("%.1f%%", percent)), hjust = -0.1, size = 3) +
    scale_fill_manual(values = c("FALSE" = bar_col, "TRUE" = resid_col),
                      guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    coord_flip() +
    labs(x = "", y = "Percent of variance in focal coefficient",
         title = title) +
    theme_sca()
}
