# Vibration of effects----------------------------------------------------------

#' Summarise the vibration of effects across a specification curve
#'
#' @description
#' `sca_voe()` computes the vibration-of-effects summary of Patel, Burford, and
#' Ioannidis (2015) over the specifications estimated by [sca()]: how much the
#' focal estimate and its p-value "vibrate" across the modelling choices. It
#' reports the spread of the estimate between two reference percentiles
#' (by default the 1st and 99th), the relative p-value `RP` (the difference
#' between the same percentiles of -log10(p)), and whether the curve shows a
#' Janus effect -- the estimate pointing in opposite directions at the two
#' percentiles, so the sign of the result depends on the controls chosen.
#'
#' For Cox models (`sca(family = "cox")`) the estimate spread is reported as
#' the relative hazard ratio (the ratio of the hazard ratios at the two
#' percentiles), matching the original paper; for other models it is the
#' difference between the percentile estimates.
#'
#' @param sca_data A data frame returned by [sca()].
#' @param probs A length-two numeric vector of increasing percentiles (strictly
#'              between 0 and 1) at which the spread is measured. Defaults to
#'              `c(0.01, 0.99)`, the 1st and 99th percentiles used by Patel,
#'              Burford, and Ioannidis (2015). With small specification curves
#'              these percentiles are effectively the minimum and maximum.
#'
#' @return A one-row data frame with columns `n_specs`, `estimate_lo` and
#'         `estimate_hi` (the percentile estimates, on the coefficient scale),
#'         `relative_effect` and `relative_effect_type` (the hazard-ratio ratio
#'         for Cox models, otherwise the difference), `RP` (the -log10(p)
#'         spread), `janus` (`TRUE` when the percentile estimates have opposite
#'         signs), `median_estimate`, and `median_p`.
#'
#' @references Patel, C.J., Burford, B., & Ioannidis, J.P.A. (2015). Assessment
#'   of vibration of effects due to model specification can demonstrate the
#'   instability of observational associations. *Journal of Clinical
#'   Epidemiology*, 68(9), 1046-1058.
#'
#' @seealso [plot_voe()] to visualise the vibration as a volcano plot.
#'
#' @export
#'
#' @examples
#' s <- sca(y = "Salnty", x = "T_degC", controls = c("ChlorA", "O2Sat"),
#'          data = bottles, progress_bar = FALSE)
#' sca_voe(s)
sca_voe <- function(sca_data, probs = c(0.01, 0.99)){
  if(!is.data.frame(sca_data) || !all(c("coef", "p") %in% names(sca_data))){
    stop("`sca_data` does not look like sca() output (missing one of the ",
         "`coef` or `p` columns).", call.=FALSE)
  }
  if(!is.numeric(probs) || length(probs) != 2 || any(!is.finite(probs)) ||
     any(probs <= 0) || any(probs >= 1) || probs[1] >= probs[2]){
    stop("`probs` must be two increasing probabilities strictly between 0 ",
         "and 1.", call.=FALSE)
  }

  fam <- attr(sca_data, "family", exact = TRUE)
  is_ratio <- identical(fam, "cox")

  q <- stats::quantile(sca_data$coef, probs = probs, na.rm = TRUE,
                       names = FALSE)
  # p-values of exactly zero would make -log10(p) infinite; clamp at the
  # smallest representable double so RP stays finite.
  logp <- -log10(pmax(sca_data$p, .Machine$double.xmin))
  qp <- stats::quantile(logp, probs = probs, na.rm = TRUE, names = FALSE)

  data.frame(
    n_specs = nrow(sca_data),
    estimate_lo = q[1],
    estimate_hi = q[2],
    relative_effect = if(is_ratio) exp(q[2]) / exp(q[1]) else q[2] - q[1],
    relative_effect_type = if(is_ratio) "ratio (RHR)" else "difference",
    RP = qp[2] - qp[1],
    janus = q[1] < 0 && q[2] > 0,
    median_estimate = stats::median(sca_data$coef, na.rm = TRUE),
    median_p = stats::median(sca_data$p, na.rm = TRUE),
    stringsAsFactors = FALSE)
}

#' Volcano plot of the vibration of effects
#'
#' @description
#' `plot_voe()` draws the vibration-of-effects volcano plot of Patel, Burford,
#' and Ioannidis (2015) for an [sca()] result: each specification is a point
#' at its focal estimate (x) and -log10(p) (y), coloured by significance as in
#' the other speccurvieR plots. A vertical reference line marks the null (a
#' hazard ratio of 1 for Cox models, zero otherwise) and a dashed horizontal
#' line marks p = 0.05. For Cox models the x axis shows hazard ratios on a
#' logarithmic scale. When the curve is large enough, density contours are
#' overlaid to show where specifications concentrate.
#'
#' @inheritParams sca_voe
#' @param title A string with the plot title. Defaults to an empty string.
#' @param point_size A number for the size of the points. Defaults to `NULL`,
#'                   sizing points by the number of specifications as in
#'                   [plot_curve()].
#' @param show_contours A boolean for overlaying 2-D density contours. Defaults
#'                      to `NULL`, which draws them when there are at least 50
#'                      specifications.
#'
#' @return A `ggplot` object.
#'
#' @references Patel, C.J., Burford, B., & Ioannidis, J.P.A. (2015). Assessment
#'   of vibration of effects due to model specification can demonstrate the
#'   instability of observational associations. *Journal of Clinical
#'   Epidemiology*, 68(9), 1046-1058.
#'
#' @seealso [sca_voe()] for the numeric summary.
#'
#' @export
#'
#' @examples
#' s <- sca(y = "Salnty", x = "T_degC", controls = c("ChlorA", "O2Sat"),
#'          data = bottles, progress_bar = FALSE)
#' plot_voe(s)
plot_voe <- function(sca_data, title = "", point_size = NULL,
                     show_contours = NULL){
  if(!is.data.frame(sca_data) ||
     !all(c("coef", "p", "sig.level") %in% names(sca_data))){
    stop("`sca_data` does not look like sca() output (missing one of the ",
         "`coef`, `p`, or `sig.level` columns).", call.=FALSE)
  }

  fam <- attr(sca_data, "family", exact = TRUE)
  is_ratio <- identical(fam, "cox")

  if(is.null(point_size)) point_size <- spec_point_size(sca_data)
  if(is.null(show_contours)) show_contours <- nrow(sca_data) >= 50

  df <- data.frame(
    est = if(is_ratio) exp(sca_data$coef) else sca_data$coef,
    logp = -log10(pmax(sca_data$p, .Machine$double.xmin)),
    sig.level = factor(sca_data$sig.level, levels = names(sca_sig_colors())))

  gg <- ggplot(df, aes(x = est, y = logp)) +
    geom_vline(xintercept = if(is_ratio) 1 else 0, color = "red",
               linetype = "dashed", linewidth = .6, alpha = .4) +
    geom_hline(yintercept = -log10(0.05), color = "grey55",
               linetype = "dotted", linewidth = .6) +
    {if(show_contours) geom_density_2d(color = "grey75", linewidth = .3)} +
    # Points are filled by significance with a thin white outline, matching
    # the package's other plots.
    geom_point(aes(fill = sig.level), shape = 21, color = "white",
               stroke = .4, size = point_size) +
    scale_fill_manual(values = sca_sig_colors(), drop = TRUE) +
    labs(title = title,
         x = if(is_ratio) "Hazard ratio" else "Coefficient",
         y = expression(-log[10](p))) +
    theme_sca() +
    guides(fill = guide_legend(override.aes = list(size = 2)))

  if(is_ratio) gg <- gg + scale_x_log10()

  gg
}
