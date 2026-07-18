# Publication-ready results tables----------------------------------------------

# Human-readable null-mechanism labels, shared with print.sca_test().
sca_null_labels <- function(){
  c(shuffle_x = "shuffle x",
    freedman_lane = "Freedman-Lane (control superset, common sample)",
    residual_bootstrap = "residual bootstrap (common sample)")
}

# Internal: render the assembled label/value data frame in the requested format.
# The default keeps speccurvieR dependency-free (a classed data frame with a
# print method); the richer formats are gated on Suggested packages.
sca_table_render <- function(df, format, note){
  attr(df, "note") <- note
  if(format == "data.frame"){
    class(df) <- c("sca_table", "data.frame")
    return(df)
  }
  if(format %in% c("markdown", "latex")){
    if(!requireNamespace("knitr", quietly = TRUE)){
      stop("format = \"", format, "\" requires the knitr package.",
           call. = FALSE)
    }
    kfmt <- if(format == "markdown") "pipe" else "latex"
    cap <- note
    # knitr escapes table cells but not the caption, so escape LaTeX specials
    # (the notes contain "sca_test()").
    if(kfmt == "latex") cap <- gsub("_", "\\_", cap, fixed = TRUE)
    return(knitr::kable(df, format = kfmt, booktabs = TRUE, caption = cap,
                        col.names = c("", ""), row.names = FALSE))
  }
  if(format == "gt"){
    if(!requireNamespace("gt", quietly = TRUE)){
      stop("format = \"gt\" requires the gt package; install it or use ",
           "format = \"data.frame\".", call. = FALSE)
    }
    return(gt::tab_source_note(gt::gt(df), note))
  }
  if(format == "kableExtra"){
    if(!requireNamespace("kableExtra", quietly = TRUE)){
      stop("format = \"kableExtra\" requires the kableExtra package.",
           call. = FALSE)
    }
    return(kableExtra::kable_styling(
      knitr::kable(df, caption = note, col.names = c("", ""),
                   row.names = FALSE)))
  }
  if(format == "flextable"){
    if(!requireNamespace("flextable", quietly = TRUE)){
      stop("format = \"flextable\" requires the flextable package.",
           call. = FALSE)
    }
    return(flextable::add_footer_lines(flextable::flextable(df), note))
  }
}

#' Render a specification curve result as a publication table
#'
#' @description
#' `sca_table()` turns the output of [sca()] or [sca_test()] into a compact
#' results block suitable for a manuscript: the number of specifications, the
#' median estimate, the share of significant specifications, and (for
#' [sca_test()]) the joint-inference statistics and their permutation p-values.
#'
#' The default `format = "data.frame"` returns a small two-column
#' (`label` / `value`) data frame with a `print` method and adds no
#' dependencies. The other formats render that frame with an optional package:
#' `"markdown"` / `"latex"` use `knitr`, and `"gt"`, `"kableExtra"`, and
#' `"flextable"` use the package of the same name (install it first). A
#' [tidy()]/[glance()]-aware tool such as `modelsummary` can consume the object
#' directly without `sca_table()`.
#'
#' @param x An object returned by [sca()] or [sca_test()].
#' @param format The output format: one of `"data.frame"` (default),
#'   `"markdown"`, `"latex"`, `"gt"`, `"kableExtra"`, or `"flextable"`.
#' @param digits Number of significant digits for estimates. Defaults to `3`.
#' @param ... Ignored.
#'
#' @return For `format = "data.frame"`, a data frame of class `"sca_table"`;
#'   otherwise an object of the corresponding rendering package.
#'
#' @seealso [sca_report()] for a prose summary; [tidy()]/[glance()] for tidy
#'   data frames.
#'
#' @export
#'
#' @examples
#' s <- sca(y = "Salnty", x = "T_degC", controls = c("ChlorA", "O2Sat"),
#'          data = bottles, progress_bar = FALSE, parallel = FALSE)
#' sca_table(s)
#' \donttest{
#' result <- sca_test(y = "Salnty", x = "T_degC", controls = c("ChlorA", "O2Sat"),
#'                    data = bottles, n_permutations = 50, seed = 1,
#'                    progress_bar = FALSE)
#' sca_table(result)
#' }
sca_table <- function(x, format = c("data.frame", "markdown", "latex", "gt",
                                    "kableExtra", "flextable"),
                      digits = 3, ...){
  UseMethod("sca_table")
}

#' @export
sca_table.default <- function(x, ...){
  stop("sca_table() requires an object from sca() or sca_test().", call. = FALSE)
}

#' @export
sca_table.sca_test <- function(x, format = c("data.frame", "markdown", "latex",
                                             "gt", "kableExtra", "flextable"),
                               digits = 3, ...){
  format <- match.arg(format)
  p <- x$params
  labels <- sca_test_labels()
  null_label <- sca_null_labels()
  nt <- if(is.null(p$null_type)) "shuffle_x" else p$null_type

  rows <- list(
    c("Focal variable", p$x),
    c("Specifications", as.character(p$n_specs)),
    c("Observations",
      if(isTRUE(p$common_sample)) "common sample" else "per specification"),
    c("Permutations", as.character(p$n_used)),
    c("Null hypothesis", null_label[[nt]]),
    c("Direction", p$direction),
    c("alpha", format(p$alpha)))
  for(s in p$test_stats){
    val <- x$observed[[s]]
    vstr <- if(s %in% c("share_significant", "share_sign")){
      sca_fmt_pct(val)
    } else {
      sca_fmt_num(val, digits)
    }
    rows <- c(rows, list(c(labels[[s]],
                           paste0(vstr, "  (p = ", sca_fmt_p(x$p_values[[s]]),
                                  ")"))))
  }

  # Per-specification FWER rows, when attached (keep_curves + a
  # confound-preserving null).
  fwer <- x$null_curves$fwer
  if(!is.null(fwer)){
    s <- fwer$summary
    rows <- c(rows,
              list(c("Significant after correction",
                     paste0(s$n_significant, " of ", s$n_specs_tested)),
                   c("Smallest corrected p", sca_fmt_p(s$min_p_adj))))
  }

  df <- data.frame(label = vapply(rows, `[`, character(1), 1),
                   value = vapply(rows, `[`, character(1), 2),
                   stringsAsFactors = FALSE)
  note <- paste0("p-values are permutation-based; resolution floor = ",
                 sca_fmt_p(1 / (p$n_used + 1)), ".")
  sca_table_render(df, format, note)
}

#' @export
sca_table.sca <- function(x, format = c("data.frame", "markdown", "latex", "gt",
                                        "kableExtra", "flextable"),
                          digits = 3, ...){
  format <- match.arg(format)
  g <- glance.sca(x)
  sig_pos <- mean(x$coef > 0 & x$p < 0.05, na.rm = TRUE)
  sig_neg <- mean(x$coef < 0 & x$p < 0.05, na.rm = TRUE)
  obs_str <- if(g$common_sample){
    as.character(g$n_obs_min)
  } else {
    paste0(g$n_obs_min, "-", g$n_obs_max, " (varies by specification)")
  }

  rows <- list(
    c("Focal variable", g$focal),
    c("Specifications", as.character(g$n_specs)),
    c("Observations", obs_str),
    c("Median estimate", sca_fmt_num(g$median_estimate, digits)),
    c("Estimate range",
      paste0("[", sca_fmt_num(g$estimate_min, digits), ", ",
             sca_fmt_num(g$estimate_max, digits), "]")),
    c("Share significant (p<.05)", sca_fmt_pct(g$share_significant)),
    c("Share significant, positive", sca_fmt_pct(sig_pos)),
    c("Share significant, negative", sca_fmt_pct(sig_neg)),
    c("Sign agreement", sca_fmt_pct(max(g$share_positive, g$share_negative))))

  df <- data.frame(label = vapply(rows, `[`, character(1), 1),
                   value = vapply(rows, `[`, character(1), 2),
                   stringsAsFactors = FALSE)
  note <- "Descriptive only; for joint inference run sca_test()."
  sca_table_render(df, format, note)
}

#' @export
print.sca_table <- function(x, ...){
  w <- max(nchar(x$label))
  for(i in seq_len(nrow(x))){
    cat(sprintf("%-*s  %s\n", w, x$label[i], x$value[i]))
  }
  note <- attr(x, "note")
  if(!is.null(note)) cat("\n", note, "\n", sep = "")
  invisible(x)
}
