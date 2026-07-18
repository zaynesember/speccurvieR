# Prose results paragraphs------------------------------------------------------

#' Write a results paragraph for a specification curve analysis
#'
#' @description
#' `sca_report()` returns a one-paragraph, manuscript-ready description of a
#' specification curve. For an [sca_test()] object it states the median
#' estimate, the share of significant specifications, the joint-inference
#' statistics with their permutation p-values, and (when the full set of
#' statistics is available) whether the null of no effect is rejected. For an
#' [sca()] curve it gives a purely descriptive summary.
#'
#' @param x An object returned by [sca()] or [sca_test()].
#' @param digits Number of significant digits for estimates. Defaults to `3`.
#' @param ... Ignored.
#'
#' @return A length-one character string.
#'
#' @seealso [sca_table()] for a tabular summary.
#'
#' @export
#'
#' @examples
#' s <- sca(y = "Salnty", x = "T_degC", controls = c("ChlorA", "O2Sat"),
#'          data = bottles, progress_bar = FALSE, parallel = FALSE)
#' sca_report(s)
#' \donttest{
#' result <- sca_test(y = "Salnty", x = "T_degC", controls = c("ChlorA", "O2Sat"),
#'                    data = bottles, n_permutations = 50, seed = 1,
#'                    progress_bar = FALSE)
#' sca_report(result)
#' }
sca_report <- function(x, digits = 3, ...){
  UseMethod("sca_report")
}

#' @export
sca_report.default <- function(x, ...){
  stop("sca_report() requires an object from sca() or sca_test().",
       call. = FALSE)
}

#' @export
sca_report.sca_test <- function(x, digits = 3, ...){
  p <- x$params
  obs <- x$observed
  pv <- x$p_values
  ts <- p$test_stats
  nt <- if(is.null(p$null_type)) "shuffle_x" else p$null_type
  null_word <- c(shuffle_x = "shuffle-x", freedman_lane = "Freedman-Lane",
                 residual_bootstrap = "residual-bootstrap")[[nt]]
  dir_word <- if(p$direction == "two.sided"){
    "statistically significant"
  } else {
    paste0("significant in the predicted ", p$direction, " direction")
  }

  obs_clause <- if(isTRUE(p$common_sample)) " on a common sample" else ""
  parts <- sprintf(
    "A specification curve analysis estimated the effect of %s across %d specifications%s.",
    p$x, p$n_specs, obs_clause)

  if("median" %in% ts){
    med <- sprintf("The median estimate was %s",
                   sca_fmt_num(obs[["median"]], digits))
    if("share_significant" %in% ts){
      med <- paste0(med, sprintf(" (%s of specifications were %s at the %s level)",
                                 sca_fmt_pct(obs[["share_significant"]]), dir_word,
                                 format(p$alpha)))
    }
    parts <- c(parts, paste0(med, "."))
  }

  # Joint-inference sentence: list the computed statistics with their p-values.
  clauses <- character(0)
  if("median" %in% ts){
    clauses <- c(clauses, sprintf("the median estimate (p = %s)",
                                  sca_fmt_p(pv[["median"]])))
  }
  if("share_significant" %in% ts){
    clauses <- c(clauses,
                 sprintf("the share of significant specifications (p = %s)",
                         sca_fmt_p(pv[["share_significant"]])))
  }
  if("stouffer" %in% ts){
    clauses <- c(clauses, sprintf("Stouffer's combined test (Z = %s, p = %s)",
                                  sca_fmt_num(obs[["stouffer"]], digits),
                                  sca_fmt_p(pv[["stouffer"]])))
  }
  if("share_sign" %in% ts){
    clauses <- c(clauses,
                 sprintf("the share with the dominant sign (p = %s)",
                         sca_fmt_p(pv[["share_sign"]])))
  }

  # Join the statistic clauses, with a serial "and" only for three or more.
  clause_str <- if(length(clauses) > 2){
    paste0(paste(utils::head(clauses, -1), collapse = ", "), ", and ",
           utils::tail(clauses, 1))
  } else {
    paste(clauses, collapse = " and ")
  }

  # Only assert a verdict when the full canonical trio was computed and none of
  # its p-values is missing; otherwise report the statistics without a decision.
  trio <- all(c("median", "share_significant", "stouffer") %in% ts) &&
    !anyNA(pv[c("median", "share_significant", "stouffer")])
  if(trio){
    rejects <- pv[["median"]] < p$alpha &&
      (pv[["share_significant"]] < p$alpha || pv[["stouffer"]] < p$alpha)
    verdict <- if(rejects) "rejected the null of no effect" else
      "did not reject the null of no effect"
    ji <- sprintf("Joint inference via %s permutation tests (%d permutations) %s: %s.",
                  null_word, p$n_used, verdict, clause_str)
  } else {
    ji <- sprintf("Joint inference via %s permutation tests (%d permutations) yielded %s.",
                  null_word, p$n_used, clause_str)
  }

  # Per-specification FWER sentence, when attached (keep_curves + a
  # confound-preserving null). Reports the K-of-N count and the smallest
  # adjusted p-value in plain language.
  fwer_sentence <- NULL
  fwer <- x$null_curves$fwer
  if(!is.null(fwer)){
    s <- fwer$summary
    fwer_sentence <- if(s$n_significant > 0){
      sprintf(paste0("After family-wise error-rate correction for the %d ",
                     "specifications searched, %d remained significant ",
                     "(smallest adjusted p = %s)."),
              s$n_specs_tested, s$n_significant, sca_fmt_p(s$min_p_adj))
    } else {
      sprintf(paste0("After family-wise error-rate correction for the %d ",
                     "specifications searched, no individual specification ",
                     "remained significant (smallest adjusted p = %s)."),
              s$n_specs_tested, sca_fmt_p(s$min_p_adj))
    }
  }

  paste(c(parts, ji, fwer_sentence), collapse = " ")
}

#' @export
sca_report.sca <- function(x, digits = 3, ...){
  g <- glance.sca(x)
  n_clause <- if(g$common_sample){
    sprintf("N = %d", g$n_obs_min)
  } else {
    sprintf("N ranged from %d to %d across specifications", g$n_obs_min,
            g$n_obs_max)
  }
  sprintf(paste0("A specification curve analysis estimated the effect of %s ",
                 "across %d specifications (%s). The median estimate was %s ",
                 "(%s significant at p < .05), with estimates from %s to %s and ",
                 "%s sign agreement. This curve is descriptive; run sca_test() ",
                 "for joint inference."),
          g$focal, g$n_specs, n_clause, sca_fmt_num(g$median_estimate, digits),
          sca_fmt_pct(g$share_significant), sca_fmt_num(g$estimate_min, digits),
          sca_fmt_num(g$estimate_max, digits),
          sca_fmt_pct(max(g$share_positive, g$share_negative)))
}
