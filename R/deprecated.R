# Deprecated camelCase aliases -------------------------------------------------

# Internal: rename any deprecated argument names in `args` to their snake_case
# equivalents. Silent -- the calling alias already warns at the function level.
remap_args <- function(args, mapping){
  nm <- names(args)
  if(is.null(nm)) return(args)
  for(i in seq_along(nm)){
    if(nzchar(nm[i]) && nm[i] %in% names(mapping)){
      nm[i] <- mapping[[nm[i]]]
    }
  }
  names(args) <- nm
  args
}

#' Deprecated functions in speccurvieR
#'
#' @description
#' As of version 0.6.0 the package uses snake_case names throughout. These
#' camelCase aliases are retained for backward compatibility: they emit a
#' deprecation warning and forward to their snake_case replacements (translating
#' any deprecated camelCase argument names automatically). They will be removed
#' in a future release.
#'
#' @param ... Arguments passed on to the replacement function.
#'
#' @return The value returned by the corresponding replacement function.
#'
#' @name speccurvieR-deprecated
#' @keywords internal
NULL

#' @rdname speccurvieR-deprecated
#' @export
plotCurve <- function(...){
  .Deprecated("plot_curve")
  do.call(plot_curve, remap_args(list(...),
    c(showIndex="show_index", plotVars="plot_vars", plotSE="plot_se",
      medianLine="median_line", pointSize="point_size")))
}

#' @rdname speccurvieR-deprecated
#' @export
plotVars <- function(...){
  .Deprecated("plot_vars")
  do.call(plot_vars, remap_args(list(...), c(colorControls="color_controls")))
}

#' @rdname speccurvieR-deprecated
#' @export
plotRMSE <- function(...){
  .Deprecated("plot_rmse")
  do.call(plot_rmse, remap_args(list(...),
    c(showIndex="show_index", plotVars="plot_vars")))
}

#' @rdname speccurvieR-deprecated
#' @export
plotR2Adj <- function(...){
  .Deprecated("plot_r2_adj")
  do.call(plot_r2_adj, remap_args(list(...),
    c(showIndex="show_index", plotVars="plot_vars")))
}

#' @rdname speccurvieR-deprecated
#' @export
plotAIC <- function(...){
  .Deprecated("plot_aic")
  do.call(plot_aic, remap_args(list(...),
    c(showIndex="show_index", plotVars="plot_vars")))
}

#' @rdname speccurvieR-deprecated
#' @export
plotDeviance <- function(...){
  .Deprecated("plot_deviance")
  do.call(plot_deviance, remap_args(list(...),
    c(showIndex="show_index", plotVars="plot_vars")))
}

#' @rdname speccurvieR-deprecated
#' @export
plotControlDistributions <- function(...){
  .Deprecated("plot_control_distributions")
  do.call(plot_control_distributions,
          remap_args(list(...), c(zeroLine="zero_line")))
}

#' @rdname speccurvieR-deprecated
#' @export
plotSE <- function(...){
  .Deprecated("plot_se")
  do.call(plot_se, list(...))
}

#' @rdname speccurvieR-deprecated
#' @export
plotInfluence <- function(...){
  .Deprecated("plot_influence")
  do.call(plot_influence, list(...))
}

#' @rdname speccurvieR-deprecated
#' @export
plotCoefFit <- function(...){
  .Deprecated("plot_coef_fit")
  do.call(plot_coef_fit, list(...))
}

#' @rdname speccurvieR-deprecated
#' @export
plotMultiSE <- function(...){
  .Deprecated("plot_multi_se")
  do.call(plot_multi_se, remap_args(list(...), c(fixedEffects="fixed_effects")))
}

#' @rdname speccurvieR-deprecated
#' @export
controlExtractor <- function(...){
  .Deprecated("control_extractor")
  do.call(control_extractor, list(...))
}

#' @rdname speccurvieR-deprecated
#' @export
unAsIs <- function(...){
  .Deprecated("un_as_is")
  do.call(un_as_is, list(...))
}
