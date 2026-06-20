# User-facing functions---------------------------------------------------------

#' Perform specification curve analysis
#'
#' @description
#' sca() is the workhorse function of the package--this estimates models with every
#' possible combination of the controls supplied and returns a data frame
#' where each row contains the pertinent information and parameters for a
#' given model by default. This data frame can then be input to plot_curve()
#' or any other plotting function in the package. Alternatively, if
#' `return_formulae = TRUE`, it returns a list of formula objects with every
#' possible combination of controls.
#'
#' @param y A string containing the column name of the dependent variable in
#'          data. Alternatively, a two-sided formula specifying the whole model,
#'          e.g. `y ~ x + control1 + control2 | fixedEffect`. When a formula is
#'          supplied the first right-hand-side term is taken as the independent
#'          variable `x`, the remaining terms as `controls`, and anything after
#'          `|` as `fixed_effects`; the `x`, `controls`, and `fixed_effects`
#'          arguments are then taken from the formula. `data` may be passed
#'          positionally in this case, e.g. `sca(y ~ x + z, data)`.
#' @param x A string containing the column name of the independent variable in
#'          data.
#' @param controls A vector of strings containing the column names of the
#'                 control variables in data.
#' @param data A dataframe containing y, x, controls, and (optionally) the
#'             variables to be used for fixed effects or clustering.
#' @param weights Optional string with the column name in `data` that contains
#'                weights.
#' @param family A string indicating the family of models to be used. Defaults
#'               to "linear" for OLS regression but supports all families
#'               supported by `glm()`.
#' @param link A string specifying the link function to be used for the model.
#'             Defaults to `NULL` for OLS regression using `lm()` or
#'             `fixest::feols()` depending on whether fixed effects are supplied.
#'             Supports all link functions supported by the family parameter of
#'             `glm()`.
#' @param fixed_effects A string containing the column name of the variable
#'                     in data desired for fixed effects. Defaults to NULL in
#'                     which case no fixed effects are included.
#' @param return_formulae A boolean. When `TRUE` a list of model formula objects
#'                       is returned but the models are not estimated. Defaults
#'                       to `FALSE` in which case a dataframe of model results
#'                       is returned.
#' @param progress_bar A boolean indicating whether the user wants a progress bar
#'                    for model estimation. Defaults to `TRUE`.
#' @param parallel A boolean indicating whether to parallelize model estimation.
#'                 Parallelization only offers a speed advantage when a large
#'                 (> 1000) number of models is being estimated. Defaults to
#'                 `FALSE`.
#' @param workers An integer indicating the number of workers to use for
#'                parallelization. Defaults to 2.
#' @param ... Deprecated camelCase arguments (`fixedEffects`, `returnFormulae`,
#'            `progressBar`); use the snake_case equivalents instead.
#'
#' @return When `return_formulae` is `FALSE`, a dataframe where each row contains
#'         the independent variable coefficient estimate, standard error,
#'         test statistic, p-value, model specification, and measures of model
#'         fit.
#'
#' @export
#'
#' @examples
#' sca(y = "Salnty", x = "T_degC", controls = c("ChlorA", "O2Sat"),
#'     data = bottles, progress_bar = TRUE, parallel = FALSE);
#' # Equivalent call using the formula interface:
#' sca(Salnty ~ T_degC + ChlorA + O2Sat, data = bottles, progress_bar = FALSE);
#' # Formula interface with an interaction control and fixed effects:
#' sca(Salnty ~ T_degC + ChlorA + ChlorA*O2Sat | Sta_ID, data = bottles,
#'     progress_bar = FALSE);
#' \donttest{
#' sca(y = "Salnty", x = "T_degC", controls = c("ChlorA*NO3uM", "O2Sat*NO3uM"),
#'     data = bottles, progress_bar = TRUE, parallel = TRUE, workers = 2);
#' }
#' sca(y = "Salnty", x = "T_degC", controls = c("ChlorA", "O2Sat*NO3uM"),
#'     data = bottles, progress_bar = TRUE, parallel = FALSE,
#'     return_formulae = TRUE);
sca <- function(y, x, controls, data, weights=NULL,
                family="linear", link=NULL,
                fixed_effects=NULL, return_formulae=FALSE,
                progress_bar=TRUE, parallel=FALSE, workers=2, ...){

  # Backward compatibility: translate deprecated camelCase argument names.
  .dots <- list(...)
  if("fixedEffects" %in% names(.dots)){
    .Deprecated(msg="`fixedEffects` is deprecated; use `fixed_effects`.")
    fixed_effects <- .dots$fixedEffects
  }
  if("returnFormulae" %in% names(.dots)){
    .Deprecated(msg="`returnFormulae` is deprecated; use `return_formulae`.")
    return_formulae <- .dots$returnFormulae
  }
  if("progressBar" %in% names(.dots)){
    .Deprecated(msg="`progressBar` is deprecated; use `progress_bar`.")
    progress_bar <- .dots$progressBar
  }

  # Formula interface: sca(y ~ x + control1 + control2 | fe, data = ...).
  # The response is y, the first right-hand-side term is the focal independent
  # variable, the remaining terms are controls, and anything after `|` is
  # treated as fixed effects.
  if(inherits(y, "formula")){
    # Allow data to be passed positionally, i.e. sca(y ~ ..., data).
    if(missing(data) && !missing(x) && is.data.frame(x)){
      data <- x
    }
    else if(!missing(x) || !missing(controls)){
      warning("`x` and `controls` are ignored when `y` is a formula.")
    }
    parsed <- formula_to_args(y)
    y <- parsed$y
    x <- parsed$x
    controls <- parsed$controls
    if(!is.null(parsed$fixed_effects)) fixed_effects <- parsed$fixed_effects
  }

  # Treat the common alias "gaussian" as ordinary least squares.
  if(family=="gaussian") family <- "linear"

  if(family!="linear" & !is.null(fixed_effects))
  {
    warning(paste0("Fixed effects unsupported for models other than OLS ",
                   "regression. Ignoring fixed effects."))
    # Actually drop the fixed effects so downstream estimation/extraction uses
    # the glm path, as the warning promises.
    fixed_effects <- NULL
  }

  # Build the glm family object, defaulting to the family's canonical link when
  # `link` is NULL, with a clear error for an unrecognised family.
  if(family!="linear"){
    fam_fun <- tryCatch(match.fun(family),
                        error=function(e)
                          stop("'", family,
                               "' is not a recognised model family.",
                               call.=FALSE))
    fam_obj <- if(is.null(link)) fam_fun() else fam_fun(link=link)
  }

  # Just generate the formulae and return if desired
  if(return_formulae){
    if(!is.null(fixed_effects)){
      return(formula_builder(y=y, x=x, controls=controls,
                             fixed_effects=fixed_effects))
    }
    else{
      return(formula_builder(y=y, x=x, controls=controls))
    }
  }

  # Validate that the requested columns exist in `data` (interaction syntax in
  # x/controls is split so each underlying variable is checked).
  vars <- unique(trimws(unlist(strsplit(c(y, x, controls), "[*:]"))))
  check_columns(data, vars, "Variable(s)")
  if(!is.null(fixed_effects)){
    check_columns(data, fixed_effects, "Fixed-effects variable(s)")
  }
  if(!is.null(weights)) check_columns(data, weights, "Weights variable")

  # Build the model formulae (with or without fixed effects)
  if(is.null(fixed_effects)){
    formulae <- formula_builder(y=y, x=x, controls=controls)
  }
  else{
    formulae <- formula_builder(y=y, x=x, controls=controls,
                                fixed_effects=fixed_effects)
  }

  # Estimator for a single specification. Dispatches on fixed effects, family,
  # and weights and returns a model summary. Defined as a closure so that, when
  # estimation is parallelised, it carries `data`, `weights`, `family`,
  # `fam_obj`, and `fixed_effects` to the workers along with the function.
  estimate_one <- function(f){
    if(!is.null(fixed_effects)){
      if(is.null(weights)){
        summary(feols(f, data=data))
      }
      else{
        summary(feols(f, data=data, weights=data[[weights]]))
      }
    }
    else if(family=="linear"){
      if(is.null(weights)){
        summary(lm(f, data=data))
      }
      else{
        environment(f) <- environment()
        summary(lm(f, data=data, weights=get(weights)))
      }
    }
    else{
      if(is.null(weights)){
        summary(glm(f, data=data, family=fam_obj))
      }
      else{
        environment(f) <- environment()
        summary(glm(f, data=data, weights=get(weights), family=fam_obj))
      }
    }
  }

  # Estimate every specification, parallelising and/or showing a progress bar
  # as requested.
  if(parallel){
    cl <- makePSOCKcluster(rep("localhost", workers))
    clusterEvalQ(cl, library(fixest))
    clusterExport(cl, "data", envir=environment())

    if(progress_bar){
      print.noquote(paste("Estimating", length(formulae),
                          "models in parallel with", workers, "workers"))
      models <- pblapply(formulae, estimate_one, cl=cl)
    }
    else{
      models <- parLapply(cl, formulae, estimate_one)
    }

    stopCluster(cl)
  }
  else{
    if(progress_bar){
      print.noquote(paste("Estimating", length(formulae), "models"))
      models <- pblapply(formulae, estimate_one)
    }
    else{
      models <- lapply(formulae, estimate_one)
    }
  }

  # OLS models
  if(family=="linear"){

    # No fixed effects
    if(is.null(fixed_effects)){

      # Get each value of interest across models
      coef <- lapply(X=models, function(x2) x2$coefficients[x,1])
      se <- lapply(X=models, function(x2) x2$coefficients[x,2])
      statistic <- lapply(X=models, function(x2) x2$coefficients[x,3])
      p <- lapply(X=models, function(x2) x2$coefficients[x,4])
      terms <- lapply(X=models, FUN=function(x2) row.names(x2$coefficients))
      RMSE <- lapply(X=models, FUN=function(x2) sqrt(mean(x2$residuals^2)))
      adjR <- lapply(X=models, function(x2) x2$adj.r.squared)
      control_coefs <- lapply(X=models,
                              FUN=function(x2, x3) control_extractor(x2,x3),
                              x3=x)

    }
    # Fixed effects
    else{
      # Get each value of interest across models
      coef <- lapply(X=models, function(x2) x2$coeftable[x,1])
      se <- lapply(X=models, function(x2) x2$coeftable[x,2])
      statistic <- lapply(X=models, function(x2) x2$coeftable[x,3])
      p <- lapply(X=models, function(x2) x2$coeftable[x,4])
      terms <- lapply(X=models, FUN=function(x2) row.names(x2$coeftable))
      RMSE <- lapply(X=models, FUN=function(x2) fitstat(x2, type="rmse",
                                                        verbose=FALSE)[[1]])
      adjR <- lapply(X=models, FUN=function(x2) fitstat(x2, type="war2",
                                                        verbose=FALSE)[[1]])
      control_coefs <- lapply(X=models,
                              FUN=function(x2,x3)
                                control_extractor(x2, x3, feols_model=TRUE),
                              x3=x)
    }


    # Store values in a data frame to be returned
    retVal <- data.frame(coef=unlist(coef), se=unlist(se),
                         statistic=unlist(statistic),
                         p=unlist(p), RMSE=unlist(RMSE), adjR=unlist(adjR))

    # R doesn't like it when these kinds of objects are assigned above
    retVal$terms <- terms
    retVal$control_coefs <- control_coefs

    retVal <- retVal %>%
      mutate(
        sig.level=case_when(
          p < .005 ~ "p < .005",
          p < .05 ~ "p < .05",
          p < .1 ~ "p < .1",
          p >= .1 ~ "p >= .1",
          TRUE ~ NA_character_
        )) %>%
      arrange(coef) %>%
      mutate(index=row_number())

  }
  # glm models
  else{
    # Get each value of interest across models
    coef <- lapply(X=models, function(x2) x2$coefficients[x,1])
    se <- lapply(X=models, function(x2) x2$coefficients[x,2])
    statistic <- lapply(X=models, function(x2) x2$coefficients[x,3])
    p <- lapply(X=models, function(x2) x2$coefficients[x,4])
    terms <- lapply(X=models, FUN=function(x2) row.names(x2$coefficients))
    AIC <- lapply(X=models, FUN=function(x2) x2$aic)
    deviance <- lapply(X=models, FUN=function(x2) x2$deviance)
    control_coefs <- lapply(X=models,
                            FUN=function(x2,x3,
                                         x4) control_extractor(x2,x3),x3=x)


    # Store values in a data frame to be returned
    retVal <- data.frame(coef=unlist(coef), se=unlist(se),
                         statistic=unlist(statistic),
                         p=unlist(p), AIC=unlist(AIC),
                         deviance=unlist(deviance))

    # R doesn't like it when these kinds of objects are assigned above
    retVal$terms <- terms
    retVal$control_coefs <- control_coefs

    retVal <- retVal %>%
      mutate(
        sig.level=case_when(
          p < .005 ~ "p < .005",
          p < .05 ~ "p < .05",
          p < .1 ~ "p < .1",
          p >= .1 ~ "p >= .1",
          TRUE ~ NA_character_
        )) %>%
      arrange(coef) %>%
      mutate(index=row_number())
  }

  # Build dummy columns for terms present in each model for visualization
  temp <- data.frame(matrix(ncol = length(controls), nrow = nrow(retVal)))

  control_names <- str_replace(controls, fixed("*"), fixed(":"))

  colnames(temp) <- control_names

  retVal <- cbind(retVal, temp)

  for(c in control_names){
    # Hiding the following warning:
    # In stri_detect_fixed(string, pattern, negate = negate,
    # opts_fixed = opts(pattern)): argument is not an atomic vector; coercing
    suppressWarnings(retVal[c] <- ifelse(str_detect(retVal$terms, fixed(c)),
                                         1, 0))
  }

  # Remove duplicate columns
  retVal <- retVal %>% select(where(~!all(is.na(.x))))

  return(retVal)

}

#' Plots a specification curve.
#'
#' @description
#' plot_curve() takes the data frame output of sca() and produces a ggplot of
#' the independent variable's coefficient (as indicated in the call to sca())
#' across model specifications. By default a panel is added showing which
#' control variables are present in each model. The combined plot is returned as
#' a `patchwork` object, so it can be further customised with ggplot2 and
#' patchwork operators (e.g. `& theme_sca(base_size = 14)`).
#'
#' @param sca_data A data frame returned by `sca()` containing model estimates
#'                 from the specification curve analysis.
#' @param title A string to use as the plot title. Defaults to an empty string,
#'              `""`.
#' @param show_index A boolean indicating whether to label the model index on the
#'                  the x-axis. Defaults to `TRUE`.
#' @param plot_vars A boolean indicating whether to include a panel on the plot
#'                 showing which variables are present in each model. Defaults
#'                 to `TRUE`.
#' @param ylab A string to be used as the y-axis label. Defaults to
#'             `"Coefficient"`.
#' @param plot_se A string indicating whether to display standard errors as
#'               bars or plots. For bars `plot_se = "bar"`, for ribbons
#'               `plot_se = "ribbon"`. If any other value is supplied then no
#'               standard errors are included. Defaults to `"bar"`.
#' @param median_line A boolean indicating whether to add a dotted line at the
#'                   median coefficient across specifications. Defaults to
#'                   `FALSE`.
#' @param point_size A number giving the size of the plotted points. Defaults to
#'                  `NULL`, in which case a size is chosen automatically from the
#'                  number of specifications.
#'
#' @return If `plot_vars = TRUE` a `patchwork` object combining the curve and the
#'         variable panel; if `plot_vars = FALSE` a ggplot object. Both can be
#'         further customised with ggplot2 (and patchwork) operators.
#'
#' @export
#'
#' @examples
#' plot_curve(sca_data = sca(y="Salnty", x="T_degC", c("ChlorA", "O2Sat"),
#'                          data=bottles, progress_bar=TRUE, parallel=FALSE),
#'                      title = "Salinity and Temperature Models",
#'                      show_index = TRUE, plot_vars = TRUE,
#'                      ylab = "Coefficient value", plot_se = "ribbon");
#' plot_curve(sca_data = sca(y="Salnty", x="T_degC",
#'                          c("ChlorA*O2Sat", "ChlorA", "O2Sat"),
#'                          data=bottles, progress_bar=FALSE, parallel=FALSE),
#'                      show_index = TRUE, plot_vars = TRUE,
#'                      plot_se = "ribbon");
#' \donttest{
#' plot_curve(sca_data = sca(y="Salnty", x="T_degC",
#'                          c("ChlorA*NO3uM", "O2Sat", "ChlorA", "NO3uM"),
#'                          data=bottles,
#'                          progress_bar = TRUE, parallel = TRUE, workers=2),
#'           plot_se="");
#' }
plot_curve <- function(sca_data, title="", show_index=TRUE, plot_vars=TRUE,
                         ylab="Coefficient", plot_se="bar", median_line=FALSE,
                         point_size=NULL){

  if("control_coefs" %in% names(sca_data)){
    sca_data <- sca_data %>% select(-control_coefs)
  }

  if(is.null(point_size)) point_size <- spec_point_size(sca_data)

  # Order the significance bins so the colour scale is consistent across plots.
  sca_data <- sca_data %>%
    mutate(sig.level = factor(sig.level, levels = names(sca_sig_colors())))

  if(tolower(plot_se)=="ribbon"){
    sca_data <- sca_data %>%
      mutate(ribbon.group = cumsum(sig.level != stats::lag(sig.level,
                                                    def = first(sig.level))))
  }

  sc1 <- ggplot(data=sca_data, aes(y=coef, x=index)) +
    geom_hline(yintercept = 0, color="red", linetype="dashed", linewidth=.6) +
    {if(median_line) geom_hline(yintercept = stats::median(sca_data$coef),
                               color="grey40", linetype="dotted",
                               linewidth=.6)} +
    {if(tolower(plot_se)=="ribbon") geom_ribbon(aes(ymin=coef-se, ymax=coef+se,
                                           group=factor(ribbon.group),
                                           fill=sig.level),
                                       alpha=.4)} +
    {if(tolower(plot_se)=="bar") geom_errorbar(aes(ymin=coef-se, ymax=coef+se,
                                          color=sig.level),
                                      width=0.25)} +
    {if(!tolower(plot_se) %in% c("ribbon", "bar"))
      geom_point(aes(color=sig.level), size=point_size)} +
    {if(tolower(plot_se) %in% c("ribbon", "bar"))
      geom_point(color="black", size=point_size)} +
    {if(tolower(plot_se)!="ribbon")
      scale_color_manual(values=sca_sig_colors(), drop=TRUE)} +
    {if(tolower(plot_se)=="ribbon")
      scale_fill_manual(values=sca_sig_colors(), drop=TRUE)} +
    labs(title=title, x="", y=ylab) +
    theme_sca() +
    theme(axis.text.x = {if(show_index) element_text() else element_blank()}) +
    guides(color = guide_legend(override.aes = list(size=2)),
           fill  = guide_legend(override.aes = list(size=2)))

  if(plot_vars){
    sc2 <- plot_vars(sca_data)
    return(patchwork::wrap_plots(sc1, sc2, ncol=1, heights=c(3, 1)))
  }
  else{
    return(sc1)
  }
}

#' Plots the variables in each model.
#'
#' @description
#' plot_vars() plots the variables included in each model specification in order
#' of model index. Returns a ggplot object that can then be combined with the
#' output of other functions like plot_rmse() if further customization of each
#' plot is desired.
#'
#' @inheritParams plot_curve
#' @param color_controls A boolean indicating whether to give each variable a
#'                      color to improve readability. Defaults to `FALSE`.
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @examples
#' plot_vars(sca_data = sca(y = "Salnty", x = "T_degC",
#'                         controls = c("ChlorA", "O2Sat"),
#'                         data = bottles, progress_bar = TRUE,
#'                         parallel = FALSE),
#'                      title = "Model Variable Specifications");
#' plot_vars(sca_data = sca(y = "Salnty", x = "T_degC",
#'                         controls = c("ChlorA*O2Sat"),
#'                         data = bottles, progress_bar = FALSE,
#'                         parallel = FALSE),
#'                      color_controls = TRUE);
#' \donttest{
#' plot_vars(sca_data = sca(y = "Salnty", x = "T_degC",
#'                         controls = c("ChlorA*NO3uM", "O2Sat*NO3uM"),
#'                         data = bottles,
#'                         progress_bar = TRUE, parallel = TRUE, workers = 2));
#' }
plot_vars <- function(sca_data, title="", color_controls=FALSE){

  if("control_coefs" %in% names(sca_data)){
    sca_data <- sca_data %>% select(-control_coefs)
  }

  scp_data <- scp(sca_data)

  markSize <- 10/length(scp_data[[2]])

  margin <- {if(title=="") unit(c(-5,2,-5,2), "points")
             else unit(c(5,2,-5,2), "points")}

  if(color_controls){
    sc <- ggplot(data=scp_data[[1]],
                  aes(x=index,y=factor(controlID), color=factor(controlID))
    ) +
      geom_point(shape="|", size=markSize) +
      labs(y="", x="", title=title) +
      scale_y_discrete(labels=scp_data[[2]], expand=c(.25,.25)) +
      theme_void() +
      theme(
        legend.position = "none",
        axis.text.y = element_text(size=6, hjust=0),
        axis.text.x = element_blank(),
        plot.margin = margin
      )
  }
  else{
      sc <- ggplot(data=scp_data[[1]],
                   aes(x=index,y=factor(controlID))
      ) +
        geom_point(shape="|", size=markSize) +
        labs(y="", x="", title=title) +
        scale_y_discrete(labels=scp_data[[2]], expand=c(.25,.25)) +
        theme_void() +
        theme(
          legend.position = "none",
          axis.text.y = element_text(size=6, hjust=0),
          axis.text.x = element_blank(),
          plot.margin = margin
        )
    }
  return(sc)
}

# Internal: point size for the specification-curve scatter plots, scaled down
# as the number of control-indicator columns grows so dense plots stay legible.
spec_point_size <- function(sca_data){
  -.25 * (ncol(sca_data) - 7) + (13 / 4)
}

# Internal: shared implementation behind plot_rmse()/plot_r2_adj()/plot_aic()/
# plot_deviance(), which differ only in the metric column, the axis label, and
# the message shown when that metric is absent from the sca() output.
plot_metric <- function(sca_data, metric, ylab, missing_message,
                        title="", show_index=TRUE, plot_vars=TRUE){

  if(!metric %in% colnames(sca_data)){
    message(missing_message)
    return(invisible(NULL))
  }

  sca_data <- sca_data %>% select(-control_coefs)

  point_size <- spec_point_size(sca_data)

  sc1 <- ggplot(data=sca_data, aes(x=.data$index, y=.data[[metric]])) +
    geom_point(size=point_size) +
    labs(title=title, x="", y=ylab) +
    theme_sca() +
    theme(axis.text.x = {if(show_index) element_text() else element_blank()})

  if(plot_vars){
    sc2 <- plot_vars(sca_data)
    return(patchwork::wrap_plots(sc1, sc2, ncol=1, heights=c(3, 1)))
  }
  else{
    return(sc1)
  }
}

#' Plots RMSE across model specifications.
#'
#' @description
#' plot_rmse() plots the root mean square error across model specifications. Only
#' available for linear regression models.
#'
#' @inheritParams plot_curve
#' @param show_index A boolean indicating whether to label the model index on the
#'                  the x-axis. Defaults to `TRUE`.
#' @param plot_vars A boolean indicating whether to include a panel on the plot
#'                 showing which variables are present in each model. Defaults
#'                 to `TRUE`.
#'
#' @return If `plot_vars = TRUE` a `patchwork` object combining the plot and the
#'         variable panel; if `plot_vars = FALSE` a ggplot object.
#'
#' @export
#'
#' @examples
#' plot_rmse(sca_data = sca(y="Salnty", x="T_degC", c("ChlorA", "O2Sat"),
#'                          data=bottles, progress_bar=TRUE, parallel=FALSE),
#'                      title = "RMSE");
#' plot_rmse(sca_data = sca(y="Salnty", x="T_degC", c("ChlorA*O2Sat"),
#'                          data=bottles, progress_bar=FALSE, parallel=FALSE),
#'                      show_index = FALSE, plot_vars = FALSE);
#' \donttest{
#' plot_rmse(sca_data = sca(y="Salnty", x="T_degC",
#'                          c("ChlorA*NO3uM", "O2Sat*NO3uM"), data=bottles,
#'                          progress_bar = TRUE, parallel=TRUE, workers=2));
#' }
plot_rmse <- function(sca_data, title="", show_index=TRUE, plot_vars=TRUE){
  plot_metric(sca_data, metric="RMSE", ylab="RMSE",
              missing_message=paste0("RMSE not found. Are your models nonlinear? ",
                                     "Try plot_aic() or plot_deviance() instead."),
              title=title, show_index=show_index, plot_vars=plot_vars)
}

#' Plots the adj. R-squared across model specifications.
#'
#' @description
#' plot_r2_adj() plots the adjusted R-squared across model specifications. Only
#' available for linear regression models. Note when fixed effects are
#' are specified the within adjusted R-squared is used (i.e. `fixest::r2()`
#' with `type="war2"`).
#'
#' @inheritParams plot_rmse
#'
#' @return If `plot_vars = TRUE` a `patchwork` object combining the plot and the
#'         variable panel; if `plot_vars = FALSE` a ggplot object.
#'
#' @export
#'
#' @examples
#' plot_r2_adj(sca_data = sca(y = "Salnty", x = "T_degC",
#'                          controls = c("ChlorA", "O2Sat"),
#'                          data = bottles, progress_bar = TRUE,
#'                          parallel = FALSE),
#'                      title = "Adjusted R^2");
#' plot_r2_adj(sca_data = sca(y="Salnty", x="T_degC",
#'                          controls = c("ChlorA*O2Sat"),
#'                          data = bottles, progress_bar = FALSE,
#'                          parallel = FALSE),
#'                      show_index = FALSE, plot_vars = FALSE);
#' \donttest{
#' plot_r2_adj(sca_data = sca(y = "Salnty", x = "T_degC",
#'                          controls = c("ChlorA*NO3uM", "O2Sat*NO3uM"),
#'                          data = bottles,
#'                          progress_bar = TRUE, parallel = TRUE, workers = 2));
#' }
plot_r2_adj <- function(sca_data, title="", show_index=TRUE, plot_vars=TRUE){
  plot_metric(sca_data, metric="adjR", ylab=bquote('Adj. R'^2),
              missing_message=paste0("Adj. R^2 not found. Are your models nonlinear? ",
                                     "Try plot_aic() or plot_deviance() instead."),
              title=title, show_index=show_index, plot_vars=plot_vars)
}

#' Plots the AIC across model specifications.
#'
#' @description
#' plot_aic() plots the Akaike information criterion across model specifications.
#' Only available for nonlinear regression models.
#'
#' @inheritParams plot_rmse
#'
#' @return If `plot_vars = TRUE` a `patchwork` object combining the plot and the
#'         variable panel; if `plot_vars = FALSE` a ggplot object.
#'
#' @export
#'
#' @examples
#' plot_aic(sca_data = sca(y = "Salnty", x = "T_degC",
#'                        controls = c("ChlorA", "O2Sat"),
#'                        data = bottles, progress_bar = TRUE, parallel = FALSE),
#'                      title = "AIC");
#' plot_aic(sca_data = sca(y = "Salnty", x = "T_degC",
#'                        controls = c("ChlorA*O2Sat"),
#'                        data = bottles, progress_bar = FALSE,
#'                        parallel = FALSE),
#'                        show_index = FALSE, plot_vars = FALSE);
#' \donttest{
#' plot_aic(sca_data = sca(y = "Salnty", x = "T_degC",
#'                          controls = c("ChlorA*NO3uM", "O2Sat*NO3uM"),
#'                          data = bottles,
#'                          progress_bar = TRUE, parallel = TRUE, workers = 2));
#' }
plot_aic <- function(sca_data, title="", show_index=TRUE, plot_vars=TRUE){
  plot_metric(sca_data, metric="AIC", ylab="AIC",
              missing_message=paste0("AIC not found. Are your models linear? ",
                                     "Try plot_r2_adj() or plot_rmse instead."),
              title=title, show_index=show_index, plot_vars=plot_vars)
}

#' Plots the deviance of residuals across model specifications.
#'
#' @description
#' plot_deviance() plots the deviance of residuals across model specifications.
#' Only available for linear regression models.
#'
#' @inheritParams plot_rmse
#'
#' @return If `plot_vars = TRUE` a `patchwork` object combining the plot and the
#'         variable panel; if `plot_vars = FALSE` a ggplot object.
#'
#' @export
#'
#' @examples
#' plot_deviance(sca_data = sca(y = "Salnty", x = "T_degC",
#'                             controls = c("ChlorA", "O2Sat"),
#'                             data = bottles, progress_bar = TRUE,
#'                             parallel = FALSE),
#'                      title = "Model Deviance");
#' plot_deviance(sca_data = sca(y = "Salnty", x = "T_degC",
#'                             controls = c("ChlorA*O2Sat"),
#'                             data = bottles, progress_bar = FALSE,
#'                             parallel = FALSE),
#'                      show_index = FALSE, plot_vars = FALSE);
#' \donttest{
#' plot_deviance(sca_data = sca(y = "Salnty", x="T_degC",
#'                          controls = c("ChlorA*NO3uM", "O2Sat*NO3uM"),
#'                          data = bottles, progress_bar = TRUE, parallel = TRUE,
#'                          workers = 2));
#' }
plot_deviance <- function(sca_data, title="", show_index=TRUE, plot_vars=TRUE){
  plot_metric(sca_data, metric="deviance", ylab="Deviance",
              missing_message=paste0("Deviance of residuals not found. ",
                                     "Are your models linear? Try plot_aic(), ",
                                     "plot_r2_adj(), or plot_rmse() instead."),
              title=title, show_index=show_index, plot_vars=plot_vars)
}

#' Plots control variable distributions.
#'
#' @description
#' plot_control_distributions() plots the distribution of coefficients for each
#' control variable included in the model specifications.
#'
#' @inheritParams plot_rmse
#' @param type A string indicating what type of distribution plot to produce.
#'             When `type = "density"` density plots are produced. When
#'             `type = "hist"` or `type = "histogram"` histograms are produced.
#'             Defaults to `"density"`.
#' @param zero_line A boolean indicating whether to draw a dashed reference line
#'                 at zero, making it easy to see whether a control's effect
#'                 crosses zero. Defaults to `TRUE`.
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @examples
#' plot_control_distributions(sca_data = sca(y="Salnty", x="T_degC",
#'                                     controls = c("ChlorA", "O2Sat"),
#'                                     data = bottles,
#'                                     progress_bar = TRUE, parallel = FALSE),
#'                          title = "Control Variable Distributions")
#' plot_control_distributions(sca_data = sca(y = "Salnty", x="T_degC",
#'                                     controls = c("ChlorA*O2Sat"),
#'                                     data = bottles,
#'                                     progress_bar = FALSE, parallel = FALSE),
#'                          type = "hist")
#' \donttest{
#' plot_control_distributions(sca_data = sca(y = "Salnty", x = "T_degC",
#'                                     controls = c("ChlorA*NO3uM",
#'                                                  "O2Sat*NO3uM"),
#'                                     data = bottles, progress_bar = TRUE,
#'                                     parallel = TRUE, workers = 2),
#'                          type = "density")
#' }
plot_control_distributions <- function(sca_data, title="", type="density",
                                     zero_line=TRUE){

  histData <- bind_rows(un_as_is(sca_data$control_coefs))

  rownames(histData) <- NULL

  histData$term <- sapply(sapply(histData$term, str_split, pattern=":"),
                          paste0, collapse=" %*% ")

  n_facets <- length(unique(histData$term))

  # A single cohesive fill from the package palette rather than a clashing
  # colour per facet (the colour carried no information).
  fillColor <- sca_sig_colors()[["p < .005"]]

  sc1 <- histData %>%
    ggplot(aes(x=coef)) +
      {if(tolower(type)=="hist" | tolower(type)=="histogram")
         geom_histogram(fill=fillColor, color="white")
       else if (tolower(type)=="density")
         geom_density(fill=fillColor, color="grey20", alpha=.85)} +
      {if(zero_line) geom_vline(xintercept=0, color="red", linetype="dashed",
                               linewidth=.5)} +
      labs(x="", y="", title=title) +
      theme_sca() +
      theme(legend.position="none") +
    {if(n_facets>16) theme(axis.text.x=element_text(size=4),
                           axis.text.y=element_text(size=4),
                           strip.text=element_text(size=6))
      else theme()} +
      facet_wrap(~factor(term), scales="free", labeller = label_parsed)
  return(sc1)
}


# Internal: bootstrapped SEs for se_compare(). Estimates one column per
# (boot_samples, boot_sample_size) combination and returns a named matrix (rows =
# coefficients), or NULL if no estimates were produced. `fe_suffix` tags the
# column names for fixed-effects models.
boot_ses <- function(data, formula, n_x, boot_samples, boot_sample_size,
                     weights=NULL, fe_suffix=""){
  if(length(boot_samples)==1 & length(boot_sample_size)==1){
    samples <- boot_samples
    sample_sizes <- boot_sample_size

    if(is.null(weights)){
      boot <- se_boot(data=data, formula=formula, n_x=n_x,
                      n_samples=boot_samples[[1]],
                      sample_size=boot_sample_size[[1]])
    }
    else{
      boot <- se_boot(data=data, formula=formula, n_x=n_x,
                      n_samples=boot_samples[[1]],
                      sample_size=boot_sample_size[[1]],
                      weights=weights)
    }

    if(is.null(boot)) return(NULL)
    boot <- matrix(boot, ncol=1, dimnames=list(names(boot), NULL))
  }
  else{
    samples <- rep(boot_samples, length(boot_sample_size))
    sample_sizes <- sort(rep(boot_sample_size, length(boot_samples)))

    if(is.null(weights)){
      boot <- mapply(FUN=se_boot, n_samples=samples,
                     sample_size=sample_sizes,
                     MoreArgs=list(data=data, formula=formula, n_x=n_x))
    }
    else{
      boot <- mapply(FUN=se_boot, n_samples=samples,
                     sample_size=sample_sizes,
                     MoreArgs=list(data=data, formula=formula, n_x=n_x,
                                   weights=weights))
    }

    if(is.null(boot)) return(NULL)
  }

  colnames(boot) <- paste("bootstrap_", "k", samples, "n", sample_sizes,
                          fe_suffix, sep="")
  boot
}

#' Compare different kinds of standard errors
#'
#' @description
#' se_compare() takes in a regression formula (with or without fixed effects),
#' data, and the types of standard errors desired, including clustered,
#' heteroskedasticity-consistent, and bootstrapped. It then returns a data
#' frame with coefficient and standard error estimates for easy comparison and
#' plotting.
#'
#' @param formula A string containing a regression formula, with or without
#'                fixed effects.
#' @param data A data frame containing the variables provided in `formula` and
#'             any clustering variables passed to `cluster`.
#' @param weights Optional string with the column name in `data` that contains
#'                weights.
#' @param types A string or vector of strings specifying what types of
#'              standard errors are desired. Defaults to "all".
#'
#'              The following types are supported for non-fixed effects models:
#'
#'                With clustering: "HC0, "HC1", "HC2", "HC3".
#'
#'                Without clustering: "iid" (i.e. normal standard errors),
#'                                    "HC0, "HC1", "HC2", "HC3",
#'                                    "HC4", "HC4m", "HC5",
#'                                    "bootstrapped".
#'
#'              The following types are supported for fixed effects models:
#'
#'                With clustering: "CL_FE" (clustered by fixed effects, i.e.
#'                                 the default standard errors reported by
#'                                 `feols()` if no clusters are supplied), if
#'                                 clusters are supplied then the conventional
#'                                 clustered standard errors from `feols()` are
#'                                 estimated for each clustering variable. Two-
#'                                 way clustered standard errors are not
#'                                 supported at this time.
#'
#'                Without clustering: "HC0, "HC1", "HC2", "HC3",
#'                                    "HC4", "HC4m", "HC5",
#'                                    "bootstrapped".
#' @param cluster A string or vector of strings specifying variables present in
#'                `data` to be used for clustering standard errors.
#' @param clustered_only A boolean indicating whether only standard errors with
#'                      clustering should be estimated, defaults to `FALSE`.
#' @param fixed_effects_only A boolean indicating whether only standard errors for
#'                         fixed effects models should be estimated, defaults to
#'                         `FALSE`.
#' @param boot_samples An integer or vector of integers indicating how many times
#'                    the model should be estimated with a random subset of the
#'                    data. If a vector then every combination of `boot_samples`
#'                    and `boot_sample_size` are estimated.
#' @param boot_sample_size An integer or vector of integers indicating how many
#'                       observations are in each random subset of the data.
#'                       If a vector then every combination of `boot_samples`
#'                       and `boot_sample_size` are estimated.
#' @param ... Deprecated camelCase arguments (`clusteredOnly`,
#'            `fixedEffectsOnly`, `bootSamples`, `bootSampleSize`); use the
#'            snake_case equivalents instead.
#'
#' @return A data frame where row represents an independent variable in the
#'         model and each column a type of standard error. Coefficient estimates
#'         for each variable are also included (column `"estimate"` for
#'         non-fixed effects model and column `"estimate_FE"` for fixed effects
#'         models). Columns are automatically named to specify the standard
#'         error type.
#'
#'         Some examples:
#'
#'          "iid" = normal standard errors, i.e. assuming homoskedasticity
#'
#'          "CL_FE" = standard errors clustered by fixed effects
#'
#'          "bootstrap_k8n300_FE" =  bootstrapped standard errors for a fixed
#'                                   effects model where `boot_samples = 8` and
#'                                   `boot_sample_size = 300`
#'
#'          "CL_Depth_ID_FE" = standard errors clustered by the variable
#'                               "Depth_ID" for a model with fixed effects
#'
#'          "HC0_Sta_ID" = HC0 standard errors clustered by the variable
#'                           "Sta_ID"
#'
#'          Note: for fixed effects models the "(Intercept)" row will be all
#'          `NA` because the intercept is not reported by `feols()` when fixed
#'          effects are present.
#'
#' @export
#'
#' @examples
#'
#' se_compare(formula = "Salnty ~ T_degC + ChlorA + O2Sat | Sta_ID",
#'            data = bottles, types = "all", cluster = c("Depth_ID", "Sta_ID"),
#'            fixed_effects_only = FALSE, boot_samples=c(4, 8, 10),
#'            boot_sample_size=c(300, 500))
#'
#' se_compare(formula = "Salnty ~ T_degC + ChlorA + O2Sat", data = bottles,
#'            types = "bootstrapped", boot_samples = c(8, 10),
#'            boot_sample_size = c(300, 500))
#'
#' se_compare(formula = "Salnty ~ T_degC + ChlorA", data = bottles,
#'            types = c("HC0", "HC1", "HC3"))
#'
se_compare <- function(formula, data, weights=NULL,
                       types="all", cluster=NULL,
                       clustered_only=FALSE, fixed_effects_only=FALSE,
                       boot_samples=NULL, boot_sample_size=NULL, ...){

  # Backward compatibility: translate deprecated camelCase argument names.
  .dots <- list(...)
  if("clusteredOnly" %in% names(.dots)){
    .Deprecated(msg="`clusteredOnly` is deprecated; use `clustered_only`.")
    clustered_only <- .dots$clusteredOnly
  }
  if("fixedEffectsOnly" %in% names(.dots)){
    .Deprecated(msg="`fixedEffectsOnly` is deprecated; use `fixed_effects_only`.")
    fixed_effects_only <- .dots$fixedEffectsOnly
  }
  if("bootSamples" %in% names(.dots)){
    .Deprecated(msg="`bootSamples` is deprecated; use `boot_samples`.")
    boot_samples <- .dots$bootSamples
  }
  if("bootSampleSize" %in% names(.dots)){
    .Deprecated(msg="`bootSampleSize` is deprecated; use `boot_sample_size`.")
    boot_sample_size <- .dots$bootSampleSize
  }

  # Create objects that will store the standard errors
  ses_CL <- NULL
  ses_HC <- NULL
  ses_other <- NULL

  # Create the object we will eventually return
  ses <- NULL

  # Whether the formula specifies fixed effects (a pipe). Captured up front
  # because `formula` is later stripped of its fixed effects for the non-FE
  # model. "CL_FE" is a fixed-effects-only type, so the non-FE branch should
  # not flag it as invalid when fixed effects are present.
  has_fe <- grepl("|", formula, fixed=TRUE)

  # Validate the columns referenced by the formula and (if supplied) the
  # weights. Clustering variables are validated where they are used, with a
  # warning rather than an error.
  formula_vars <- setdiff(all.vars(stats::as.formula(formula)), ".")
  check_columns(data, formula_vars, "Variable(s)")
  if(!is.null(weights)) check_columns(data, weights, "Weights variable")

  # If the formula contains a pipe then fixed effects are assumed to be
  # present and models are estimated with feols() rather than lm()
  if(grepl("|", formula, fixed=TRUE)){

    if(is.null(weights)){
      model_fe <- tryCatch(feols(as.formula(formula), data=data),
                           error=function(cond){
                             message("Fixed effects model estimation failed.",
                                     cond)
                             return(NULL)
                           })
    }
    else{
      model_fe <- tryCatch(
        feols(as.formula(formula), data=data, weights=data[[weights]]),
                           error=function(cond){
                             message("Fixed effects model estimation failed.",
                                     cond)
                             return(NULL)
                           })
    }

      if(is.null(model_fe)){
        message("Fixed effects model estimation failed.")
      }
      else{
      # Add FE model coefficients to the matrix
      ses <- cbind(ses, matrix(c("(Intercept)"=NA, model_fe$coefficients), ncol=1,
                               dimnames=list(c("(Intercept)",
                                               names(model_fe$coefficients)),
                                             c("estimate_FE"))))

      # Case when user wants to cluster by FEs (i.e. the default SEs reported
      # by feols()) or bootstrap
      if(!clustered_only){
        types_other <- c("CL_FE","bootstrapped")

        if(!"all" %in% types){

          if(length(setdiff(types, types_other))!=0){
            warning(paste0(setdiff(types, types_other),
                           " not a valid type for SEs in FE model, ignoring.",
                           collapse="\n"))
          }

          types_other <- types[types %in% types_other]

        }

        # Get the default standard errors from feols() output
        if("CL_FE" %in% types_other){
          ses_other <- cbind(ses_other, "CL_FE"=c("(Intercept)"=NA,
                                                  coeftable(model_fe)[,2]))
        }

        # Get bootstrapped SEs
        if("bootstrapped" %in% types_other & !is.null(boot_samples) &
           !is.null(boot_sample_size)){

          boot <- boot_ses(data=data, formula=formula,
                           n_x=length(model_fe$coefficients),
                           boot_samples=boot_samples, boot_sample_size=boot_sample_size,
                           weights=weights, fe_suffix="_FE")

          if(!is.null(boot)){
            ses_other <- cbind(ses_other, boot)
          }
        }
        # Attach bootstapped/default SEs to the object to be returned
        ses <- cbind(ses, ses_other)
      }

      # Estimate clustered standard errors for FE model for variables other than
      # the FEs
      if(!is.null(cluster)){

        if(length(setdiff(cluster, colnames(data)))!=0){
          warning(paste0(setdiff(cluster, colnames(data)),
                         " not a valid clustering variable, ignoring.",
                         collapse="\n"))

          cluster <- cluster[cluster %in% colnames(data)]
        }

        # NB: for fixed-effects models the clustered SEs are feols' default
        # cluster-robust SEs estimated per clustering variable; they do not
        # depend on `types`, so there is no type validation here.

        # Estimate standard errors clustered by each desired variable
        ses_CL <- sapply(cluster, FUN=function(c){
          (feols(as.formula(formula), data=data,
                 cluster=data[c]))$coeftable[,2]})

        # Label them nicely
        labs <- c()
        for(c in cluster){
          labs <- c(labs, paste0("CL", "_", c, "_FE"))
        }

        colnames(ses_CL) <- labs

        ses_CL <- rbind("(Intercept)"=NA, ses_CL)

        ses <- cbind(ses, ses_CL)
      }
    }

  }
  # Case when a non-FE model is desired
  if(!fixed_effects_only){

    # Allocate objects to hold SEs
    ses_other <- NULL
    ses_HC <- NULL
    ses_CL <- NULL

    # If the formula has FEs remove them
    if(grepl("|",formula, fixed=TRUE)){
      formula <- str_trim(str_split(formula, fixed("|"))[[1]][[1]])
    }

    # Estimate the non-FE model and get the coefficients
    if(is.null(weights)){
      model <- lm(formula=as.formula(formula), data=data)
    }
    else{
      fmla <- as.formula(formula)
      environment(fmla) <- environment()
      model <- lm(formula=fmla, data=data, weights=get(weights))
    }

    ses <- cbind(ses, matrix(model$coefficients, ncol=1,
                             dimnames=list(c(names(model$coefficients)),
                                           c("estimate"))))

    # Parse the user's desired SE types
    if(!clustered_only){
      types_HC <- c("HC0", "HC1", "HC2", "HC3", "HC4", "HC4m", "HC5")
      types_other <- c("iid", "bootstrapped")

      if(!"all" %in% types){

        # "CL_FE" is handled by the FE branch when fixed effects are present,
        # so it is not an invalid type here in that case.
        invalid <- setdiff(types, c(types_HC, types_other,
                                    if(has_fe) "CL_FE"))
        if(length(invalid)!=0){
          warning(paste0(invalid,
                         " not a valid type for SEs, ignoring.", collapse="\n"))
        }

        types_HC <- types[types %in% types_HC]
        types_other <- types[types %in% types_other]

      }
      # Get the normal iid standard errors
      if("iid" %in% types_other){
        ses_other <- cbind(ses_other, "iid"=summary(model)$coefficients[,2])
      }

      # Get HC standard errors
      ses_HC <- sapply(types_HC,
                       function(x) coeftest(model, vcov.=vcovHC, type=x)[,2])

      # Get bootstrapped standard errors
      if("bootstrapped" %in% types_other & !is.null(boot_samples) &
         !is.null(boot_sample_size)){

        boot <- boot_ses(data=data, formula=formula,
                         n_x=length(model$coefficients)-1,
                         boot_samples=boot_samples, boot_sample_size=boot_sample_size,
                         weights=weights, fe_suffix="")

        if(!is.null(boot)){
          ses_other <- cbind(ses_other, boot)
        }
      }

      # Attach SEs to the return object
      ses <- cbind(ses, ses_other, ses_HC)
    }

    # Case when clustered SEs for non-FE model are desired
    if(!is.null(cluster)){

      if(length(setdiff(cluster, colnames(data)))!=0){
        warning(paste0(setdiff(cluster, colnames(data)),
                       " not a valid clustering variable, ignoring.",
                       collapse="\n"))

        cluster <- cluster[cluster %in% colnames(data)]
      }

      types_CL <- c("HC0", "HC1", "HC2", "HC3")

      if(!"all" %in% types){

        # Bootstrapped/iid SEs are not clustered, and "CL_FE" is a
        # fixed-effects-only type, so none of them are invalid here.
        not_clustered <- c("bootstrapped", "iid", if(has_fe) "CL_FE")
        invalid <- setdiff(types[!types %in% not_clustered], types_CL)
        if(length(invalid)!=0){
          warning(paste0(invalid,
                         " not a valid type for clustered SEs, ignoring.",
                         collapse="\n"))
        }

        types_CL <- types[types %in% types_CL]
      }

      # Estimate and extract clustered SEs
      if(length(types_CL)>0){

        ses_CL <- sapply(cluster, FUN=function(c, types){
          sapply(types, function(x){
            coeftest(model, vcov.=vcovCL, type=x, cluster=data[c])[,2]
          })
        }, types=types_CL, simplify=FALSE)

        ses_CL <- do.call(cbind, ses_CL)

        # Label with type and clustering variable
        labs <- c()
        for(c in cluster){
          for(t in types_CL){
            labs <- c(labs, paste0(t, "_", c))
          }
        }

        colnames(ses_CL) <- labs
      }

      ses <- cbind(ses, ses_CL)
    }
  }

  # Coerce matrix to a data frame and return
  return(as.data.frame(apply(ses, FUN=unlist, MARGIN=2)))
}

#' Plots standard error estimates across types.
#'
#' @description
#' plot_se() takes the data frame output of `se_compare()` and plots, for each
#' coefficient, the estimate together with a confidence interval derived from
#' every available type of standard error. This makes it easy to see how
#' inference about a coefficient changes with the choice of standard error.
#'
#' @param se_data A data frame returned by `se_compare()`.
#' @param level The confidence level used for the intervals. Defaults to `0.95`.
#' @param intercept A boolean indicating whether to include the `"(Intercept)"`
#'                  coefficient. Defaults to `FALSE`.
#' @param title A string to use as the plot title. Defaults to an empty string,
#'              `""`.
#'
#' @return A ggplot object with one facet per coefficient; within each facet the
#'         estimate is plotted against each standard error type with a
#'         confidence interval, and the colour indicates whether that interval
#'         excludes zero.
#'
#' @export
#'
#' @examples
#' plot_se(se_compare(formula = "Salnty ~ T_degC + ChlorA", data = bottles,
#'                   types = c("iid", "HC0", "HC3")));
#' plot_se(se_compare(formula = "Salnty ~ T_degC + ChlorA", data = bottles,
#'                   types = "HC1", cluster = c("Sta_ID", "Depth_ID")),
#'        level = 0.9);
plot_se <- function(se_data, level=0.95, intercept=FALSE, title=""){

  df <- as.data.frame(se_data)
  df$term <- rownames(df)
  rownames(df) <- NULL

  has_fe  <- "estimate_FE" %in% names(df)
  has_ols <- "estimate" %in% names(df)
  est_cols <- intersect(c("estimate", "estimate_FE"), names(df))
  se_cols <- setdiff(names(df), c(est_cols, "term"))

  if(length(se_cols) == 0){
    message("No standard error columns found in se_data.")
    return(invisible(NULL))
  }

  # Half-width multiplier for the requested confidence level.
  z <- stats::qnorm(1 - (1 - level) / 2)

  long <- se_data %>%
    as.data.frame() %>%
    mutate(term = rownames(.)) %>%
    pivot_longer(cols=all_of(se_cols), names_to="se_type", values_to="se") %>%
    # Pair each SE type with its estimate (fixed-effects types end in "_FE").
    mutate(is_fe = grepl("_FE$", se_type))

  if(has_fe & has_ols){
    long <- long %>% mutate(estimate = ifelse(is_fe, estimate_FE, estimate))
  }
  else if(has_fe){
    long <- long %>% mutate(estimate = estimate_FE)
  }

  if(!intercept){
    long <- long %>% filter(term != "(Intercept)")
  }

  long <- long %>%
    filter(!is.na(se), !is.na(estimate)) %>%
    mutate(lower = estimate - z * se,
           upper = estimate + z * se,
           sig   = (lower > 0) | (upper < 0))

  if(nrow(long) == 0){
    message("Nothing to plot after removing missing values.")
    return(invisible(NULL))
  }

  ggplot(data=long, aes(x=se_type, y=estimate, color=sig)) +
    geom_hline(yintercept=0, color="red", linetype="dashed") +
    geom_pointrange(aes(ymin=lower, ymax=upper)) +
    facet_wrap(~term, scales="free_y") +
    scale_color_manual(values=c("FALSE"="#999999", "TRUE"="#0072B2"),
                       labels=c("FALSE"="CI includes 0",
                                "TRUE"="CI excludes 0"),
                       drop=FALSE) +
    labs(title=title, x="Standard error type", y="Estimate", color="") +
    theme_sca() +
    theme(axis.text.x=element_text(angle=45, hjust=1))
}

#' Plots how each control influences the independent variable's coefficient.
#'
#' @description
#' plot_influence() shows, for every control variable, the distribution of the
#' independent variable's coefficient across the specifications that include
#' versus exclude that control. It makes clear which modelling choices move the
#' estimate, and by how much.
#'
#' @param sca_data A data frame returned by `sca()`.
#' @param title A string to use as the plot title. Defaults to `""`.
#'
#' @return A ggplot object with one facet per control comparing the coefficient
#'         when that control is excluded versus included.
#'
#' @export
#'
#' @examples
#' plot_influence(sca(y = "Salnty", x = "T_degC",
#'                   controls = c("ChlorA", "O2Sat", "NO2uM"),
#'                   data = bottles, progress_bar = FALSE));
plot_influence <- function(sca_data, title=""){

  controls <- sca_control_cols(sca_data)
  if(length(controls) == 0){
    message("No control indicator columns found in sca_data.")
    return(invisible(NULL))
  }

  long <- sca_data %>%
    select(all_of(c("coef", controls))) %>%
    pivot_longer(all_of(controls), names_to="control", values_to="included") %>%
    mutate(included = factor(ifelse(included == 1, "Included", "Excluded"),
                             levels = c("Excluded", "Included")))

  ggplot(long, aes(x=included, y=coef, fill=included)) +
    geom_hline(yintercept=0, color="red", linetype="dashed", linewidth=.5) +
    geom_boxplot(outlier.size=.6, alpha=.9) +
    facet_wrap(~control) +
    scale_fill_manual(values=c("Excluded"="#9E9E9E", "Included"="#3182BD")) +
    labs(title=title, x="", y="Coefficient") +
    theme_sca() +
    theme(legend.position="none")
}

#' Plots the independent variable's coefficient against model fit.
#'
#' @description
#' plot_coef_fit() plots the independent variable's coefficient against a measure
#' of model fit across specifications, revealing whether better-fitting models
#' tend to produce systematically different estimates (i.e. whether your
#' best-fitting specifications are outliers).
#'
#' @param sca_data A data frame returned by `sca()`.
#' @param metric A string naming the fit measure to plot against, one of
#'               `"RMSE"`, `"adjR"`, `"AIC"`, or `"deviance"`. Defaults to
#'               `NULL`, in which case the first measure available in `sca_data`
#'               is used.
#' @param title A string to use as the plot title. Defaults to `""`.
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @examples
#' plot_coef_fit(sca(y = "Salnty", x = "T_degC",
#'                 controls = c("ChlorA", "O2Sat", "NO2uM"),
#'                 data = bottles, progress_bar = FALSE));
plot_coef_fit <- function(sca_data, metric=NULL, title=""){

  available <- intersect(c("RMSE", "adjR", "AIC", "deviance"), names(sca_data))
  if(length(available) == 0){
    message("No model-fit columns found in sca_data.")
    return(invisible(NULL))
  }
  if(is.null(metric)) metric <- available[1]
  if(!metric %in% available){
    stop("`metric` must be one of: ", paste(available, collapse=", "),
         call.=FALSE)
  }

  axis_labels <- c(RMSE="RMSE", adjR="Adjusted R-squared", AIC="AIC",
                   deviance="Deviance")

  sca_data <- sca_data %>%
    mutate(sig.level = factor(sig.level, levels = names(sca_sig_colors())))

  ggplot(sca_data, aes(x=.data[[metric]], y=coef, color=sig.level)) +
    geom_hline(yintercept=0, color="red", linetype="dashed", linewidth=.5) +
    geom_point(size=2) +
    scale_color_manual(values=sca_sig_colors(), drop=TRUE) +
    labs(title=title, x=axis_labels[[metric]], y="Coefficient") +
    theme_sca()
}

#' Plots a specification curve under multiple standard error types.
#'
#' @description
#' plot_multi_se() estimates every specification (as `sca()` does) and, for each,
#' computes the independent variable's standard error under several types via
#' `se_compare()`. It then plots the specification curve faceted by standard
#' error type: the coefficient estimates are identical across facets, but the
#' confidence intervals -- and hence which specifications are "significant" --
#' change with the choice of standard error, showing how sensitive your
#' conclusions are to that choice.
#'
#' @inheritParams sca
#' @param types A vector of standard error types to compare, passed to
#'              `se_compare()`: `"iid"`, the `"HC*"` types, or (with `cluster`)
#'              clustered types. Defaults to `c("iid", "HC3")`. Bootstrapped
#'              standard errors are not supported here; use `se_compare()`
#'              directly for those.
#' @param cluster Optional clustering variable(s) passed to `se_compare()`.
#' @param level The confidence level used for the intervals. Defaults to `0.95`.
#' @param title A string to use as the plot title. Defaults to `""`.
#'
#' @return A ggplot object: the specification curve faceted by standard error
#'         type, points coloured by significance under each type.
#'
#' @export
#'
#' @examples
#' plot_multi_se(y = "Salnty", x = "T_degC", controls = c("ChlorA", "O2Sat"),
#'             data = bottles, types = c("iid", "HC1", "HC3"));
plot_multi_se <- function(y, x, controls, data, types=c("iid", "HC3"),
                        cluster=NULL, fixed_effects=NULL, level=0.95, title=""){

  # Reuse sca()'s machinery to build the specification formulae.
  formulae <- sca(y=y, x=x, controls=controls, data=data,
                  fixed_effects=fixed_effects, return_formulae=TRUE)

  est_col <- if(!is.null(fixed_effects)) "estimate_FE" else "estimate"

  # For each specification, pull the focal variable's coefficient and its SE
  # under every requested type from se_compare().
  rows <- lapply(formulae, function(f){
    fstr <- paste(deparse(f), collapse=" ")
    res <- tryCatch(
      suppressWarnings(suppressMessages(
        se_compare(formula=fstr, data=data, types=types, cluster=cluster,
                   fixed_effects_only=!is.null(fixed_effects)))),
      error=function(e) NULL)
    if(is.null(res) || !x %in% rownames(res) || !est_col %in% colnames(res)){
      return(NULL)
    }
    se_cols <- setdiff(colnames(res), c("estimate", "estimate_FE"))
    data.frame(coef=res[x, est_col], se_type=se_cols,
               se=as.numeric(res[x, se_cols]), stringsAsFactors=FALSE)
  })

  long <- bind_rows(rows)
  if(nrow(long) == 0){
    message("No standard errors could be computed for the focal variable.")
    return(invisible(NULL))
  }

  z <- stats::qnorm(1 - (1 - level) / 2)

  long <- long %>%
    filter(!is.na(se)) %>%
    arrange(se_type, coef) %>%
    group_by(se_type) %>%
    mutate(index = row_number()) %>%
    ungroup() %>%
    mutate(
      p = 2 * stats::pnorm(-abs(coef / se)),
      sig.level = factor(case_when(
        p < .005 ~ "p < .005",
        p < .05  ~ "p < .05",
        p < .1   ~ "p < .1",
        TRUE     ~ "p >= .1"
      ), levels = names(sca_sig_colors())),
      lower = coef - z * se,
      upper = coef + z * se
    )

  ggplot(long, aes(x=index, y=coef)) +
    geom_hline(yintercept=0, color="red", linetype="dashed", linewidth=.5) +
    geom_errorbar(aes(ymin=lower, ymax=upper, color=sig.level), width=.25) +
    geom_point(color="black", size=.9) +
    facet_wrap(~se_type) +
    scale_color_manual(values=sca_sig_colors(), drop=TRUE) +
    labs(title=title, x="", y="Coefficient") +
    theme_sca()
}
