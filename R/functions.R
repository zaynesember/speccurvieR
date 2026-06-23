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
#' @param common_sample A boolean. When `TRUE`, every specification is fit on the
#'                      same sample: the rows that are complete across all model
#'                      variables (the dependent, independent, control, fixed-
#'                      effects, and weights variables). When `FALSE` (the
#'                      default) each specification uses its own complete cases,
#'                      so specifications with different controls may be fit on
#'                      different samples; the returned `n_obs` column reveals
#'                      this. See also [plot_samplesizes()].
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
#'         test statistic, p-value, model specification, measures of model fit,
#'         and `n_obs`, the number of observations the specification was fit on.
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
                fixed_effects=NULL, common_sample=FALSE, return_formulae=FALSE,
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

  # Resolve family/link into a normalised family string and (for glm families)
  # a family object. Treats "gaussian" as OLS and errors on an unknown family.
  resolved <- resolve_family(family, link)
  family <- resolved$family
  fam_obj <- resolved$fam_obj

  if(family!="linear" & !is.null(fixed_effects))
  {
    warning(paste0("Fixed effects unsupported for models other than OLS ",
                   "regression. Ignoring fixed effects."))
    # Actually drop the fixed effects so downstream estimation/extraction uses
    # the glm path, as the warning promises.
    fixed_effects <- NULL
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

  # Common-sample mode: by default each specification is fit on its own
  # complete cases, so specifications with different control sets can be fit on
  # different samples (and the curve then conflates control effects with sample
  # changes). When `common_sample = TRUE`, restrict `data` to the rows that are
  # complete across every model variable so all specifications share one sample.
  if(common_sample){
    union_vars <- unique(c(vars, fixed_effects, weights))
    cc <- stats::complete.cases(data[, union_vars, drop=FALSE])
    data <- data[cc, , drop=FALSE]
  }

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

  # Guard: the focal variable must correspond to a single model coefficient.
  # A factor, interaction, or transformed `x` expands to differently-named
  # rows (e.g. "xLevel2", "a:b") or none at all, in which case the name-based
  # extraction below would error with "subscript out of bounds" or silently
  # return NA for every specification. Check EVERY specification, not just the
  # first: a control name that collides with x, or a collinear drop, can leave
  # the focal coefficient out of some specifications but not others.
  has_x <- vapply(models, function(m){
    tms <- if(!is.null(fixed_effects)) rownames(m$coeftable)
           else rownames(m$coefficients)
    x %in% tms
  }, logical(1))
  if(!all(has_x)){
    stop("`x` (\"", x, "\") does not correspond to a single model coefficient ",
         "in every specification. It may be a factor, interaction, or ",
         "transformed term, or its name may collide with a control variable; ",
         "specification curve analysis requires a single focal coefficient.",
         call.=FALSE)
  }

  # Number of observations each specification was actually fit on (after
  # listwise deletion and, for fixed-effects models, fixest singleton removal),
  # so users can see when specifications use different samples. Extracted from
  # each model summary: feols reports `nobs`; for glm the null model has n - 1
  # degrees of freedom; for lm the rank plus residual degrees of freedom give n.
  n_obs <- vapply(models, function(m){
    if(!is.null(fixed_effects)) m$nobs
    else if(family != "linear") m$df.null + 1
    else sum(m$df[1:2])
  }, numeric(1))

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
    retVal$n_obs <- n_obs

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
                            FUN=function(x2) control_extractor(x2, x))


    # Store values in a data frame to be returned
    retVal <- data.frame(coef=unlist(coef), se=unlist(se),
                         statistic=unlist(statistic),
                         p=unlist(p), AIC=unlist(AIC),
                         deviance=unlist(deviance))

    # R doesn't like it when these kinds of objects are assigned above
    retVal$terms <- terms
    retVal$control_coefs <- control_coefs
    retVal$n_obs <- n_obs

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
    # Exact membership against each model's term set. `terms` is a list-column
    # (one character vector of term names per specification), so test `c %in%`
    # the row's terms rather than substring-matching a deparsed string: the
    # latter set false 1s whenever a control name was a substring of another
    # term (e.g. "O2" inside "O2Sat", or a name inside an interaction term).
    retVal[[c]] <- vapply(retVal$terms,
                          function(t) as.integer(c %in% t), integer(1))
  }

  # Remove duplicate columns
  retVal <- retVal %>% select(where(~!all(is.na(.x))))

  # Tag the result so tidy()/glance()/sca_table()/sca_report() can dispatch on
  # it. It remains a data frame in every other respect (no print method, so
  # derived frames behave normally). The focal variable and (normalised) family
  # are stored as attributes so the reporting helpers need not re-derive them
  # (which is ambiguous for a single-control curve).
  attr(retVal, "x") <- x
  attr(retVal, "family") <- family
  class(retVal) <- c("sca", "data.frame")
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

  if(!all(c("coef", "index", "sig.level") %in% names(sca_data))){
    stop("`sca_data` does not look like sca() output (missing one of the ",
         "`coef`, `index`, or `sig.level` columns).", call.=FALSE)
  }

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

  if(!all(c("coef", "index", "sig.level") %in% names(sca_data))){
    stop("`sca_data` does not look like sca() output (missing one of the ",
         "`coef`, `index`, or `sig.level` columns).", call.=FALSE)
  }

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

#' Plots the number of observations across model specifications.
#'
#' @description
#' plot_samplesizes() plots `n_obs`, the number of observations each
#' specification was fit on, against the specification index. It makes visible
#' whether specifications were fit on different samples -- with default listwise
#' deletion they often are -- which is a reason to consider `sca(common_sample
#' = TRUE)`.
#'
#' @inheritParams plot_rmse
#'
#' @return A ggplot object.
#'
#' @seealso [sca()] and its `common_sample` argument.
#'
#' @export
#'
#' @examples
#' plot_samplesizes(sca(y = "Salnty", x = "T_degC",
#'                      controls = c("ChlorA", "O2Sat", "NO2uM"),
#'                      data = bottles, progress_bar = FALSE, parallel = FALSE))
plot_samplesizes <- function(sca_data, title=""){
  df <- as.data.frame(sca_data)
  if(!"n_obs" %in% names(df)){
    stop("`sca_data` has no `n_obs` column; was it produced by sca()?",
         call.=FALSE)
  }
  ggplot(df, aes(x=index, y=n_obs)) +
    geom_col(fill=sca_sig_colors()[["p < .005"]], width=0.8) +
    labs(x="Specification (ranked by estimate)", y="Observations",
         title=title) +
    theme_sca()
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
                     weights=NULL, fe_suffix="", fam_obj=NULL){
  if(length(boot_samples)==1 & length(boot_sample_size)==1){
    samples <- boot_samples
    sample_sizes <- boot_sample_size

    if(is.null(weights)){
      boot <- se_boot(data=data, formula=formula, n_x=n_x,
                      n_samples=boot_samples[[1]],
                      sample_size=boot_sample_size[[1]],
                      fam_obj=fam_obj)
    }
    else{
      boot <- se_boot(data=data, formula=formula, n_x=n_x,
                      n_samples=boot_samples[[1]],
                      sample_size=boot_sample_size[[1]],
                      weights=weights, fam_obj=fam_obj)
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
                     MoreArgs=list(data=data, formula=formula, n_x=n_x,
                                   fam_obj=fam_obj))
    }
    else{
      boot <- mapply(FUN=se_boot, n_samples=samples,
                     sample_size=sample_sizes,
                     MoreArgs=list(data=data, formula=formula, n_x=n_x,
                                   weights=weights, fam_obj=fam_obj))
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
#' @param formula A regression formula, with or without fixed effects, given
#'                either as a string (`"y ~ x | fe"`) or as a formula object
#'                (`y ~ x | fe`).
#' @param data A data frame containing the variables provided in `formula` and
#'             any clustering variables passed to `cluster`.
#' @param weights Optional string with the column name in `data` that contains
#'                weights.
#' @param family A string indicating the family of models to be used. Defaults
#'               to "linear" for OLS regression but supports all families
#'               supported by `glm()`. When a non-linear family is supplied the
#'               models are estimated with `glm()`; fixed effects are not
#'               supported in that case and are ignored with a warning.
#' @param link A string specifying the link function to be used for the model.
#'             Defaults to `NULL`, using OLS regression via `lm()` (or
#'             `fixest::feols()` when fixed effects are supplied). For a
#'             non-linear `family` the canonical link is used when `link` is
#'             `NULL`. Supports all link functions supported by the family
#'             parameter of `glm()`.
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
#'                With clustering: "CL_FE" (standard errors clustered by the
#'                                 first fixed effect), if clusters are supplied
#'                                 then the conventional clustered standard
#'                                 errors from `feols()` are estimated for each
#'                                 clustering specification. Two-way (and
#'                                 multiway) clustering is supported; see the
#'                                 `cluster` argument.
#'
#'                Without clustering: "HC0, "HC1", "HC2", "HC3",
#'                                    "HC4", "HC4m", "HC5",
#'                                    "bootstrapped".
#' @param cluster Variables in `data` to cluster the standard errors on. Either
#'                a character vector, in which case each element is used for a
#'                separate **one-way** clustering, or a list of character
#'                vectors, in which case each element is clustered on **jointly**
#'                (one-way when the element names a single variable, two-way or
#'                higher when it names several). For example
#'                `cluster = list("a", "b", c("a", "b"))` produces one-way SEs
#'                clustered by `a`, one-way by `b`, and two-way clustered by `a`
#'                and `b`. Multiway columns are labelled with the clustering
#'                dimensions joined by `_BY_` (e.g. `"HC1_a_BY_b"`, or
#'                `"CL_a_BY_b_FE"` for a fixed-effects model); the dimensions are
#'                sorted, so the label is the same regardless of the order they
#'                are listed in. Unknown variables are dropped with a warning.
#'                Defaults to `NULL` (no clustering).
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
#'          "CL_FE" = standard errors clustered by the first fixed effect
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
#'          "HC0_Depth_ID_BY_Sta_ID" = HC0 standard errors two-way clustered by
#'                                       "Depth_ID" and "Sta_ID"
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
#' # Two-way (and multiway) clustering: pass a list, where each element names
#' # the dimensions to cluster on jointly. Here: one-way by Sta_ID, one-way by
#' # Depth_ID, and two-way by both.
#' se_compare(formula = "Salnty ~ T_degC + ChlorA", data = bottles,
#'            types = "HC1",
#'            cluster = list("Sta_ID", "Depth_ID", c("Sta_ID", "Depth_ID")))
#'
#' # Logistic regression: compare standard error types for a binary outcome.
#' bottles$saline <- as.integer(bottles$Salnty >
#'                                stats::median(bottles$Salnty, na.rm = TRUE))
#' se_compare(formula = "saline ~ T_degC + ChlorA", data = bottles,
#'            family = "binomial", types = c("iid", "HC0", "HC3"))
#'
se_compare <- function(formula, data, weights=NULL,
                       family="linear", link=NULL,
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

  # `formula` may be supplied either as a string ("y ~ x | fe") or as a formula
  # object (y ~ x | fe). se_compare() works with the string form throughout --
  # grepl()/str_split() detect the fixed-effects pipe and as.formula() refits --
  # so normalise a formula object to a one-line string up front. (as.character()
  # on a formula returns a length-3 vector, "~"/lhs/rhs, which would make the
  # scalar grepl("|", formula) below length 3 and error in the if() that uses it.)
  if(inherits(formula, "formula")){
    formula <- paste(deparse(formula), collapse = " ")
  }

  # Create objects that will store the standard errors
  ses_CL <- NULL
  ses_HC <- NULL
  ses_other <- NULL

  # Create the object we will eventually return
  ses <- NULL

  # Resolve family/link into a normalised family string and (for glm families)
  # a family object. Treats "gaussian" as OLS and errors on an unknown family.
  resolved <- resolve_family(family, link)
  family <- resolved$family
  fam_obj <- resolved$fam_obj
  is_glm <- family != "linear"

  # Whether the formula specifies fixed effects (a pipe). Captured up front
  # because `formula` is later stripped of its fixed effects for the non-FE
  # model. "CL_FE" is a fixed-effects-only type, so the non-FE branch should
  # not flag it as invalid when fixed effects are present.
  has_fe <- grepl("|", formula, fixed=TRUE)

  # Fixed effects are only supported for OLS. For a glm family, mirror sca():
  # warn and drop the fixed effects so estimation falls back to the glm path.
  # Stripping the pipe here means the FE branch below is skipped and "CL_FE"
  # is (correctly) no longer treated as a valid type.
  if(is_glm & has_fe){
    warning(paste0("Fixed effects unsupported for models other than OLS ",
                   "regression. Ignoring fixed effects."))
    formula <- str_trim(str_split(formula, fixed("|"))[[1]][[1]])
    has_fe <- FALSE
  }

  # Validate the columns referenced by the formula and (if supplied) the
  # weights. Clustering variables are validated where they are used, with a
  # warning rather than an error.
  formula_vars <- setdiff(all.vars(stats::as.formula(formula)), ".")
  check_columns(data, formula_vars, "Variable(s)")
  if(!is.null(weights)) check_columns(data, weights, "Weights variable")

  # Normalise the clustering request once, before fitting either model, so any
  # warning about unknown clustering variables fires a single time. `cluster`
  # becomes a cleaned list of clustering specifications (each a character vector
  # of dimensions clustered jointly): a length-1 spec is a one-way clustering
  # (the historical behaviour, preserved byte-for-byte) and a longer spec is a
  # multiway clustering. See normalize_cluster_spec().
  cluster <- normalize_cluster_spec(cluster, colnames(data))

  # If the formula contains a pipe then fixed effects are assumed to be
  # present and models are estimated with feols() rather than lm()
  if(grepl("|", formula, fixed=TRUE)){

    if(is.null(weights)){
      model_fe <- tryCatch(feols(as.formula(formula), data=data),
                           error=function(cond){
                             message("Fixed effects model estimation failed: ",
                                     conditionMessage(cond))
                             return(NULL)
                           })
    }
    else{
      model_fe <- tryCatch(
        feols(as.formula(formula), data=data, weights=data[[weights]]),
                           error=function(cond){
                             message("Fixed effects model estimation failed: ",
                                     conditionMessage(cond))
                             return(NULL)
                           })
    }

      if(!is.null(model_fe)){
      # The first fixed effect, used below to compute "CL_FE". Historically
      # feols() clustered its default standard errors by the first fixed
      # effect, but modern fixest (>= 0.10) defaults to IID, so we cluster by
      # the first fixed effect explicitly to keep "CL_FE" cluster-robust as its
      # name and documentation promise, regardless of fixest version.
      fe_part <- str_trim(str_split(formula, fixed("|"))[[1]][[2]])
      first_fe <- str_trim(str_split(fe_part, fixed("+"))[[1]][[1]])

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

        # Standard errors clustered by the first fixed effect. When the first
        # fixed effect is a plain column it is used as the clustering variable;
        # otherwise (e.g. an interaction fixed effect such as "a^b") we fall
        # back to feols()'s default standard errors.
        if("CL_FE" %in% types_other){
          cl_fe_se <- if(first_fe %in% colnames(data)){
            summary(model_fe, cluster=data[first_fe])$coeftable[,2]
          } else {
            coeftable(model_fe)[,2]
          }
          ses_other <- cbind(ses_other, "CL_FE"=c("(Intercept)"=NA, cl_fe_se))
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

      # Estimate clustered standard errors for the FE model. `cluster` is the
      # cleaned list of specifications: a one-element spec gives feols' one-way
      # cluster-robust SE (unchanged), a multi-element spec gives joint multiway
      # clustering (Cameron-Gelbach-Miller), which fixest computes natively when
      # passed the multi-column data frame `data[dims]`. These SEs are feols'
      # cluster-robust SEs and do not depend on `types`, so there is no type
      # validation here.
      if(!is.null(cluster)){

        cl_cols <- list()
        for(dims in cluster){
          key <- paste(dims, collapse="_BY_")
          cl_cols[[paste0("CL_", key, "_FE")]] <-
            feols(as.formula(formula), data=data,
                  cluster=data[dims])$coeftable[,2]
        }

        ses_CL <- do.call(cbind, cl_cols)
        colnames(ses_CL) <- names(cl_cols)

        ses_CL <- rbind("(Intercept)"=NA, ses_CL)

        ses <- cbind(ses, ses_CL)
      }
    }

  }
  # Case when a non-FE model is desired. A glm family always estimates a non-FE
  # model (any fixed effects were dropped above), so the branch is forced for
  # glm even when `fixed_effects_only` was requested.
  if(!fixed_effects_only | is_glm){

    # Allocate objects to hold SEs
    ses_other <- NULL
    ses_HC <- NULL
    ses_CL <- NULL

    # If the formula has FEs remove them
    if(grepl("|",formula, fixed=TRUE)){
      formula <- str_trim(str_split(formula, fixed("|"))[[1]][[1]])
    }

    # Estimate the non-FE model and get the coefficients. A non-linear family
    # uses glm() with the resolved family object; otherwise lm(). The HC and
    # clustered SE machinery below (coeftest + vcovHC/vcovCL) accepts glm and
    # lm objects alike, so only the model-fitting call differs.
    if(is.null(weights)){
      model <- if(is_glm) glm(formula=as.formula(formula), data=data,
                              family=fam_obj)
               else lm(formula=as.formula(formula), data=data)
    }
    else{
      fmla <- as.formula(formula)
      environment(fmla) <- environment()
      model <- if(is_glm) glm(formula=fmla, data=data, family=fam_obj,
                              weights=get(weights))
               else lm(formula=fmla, data=data, weights=get(weights))
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
                         weights=weights, fe_suffix="", fam_obj=fam_obj)

        if(!is.null(boot)){
          ses_other <- cbind(ses_other, boot)
        }
      }

      # Attach SEs to the return object
      ses <- cbind(ses, ses_other, ses_HC)
    }

    # Case when clustered SEs for non-FE model are desired. `cluster` is the
    # cleaned list of specifications.
    if(!is.null(cluster)){

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

      # Estimate and extract clustered SEs. A one-element spec reduces exactly to
      # the historical one-way clustering; a multi-element spec gives joint
      # multiway clustering (Cameron-Gelbach-Miller), which sandwich::vcovCL
      # computes natively from the multi-column data frame. Passing data[dims]
      # as a data FRAME (never a bare vector) is load-bearing: vcovCL aligns it
      # to the rows the model actually used (via na.action), so specifications
      # fit on NA-dropped data stay correct.
      if(length(types_CL)>0){

        cl_cols <- list()
        for(dims in cluster){
          key <- paste(dims, collapse="_BY_")
          for(t in types_CL){
            # Some clustered HC variants (notably HC3, whose leverage adjustment
            # divides by 1 - h_ii) are numerically unstable and can fail with a
            # raw LAPACK "singular" error on certain designs. Catch that per
            # type/cluster combination: warn clearly and skip just that column,
            # so the other requested standard errors are still returned.
            se <- tryCatch(
              coeftest(model, vcov.=vcovCL, type=t, cluster=data[dims])[,2],
              error = function(e){
                warning(t, " standard errors clustered by ",
                        gsub("_BY_", " + ", key, fixed=TRUE),
                        " could not be computed and were skipped (",
                        conditionMessage(e), ").", call.=FALSE)
                NULL
              })
            if(!is.null(se)) cl_cols[[paste0(t, "_", key)]] <- se
          }
        }

        if(length(cl_cols) > 0){
          ses_CL <- do.call(cbind, cl_cols)
          colnames(ses_CL) <- names(cl_cols)
        }
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
