# User-facing functions---------------------------------------------------------

#' Perform specification curve analysis
#'
#' @description
#' sca() is the workhorse function of the package--this estimates models with every
#' possible combination of the controls supplied and returns a data frame
#' where each row contains the pertinent information and parameters for a
#' given model by default. This data frame can then be input to plotCurve()
#' or any other plotting function in the package. Alternatively, if
#' `returnFormulae = TRUE`, it returns a list of formula objects with every
#' possible combination of controls.
#'
#' @param y A string containing the column name of the dependent variable in
#'          data.
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
#' @param fixedEffects A string containing the column name of the variable
#'                     in data desired for fixed effects. Defaults to NULL in
#'                     which case no fixed effects are included.
#' @param returnFormulae A boolean. When `TRUE` a list of model formula objects
#'                       is returned but the models are not estimated. Defaults
#'                       to `FALSE` in which case a dataframe of model results
#'                       is returned.
#' @param progressBar A boolean indicating whether the user wants a progress bar
#'                    for model estimation. Defaults to `TRUE`.
#' @param parallel A boolean indicating whether to parallelize model estimation.
#'                 Parallelization only offers a speed advantage when a large
#'                 (> 1000) number of models is being estimated. Defaults to
#'                 `FALSE`.
#' @param workers An integer indicating the number of workers to use for
#'                parallelization. Defaults to 2.
#'
#' @return When `returnFormulae` is `FALSE`, a dataframe where each row contains
#'         the independent variable coefficient estimate, standard error,
#'         test statistic, p-value, model specification, and measures of model
#'         fit.
#'
#' @export
#'
#' @examples
#' sca(y = "Salnty", x = "T_degC", controls = c("ChlorA", "O2Sat"),
#'     data = bottles, progressBar = TRUE, parallel = FALSE);
#' \donttest{
#' sca(y = "Salnty", x = "T_degC", controls = c("ChlorA*NO3uM", "O2Sat*NO3uM"),
#'     data = bottles, progressBar = TRUE, parallel = TRUE, workers = 2);
#' }
#' sca(y = "Salnty", x = "T_degC", controls = c("ChlorA", "O2Sat*NO3uM"),
#'     data = bottles, progressBar = TRUE, parallel = FALSE,
#'     returnFormulae = TRUE);
sca <- function(y, x, controls, data, weights=NULL,
                family="linear", link=NULL,
                fixedEffects=NULL, returnFormulae=FALSE,
                progressBar=TRUE, parallel=FALSE, workers=2){

  # Treat the common alias "gaussian" as ordinary least squares.
  if(family=="gaussian") family <- "linear"

  if(family!="linear" & !is.null(fixedEffects))
  {
    warning(paste0("Fixed effects unsupported for models other than OLS ",
                   "regression. Ignoring fixed effects."))
    # Actually drop the fixed effects so downstream estimation/extraction uses
    # the glm path, as the warning promises.
    fixedEffects <- NULL
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
  if(returnFormulae){
    if(!is.null(fixedEffects)){
      return(formula_builder(y=y, x=x, controls=controls,
                             fixedEffects=fixedEffects))
    }
    else{
      return(formula_builder(y=y, x=x, controls=controls))
    }
  }

  # Validate that the requested columns exist in `data` (interaction syntax in
  # x/controls is split so each underlying variable is checked).
  vars <- unique(trimws(unlist(strsplit(c(y, x, controls), "[*:]"))))
  check_columns(data, vars, "Variable(s)")
  if(!is.null(fixedEffects)){
    check_columns(data, fixedEffects, "Fixed-effects variable(s)")
  }
  if(!is.null(weights)) check_columns(data, weights, "Weights variable")

  # Build the model formulae (with or without fixed effects)
  if(is.null(fixedEffects)){
    formulae <- formula_builder(y=y, x=x, controls=controls)
  }
  else{
    formulae <- formula_builder(y=y, x=x, controls=controls,
                                fixedEffects=fixedEffects)
  }

  # Estimator for a single specification. Dispatches on fixed effects, family,
  # and weights and returns a model summary. Defined as a closure so that, when
  # estimation is parallelised, it carries `data`, `weights`, `family`,
  # `fam_obj`, and `fixedEffects` to the workers along with the function.
  estimate_one <- function(f){
    if(!is.null(fixedEffects)){
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

    if(progressBar){
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
    if(progressBar){
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
    if(is.null(fixedEffects)){

      # Get each value of interest across models
      coef <- lapply(X=models, function(x2) x2$coefficients[x,1])
      se <- lapply(X=models, function(x2) x2$coefficients[x,2])
      statistic <- lapply(X=models, function(x2) x2$coefficients[x,3])
      p <- lapply(X=models, function(x2) x2$coefficients[x,4])
      terms <- lapply(X=models, FUN=function(x2) row.names(x2$coefficients))
      RMSE <- lapply(X=models, FUN=function(x2) sqrt(mean(x2$residuals^2)))
      adjR <- lapply(X=models, function(x2) x2$adj.r.squared)
      control_coefs <- lapply(X=models,
                              FUN=function(x2, x3) controlExtractor(x2,x3),
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
                                controlExtractor(x2, x3, feols_model=TRUE),
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
                                         x4) controlExtractor(x2,x3),x3=x)


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
#' plotCurve() takes the data frame output of sca() and produces a ggplot of
#' the independent variable's coefficient (as indicated in the call to sca())
#' across model specifications. By default a panel is added showing which
#' control variables are present in each model. Note that the ggplot output by
#' this function can only be further customized when `plotVars = FALSE`, i.e.
#' when the control variable panel is not included.
#'
#' @param sca_data A data frame returned by `sca()` containing model estimates
#'                 from the specification curve analysis.
#' @param title A string to use as the plot title. Defaults to an empty string,
#'              `""`.
#' @param showIndex A boolean indicating whether to label the model index on the
#'                  the x-axis. Defaults to `TRUE`.
#' @param plotVars A boolean indicating whether to include a panel on the plot
#'                 showing which variables are present in each model. Defaults
#'                 to `TRUE`.
#' @param ylab A string to be used as the y-axis label. Defaults to
#'             `"Coefficient"`.
#' @param plotSE A string indicating whether to display standard errors as
#'               bars or plots. For bars `plotSE = "bar"`, for ribbons
#'               `plotSE = "ribbon"`. If any other value is supplied then no
#'               standard errors are included. Defaults to `"bar"`.
#'
#' @return If `plotVars = TRUE` returns a grid grob (i.e. the output of a call
#'         to `grid.draw`). If `plotVars =  FALSE` returns a ggplot object.
#'
#' @export
#'
#' @examples
#' plotCurve(sca_data = sca(y="Salnty", x="T_degC", c("ChlorA", "O2Sat"),
#'                          data=bottles, progressBar=TRUE, parallel=FALSE),
#'                      title = "Salinity and Temperature Models",
#'                      showIndex = TRUE, plotVars = TRUE,
#'                      ylab = "Coefficient value", plotSE = "ribbon");
#' plotCurve(sca_data = sca(y="Salnty", x="T_degC",
#'                          c("ChlorA*O2Sat", "ChlorA", "O2Sat"),
#'                          data=bottles, progressBar=FALSE, parallel=FALSE),
#'                      showIndex = TRUE, plotVars = TRUE,
#'                      plotSE = "ribbon");
#' \donttest{
#' plotCurve(sca_data = sca(y="Salnty", x="T_degC",
#'                          c("ChlorA*NO3uM", "O2Sat", "ChlorA", "NO3uM"),
#'                          data=bottles,
#'                          progressBar = TRUE, parallel = TRUE, workers=2),
#'           plotSE="");
#' }
plotCurve <- function(sca_data, title="", showIndex=TRUE, plotVars=TRUE,
                         ylab="Coefficient", plotSE="bar"){

  if("control_coefs" %in% names(sca_data)){
    sca_data <- sca_data %>% select(-control_coefs)
  }

  pointSize <- spec_point_size(sca_data)

  if(tolower(plotSE)=="ribbon"){
    sca_data <- sca_data %>%
      mutate(ribbon.group = cumsum(sig.level != stats::lag(sig.level,
                                                    def = first(sig.level))))
  }

  margin <- {if(title=="") unit(c(-15,2,-5,2), "points")
             else unit(c(5,2,-5,2), "points")}

  sc1 <- ggplot(data=sca_data, aes(y=coef, x=index)) +
    geom_hline(yintercept = 0, color="red", linetype="dashed", linewidth=.75) +
    {if(plotSE=="ribbon") geom_ribbon(aes(ymin=coef-se, ymax=coef+se,
                                           group=factor(ribbon.group),
                                           fill=factor(sig.level)),
                                       alpha=.4)} +
    {if(tolower(plotSE)=="bar") geom_errorbar(aes(ymin=coef-se, ymax=coef+se,
                                          color=factor(sig.level)),
                                      width=0.25)} +
    {if(!tolower(plotSE) %in% c("ribbon",
                       "bar")) geom_point(aes(color=as.factor(sig.level)),
                                          size=pointSize)} +
    {if(tolower(plotSE) %in% c("ribbon",
                      "bar")) geom_point(color="black",size=pointSize)} +
    labs(title=title, x="", y=ylab) +
    theme_bw() +
    theme(
      axis.text.x = {if(showIndex) element_text()
                      else element_blank()},
      axis.title.y = element_text(vjust=-0.5),
      legend.position="top",
      legend.title=element_blank(),
      plot.margin = {if(title=="") unit(c(-15,1,-5,1), "points")
                     else unit(c(5,1,-5,1), "points")}
    ) +
    guides(color = guide_legend(override.aes = list(size=2))) +
    guides(fill = guide_legend(override.aes = list(size=2)))


  if(plotVars){
    sc2 <- plotVars(sca_data)

    grid::grid.newpage()

    return(grid::grid.draw(rbind(ggplotGrob(sc1), ggplotGrob(sc2))))
  }
  else{
    return(sc1)
  }
}

#' Plots the variables in each model.
#'
#' @description
#' plotVars() plots the variables included in each model specification in order
#' of model index. Returns a ggplot object that can then be combined with the
#' output of other functions like plotRMSE() if further customization of each
#' plot is desired.
#'
#' @inheritParams plotCurve
#' @param colorControls A boolean indicating whether to give each variable a
#'                      color to improve readability. Defaults to `FALSE`.
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @examples
#' plotVars(sca_data = sca(y = "Salnty", x = "T_degC",
#'                         controls = c("ChlorA", "O2Sat"),
#'                         data = bottles, progressBar = TRUE,
#'                         parallel = FALSE),
#'                      title = "Model Variable Specifications");
#' plotVars(sca_data = sca(y = "Salnty", x = "T_degC",
#'                         controls = c("ChlorA*O2Sat"),
#'                         data = bottles, progressBar = FALSE,
#'                         parallel = FALSE),
#'                      colorControls = TRUE);
#' \donttest{
#' plotVars(sca_data = sca(y = "Salnty", x = "T_degC",
#'                         controls = c("ChlorA*NO3uM", "O2Sat*NO3uM"),
#'                         data = bottles,
#'                         progressBar = TRUE, parallel = TRUE, workers = 2));
#' }
plotVars <- function(sca_data, title="", colorControls=FALSE){

  if("control_coefs" %in% names(sca_data)){
    sca_data <- sca_data %>% select(-control_coefs)
  }

  scp_data <- scp(sca_data)

  markSize <- 10/length(scp_data[[2]])

  margin <- {if(title=="") unit(c(-5,2,-5,2), "points")
             else unit(c(5,2,-5,2), "points")}

  if(colorControls){
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

# Internal: shared implementation behind plotRMSE()/plotR2Adj()/plotAIC()/
# plotDeviance(), which differ only in the metric column, the axis label, and
# the message shown when that metric is absent from the sca() output.
plot_metric <- function(sca_data, metric, ylab, missing_message,
                        title="", showIndex=TRUE, plotVars=TRUE){

  if(!metric %in% colnames(sca_data)){
    message(missing_message)
    return(invisible(NULL))
  }

  sca_data <- sca_data %>% select(-control_coefs)

  pointSize <- spec_point_size(sca_data)

  margin <- {if(title=="") unit(c(-5,2,-5,2), "points")
             else unit(c(5,2,-5,2), "points")}

  sc1 <- ggplot(data=sca_data, aes(x=.data$index, y=.data[[metric]])) +
    geom_point(size=pointSize) +
    labs(title=title, x="", y=ylab) +
    theme_bw() +
    theme(
      axis.text.x = {if(showIndex) element_text()
                     else element_blank()},
      legend.title=element_blank(),
      legend.key.size = unit(.4, 'cm'),
      plot.margin = margin
    )

  if(plotVars){
    sc2 <- plotVars(sca_data)

    grid::grid.newpage()

    return(grid::grid.draw(rbind(ggplotGrob(sc1), ggplotGrob(sc2))))
  }
  else{
    return(sc1)
  }
}

#' Plots RMSE across model specifications.
#'
#' @description
#' plotRMSE() plots the root mean square error across model specifications. Only
#' available for linear regression models.
#'
#' @inheritParams plotCurve
#' @param showIndex A boolean indicating whether to label the model index on the
#'                  the x-axis. Defaults to `TRUE`.
#' @param plotVars A boolean indicating whether to include a panel on the plot
#'                 showing which variables are present in each model. Defaults
#'                 to `TRUE`.
#'
#' @return If `plotVars = TRUE` returns a grid grob (i.e. the output of a call
#'         to `grid.draw`). If `plotVars =  FALSE` returns a ggplot object.
#'
#' @export
#'
#' @examples
#' plotRMSE(sca_data = sca(y="Salnty", x="T_degC", c("ChlorA", "O2Sat"),
#'                          data=bottles, progressBar=TRUE, parallel=FALSE),
#'                      title = "RMSE");
#' plotRMSE(sca_data = sca(y="Salnty", x="T_degC", c("ChlorA*O2Sat"),
#'                          data=bottles, progressBar=FALSE, parallel=FALSE),
#'                      showIndex = FALSE, plotVars = FALSE);
#' \donttest{
#' plotRMSE(sca_data = sca(y="Salnty", x="T_degC",
#'                          c("ChlorA*NO3uM", "O2Sat*NO3uM"), data=bottles,
#'                          progressBar = TRUE, parallel=TRUE, workers=2));
#' }
plotRMSE <- function(sca_data, title="", showIndex=TRUE, plotVars=TRUE){
  plot_metric(sca_data, metric="RMSE", ylab="RMSE",
              missing_message=paste0("RMSE not found. Are your models nonlinear? ",
                                     "Try plotAIC() or plotDeviance() instead."),
              title=title, showIndex=showIndex, plotVars=plotVars)
}

#' Plots the adj. R-squared across model specifications.
#'
#' @description
#' plotR2Adj() plots the adjusted R-squared across model specifications. Only
#' available for linear regression models. Note when fixed effects are
#' are specified the within adjusted R-squared is used (i.e. `fixest::r2()`
#' with `type="war2"`).
#'
#' @inheritParams plotRMSE
#'
#' @return If `plotVars = TRUE` returns a grid grob (i.e. the output of a call
#'         to `grid.draw`). If `plotVars =  FALSE` returns a ggplot object.
#'
#' @export
#'
#' @examples
#' plotR2Adj(sca_data = sca(y = "Salnty", x = "T_degC",
#'                          controls = c("ChlorA", "O2Sat"),
#'                          data = bottles, progressBar = TRUE,
#'                          parallel = FALSE),
#'                      title = "Adjusted R^2");
#' plotR2Adj(sca_data = sca(y="Salnty", x="T_degC",
#'                          controls = c("ChlorA*O2Sat"),
#'                          data = bottles, progressBar = FALSE,
#'                          parallel = FALSE),
#'                      showIndex = FALSE, plotVars = FALSE);
#' \donttest{
#' plotR2Adj(sca_data = sca(y = "Salnty", x = "T_degC",
#'                          controls = c("ChlorA*NO3uM", "O2Sat*NO3uM"),
#'                          data = bottles,
#'                          progressBar = TRUE, parallel = TRUE, workers = 2));
#' }
plotR2Adj <- function(sca_data, title="", showIndex=TRUE, plotVars=TRUE){
  plot_metric(sca_data, metric="adjR", ylab=bquote('Adj. R'^2),
              missing_message=paste0("Adj. R^2 not found. Are your models nonlinear? ",
                                     "Try plotAIC() or plotDeviance() instead."),
              title=title, showIndex=showIndex, plotVars=plotVars)
}

#' Plots the AIC across model specifications.
#'
#' @description
#' plotAIC() plots the Akaike information criterion across model specifications.
#' Only available for nonlinear regression models.
#'
#' @inheritParams plotRMSE
#'
#' @return If `plotVars = TRUE` returns a grid grob (i.e. the output of a call
#'         to `grid.draw`). If `plotVars =  FALSE` returns a ggplot object.
#'
#' @export
#'
#' @examples
#' plotAIC(sca_data = sca(y = "Salnty", x = "T_degC",
#'                        controls = c("ChlorA", "O2Sat"),
#'                        data = bottles, progressBar = TRUE, parallel = FALSE),
#'                      title = "AIC");
#' plotAIC(sca_data = sca(y = "Salnty", x = "T_degC",
#'                        controls = c("ChlorA*O2Sat"),
#'                        data = bottles, progressBar = FALSE,
#'                        parallel = FALSE),
#'                        showIndex = FALSE, plotVars = FALSE);
#' \donttest{
#' plotAIC(sca_data = sca(y = "Salnty", x = "T_degC",
#'                          controls = c("ChlorA*NO3uM", "O2Sat*NO3uM"),
#'                          data = bottles,
#'                          progressBar = TRUE, parallel = TRUE, workers = 2));
#' }
plotAIC <- function(sca_data, title="", showIndex=TRUE, plotVars=TRUE){
  plot_metric(sca_data, metric="AIC", ylab="AIC",
              missing_message=paste0("AIC not found. Are your models linear? ",
                                     "Try plotR2Adj() or plotRMSE instead."),
              title=title, showIndex=showIndex, plotVars=plotVars)
}

#' Plots the deviance of residuals across model specifications.
#'
#' @description
#' plotDeviance() plots the deviance of residuals across model specifications.
#' Only available for linear regression models.
#'
#' @inheritParams plotRMSE
#'
#' @return If `plotVars = TRUE` returns a grid grob (i.e. the output of a call
#'         to `grid.draw`). If `plotVars =  FALSE` returns a ggplot object.
#'
#' @export
#'
#' @examples
#' plotDeviance(sca_data = sca(y = "Salnty", x = "T_degC",
#'                             controls = c("ChlorA", "O2Sat"),
#'                             data = bottles, progressBar = TRUE,
#'                             parallel = FALSE),
#'                      title = "Model Deviance");
#' plotDeviance(sca_data = sca(y = "Salnty", x = "T_degC",
#'                             controls = c("ChlorA*O2Sat"),
#'                             data = bottles, progressBar = FALSE,
#'                             parallel = FALSE),
#'                      showIndex = FALSE, plotVars = FALSE);
#' \donttest{
#' plotDeviance(sca_data = sca(y = "Salnty", x="T_degC",
#'                          controls = c("ChlorA*NO3uM", "O2Sat*NO3uM"),
#'                          data = bottles, progressBar = TRUE, parallel = TRUE,
#'                          workers = 2));
#' }
plotDeviance <- function(sca_data, title="", showIndex=TRUE, plotVars=TRUE){
  plot_metric(sca_data, metric="deviance", ylab="Deviance",
              missing_message=paste0("Deviance of residuals not found. ",
                                     "Are your models linear? Try plotAIC(), ",
                                     "plotR2Adj(), or plotRMSE() instead."),
              title=title, showIndex=showIndex, plotVars=plotVars)
}

#' Plots control variable distributions.
#'
#' @description
#' plotControlDistributions() plots the distribution of coefficients for each
#' control variable included in the model specifications.
#'
#' @inheritParams plotRMSE
#' @param type A string indicating what type of distribution plot to produce.
#'             When `type = "density"` density plots are produced. When
#'             `type = "hist"` or `type = "histogram"` histograms are produced.
#'             Defaults to `"density"`.
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @examples
#' plotControlDistributions(sca_data = sca(y="Salnty", x="T_degC",
#'                                     controls = c("ChlorA", "O2Sat"),
#'                                     data = bottles,
#'                                     progressBar = TRUE, parallel = FALSE),
#'                          title = "Control Variable Distributions")
#' plotControlDistributions(sca_data = sca(y = "Salnty", x="T_degC",
#'                                     controls = c("ChlorA*O2Sat"),
#'                                     data = bottles,
#'                                     progressBar = FALSE, parallel = FALSE),
#'                          type = "hist")
#' \donttest{
#' plotControlDistributions(sca_data = sca(y = "Salnty", x = "T_degC",
#'                                     controls = c("ChlorA*NO3uM",
#'                                                  "O2Sat*NO3uM"),
#'                                     data = bottles, progressBar = TRUE,
#'                                     parallel = TRUE, workers = 2),
#'                          type = "density")
#' }
plotControlDistributions <- function(sca_data, title="", type="density"){

  histData <- bind_rows(unAsIs(sca_data$control_coefs))

  rownames(histData) <- NULL

  histData$term <- sapply(sapply(histData$term, str_split, pattern=":"),
                          paste0, collapse=" %*% ")

  n_facets <- length(unique(histData$term))

  sc1 <- histData %>%
    ggplot(aes(x=coef, fill=factor(term))) +
      {if(tolower(type)=="hist" | tolower(type)=="histogram") geom_histogram()
       else if (tolower(type)=="density") geom_density()} +
      labs(x="", y="", title=title) +
      theme_bw() +
      theme(
        legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
      ) +
    {if(n_facets>16) theme(axis.text.x=element_text(size=4),
                           axis.text.y=element_text(size=4),
                           strip.text=element_text(size=6))
      else theme()} +
      facet_wrap(~factor(term), scales="free", labeller = label_parsed)
  return(sc1)
}


# Internal: bootstrapped SEs for se_compare(). Estimates one column per
# (bootSamples, bootSampleSize) combination and returns a named matrix (rows =
# coefficients), or NULL if no estimates were produced. `fe_suffix` tags the
# column names for fixed-effects models.
boot_ses <- function(data, formula, n_x, bootSamples, bootSampleSize,
                     weights=NULL, fe_suffix=""){
  if(length(bootSamples)==1 & length(bootSampleSize)==1){
    samples <- bootSamples
    sample_sizes <- bootSampleSize

    if(is.null(weights)){
      boot <- se_boot(data=data, formula=formula, n_x=n_x,
                      n_samples=bootSamples[[1]],
                      sample_size=bootSampleSize[[1]])
    }
    else{
      boot <- se_boot(data=data, formula=formula, n_x=n_x,
                      n_samples=bootSamples[[1]],
                      sample_size=bootSampleSize[[1]],
                      weights=weights)
    }

    if(is.null(boot)) return(NULL)
    boot <- matrix(boot, ncol=1, dimnames=list(names(boot), NULL))
  }
  else{
    samples <- rep(bootSamples, length(bootSampleSize))
    sample_sizes <- sort(rep(bootSampleSize, length(bootSamples)))

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
#' @param clusteredOnly A boolean indicating whether only standard errors with
#'                      clustering should be estimated, defaults to `FALSE`.
#' @param fixedEffectsOnly A boolean indicating whether only standard errors for
#'                         fixed effects models should be estimated, defaults to
#'                         `FALSE`.
#' @param bootSamples An integer or vector of integers indicating how many times
#'                    the model should be estimated with a random subset of the
#'                    data. If a vector then every combination of `bootSamples`
#'                    and `bootSampleSize` are estimated.
#' @param bootSampleSize An integer or vector of integers indicating how many
#'                       observations are in each random subset of the data.
#'                       If a vector then every combination of `bootSamples`
#'                       and `bootSampleSize` are estimated.
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
#'                                   effects model where `bootSamples = 8` and
#'                                   `bootSampleSize = 300`
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
#'            fixedEffectsOnly = FALSE, bootSamples=c(4, 8, 10),
#'            bootSampleSize=c(300, 500))
#'
#' se_compare(formula = "Salnty ~ T_degC + ChlorA + O2Sat", data = bottles,
#'            types = "bootstrapped", bootSamples = c(8, 10),
#'            bootSampleSize = c(300, 500))
#'
#' se_compare(formula = "Salnty ~ T_degC + ChlorA", data = bottles,
#'            types = c("HC0", "HC1", "HC3"))
#'
se_compare <- function(formula, data, weights=NULL,
                       types="all", cluster=NULL,
                       clusteredOnly=FALSE, fixedEffectsOnly=FALSE,
                       bootSamples=NULL, bootSampleSize=NULL){

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
      if(!clusteredOnly){
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
        if("bootstrapped" %in% types_other & !is.null(bootSamples) &
           !is.null(bootSampleSize)){

          boot <- boot_ses(data=data, formula=formula,
                           n_x=length(model_fe$coefficients),
                           bootSamples=bootSamples, bootSampleSize=bootSampleSize,
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
  if(!fixedEffectsOnly){

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
    if(!clusteredOnly){
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
      if("bootstrapped" %in% types_other & !is.null(bootSamples) &
         !is.null(bootSampleSize)){

        boot <- boot_ses(data=data, formula=formula,
                         n_x=length(model$coefficients)-1,
                         bootSamples=bootSamples, bootSampleSize=bootSampleSize,
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
#' plotSE() takes the data frame output of `se_compare()` and plots, for each
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
#' plotSE(se_compare(formula = "Salnty ~ T_degC + ChlorA", data = bottles,
#'                   types = c("iid", "HC0", "HC3")));
#' plotSE(se_compare(formula = "Salnty ~ T_degC + ChlorA", data = bottles,
#'                   types = "HC1", cluster = c("Sta_ID", "Depth_ID")),
#'        level = 0.9);
plotSE <- function(se_data, level=0.95, intercept=FALSE, title=""){

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
    scale_color_manual(values=c("FALSE"="grey50", "TRUE"="black"),
                       labels=c("FALSE"="CI includes 0",
                                "TRUE"="CI excludes 0"),
                       drop=FALSE) +
    labs(title=title, x="Standard error type", y="Estimate", color="") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          legend.position="top")
}
