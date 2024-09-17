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
#' sca(y = "Salnty", x = "T_degC", controls = c("ChlorA*NO3uM", "O2Sat*NO3uM"),
#'     data = bottles, progressBar = TRUE, parallel = TRUE, workers = 2);
#' sca(y = "Salnty", x = "T_degC", controls = c("ChlorA", "O2Sat*NO3uM"),
#'     data = bottles, progressBar = TRUE, parallel = FALSE,
#'     returnFormulae = TRUE);
sca <- function(y, x, controls, data, family="linear", link=NULL,
                fixedEffects=NULL, returnFormulae=FALSE,
                progressBar=TRUE, parallel=FALSE, workers=2){

  if(family!="linear" & !is.null(fixedEffects))
    {
    warning(paste0("Fixed effects unsupported for models other than OLS ",
                   "regression. Ignoring fixed effects."))
  }

  # General family argument for glm
  if(family!="linear"){
    family_link <- paste0(family, "(link=\"", link, "\")" )
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

  # With parallel computing
  if(parallel){

    cl <- makePSOCKcluster(rep("localhost", workers))

    # Load needed package into each cluster
    clusterEvalQ(cl, library(fixest))

    # No fixed effects specified
    if(is.null(fixedEffects)){
      # Build the formulae
      formulae <- formula_builder(y=y, x=x, controls=controls)

      # Estimate the models with lm()
      clusterExport(cl, "formulae", envir=environment())
      clusterExport(cl, "data", envir=environment())

      # Show progress bar if desired
      if(progressBar){
        print.noquote(paste("Estimating", length(formulae),
                            "models in parallel with",
                    workers, "workers"))

        if(family=="linear"){
          system.time(models <- pblapply(
            formulae, function(x2) summary(lm(x2, data=data)), cl=cl))
        }
        else{
          system.time(models <- pblapply(
            formulae, function(x2) summary(
              glm(x2, data=data, family=eval(parse(text=family_link)))),
                                         cl=cl))
        }
      }
      else{
        if(family=="linear"){
          models <- parLapply(
            cl, formulae, function(x2) summary(lm(x2, data=data)))
        }
        else{
          models <- parLapply(
            cl, formulae, function(x2) summary(
              glm(x2, data=data, family=eval(parse(text=family_link)))))
        }
      }
    }
    # Fixed effects specified
    else{
      # Build the formulae
      formulae <- formula_builder(y=y, x=x, controls=controls,
                                       fixedEffects=fixedEffects)

      clusterExport(cl, "formulae", envir=environment())
      clusterExport(cl, "data", envir=environment())

      if(progressBar){
        print.noquote(paste("Estimating", length(formulae),
                            "models in parallel with",
                    workers, "workers"))
        system.time(models <- pblapply(formulae,
                                       function(x2) summary(feols(x2,data=data)),
                                       cl=cl))
      }
      else{
        models <- parLapply(cl, formulae,
                            function(x2) summary(feols(x2, data=data)))
      }
    }
  }
  # Without parallel computing
  else{
    # No fixed effects specified
    if(is.null(fixedEffects)){
      # Build the formulae
      formulae <- formula_builder(y=y, x=x, controls=controls)

      if(progressBar){
        print.noquote(paste("Estimating", length(formulae), "models"))
        if(family=="linear"){
          system.time(models <- pblapply(
            formulae, function(x2) summary(lm(x2, data=data))))
        }
        else{
          system.time(models <- pblapply(
            formulae, function(x2) summary(
              glm(x2, data=data, family=eval(parse(text=family_link))))))
        }
      }
      else{
        if(family=="linear"){
          models <- lapply(formulae, function(x2) summary(lm(x2, data=data)))
        }
        else{
          models <- lapply(formulae,
                           function(x2) summary(
                             glm(x2, data=data,
                                 family=eval(parse(text=family_link)))))
        }
      }
    }
    # Fixed effects specified
    else{
      # Build the formulae
      formulae <- formula_builder(y=y, x=x, controls=controls,
                                       fixedEffects=fixedEffects)

      if(progressBar){
        print.noquote(paste("Estimating", length(formulae), "models"))
        system.time(models <- pblapply(
          X=formulae, function(x2) summary(feols(x2, data=data))))
      }
      else{
        models <- lapply(X=formulae, function(x2) summary(feols(x2, data=data)))
      }
    }
  }

  # Garbage collection for parallel connections
  if(parallel) stopCluster(cl=cl)

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
                                                        verbose=F)[[1]])
      adjR <- lapply(X=models, FUN=function(x2) fitstat(x2, type="war2",
                                                        verbose=F)[[1]])
      control_coefs <- lapply(X=models,
                              FUN=function(x2,x3)
                                controlExtractor(x2, x3, feols_model=T),
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
          T ~ NA_character_
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
          T ~ NA_character_
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
#' plotCurve(sca_data = sca(y="Salnty", x="T_degC",
#'                          c("ChlorA*NO3uM", "O2Sat", "ChlorA", "NO3uM"),
#'                          data=bottles,
#'                          progressBar = TRUE, parallel = TRUE, workers=2),
#'           plotSE="");
plotCurve <- function(sca_data, title="", showIndex=TRUE, plotVars=TRUE,
                         ylab="Coefficient", plotSE="bar"){

  if("control_coefs" %in% names(sca_data)){
    sca_data <- sca_data %>% select(-control_coefs)
  }

  pointSize <- -.25*(ncol(sca_data)-7)+(13/4)

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
#' plotVars(sca_data = sca(y = "Salnty", x = "T_degC",
#'                         controls = c("ChlorA*NO3uM", "O2Sat*NO3uM"),
#'                         data = bottles,
#'                         progressBar = TRUE, parallel = TRUE, workers = 2));
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
#' plotRMSE(sca_data = sca(y="Salnty", x="T_degC",
#'                          c("ChlorA*NO3uM", "O2Sat*NO3uM"), data=bottles,
#'                          progressBar = TRUE, parallel=TRUE, workers=2));
plotRMSE <- function(sca_data, title="", showIndex=TRUE, plotVars=TRUE){

  if(!"RMSE" %in% colnames(sca_data)){
    message(paste0("RMSE not found. Are your models nonlinear? ",
                   "Try plotAIC() or plotDeviance() instead."))
    return(invisible(NULL))
  }

  sca_data <- sca_data %>% select(-control_coefs)

  pointSize <- -.25*(ncol(sca_data)-7)+(13/4)

  margin <- {if(title=="") unit(c(-5,2,-5,2), "points")
             else unit(c(5,2,-5,2), "points")}

  sc1 <- ggplot(data=sca_data, aes(y=RMSE, x=index)) +
    geom_point(size=pointSize) +
    labs(title=title, x="", y="RMSE") +
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
#' plotR2Adj(sca_data = sca(y = "Salnty", x = "T_degC",
#'                          controls = c("ChlorA*NO3uM", "O2Sat*NO3uM"),
#'                          data = bottles,
#'                          progressBar = TRUE, parallel = TRUE, workers = 2));
plotR2Adj <- function(sca_data, title="", showIndex=TRUE, plotVars=TRUE){

  if(!"adjR" %in% colnames(sca_data)){
    message(paste0("Adj. R^2 not found. Are your models nonlinear? ",
                   "Try plotAIC() or plotDeviance() instead."))
    return(invisible(NULL))
  }

  sca_data <- sca_data %>% select(-control_coefs)

  pointSize <- -.25*(ncol(sca_data)-7)+(13/4)

  margin <- {if(title=="") unit(c(-5,2,-5,2), "points")
    else unit(c(5,2,-5,2), "points")}

  sc1 <- ggplot(data=sca_data, aes(y=adjR, x=index)) +
    geom_point(size=pointSize) +
    labs(title=title, x="", y=bquote('Adj. R'^2)) +
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
#' plotAIC(sca_data = sca(y = "Salnty", x = "T_degC",
#'                          controls = c("ChlorA*NO3uM", "O2Sat*NO3uM"),
#'                          data = bottles,
#'                          progressBar = TRUE, parallel = TRUE, workers = 2));
plotAIC <- function(sca_data, title="", showIndex=TRUE, plotVars=TRUE){

  if(!"AIC" %in% colnames(sca_data)){
    message(paste0("AIC not found. Are your models linear? ",
                   "Try plotR2Adj() or plotRMSE instead."))
    return(invisible(NULL))
  }

  sca_data <- sca_data %>% select(-control_coefs)

  pointSize <- -.25*(ncol(sca_data)-7)+(13/4)

  margin <- {if(title=="") unit(c(-5,2,-5,2), "points")
    else unit(c(5,2,-5,2), "points")}

  sc1 <- ggplot(data=sca_data, aes(y=AIC, x=index)) +
    geom_point(size=pointSize) +
    labs(title=title, x="", y="AIC") +
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
#' plotDeviance(sca_data = sca(y = "Salnty", x="T_degC",
#'                          controls = c("ChlorA*NO3uM", "O2Sat*NO3uM"),
#'                          data = bottles, progressBar = TRUE, parallel = TRUE,
#'                          workers = 2));
plotDeviance <- function(sca_data, title="", showIndex=TRUE, plotVars=TRUE){

  if(!"deviance" %in% colnames(sca_data)){
    message(paste0("Deviance of residuals not found. ",
                   "Are your models linear? Try plotAIC(), ",
                   "plotR2Adj(), or plotRMSE() instead."))
    return(invisible(NULL))
  }

  sca_data <- sca_data %>% select(-control_coefs)

  pointSize <- -.25*(ncol(sca_data)-7)+(13/4)

  margin <- {if(title=="") unit(c(-5,2,-5,2), "points")
    else unit(c(5,2,-5,2), "points")}

  sc1 <- ggplot(data=sca_data, aes(y=deviance, x=index)) +
    geom_point(size=pointSize) +
    labs(title=title, x="", y="Deviance") +
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
#' plotControlDistributions(sca_data = sca(y = "Salnty", x = "T_degC",
#'                                     controls = c("ChlorA*NO3uM",
#'                                                  "O2Sat*NO3uM"),
#'                                     data = bottles, progressBar = TRUE,
#'                                     parallel = TRUE, workers = 2),
#'                          type = "density")
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
#' @param types A string or vector of strings specifying what types of
#'              standard errors are desired. Defaults to `"all"`.
#'              The following types are supported for non-fixed effects models:
#'                With clustering: `"HC0`, `"HC1"`, `"HC2"`, `"HC3"`.
#'                Without clustering: `"iid"` (i.e. normal standard errors),
#'                                    `"HC0`, `"HC1"`, `"HC2"`, `"HC3"`,
#'                                    `"HC4"`, `"HC4m"`, `"HC5"`,
#'                                    `"bootstrapped"`.
#'              The following types are supported for fixed effects models:
#'                With clustering: `"CL_FE"` (clustered by fixed effects, i.e.
#'                                 the default standard errors reported by
#'                                 `feols()` if no clusters are supplied), if
#'                                 clusters are supplied then the conventional
#'                                 clustered standard errors from `feols()` are
#'                                 estimated for each clustering variable. Two-
#'                                 way clustered standard errors are not
#'                                 supported at this time.
#'                Without clustering: `"HC0`, `"HC1"`, `"HC2"`, `"HC3"`,
#'                                    `"HC4"`, `"HC4m"`, `"HC5"`,
#'                                    `"bootstrapped"`.
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
#'         Some examples:
#'          `"iid"` = normal standard errors, i.e. assuming homoskedasticity
#'          `"CL_FE"` = standard errors clustered by fixed effects
#'          `"bootstrap_k8n300_FE"` =  bootstrapped standard errors for a fixed
#'                                   effects model where `bootSamples = 8` and
#'                                   `bootSampleSize = 300`
#'          `"CL_Depth_ID_FE"` = standard errors clustered by the variable
#'                               `"Depth_ID"` for a model with fixed effects
#'          `"HC0_Sta_ID`" = HC0 standard errors clustered by the variable
#'                           `"Sta_ID"`
#'          Note: for fixed effects models the `"(Intercept)"` row will be all
#'          `NA` because the intercept is not reported by `feols()` when fixed
#'          effects are present.
#'
#' @export
#'
#' @examples
#'
#' se_compare(formula = "Salnty ~ T_degC + ChlorA + O2Sat | Sta_ID",
#'            data = bottles, types = "all", cluster = c("Depth_ID", "Sta_ID"),
#'            fixedEffectsOnly = F, bootSamples=c(4, 8, 10),
#'            bootSampleSize=c(300, 500))
#'
#' se_compare(formula = "Salnty ~ T_degC + ChlorA + O2Sat", data = bottles,
#'            types = "bootstrapped", bootSamples = c(8, 10),
#'            bootSampleSize = c(300, 500))
#'
#' se_compare(formula = "Salnty ~ T_degC + ChlorA", data = bottles,
#'            types = c("HC0", "HC1", "HC3"))
#'
se_compare <- function(formula, data, types="all", cluster=NULL,
                       clusteredOnly=FALSE, fixedEffectsOnly=FALSE,
                       bootSamples=NULL, bootSampleSize=NULL){

  # Create objects that will store the standard errors
  ses_CL <- NULL
  ses_HC <- NULL
  ses_other <- NULL

  # Create the object we will eventually return
  ses <- NULL

  # If the formula contains a pipe then fixed effects are assumed to be
  # present and models are estimated with feols() rather than lm()
  if(grepl("|", formula, fixed=T)){

    model_fe <- tryCatch(feols(as.formula(formula), data=data),
                         error=function(cond){
                           message("Fixed effects model estimation failed.",
                                   cond)
                           return(NULL)
                         })

    if(is.null(model_fe)){
      break;
    }

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

        if(length(setdiff(types, types_other)!=0)){
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
      if("bootstrapped" %in% types_other){
        n_x <- length(model_fe$coefficients)

        if(length(bootSamples)==1 & length(bootSampleSize==1)){

          boot <- se_boot(data=data, formula=formula, n_x=n_x,
                          n_samples=bootSamples[[1]],
                          sample_size=bootSampleSize[[1]])

          if(!is.null(boot)){
            ses_other <- cbind(ses_other, boot)

            colnames(ses_other)[ncol(ses_other)] <- paste("bootstrap_", "k",
                                                          samples, "n",
                                                          sample_sizes,
                                                          "_FE", sep="")
          }
        }
        else{
          samples <- rep(bootSamples, length(bootSampleSize))
          sample_sizes <- sort(rep(bootSampleSize, length(bootSamples)))

          boot <- mapply(FUN=se_boot, n_samples=samples,
                         sample_size=sample_sizes,
                         MoreArgs=list(data=data, formula=formula, n_x=n_x))

          if(!is.null(boot)){

            colnames(boot) <- paste("bootstrap_", "k", samples, "n",
                                    sample_sizes, "_FE", sep="")

            ses_other <- cbind(ses_other, boot)
          }
        }
      }
      # Attach bootstapped/default SEs to the object to be returned
      ses <- cbind(ses, ses_other)
    }

    # Estimate clustered standard errors for FE model for variables other than
    # the FEs
    if(!is.null(cluster)){

      if(length(setdiff(cluster, colnames(data))!=0)){
        warning(paste0(setdiff(cluster, colnames(data)),
                       " not a valid clustering variable, ignoring.",
                       collapse="\n"))

        cluster <- cluster[cluster %in% colnames(data)]
      }

      if(!"all" %in% types){
        if(length(setdiff(types[!types %in% c("bootstrapped", "iid")],
                          types_CL))!=0){
          warning(paste0(setdiff(setdiff(types, types_CL),
                                 c("bootstrapped", "iid")),
                         " not a valid type for clustered SEs, ignoring.",
                         collapse="\n"))
        }
      }

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
  # Case when a non-FE model is desired
  if(!fixedEffectsOnly){

    # Allocate objects to hold SEs
    ses_other <- NULL
    ses_HC <- NULL
    ses_CL <- NULL

    # If the formula has FEs remove them
    if(grepl("|",formula, fixed=T)){
      formula <- str_trim(str_split(formula, fixed("|"))[[1]][[1]])
    }

    # Estimate the non-FE model and get the coefficients
    model <- lm(formula=as.formula(formula), data=data)

    ses <- cbind(ses, matrix(model$coefficients, ncol=1,
                             dimnames=list(c(names(model$coefficients)), c("estimate"))))

    # Parse the user's desired SE types
    if(!clusteredOnly){
      types_HC <- c("HC0", "HC1", "HC2", "HC3", "HC4", "HC4m", "HC5")
      types_other <- c("iid", "bootstrapped")

      if(!"all" %in% types){

        if(length(setdiff(types, c(types_HC, types_other))!=0)){
          warning(paste0(setdiff(types, c(types_HC, types_other)),
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
                       function(x) coeftest(model, vcov=vcovHC, type=x)[,2])

      # Get bootstrapped standard errors
      if("bootstrapped" %in% types_other){
        n_x <- length(model$coefficients)-1

        if(length(bootSamples)==1 & length(bootSampleSize==1)){

          boot <- se_boot(data=data, formula=formula, n_x=n_x,
                          n_samples=bootSamples[[1]],
                          sample_size=bootSampleSize[[1]])
        }
        else{
          samples <- rep(bootSamples, length(bootSampleSize))
          sample_sizes <- sort(rep(bootSampleSize, length(bootSamples)))

          boot <- mapply(FUN=se_boot, n_samples=samples,
                         sample_size=sample_sizes,
                         MoreArgs=list(data=data, formula=formula, n_x=n_x))
        }

        if(!is.null(boot)){
          colnames(boot) <- paste("bootstrap_", "k", samples, "n",
                                  sample_sizes, sep="")

          ses_other <- cbind(ses_other, boot)
        }
      }

      # Attach SEs to the return object
      ses <- cbind(ses, ses_other, ses_HC)
    }

    # Case when clustered SEs for non-FE model are desired
    if(!is.null(cluster)){

      if(length(setdiff(cluster, colnames(data))!=0)){
        warning(paste0(setdiff(cluster, colnames(data)),
                       " not a valid clustering variable, ignoring.",
                       collapse="\n"))

        cluster <- cluster[cluster %in% colnames(data)]
      }

      types_CL <- c("HC0", "HC1", "HC2", "HC3")

      if(!"all" %in% types){

        if(length(setdiff(types[!types %in% c("bootstrapped", "iid")],
                          types_CL))!=0){
          warning(paste0(setdiff(setdiff(types, types_CL),
                                 c("bootstrapped", "iid")),
                         " not a valid type for clustered SEs, ignoring.",
                         collapse="\n"))
        }

        types_CL <- types[types %in% types_CL]
      }

      # Estimate and extract clustered SEs
      if(length(types_CL)>0){

        ses_CL <- sapply(cluster, FUN=function(c, types){
          sapply(types, function(x){
            coeftest(model, vcov=vcovCL, type=x, cluster=data[c])[,2]
          })
        }, types=types_CL, simplify=F)

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
  return(as.data.frame(ses))
}
