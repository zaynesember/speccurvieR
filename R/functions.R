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

  if(family=="linear"){

    # No fixed effects
    if(is.null(fixedEffects)){

      # Get each value of interest across models
      coef <- lapply(X=models, function(x2) x2$coefficients[,1])
      se <- lapply(X=models, function(x2) x2$coefficients[,2])
      statistic <- lapply(X=models, function(x2) x2$coefficients[,3])
      p <- lapply(X=models, function(x2) x2$coefficients[,4])
      terms <- lapply(X=models, FUN=function(x2) row.names(x2$coefficients))
      RMSE <- lapply(X=models, FUN=function(x2) sqrt(mean(x2$residuals^2)))
      adjR <- lapply(X=models, function(x2) x2$adj.r.squared)
      control_coefs <- lapply(X=models,
                              FUN=function(x2, x3) controlExtractor(x2,x3)$terms,
                              x3=x)

    }
    # Fixed effects
    else{
      # Get each value of interest across models
      coef <- lapply(X=models, function(x2) x2$coeftable[,1])
      se <- lapply(X=models, function(x2) x2$coeftable[,2])
      statistic <- lapply(X=models, function(x2) x2$coeftable[,3])
      p <- lapply(X=models, function(x2) x2$coeftable[,4])
      terms <- lapply(X=models, FUN=function(x2) row.names(x2$coeftable))
      RMSE <- lapply(X=models, FUN=function(x2) fitstat(x2, type="rmse",
                                                        verbose=F)[[1]])
      adjR <- lapply(X=models, FUN=function(x2) fitstat(x2, type="war2",
                                                        verbose=F)[[1]])
      control_coefs <- lapply(X=models,
                              FUN=function(x2, x3) controlExtractor(x2,x3)$term,
                              x3=x)
    }

    # Get the number of rows needed for each model
    list_lengths <- lapply(coef, length)

    # Store values in a data frame to be returned
    retVal <- data.frame(terms=unlist(terms),
                         coef=unlist(coef), se=unlist(se),
                         statistic=unlist(statistic), p=unlist(p),
                         RMSE=rep(unlist(RMSE), times=list_lengths),
                         adjR=rep(unlist(adjR), times=list_lengths)) %>%
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

    retVal$control_coefs <- rep(control_coefs, times=list_lengths)

  }
  # glm models
  else{
    # Get each value of interest across models
    coef <- lapply(X=models, function(x2) x2$coefficients[,1])
    se <- lapply(X=models, function(x2) x2$coefficients[,2])
    statistic <- lapply(X=models, function(x2) x2$coefficients[,3])
    p <- lapply(X=models, function(x2) x2$coefficients[,4])
    terms <- lapply(X=models, FUN=function(x2) row.names(x2$coefficients))
    AIC <- lapply(X=models, FUN=function(x2) x2$aic)
    deviance <- lapply(X=models, FUN=function(x2) x2$deviance)
    control_coefs <- lapply(X=models,
                      FUN=function(x2,x3,
                           x4) controlExtractor(x2,x3)$term,x3=x)


    # Get the number of rows needed for each model
    list_lengths <- lapply(coef, length)

    # Store values in a data frame to be returned
    retVal <- data.frame(terms=unlist(terms),
                         coef=unlist(coef), se=unlist(se),
                         statistic=unlist(statistic), p=unlist(p),
                         AIC=rep(unlist(AIC), times=list_lengths),
                         deviance=rep(unlist(deviance), times=list_lengths)) %>%
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

    retVal$control_coefs <- rep(control_coefs, times=list_lengths)
  }

  # Build dummy columns for terms present in each model for visualization
  temp <- data.frame(matrix(ncol = length(controls), nrow = nrow(retVal)))

  colnames(temp) <- controls

  retVal <- cbind(retVal, temp)

  for(c in controls){
    retVal[c] <- ifelse(str_detect(retVal$terms, fixed(c)), 1, 0)
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
#' plotCurve(sca_data = sca(y="Salnty", x="T_degC", c("ChlorA*O2Sat"),
#'                          data=bottles, progressBar=FALSE, parallel=FALSE),
#'                      showIndex = TRUE, plotVars = TRUE,
#'                      plotSE = "ribbon");
#' plotCurve(sca_data = sca(y="Salnty", x="T_degC",
#'                          c("ChlorA*NO3uM", "O2Sat*NO3uM"), data=bottles,
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

  histData <- bind_rows(unAsIs(sca_data$control_coefs)) %>%
    rename(term=coef, coef=term)

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
