# Helper functions--------------------------------------------------------------

#' Paste together controls and independent variable
#'
#' @param controls A vector of strings containing control variable names.
#' @param x A string containing the independent variable name.
#'
#' @returns A string concatenating independent and control variables separated
#'          by '+'.
#'
#' @example
#' paste_factory(c("control1", "control2"), "independentVariable");
paste_factory <- function(controls, x){
  if(T %in% str_detect(controls, x)){
    return(paste(controls, collapse=" + "))
  }
  else return(paste(x, paste(controls, collapse=" + "), sep=" + "))
}


#' Removes duplicate control variables from user input
#'
#' @param controls A vector of strings containing control variable names.
#' @param x A string containing the independent variable name.
#'
#' @return A vector of strings containing control variable names
#'
#' @example
#' duplicate_remover(c("control1", "control2*control3"), "independentVariable");
duplicate_remover <- function(controls, x){
  # Check for interactions
  if(T %in% str_detect(controls, "\\*")){
    # Find interaction terms
    indices <- which(T==str_detect(controls, "\\*"))
    # Find controls that are in interaction terms
    extraTerms <- str_replace(str_replace(controls[indices],
                                          pattern=x,
                                          replacement=""),
                              pattern="\\*",
                              replacement="")
    # Remove controls that are already present in interaction
    return(controls[!controls %in% extraTerms])
  }
  else return(controls)
}


#' Builds models formulae with every combination of control variables possible.
#'
#' @param y A string containing the dependent variable name.
#' @param x A string containing the independent variable name.
#' @param controls A vector of strings containing control variable names.
#' @param fixedEffects A string containing the name of a variable to use for
#'                     fixed effects, defaults to `NA` indicating no fixed
#'                     effects desired.
#' @param cluster A string containing the name of a variable to use for
#'                clustered standard errors, defaults to `NA` indicating no
#'                clustering.
#'
#' @return A vector of formula objects using every possible combination of
#'         controls.
#'
#' @examples
#' formula_builder("dependentVariable", "independentVariable",
#'                 c("control1", "control2"));
#' formula_builder("dependentVariable", "independentVariable",
#'                 c("control1*control2"), fixedEffects="month");
#' formula_builder("dependentVariable", "independentVariable",
#'                 c("control1*control2", "control3"), cluster="village");
formula_builder <- function(y, x, controls, fixedEffects=NA, cluster=NA){

  # Get all combinations of controls
  powerset <- unlist(lapply(1:length(controls),
                            combinat::combn,
                            x = controls,
                            simplify = FALSE),
                     recursive=F)

  # Remove duplicate controls that are already in the interaction
  powerset <- unique(sapply(X=powerset, FUN=duplicate_remover, x=x))

  # Build right hand side of the formulae
  if(is.na(fixedEffects)){
    if(is.na(cluster)) RHS <- unique(sapply(powerset, paste_factory, x))
    else RHS <- paste(unique(sapply(powerset, paste_factory, x)), "",
                      cluster, sep=" | ")
  }
  else{
    if(is.na(cluster)) RHS <- paste(unique(sapply(powerset, paste_factory, x)),
                                   fixedEffects, sep=" | ")
    else RHS <- paste(unique(sapply(powerset, paste_factory, x)), fixedEffects,
                      cluster, sep=" | ")
  }
  # Build formulae
  formulae <- sapply(paste(y, RHS, sep=" ~ "), formula)

  return(formulae)
}


#' Extracts the control variable names and coefficients from a model summary.
#'
#' @param model A model summary object.
#' @param x A string containing the independent variable name.
#'
#' @return A dataframe with two columns, `term` contains the name of the control
#'         and `coef` contains the coefficient estimate.
#'
#' @examples
#' controlExtractor(summary(lm(y ~ x + control1 + control2, data)), "x");
controlExtractor <- function(model, x){

  r <- as.data.frame(model$coefficients[,1]) %>%
    mutate(term=row.names(.)) %>%
    filter(!row.names(.) %in% c("(Intercept)", x))

  names(r) <- c("term", "coef")

  return(r)
}


unAsIs <- function(x) {
  if("AsIs" %in% class(x)) {
    class(x) <- class(x)[-match("AsIs", class(x))]
  }
  return(x)
}


# Takes in the output of sca() and returns a list with the dataframe and
# labels to make a plot to visualize the controls included in each spec curve
# model
# Arguments:
#   spec_data = dataframe object with output from `sca()`
# Returns: list containing dataframe, controls, and control IDs
scp <- function(spec_data){
  df <- spec_data %>%
    select(-terms, -coef, -se, -statistic, -p, -sig.level) %>%
    pivot_longer(!index, names_to="control", values_to="value") %>%
    filter(value==1) %>%
    mutate(controlID = with(.,match(control, unique(control)))) %>%
    select(-value)

  df_labels <- df %>% select(control, controlID) %>% unique()

  return(list(df, setNames(as.character(df_labels$control),
                           df_labels$controlID)))
}

# User-facing functions---------------------------------------------------------

# TODO: UPDATE DOCUMENTATION
# Runs every possible combination of regression models using lm() or felm()
# Arguments:
#   y = string name of dependent variable column
#   x = string name of independent variable column
#   controls = vector of strings containing desired controls and interactions
#   data = dataframe object containing y, x, controls, and fixed effects
#          variables
#   fixedEffects = string adding fixed effects variables if desired
# Returns: dataframe object containing the coefficient, standard error, t-value,
#          p-value, list of terms in the regression, and significance level


#' The workhorse function of the package--this estimates models with every
#' possible combination of the controls supplied and returns a dataframe
#' where each row contains the pertinent information and parameters for a
#' given model by default. Alternatively, if `returnFormulae` is `TRUE`, it
#' returns a list of formula objects with every possible combination of
#' controls.
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
#'             `lfe::felm()` depending on whether fixed effects are supplied.
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
#' @examples
#' TODO
sca <- function(y, x, controls, data, family="linear", link=NULL,
                fixedEffects=NULL, returnFormulae=F,
                progressBar=T, parallel=FALSE, workers=2){

  if(family!="linear" & !is.null(fixedEffects))
    {
    warning(paste0("Fixed effects unsupported for models other than OLS ",
                   "regression. Ignoring fixed effects"))
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
    clusterEvalQ(cl, library(lfe))

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
          system.time(models <- pbsapply(
            formulae, function(x2) summary(lm(x2, data=data)), cl=cl))
        }
        else{
          system.time(models <- pbsapply(
            formulae, function(x2) summary(
              glm(x2, data=data, family=eval(parse(text=family_link)))),
                                         cl=cl))
        }
      }
      else{
        if(family=="linear"){
          models <- parSapply(
            cl, formulae, function(x2) summary(lm(x2, data=data)))
        }
        else{
          models <- parSapply(
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
        system.time(models <- pbsapply(formulae,
                                       function(x2) summary(felm(x2,data=data)),
                                       cl=cl))
      }
      else{
        models <- parSapply(cl, formulae,
                            function(x2) summary(felm(x2, data=data)))
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
          system.time(models <- pbsapply(
            formulae, function(x2) summary(lm(x2, data=data))))
        }
        else{
          system.time(models <- pbsapply(
            formulae, function(x2) summary(
              glm(x2, data=data, family=eval(parse(text=family_link))))))
        }
      }
      else{
        if(family=="linear"){
          models <- sapply(formulae, function(x2) summary(lm(x2, data=data)))
        }
        else{
          models <- sapply(formulae,
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
        system.time(models <- pbsapply(
          X=formulae, function(x2) summary(felm(x2, data=data))))
      }
      else{
        models <- sapply(X=formulae, function(x2) summary(felm(x2, data=data)))
      }
    }
  }

  # Garbage collection for parallel connections
  if(parallel) stopCluster(cl=cl)

  if(family=="linear"){
    # Extract results for IV
    vals <- apply(X=models, MARGIN=2,
                   FUN=function(x2) list(x2$coefficients[x,],
                                         sqrt(mean(x2$residuals^2)),
                                         x2$adj.r.squared,
                                         controlExtractor(x2, x)))
    # Get each value of interest across models
    coef <- unlist(lapply(lapply(vals, `[[`, 1), `[[`, 1))
    se <- unlist(lapply(lapply(vals, `[[`, 1), `[[`, 2))
    statistic <- unlist(lapply(lapply(vals, `[[`, 1), `[[`, 3))
    p <- unlist(lapply(lapply(vals, `[[`, 1), `[[`, 4))
    terms <- names(apply(X=models, MARGIN=2, FUN=function(x2) x2$terms[[1]]))
    RMSE <- unlist(lapply(lapply(vals, `[[`, 2), `[[`, 1))
    adjR <- unlist(lapply(lapply(vals, `[[`, 3), `[[`, 1))
    control_coefs <- lapply(vals, `[[`, 4)


    # Put into a dataframe
    retVal <- data.frame(terms, coef, se, statistic, p, RMSE,
                         adjR, control_coefs=I(control_coefs)) %>%
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
  else{
    # Extract results for IV
    vals <- apply(X=models, MARGIN=2,
                  FUN=function(x2) list(x2$coefficients[x,],
                                        x2$aic,
                                        x2$deviance,
                                        controlExtractor(x2, x)))
    # Get each value of interest across models
    coef <- unlist(lapply(lapply(vals, `[[`, 1), `[[`, 1))
    se <- unlist(lapply(lapply(vals, `[[`, 1), `[[`, 2))
    statistic <- unlist(lapply(lapply(vals, `[[`, 1), `[[`, 3))
    p <- unlist(lapply(lapply(vals, `[[`, 1), `[[`, 4))
    terms <- names(apply(X=models, MARGIN=2, FUN=function(x2) x2$terms[[1]]))
    AIC <- unlist(lapply(lapply(vals, `[[`, 2), `[[`, 1))
    deviance <- unlist(lapply(lapply(vals, `[[`, 3), `[[`, 1))
    control_coefs <- lapply(vals, `[[`, 4)


    # Put into a dataframe
    retVal <- data.frame(terms, coef, se, statistic, p, AIC, deviance,
                         control_coefs=I(control_coefs)) %>%
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

  colnames(temp) <- controls

  retVal <- cbind(retVal, temp)

  for(c in controls){
    retVal[c] <- ifelse(str_detect(retVal$terms, fixed(c)), 1, 0)
  }

  # Remove duplicate columns
  retVal <- retVal %>% select(where(~!all(is.na(.x))))

  return(retVal)

}

# TODO: Documentation and testing
plotCurve <- function(sca_data, title="", showIndex=T, plotVars=T,
                         ylab="Coefficient", plotSE="bar"){

  if("control_coefs" %in% names(sca_data)){
    sca_data <- sca_data %>% select(-control_coefs)
  }

  pointSize <- -.25*(ncol(sca_data)-7)+(13/4)

  if(plotSE=="ribbon"){
    sca_data <- sca_data %>%
      mutate(ribbon.group = cumsum(sig.level != lag(sig.level,
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
    {if(plotSE=="bar") geom_errorbar(aes(ymin=coef-se, ymax=coef+se,
                                          color=factor(sig.level)),
                                      width=0.25)} +
    {if(!plotSE %in% c("ribbon",
                       "bar")) geom_point(aes(color=as.factor(sig.level)),
                                          size=pointSize)} +
    {if(plotSE %in% c("ribbon",
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

# TODO: Documentation and testing
plotVars <- function(sca_data, title="", colorControls=F){

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

# TODO: Documentation and testing
plotRMSE <- function(sca_data, title="", showIndex=T, plotVars=T){

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

# TODO: Documentation and testing
plotR2Adj <- function(sca_data, title="", showIndex=T, plotVars=T){

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

# TODO: Documentation and testing
plotAIC <- function(sca_data, title="", showIndex=T, plotVars=T){

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

# TODO: Documentation and testing
plotDeviance <- function(sca_data, title="", showIndex=T, plotVars=T){

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

# TODO: Documentation
plotControlDistributions <- function(sca_data, title="", type="density"){

  histData <- bind_rows(unAsIs(sca_data$control_coefs)) %>%
    rename(term=coef, coef=term)

  rownames(histData) <- NULL

  histData$term <- sapply(sapply(histData$term, str_split, pattern=":"),
                          paste0, collapse=" %*% ")

  n_facets <- length(unique(histData$term))

  sc1 <- histData %>%
    ggplot(aes(x=coef, fill=factor(term))) +
      {if(type=="hist" | type=="histogram") geom_histogram()
       else if (type=="density") geom_density()} +
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

# TODO: AIC and Deviance plots
# TODO: Add support for grouping models
# TODO: Swap lm() for glm() and ensure support for all glm() models
# TODO: Add other measures of model fit (AIC, BIC, Log. Likelihood)
# TODO: Add support for IV via felm
# TODO: Add support random effects
# TODO: Add support for clustered SEs
# TODO: Add support for robust SEs
# TODO: Add support for comparison between SE types
# TODO: Add plot of number of observations across models
#       Maybe add some kind of power analysis?
# TODO: Compare model estimation speed to specR
