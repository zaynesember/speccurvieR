load("data/bottles.RData")
library(tidyverse)
library(fixest)
library(ggplot2)
library(parallel)
library(pbapply)
library(stringr)

# new libraries needed
library(lmtest)
library(sandwich)

ex1 <- sca(y = "Salnty", x = "T_degC", controls = c("ChlorA", "O2Sat"),
    data = bottles, progressBar = TRUE, parallel = FALSE)


# This function takes the following arguments:
#   data = a dataframe with our data
#   formula = a formula object with our regression formula
#   n_x = the number of x variables we have in the model
#   n_samples = the number of times to estimate the model with a random subset
#               of the data
#   sample_size = the number of observations to include in the subset of data
# It returns a list of bootstrapped standard errors
se_boot <- function(data, formula, n_x, n_samples, sample_size){
  # Create a matrix to store the coefficient estimates, each row contains
  # coefficients estimated from a different subset of the data
  # (ncol=n_x+1 because we are also storing the intercept estimate)
  coefs <- matrix(nrow=n_samples, ncol=n_x+1)

  # TODO: vectorize
  # Loop n_samples times, i.e. how many times we want to re-estimate the model
  for(i in 1:n_samples){
    # Estimate the model with a random subset of the data
    # sample_n is a function from the dplyr package that gives us random
    # rows from a dataframe
    model <- lm(formula, sample_n(data, sample_size))

    # Store the coefficient estimates in row i of our matrix
    coefs[i,] <- model$coefficient
  }

  # Use apply to get the std dev of each column in the matrix, i.e. our
  # bootstrapped standard errors
  retVal <- apply(coefs, FUN=sd, MARGIN=2)
  names(retVal) <- names(model$coefficient)

  return(retVal)
}


# maybe add column with x var name to sca output
# add column to sca with residuals
# https://jslsoc.sitehost.iu.edu/files_research/testing_tests/hccm/00TAS.pdf
# newey_west and driscoll_kraay require time series
# conley requires latitude
se_compare <- function(formula, data, types="all", cluster=NULL,
                       clusteredOnly=FALSE,
                       bootSamples=NULL, bootSampleSize=NULL,
                       timeSeries=FALSE, spatial=FALSE){

  ses_CL <- NULL
  ses_HC <- NULL
  ses_other <- NULL
  ses <- NULL

  # If the formula contains a pipe then fixed effects are assumed to be
  # present and models are estimated with feols() rather than lm
  # Same goes for when cluster != NULL
  if(grepl(formula, "|")){
    if(is.null(cluster & types!="all" & "cluster" %in% types)){
      warning(
        paste0("Clustered SEs included but no clustering variable provided.",
                     "Try again with a clustering variable specified")
        )
      return(NA)
    }

    if(types=="all" & timeSeries & spatial){
      if(!is.null(cluster)){
        types <- c("normal", "HC1", "cluster", "newey_west",
                   "driscoll_kraay", "conley")
      }
      else{types <- c("normal", "HC1", "newey_west",
                      "driscoll_kraay", "conley")}
    }

    else if(types=="all" & timeSeries & !spatial){
      if(!is.null(cluster)){
        types <- c("normal", "HC1", "cluster", "newey_west",
                   "driscoll_kraay")
      }
      else{types <- c("normal", "HC1", "newey_west",
                      "driscoll_kraay")}
    }

    else if(types=="all" & !timeSeries & spatial){
      if(!is.null(cluster)){
        types <- c("normal", "HC1", "cluster", "conley")
      }
      else{types <- c("normal", "HC1", "conley")}
    }

    else if(types=="all" & !timeSeries & !spatial){
      if(!is.null(cluster)){
        types <- c("normal", "HC1", "cluster")
      }
      else{types <- c("normal", "HC1")}
    }

    if(is.null(cluster)){
      models <- sapply(types, function(x) feols(fml=as.formula(formula),
                                                data=data,
                                                vcov=x))
    }
  }

  else{
    model <- lm(formula=as.formula(formula), data=data)

    ses <- matrix(model$coefficients, ncol=1,
                  dimnames=list(c(names(model$coefficients)), c("estimate")))

    if(!clusteredOnly){
      types_HC <- c("HC0", "HC1", "HC2", "HC3", "HC4", "HC5")
      types_other <- c("normal", "bootstrapped")

      if(!"all" %in% types){

        if(length(setdiff(types, c(types_HC, types_other))!=0)){
          warning(paste0(setdiff(types, c(types_HC, types_other)),
                         " not a valid type for SEs, ignoring.", collapse="\n"))
        }

        types_HC <- types[types %in% types_HC]
        types_other <- types[types %in% types_other]

      }

      if("normal" %in% types_other){
        ses_other <- cbind(ses_other, "normal"=summary(model)$coefficients[,2])
      }

      ses_HC <- sapply(types_HC,
                       function(x) coeftest(model, vcov=vcovHC, type=x)[,2])


      if("bootstrapped" %in% types_other){
        n_x <- length(model$coefficients)-1

        if(length(bootSamples)==1 & length(bootSampleSize==1)){

          boot <- se_boot(data=data, formula=formula, n_x=n_x,
                          n_samples=bootSamples[[1]], sample_size=bootSampleSize[[1]])
        }
        else{
          samples <- rep(bootSamples, length(bootSampleSize))
          sample_sizes <- sort(rep(bootSampleSize, length(bootSamples)))

          boot <- mapply(FUN=se_boot, n_samples=samples, sample_size=sample_sizes,
                          MoreArgs=list(data=data, formula=formula, n_x=n_x))
          colnames(boot) <- paste("bootstrap_", "k", samples, "n", sample_sizes, sep="")

          ses_other <- cbind(ses_other, boot)
        }
      }

      ses <- cbind(ses, ses_other, ses_HC)
    }

    # clustered no FE case
    if(!is.null(cluster)){

      if(length(setdiff(cluster, colnames(data))!=0)){
        warning(paste0(setdiff(cluster, colnames(data)),
                       " not a valid clustering variable, ignoring.", collapse="\n"))

        cluster <- cluster[cluster %in% colnames(data)]
      }

      types_CL <- c("HC0", "HC1", "HC2", "HC3")

      if(!"all" %in% types){

        if(length(setdiff(types[!types %in% c("bootstrapped", "normal")], types_CL))!=0){
          warning(paste0(setdiff(setdiff(types, types_CL),c("bootstrapped", "normal")),
                         " not a valid type for clustered SEs, ignoring.", collapse="\n"))
        }

        types_CL <- types[types %in% types_CL]
      }

      if(length(types_CL)>0){

        ses_CL <- sapply(cluster, FUN=function(c, types){
          sapply(types, function(x){
            coeftest(model, vcov=vcovCL, type=x, cluster=data[c])[,2]
          })
        }, types=types_CL, simplify=F)

        ses_CL <- do.call(cbind, ses_CL)

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

  return(ses)
}

# Figure out FE models
se_compare(formula="Salnty ~ T_degC + ChlorA + O2Sat", data=bottles, types=c("bootstrapped", "HC4", "white"),
           cluster=c("Depth_ID", "Sta_ID"),
           bootSamples=c(4, 8), bootSampleSize=c(100))



temp <- coeftest(lm(Salnty ~ T_degC + ChlorA + O2Sat, bottles), vcov=vcovHC, type="HC1")

test_model <- lm(Salnty ~ T_degC + ChlorA + O2Sat, bottles)

coeftest(test_model, vcov=vcovCL, cluster=bottles["Depth_ID"])

feols(Salnty ~ T_degC + ChlorA + O2Sat, bottles, vcov="cluster")
feols(Salnty ~ T_degC + ChlorA + O2Sat, bottles, vcov="normal") # same as iid
feols(Salnty ~ T_degC + ChlorA + O2Sat, bottles, vcov="white")
feols(Salnty ~ T_degC + ChlorA + O2Sat, bottles, vcov="HC1") # same as white


lm.boot(feols(Salnty ~ T_degC + ChlorA + O2Sat, bottles, vcov="cluster"), R=10)
# NOTES


test <- lm(Salnty ~ T_degC + ChlorA + O2Sat, bottles)
# TODO:
# SE comparison

