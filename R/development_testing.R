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

    if(!clusteredOnly){
      types_HC <- c("HC0", "HC1", "HC2", "HC3", "HC4", "HC5")
      types_other <- c("normal", "bootstrapped")

      if(types!=c("all") | types != "all"){

        types_HC <- types[types %in% types_HC]
        types_other <- types[types %in% types_other]

        if(length(setdiff(types, c(types_HC, types_other))!=0)){
          warning(paste0(setdiff(types, c(types_HC, types_other)),
                         " not valid type(s), ignoring."))
        }
      }


      ses_other <- matrix(nrow=length(model$coefficients), ncol=0)

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

        ses <- cbind(ses_other, ses_HC)
      }

    }
    # TODO: get clustered (no FE) case working
    if(!is.null(cluster)){
      print("here")
      types_CL <- c("HC0", "HC1", "HC2", "HC3")

      if(types != "all" | types != c("all")){
        types_CL <- types[types %in% types_CL]

        if(length(setdiff(types, types_CL))!=0){
          warning(paste0(setdiff(types, types_CL),
                         " not valid type(s), ignoring."))
        }
      }

      ses_CL <- sapply(types,
                       function(x){
                         coeftest(model, vcov=vcovCL, type=x, cluster=cluster)[,2]
                       }
      )
    }
  }
  print(ses_CL)
  if(!is.null(ses_CL)) ses <- cbind(ses, ses_CL)

  return(cbind(Estimate=model$coefficients, ses))
}

# FIgure out
se_compare(formula="Salnty ~ T_degC + ChlorA + O2Sat", data=bottles, types="all",
           cluster="Sta_id",
           bootSamples=c(4, 8), bootSampleSize=c(100))


test <- c(50, 150)
test2 <- c(100, 200, 300)
samples <- rep(test, length(test2))
sample_sizes <- sort(rep(test2, length(test)))

paste(samples, "samples", sample_sizes, "samplesize", sep="_")

temp <- coeftest(lm(Salnty ~ T_degC + ChlorA + O2Sat, bottles), vcov=vcovHC, type="HC1")

feols(Salnty ~ T_degC + ChlorA + O2Sat, bottles, vcov="cluster")
feols(Salnty ~ T_degC + ChlorA + O2Sat, bottles, vcov="normal") # same as iid
feols(Salnty ~ T_degC + ChlorA + O2Sat, bottles, vcov="white")
feols(Salnty ~ T_degC + ChlorA + O2Sat, bottles, vcov="HC1") # same as white


lm.boot(feols(Salnty ~ T_degC + ChlorA + O2Sat, bottles, vcov="cluster"), R=10)
# NOTES


test <- lm(Salnty ~ T_degC + ChlorA + O2Sat, bottles)
# TODO:
# SE comparison

