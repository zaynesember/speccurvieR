load("data/bottles.RData")
test_data <- data.frame(x1=rnorm(50000, mean=4, sd=10), x2=rnorm(50000, sd=50),
                        ID=rep(1:100, 500), area=rep(1:50, 1000), y=rnorm(50000))
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

  FE <- ifelse(grepl("|", formula, fixed=T), T, F)

  coefs <- matrix(nrow=n_samples, ncol=ifelse(FE, n_x, n_x+1))

  # TODO: Deal with case where subset of coefficients omitted bc of collinearity
  # TODO: vectorize
  # Loop n_samples times, i.e. how many times we want to re-estimate the model
  for(i in 1:n_samples){
    # Estimate the model with a random subset of the data
    # sample_n is a function from the dplyr package that gives us random
    # rows from a dataframe
    model <- tryCatch(
      {
        if(FE){
          suppressMessages(feols(as.formula(formula), sample_n(data, sample_size)))
        }
        else{
          suppressMessages(lm(as.formula(formula), sample_n(data, sample_size)))
        }
      },
      error=function(cond){

        message(paste0("Estimation failed during bootstrap with n_samples=",
                      n_samples, " and sample_size=", sample_size,
                      ".\nConsider respecifying bootstrap parameters or model.\n"))

        return(NULL)
      }
    )

    if(FE & length(model$coefficients) != n_x){

      message(paste0("Estimation failed due to collinearity for ",
                     paste(model$collin.var, collapse=", "),
                     " during bootstrap with n_samples=",
                     n_samples, " and sample_size=", sample_size,
                     ".\nConsider respecifying bootstrap parameters.\n"))
      return(NULL)
    }
    else if(!FE & length(model$coefficients) != n_x+1){

      message(paste0("Estimation failed due to collinearity for ",
                     paste(names(model$coefficients[is.na(model$coefficients)]), collapse=", "),
                     " during bootstrap with n_samples=",
                     n_samples, " and sample_size=", sample_size,
                     ".\nConsider respecifying bootstrap parameters.\n"))
      return(NULL)
    }

    if(is.null(model)) return(NULL)
    # Store the coefficient estimates in row i of our matrix
    coefs[i,] <- model$coefficients
  }

  # Use apply to get the std dev of each column in the matrix, i.e. our
  # bootstrapped standard errors
  retVal <- apply(coefs, FUN=sd, MARGIN=2)

  names(retVal) <- names(model$coefficients)

  if(FE) retVal <- c("(Intercept)"=NA, unlist(retVal))

  return(retVal)

}

se_boot(data=bottles, formula="Salnty ~ T_degC + ChlorA + O2Sat | Depth_ID",
        n_x=3, n_samples=4, sample_size=300)

feols(Salnty ~ T_degC + ChlorA + O2Sat | Sta_ID, sample_n(bottles, 300))


se_boot(data=test_data, formula="y ~ x1 + x2 | ID",
        n_x=2, n_samples=10, sample_size=1000)






# maybe add column with x var name to sca output
# add column to sca with residuals
# https://jslsoc.sitehost.iu.edu/files_research/testing_tests/hccm/00TAS.pdf
# newey_west and driscoll_kraay require time series
# conley requires latitude
se_compare <- function(formula, data, types="all", cluster=NULL,
                       clusteredOnly=FALSE, fixedEffectsOnly=FALSE,
                       bootSamples=NULL, bootSampleSize=NULL){

  ses_CL <- NULL
  ses_HC <- NULL
  ses_other <- NULL
  ses <- NULL

  # If the formula contains a pipe then fixed effects are assumed to be
  # present and models are estimated with feols() rather than lm
  if(grepl("|", formula, fixed=T)){

    model_fe <- tryCatch(feols(as.formula(formula), data=data),
                         error=function(cond){
                            message("Fixed effects model estimation failed.", cond)
                            return(NULL)
                         })

    if(is.null(model_fe)){
      break;
    }


    ses <- cbind(ses, matrix(c("(Intercept)"=NA, model_fe$coefficients), ncol=1,
                             dimnames=list(c("(Intercept)", names(model_fe$coefficients)),
                                           c("estimate_FE"))))
    # Unclustered FE
    if(!clusteredOnly){
      types_other <- c("iid","bootstrapped")

      if(!"all" %in% types){

        if(length(setdiff(types, types_other)!=0)){
          warning(paste0(setdiff(types, types_other),
                         " not a valid type for SEs in FE model, ignoring.",
                         collapse="\n"))
        }

        types_other <- types[types %in% types_other]

      }

      if("iid" %in% types_other){
        ses_other <- cbind(ses_other, "CL_FE"=c("(Intercept)"=NA,
                                              coeftable(model_fe)[,2]))
      }

      if("bootstrapped" %in% types_other){
        n_x <- length(model_fe$coefficients)

        if(length(bootSamples)==1 & length(bootSampleSize==1)){

          boot <- se_boot(data=data, formula=formula, n_x=n_x,
                          n_samples=bootSamples[[1]], sample_size=bootSampleSize[[1]])

          if(!is.null(boot)){
            ses_other <- cbind(ses_other, boot)

            colnames(ses_other)[ncol(ses_other)] <- paste("bootstrap_", "k", samples, "n", sample_sizes,
                                                          "_FE", sep="")
          }
        }
        else{
          samples <- rep(bootSamples, length(bootSampleSize))
          sample_sizes <- sort(rep(bootSampleSize, length(bootSamples)))

          boot <- mapply(FUN=se_boot, n_samples=samples, sample_size=sample_sizes,
                         MoreArgs=list(data=data, formula=formula, n_x=n_x))

          if(!is.null(boot)){

            colnames(boot) <- paste("bootstrap_", "k", samples, "n", sample_sizes,
                                    "_FE", sep="")

            ses_other <- cbind(ses_other, boot)
          }
        }
      }
      ses <- cbind(ses, ses_other)
    }

    # clustered case w fe
    # FIGURE OUT
    if(!is.null(cluster)){

      if(length(setdiff(cluster, colnames(data))!=0)){
        warning(paste0(setdiff(cluster, colnames(data)),
                       " not a valid clustering variable, ignoring.", collapse="\n"))

        cluster <- cluster[cluster %in% colnames(data)]
      }

      types_CL <- c("HC1")

      if(!"all" %in% types){

        if(length(setdiff(types[!types %in% c("bootstrapped", "iid")], types_CL))!=0){
          warning(paste0(setdiff(setdiff(types, types_CL),c("bootstrapped", "iid")),
                         " not a valid type for clustered SEs, ignoring.", collapse="\n"))
        }

        types_CL <- types[types %in% types_CL]
      }

      if(length(types_CL)>0){

        ses_CL <- sapply(cluster, FUN=function(c, types){
          sapply(types, function(x){
            coeftest(model_fe, vcov=vcovCL, type=x, cluster=data[c])[,2]
          })
        }, types=types_CL, simplify=F)

        # modify above

        ses_CL <- sapply(cluster, FUN=function(c, types){
          sapply(types, function(x){
            (feols(as.formula(formula), data=data, cluster=data[c]))$coeftable[,2]
          })
        }, types=types_CL, simplify=F)

        ses_CL <- do.call(cbind, ses_CL)

        labs <- c()
        for(c in cluster){
          for(t in types_CL){
            labs <- c(labs, paste0(t, "_", c, "_FE"))
          }
        }

        colnames(ses_CL) <- labs

        ses_CL <- rbind("(Intecept)"=NA, ses_CL)
      }

      ses <- cbind(ses, ses_CL)

    }
  }
  # NON-FE CASE
  if(!fixedEffectsOnly){

    if(grepl("|",formula, fixed=T)) formula <- str_trim(str_split(formula, fixed("|"))[[1]][[1]])

    model <- lm(formula=as.formula(formula), data=data)

    ses <- cbind(ses, matrix(model$coefficients, ncol=1,
                  dimnames=list(c(names(model$coefficients)), c("estimate"))))

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

      if("iid" %in% types_other){
        ses_other <- cbind(ses_other, "iid"=summary(model)$coefficients[,2])
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
        }

        if(!is.null(boot)){
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

        if(length(setdiff(types[!types %in% c("bootstrapped", "iid")], types_CL))!=0){
          warning(paste0(setdiff(setdiff(types, types_CL),c("bootstrapped", "iid")),
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

  return(as.data.frame(ses))
}


se_compare(formula="Salnty ~ T_degC + ChlorA + O2Sat | Sta_ID", data=bottles, types="all",
           cluster=c("Depth_ID", "Sta_ID"), fixedEffectsOnly=T,
           bootSamples=c(4, 8, 10), bootSampleSize=c(300))

se_compare(formula="y ~ x1 + x2 | ID", data=test_data, types=c("all"),
           cluster=c("area"), fixedEffectsOnly=F,
           bootSamples=c(10, 20), bootSampleSize=c(1000))




# broken
se_compare(formula="Salnty ~ T_degC + ChlorA + O2Sat | Depth_ID", data=bottles, types=c("all"),
           cluster=c("Depth_ID", "Sta_ID"), fixedEffectsOnly=F,
           bootSamples=c(4, 8, 10), bootSampleSize=c(300))

feols(Salnty ~ T_degC + ChlorA + O2Sat | Sta_ID, data=bottles)

grepl("|", "Salnty ~ T_degC + ChlorA + O2Sat | x1 + x2", fixed=T)


temp <- coeftest(lm(Salnty ~ T_degC + ChlorA + O2Sat, bottles), vcov=vcovHC, type="HC1")

test_model <- lm(Salnty ~ T_degC + ChlorA + O2Sat, bottles)

test_model_fe <- feols(Salnty ~ T_degC + ChlorA + O2Sat | Sta_ID, bottles)
test_model_fe_nocl <- feols(Salnty ~ T_degC + ChlorA + O2Sat, bottles)

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

