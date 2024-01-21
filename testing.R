library(tidyverse)
library(pbapply)
library(parallel)

library(lfe)
felm_model <- felm(T_degC ~ Salnty + Depthm | O2Sat, data=bottles)

felm_summary <- summary(felm_model)

library(fixest)
#load("data/bottles.Rdata")
feols_model <- feols(T_degC ~ Salnty + ChlorA + STheta | Depthm + RecInd, data=bottles)


feols_summary <- summary(feols_model)

s <- sca(y = "T_degC", x = "Salnty",
         controls = c("ChlorA", "STheta"),
         fixedEffects = "Depthm",
         data = bottles)

felm_summary

feols_summary$coefficients["Salnty"]

models_test <- readRDS("models_test.RDS")
vals_test <- readRDS("vals_test.RDS")

# to get adj r2
r2(unlist(models_test[,1]), type="ar2")


# custom function
# Not returning correct value
adj_r2_fixest <- function(x){
  y_actual <- x$residuals + x$fitted.values
  tss <- sum((y_actual-mean(y_actual))^2)

  r2 <- 1 - (x$ssr/tss)

  N <- x$nobs

  L <- length(x$fixef_vars)

  K <- length(x$coefficients)

  #adj_r2 <- 1 - (1-r2) * (N - L)/(N - L - K)

  # from bard
  # mse <- mean((y_actual-x$fitted.values)^2)
  # adj_r2 <- 1 - ((N-1-K)*mse/tss)

  # from bing for within r2 (NOT ADJUSTED)
  tss_demeaned <- sum(((y_actual-mean(y_actual))-mean(y_actual-mean(y_actual)))^2)
  adj_r2 <- 1 - (x$ssr/tss_demeaned)

  return(adj_r2)
}

adj_r2_fixest(feols_summary)

r2(feols_model, type="all")



#' Extracts adjusted R^2 from `fixest::feols` model summaries.
#'
#' @description
#' `fixest` model summaries once unlisted do not contain model fit measures.
#' This function captures the output
#'
#' @param fixest_model_summary
#'
#' @return
#' @export
#'
#' @examples
fixestAdjR2Extractor <- function(fixest_model_summary){
  temp <- capture.output(fixest_model_summary)[[11]]
  adjR2 <- as.numeric(
    strsplit(
      unlist(
        regmatches(temp,
                   gregexpr("[ ][-]{0,1}[[:digit:]]+\\.{0,1}[[:digit:]]*[ ]",
                            capture.output(feols_summary)[[11]]))), " ")[[2]][[2]])

  return(adjR2)
}


# Cases to fix:
# no parallel w/ FE
# parallel w/ FE

# TODO:
# Add support for TWFE
# Add iteration of fixed effects, i.e. show with no FE, one FE var, both FE var
