library(tidyverse)
library(pbapply)
library(parallel)

library(lfe)
felm_model <- felm(T_degC ~ Salnty + Depthm | O2Sat, data=bottles)

felm_summary <- summary(felm_model)

library(fixest)
feols_model <- feols(T_degC ~ Salnty + ChlorA + STheta | Depthm, data=bottles)


feols_summary <- summary(feols_model)

s <- sca(y = "T_degC", x = "Salnty",
         controls = c("ChlorA", "STheta"),
         fixedEffects = "Depthm",
         data = bottles)

felm_summary

feols_summary$coefficients["Salnty"]

models_test <- readRDS("testing_models.RDS")


# to get adj r2
r2(unlist(models_test[,1]), type="ar2")


# custom function
# Not returning correct value
adj_r2_fixest <- function(x){
  y_actual <- x$residuals + x$fitted.values
  tss <- sum((y_actual-mean(y_actual))^2)

  r2 <- 1 - (x$ssr/tss)

  adj_r2 <- 1 - (((1-r2)*(x$nobs-1))/(x$nobs-length(x$coefficients)))

  return(adj_r2)
}

adj_r2_fixest(feols_summary)


# Cases to fix:
# no parallel w/ FE
# parallel w/ FE

# TODO:
# Add support for TWFE
# Add iteration of fixed effects, i.e. show with no FE, one FE var, both FE var
