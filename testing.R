library(tidyverse)
library(pbapply)
library(parallel)
library(fixest)
#load("data/bottles.Rdata")

feols_model <- feols(T_degC ~ Salnty + ChlorA + STheta | Depthm + RecInd, data=bottles)


feols_summary <- summary(feols_model)

s <- sca(y = "T_degC", x = "Salnty",
         controls = c("ChlorA", "STheta"),
         fixedEffects = "Depthm",
         data = bottles)

felm_summary

lm_formulae <- list(as.formula("T_degC ~ Salnty + ChlorA + STheta"),
                    as.formula("T_degC ~ Salnty + STheta"))

lm_models <- lapply(X=lm_formulae, lm, data=bottles)

lm_summaries <- lapply(X=lm_models, summary)

lm_summaries[[1]]$coefficients

glm_model <- glm(Depthm ~ T_degC + Salnty + ChlorA, data=bottles,
                 family=poisson(link="log"))

glm_summary <- summary(glm_model)

s2 <- sca(y = "Depthm", x = "T_degC",
          controls = c("ChlorA", "STheta"),
          data = bottles, family="poisson", link="log")

# TODO:
# Add support for TWFE
# Add iteration of fixed effects, i.e. show with no FE, one FE var, both FE var
# Add user control for type of r2 estimated
