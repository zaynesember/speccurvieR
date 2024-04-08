load("data/bottles.RData")
library(tidyverse)
library(fixest)
library(ggplot2)
library(parallel)
library(pbapply)
library(stringr)

ex1 <- sca(y = "Salnty", x = "T_degC", controls = c("ChlorA", "O2Sat"),
    data = bottles, progressBar = TRUE, parallel = FALSE)

# maybe add column with x var name to sca output
# add columnn to sca with residuals
# https://jslsoc.sitehost.iu.edu/files_research/testing_tests/hccm/00TAS.pdf
se_compare <- function(x, sca, types="all", cluster=NULL,){
  X <-
  resids <-
  # for HC2
  h <- # figure out
  N <- # nobs
  K <- # def freedom
  if(types=="all"){
    # Calculate robust HC0
    HC0 <- inv(T(X) %*% X) %*% t(X) %*% diag(resids^2) %*% X %*% inv(T(X) %*% X)
    # Calculate robust HC1
    HC1 <- (N/N-K) * HC0
    # Calculate robust HC2
    HC2 <- inv(T(X) %*% X) %*% t(X) %*% diag(resids^2 / (1-h)) %*% X %*% inv(T(X) %*% X)
    # Calculate robust HC3
    HC3 <- inv(T(X) %*% X) %*% t(X) %*% diag(resids^2 / (1-h)^2) %*% X %*% inv(T(X) %*% X)
    # Calculate clustered
    # Calculate clustered robust HC1(?)
    # Calculate bootstrap
  }
}


# TODO:
# SE comparison

