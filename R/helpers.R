# Helper functions--------------------------------------------------------------

#' Builds models formulae with every combination of control variables possible.
#'
#' @param y A string containing the dependent variable name.
#' @param x A string containing the independent variable name.
#' @param controls A vector of strings containing control variable names.
#' @param fixedEffects A string containing the name of a variable to use for
#'                     fixed effects, defaults to `NA` indicating no fixed
#'                     effects desired.
#'
#' @return A vector of formula objects using every possible combination of
#'         controls.
#'
#' @export
#'
#' @examples
#' formula_builder("dependentVariable", "independentVariable",
#'                 c("control1", "control2"));
#' formula_builder("dependentVariable", "independentVariable",
#'                 c("control1*control2"), fixedEffects="month");
formula_builder <- function(y, x, controls, fixedEffects=NA){

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
    RHS <- unique(sapply(powerset, paste_factory, x))
  }
  else{
    RHS <- paste(unique(sapply(powerset, paste_factory, x)), fixedEffects,
                 sep=" | ")
  }
  # Build formulae
  formulae <- sapply(paste(y, RHS, sep=" ~ "), formula)

  return(formulae)
}

#' Paste together controls and independent variable
#'
#' @description
#' `paste_factory()` constructs the right hand side of the regression as a
#' a string i.e. "x + control1 + control2".
#'
#'
#' @inheritParams formula_builder
#'
#' @returns A string concatenating independent and control variables separated
#'          by '+'.
#'
#' @export
#'
#' @examples
#' paste_factory(controls = c("control1", "control2"),
#'               x = "independentVariable");
paste_factory <- function(controls, x){
  if(T %in% str_detect(controls, x)){
    return(paste(controls, collapse=" + "))
  }
  else return(paste(x, paste(controls, collapse=" + "), sep=" + "))
}


#' Removes duplicate control variables
#'
#' @description
#' Removes duplicate control variables from user input.
#'
#' @inheritParams formula_builder
#'
#' @return A vector of strings containing control variable names
#'
#' @export
#'
#' @examples
#' duplicate_remover(controls = c("control1", "control2*control3"),
#'                   x = "independentVariable");
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


#' Extracts the control variable names and coefficients from an lm model
#' summary.
#'
#' @description
#' Extracts the control variable names and coefficients from a model summary.
#'
#'
#' @param model A model summary object.
#' @param feols_model An indicator for whether `model` is a `fixest::feols()`
#'        model. Defaults to `FALSE`.
#' @inheritParams formula_builder
#'
#' @return A dataframe with two columns, `term` contains the name of the control
#'         and `coef` contains the coefficient estimate.
#'
#' @export
#'
#' @examples
#' m <- summary(lm(Salnty ~ STheta + T_degC, bottles))
#' controlExtractor(model = m, x = "STheta");
#'
#' m <- summary(lm(Salnty ~ STheta*T_degC + O2Sat, bottles))
#' controlExtractor(model = m, x = "STheta");
controlExtractor <- function(model, x, feols_model=F){
  if(feols_model){
    input <- model$coeftable[,1]
  }
  else{
    input <- model$coefficients[,1]
  }

  r <- as.data.frame(input) %>%
    mutate(term=row.names(.)) %>%
    filter(!row.names(.) %in% c("(Intercept)", x))

  names(r) <- c("coef", "term")

  return(r)
}


#' Removes the `AsIs` class attribute from the input.
#'
#' @description
#' Removes the `AsIs` class attribute from the input. Taken from:
#' <https://stackoverflow.com/a/12866609>
#'
#' @param x An object with the `AsIs` class attribute.
#'
#' @return An object without the `AsIs` class attribute.
#'
#' @export
#'
#' @examples
#' unAsIs(x = I(c(1:4)));
unAsIs <- function(x) {
  if("AsIs" %in% class(x)) {
    class(x) <- class(x)[-match("AsIs", class(x))]
  }
  return(x)
}

#' Prepares the output of `sca()` for plotting.
#'
#' @description
#' Takes in the data frame output by `sca()` and returns a list with the data
#' frame and labels to make a plot to visualize the controls included in each
#' spec curve model.
#'
#' @param sca_data A data frame output by `sca`.
#'
#' @return A list containing a data frame, control coefficients, and control
#'         names.
#'
#' @export
#'
#' @examples
#' scp(sca(y = "Salnty", x = "T_degC", controls = c("ChlorA", "O2Sat"),
#'         data = bottles, progressBar=TRUE, parallel=FALSE));
scp <- function(sca_data){
  if("control_coefs" %in% names(sca_data)){
    df <- sca_data %>%
      select(-terms, -coef, -se, -statistic, -p, -sig.level) %>%
      pivot_longer(-c(index, control_coefs),
                   names_to="control", values_to="value") %>%
      filter(value==1) %>%
      mutate(controlID = with(.,match(control, unique(control)))) %>%
      select(-value)
  }
  else{
    df <- sca_data %>%
      select(-terms, -coef, -se, -statistic, -p, -sig.level) %>%
      pivot_longer(-index,
                   names_to="control", values_to="value") %>%
      filter(value==1) %>%
      mutate(controlID = with(.,match(control, unique(control)))) %>%
      select(-value)
  }

  df_labels <- df %>% select(control, controlID) %>% unique()

  return(list(df, setNames(as.character(df_labels$control),
                           df_labels$controlID)))
}



# This function takes the following arguments:
#   data = a dataframe with our data
#   formula = a formula object with our regression formula
#   n_x = the number of x variables we have in the model
#   n_samples = the number of times to estimate the model with a random subset
#               of the data
#   sample_size = the number of observations to include in the subset of data
# It returns a list of bootstrapped standard errors


#' Estimates bootstrapped standard errors for regression models
#'
#' @description
#' Takes in a data frame, regression formula, and bootstrapping parameters and
#' estimates bootstrapped standard errors for models with and without fixed
#' effects.
#'
#'
#' @param data A data frame containing the variables provided in `formula`.
#' @param formula A string containing a regression formula, with or without
#'                fixed effects.
#' @param n_x An integer representing the number of independent variables in
#'            the regression.
#' @param n_samples An integer indicating how many times the model should be
#'                  estimated with a random subset of the data.
#' @param sample_size An integer indicating how many observations are in each
#'                    random subset of the data
#'
#' @return A named list containing bootstrapped standard errors for each
#'         coefficient.
#' @export
#'
#' @examples
#'
#' se_boot(data = bottles, formula = "Salnty ~ T_degC + ChlorA + O2Sat",
#'         n_x = 3, n_samples = 4, sample_size = 300)
#'
#' se_boot(data = data.frame(x1 = rnorm(50000, mean=4, sd=10),
#'                           x2 = rnorm(50000, sd=50),
#'                           ID = rep(1:100, 500),
#'                           area = rep(1:50, 1000),
#'                           y = rnorm(50000)),
#'         formula = "y ~ x1 + x2 | ID",
#'         n_x = 2, n_samples = 10, sample_size = 1000)
#'
se_boot <- function(data, formula, n_x, n_samples, sample_size){

  # Check for fixed effects in the formula
  FE <- ifelse(grepl("|", formula, fixed=T), T, F)

  # Create a matrix to store the coefficient estimates, each row contains
  # coefficients estimated from a different subset of the data
  # ncol=n_x+1 when fixed FE aren't present because
  # we are also storing the intercept estimate
  # ncol=n_x when FE are present because feols() does not report the
  # intercept
  coefs <- matrix(nrow=n_samples, ncol=ifelse(FE, n_x, n_x+1))

  # Loop n_samples times, i.e. how many times we want to re-estimate the model.
  # In the future this should be vectorized.
  for(i in 1:n_samples){
    # Estimate the model with a random subset of the data
    # sample_n is a function from the dplyr package that gives us random
    # rows from a dataframe
    model <- tryCatch(
      {
        if(FE){
          suppressMessages(feols(as.formula(formula),
                                 sample_n(data, sample_size)))
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
    # If the model can't be estimated then return NULL
    if(is.null(model)) return(NULL)
    # Catch instances where coefficients aren't estimated due to
    # collinearity. We can't in good conscience estimate bootstrapped SEs
    # from different numbers of coefficients, it could bias the SEs.
    else if(FE & length(model$coefficients) != n_x){

      message(paste0("Estimation failed due to collinearity for ",
                     paste(model$collin.var, collapse=", "),
                     " during bootstrap with n_samples=",
                     n_samples, " and sample_size=", sample_size,
                     ".\nConsider respecifying bootstrap parameters.\n"))
      return(NULL)
    }
    else if(!FE & length(model$coefficients) != n_x+1){

      message(paste0("Estimation failed due to collinearity for ",
                     paste(names(model$coefficients[is.na(model$coefficients)]),
                           collapse=", "),
                     " during bootstrap with n_samples=",
                     n_samples, " and sample_size=", sample_size,
                     ".\nConsider respecifying bootstrap parameters.\n"))
      return(NULL)
    }

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
