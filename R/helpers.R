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


# Takes in the output of `sca()` and returns a list with the dataframe and
# labels to make a plot to visualize the controls included in each spec curve
# model
# Arguments:
#   spec_data = dataframe object with output from `sca()`
# Returns: list containing dataframe, controls, and control IDs

#' Prepares the output of `sca()` for plotting.
#'
#' @description
#' Takes in the data frame output by `sca()` and returns a list with the data
#' frame and labels to make a plot to visualize the controls included in each
#' spec curve model.
#'
#' @param spec_data A data frame output by `sca`.
#'
#' @return A list containing a data frame, control coefficients, and control
#'         names.
#'
#' @export
#'
#' @examples
#' scp(sca(y = "Salnty", x = "T_degC", controls = c("ChlorA", "O2Sat"),
#'         data = bottles, progressBar=TRUE, parallel=FALSE));
scp <- function(spec_data){
  if("control_coefs" %in% names(spec_data)){
    df <- spec_data %>%
      select(-terms, -coef, -se, -statistic, -p, -sig.level) %>%
      pivot_longer(-c(index, control_coefs),
                   names_to="control", values_to="value") %>%
      filter(value==1) %>%
      mutate(controlID = with(.,match(control, unique(control)))) %>%
      select(-value)
  }
  else{
    df <- spec_data %>%
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
