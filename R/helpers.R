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


#' Removes the `AsIs` class attribute from the input. Taken from:
#' https://stackoverflow.com/a/12866609
#'
#' @param x An object with the `AsIs` class attribute.
#'
#' @return An object without the `AsIs` class attribute.
#'
#' @examples
#' unAsIs(I(c(1:4)))
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

#' Takes in the output of `sca` and returns a list with the data frame and
#' labels to make a plot to visualize the controls included in each spec curve
#' model.
#'
#' @param spec_data A data frame output by `sca`.
#'
#' @return A list containing a data frame, control coefficients, and control
#'         names.
#'
#' @examples
#' scp(sca(y="Salnty", x="T_degC", c("ChlorA", "O2Sat"), data=bottles,
#'     progressBar=T, parallel=FALSE))
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
