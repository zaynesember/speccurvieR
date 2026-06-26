# Helper functions--------------------------------------------------------------

# Internal: the control-indicator columns of an sca() data frame, i.e. the 0/1
# columns naming each control, found by removing the known result columns.
sca_control_cols <- function(sca_data){
  meta <- c("coef", "se", "statistic", "p", "RMSE", "adjR", "AIC", "deviance",
            "terms", "control_coefs", "sig.level", "index", "n_obs")
  setdiff(names(sca_data), meta)
}

# Internal: recover the focal independent variable from an sca() result. Every
# specification contains the focal variable but each control is absent from some
# specifications, so the focal variable is the single non-intercept term common
# to every specification's `terms`. Returns NA (with a warning) if it cannot be
# uniquely resolved.
sca_focal_var <- function(sca_data){
  # Prefer the focal variable stored by sca() (unambiguous even for a
  # single-control curve); fall back to the terms heuristic for older objects.
  ax <- attr(sca_data, "x", exact = TRUE)
  if(!is.null(ax) && length(ax) == 1L){
    return(ax)
  }
  if(is.null(sca_data$terms)){
    return(NA_character_)
  }
  common <- setdiff(Reduce(intersect, sca_data$terms), "(Intercept)")
  if(length(common) != 1L){
    warning("Could not uniquely identify the focal variable from `sca_data`.",
            call. = FALSE)
    return(NA_character_)
  }
  common
}

# Internal: stop with an informative message if any of `cols` are absent from
# `data`. `what` labels the offending argument in the error message.
check_columns <- function(data, cols, what){
  missing <- setdiff(cols, colnames(data))
  if(length(missing) > 0){
    stop(what, " not found in data: ", paste(missing, collapse=", "),
         call.=FALSE)
  }
}

# Internal: normalise the `cluster` argument of se_compare() into a cleaned list
# of clustering specifications, one per requested set of clustering dimensions.
#
# `cluster` may be NULL, a character vector (each element a SEPARATE one-way
# clustering, the historical behaviour), or a list of character vectors (each
# element clustered JOINTLY -- one-way when length 1, multiway when longer). A
# character vector is normalised with as.list(), so a one-way request flows
# through exactly the same path whether written `"a"`, `c("a", "b")`, or
# `list("a", "b")`.
#
# Within each specification the dimensions are de-duplicated and sorted: a
# multiway standard error is invariant to the order of its clustering
# dimensions, so sorting makes the column label canonical and lets duplicate
# specifications (e.g. `c("a", "b")` and `c("b", "a")`) collapse to one. Unknown
# columns are dropped with a single warning (matching the historical message),
# empty specifications are dropped, and duplicate specifications are removed.
# Returns NULL if nothing usable remains. The order of distinct specifications
# follows the input, so existing one-way column order is preserved.
normalize_cluster_spec <- function(cluster, data_cols){
  if(is.null(cluster)) return(NULL)
  if(is.character(cluster)) cluster <- as.list(cluster)
  if(!is.list(cluster) ||
     !all(vapply(cluster, is.character, logical(1)))){
    stop("`cluster` must be NULL, a character vector, or a list of character ",
         "vectors of column names.", call. = FALSE)
  }

  # One warning for every unknown clustering variable across all specifications.
  unknown <- setdiff(unique(unlist(cluster)), data_cols)
  if(length(unknown) > 0){
    warning(paste0(unknown, " not a valid clustering variable, ignoring.",
                   collapse = "\n"), call. = FALSE)
  }

  specs <- lapply(cluster, function(d){
    d <- sort(unique(d[d %in% data_cols]))
    if(length(d) == 0L) NULL else d
  })
  specs <- specs[!vapply(specs, is.null, logical(1))]
  if(length(specs) == 0L) return(NULL)
  # Drop duplicate specifications by their canonical (sorted) dimension key.
  keys <- vapply(specs, paste, character(1), collapse = "_BY_")
  specs[!duplicated(keys)]
}

# Internal: resolve the `family`/`link` arguments shared by sca() and
# se_compare() into a normalised family string and (for glm families) a family
# object. Treats the common alias "gaussian" as ordinary least squares
# ("linear"); for any other family it builds the family object, defaulting to
# the family's canonical link when `link` is NULL and erroring clearly on an
# unrecognised family. Returns a list with `family` (the normalised string) and
# `fam_obj` (the family object, or NULL for the linear case).
resolve_family <- function(family, link){
  if(family=="gaussian") family <- "linear"

  fam_obj <- NULL
  if(family!="linear"){
    fam_fun <- tryCatch(match.fun(family),
                        error=function(e)
                          stop("'", family,
                               "' is not a recognised model family.",
                               call.=FALSE))
    fam_obj <- if(is.null(link)) fam_fun() else fam_fun(link=link)
  }

  list(family=family, fam_obj=fam_obj)
}

# Internal: decompose a model formula into the string components sca() uses.
# The response is y; the FIRST right-hand-side term is the focal independent
# variable x; remaining terms are controls; anything after a `|` is treated as
# fixed effects (as in fixest). Right-hand-side terms are split on top-level `+`
# so interaction units (e.g. a:b or a*b) are kept whole, and spaces around the
# interaction operators are removed so `a * b` matches the vector interface's
# "a*b".
formula_to_args <- function(formula){
  if(!inherits(formula, "formula") || length(formula) != 3){
    stop("`formula` must be a two-sided formula, e.g. y ~ x + control1.",
         call.=FALSE)
  }

  norm <- function(expr){
    gsub("\\s*([*:^])\\s*", "\\1", paste(deparse(expr), collapse=""))
  }

  split_plus <- function(expr){
    if(length(expr) == 3 && identical(expr[[1]], as.name("+"))){
      c(split_plus(expr[[2]]), split_plus(expr[[3]]))
    }
    else{
      norm(expr)
    }
  }

  y <- norm(formula[[2]])

  rhs <- formula[[3]]
  fixed_effects <- NULL
  if(length(rhs) == 3 && identical(rhs[[1]], as.name("|"))){
    fixed_effects <- norm(rhs[[3]])
    rhs <- rhs[[2]]
  }

  terms <- split_plus(rhs)
  if(length(terms) < 2){
    stop("The formula must supply a focal independent variable and at least ",
         "one control, e.g. y ~ x + control1.", call.=FALSE)
  }

  list(y = y, x = terms[1], controls = terms[-1], fixed_effects = fixed_effects)
}

#' Builds models formulae with every combination of control variables possible.
#'
#' @param y A string containing the dependent variable name.
#' @param x A string containing the independent variable name.
#' @param controls A vector of strings containing control variable names.
#' @param fixed_effects A string containing the name of a variable to use for
#'                     fixed effects, defaults to `NA` indicating no fixed
#'                     effects desired.
#' @param ... Deprecated camelCase arguments (`fixedEffects`); use
#'            `fixed_effects` instead.
#'
#' @return A vector of formula objects using every possible combination of
#'         controls.
#'
#' @keywords internal
#' @export
#'
#' @examples
#' formula_builder("dependentVariable", "independentVariable",
#'                 c("control1", "control2"));
#' formula_builder("dependentVariable", "independentVariable",
#'                 c("control1*control2"), fixed_effects="month");
formula_builder <- function(y, x, controls, fixed_effects=NA, ...){

  # Backward compatibility: translate the deprecated `fixedEffects` argument.
  .dots <- list(...)
  if("fixedEffects" %in% names(.dots)){
    .Deprecated(msg="`fixedEffects` is deprecated; use `fixed_effects`.")
    fixed_effects <- .dots$fixedEffects
  }

  # Get all combinations of controls
  powerset <- unlist(lapply(1:length(controls),
                            combinat::combn,
                            x = controls,
                            simplify = FALSE),
                     recursive=FALSE)

  # Remove duplicate controls that are already in the interaction
  powerset <- unique(sapply(X=powerset, FUN=duplicate_remover, x=x))

  # Build right hand side of the formulae
  if(is.na(fixed_effects)){
    RHS <- unique(sapply(powerset, paste_factory, x))
  }
  else{
    RHS <- paste(unique(sapply(powerset, paste_factory, x)), fixed_effects,
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
#' @keywords internal
#' @export
#'
#' @examples
#' paste_factory(controls = c("control1", "control2"),
#'               x = "independentVariable");
paste_factory <- function(controls, x){
  # Prepend the focal variable x to the right-hand side unless it already
  # appears among the controls (as a standalone term, or as a component of an
  # interaction such as "x*a"). Membership is tested by exact term equality,
  # splitting each control on its interaction operators, NOT by str_detect()'s
  # regex substring test: the old test treated x as already present whenever a
  # control name merely contained x as a substring (e.g. x = "Temp" with a
  # control "TempX") or x held regex metacharacters, which silently dropped the
  # focal term from the formula and crashed sca() with "subscript out of bounds"
  # during coefficient extraction.
  control_terms <- unique(trimws(unlist(strsplit(controls, "[*:]"))))
  if(x %in% control_terms){
    paste(controls, collapse=" + ")
  } else {
    paste(x, paste(controls, collapse=" + "), sep=" + ")
  }
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
#' @keywords internal
#' @export
#'
#' @examples
#' duplicate_remover(controls = c("control1", "control2*control3"),
#'                   x = "independentVariable");
duplicate_remover <- function(controls, x){
  # When the focal variable x appears inside an interaction control (e.g.
  # "x*b"), its interaction partner ("b") is already implied by that term, so a
  # redundant standalone copy of the partner is dropped. Matching is by exact
  # term -- splitting each interaction on its "*" operator and comparing parts
  # for equality -- never by str_replace()'s regex substitution of x, which
  # mishandled a focal name that was a substring of, or shared regex
  # metacharacters with, another term.
  has_interaction <- str_detect(controls, fixed("*"))
  if(any(has_interaction)){
    partners <- unlist(lapply(controls[has_interaction], function(term){
      parts <- trimws(strsplit(term, "*", fixed=TRUE)[[1]])
      if(x %in% parts) setdiff(parts, x) else character(0)
    }))
    return(controls[!controls %in% partners])
  }
  controls
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
#' @keywords internal
#' @export
#'
#' @examples
#' m <- summary(lm(Salnty ~ STheta + T_degC, bottles))
#' control_extractor(model = m, x = "STheta");
#'
#' m <- summary(lm(Salnty ~ STheta*T_degC + O2Sat, bottles))
#' control_extractor(model = m, x = "STheta");
control_extractor <- function(model, x, feols_model=FALSE){
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
#' @keywords internal
#' @export
#'
#' @examples
#' un_as_is(x = I(c(1:4)));
un_as_is <- function(x) {
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
#' @keywords internal
#' @export
#'
#' @examples
#' scp(sca(y = "Salnty", x = "T_degC", controls = c("ChlorA", "O2Sat"),
#'         data = bottles, progress_bar=TRUE, parallel=FALSE));
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
#' effects. The model is refit on `n_samples` resamples of `sample_size`
#' observations drawn with replacement, and the standard deviation of the
#' resampled coefficients is rescaled by `sqrt(sample_size / nrow(data))` to
#' estimate the standard error of the full-sample estimator (an m-out-of-n
#' bootstrap; when `sample_size` equals `nrow(data)` this is the ordinary
#' nonparametric bootstrap).
#'
#'
#' @param data A data frame containing the variables provided in `formula`.
#' @param formula A string containing a regression formula, with or without
#'                fixed effects.
#' @param n_x An integer representing the number of independent variables in
#'            the regression.
#' @param n_samples An integer indicating how many bootstrap resamples to draw,
#'                  i.e. how many times the model is refit.
#' @param sample_size An integer indicating how many observations are drawn
#'                    (with replacement) in each bootstrap resample.
#' @param weights Optional string with the column name in `data` that contains
#'                weights.
#' @param fam_obj Optional `family` object (as returned by e.g. `binomial()`)
#'                used to fit each resampled model with `glm()`. Defaults to
#'                `NULL`, fitting linear models with `lm()` (or `feols()` when
#'                the formula contains fixed effects). Fixed effects are not
#'                supported together with a `family` object.
#'
#' @return A named list containing bootstrapped standard errors for each
#'         coefficient.
#' @keywords internal
#' @export
#'
#' @examples
#'
#' se_boot(data = bottles, formula = "Salnty ~ T_degC + ChlorA + O2Sat",
#'         n_x = 3, n_samples = 4, sample_size = 300)
#'
#' \donttest{
#' se_boot(data = data.frame(x1 = rnorm(50000, mean=4, sd=10),
#'                           x2 = rnorm(50000, sd=50),
#'                           ID = rep(1:100, 500),
#'                           area = rep(1:50, 1000),
#'                           y = rnorm(50000)),
#'         formula = "y ~ x1 + x2 | ID",
#'         n_x = 2, n_samples = 10, sample_size = 1000)
#' }
#'
se_boot <- function(data, formula, n_x, n_samples, sample_size, weights=NULL,
                    fam_obj=NULL){

  # Check for fixed effects in the formula
  FE <- ifelse(grepl("|", formula, fixed=TRUE), TRUE, FALSE)

  # A non-NULL family object selects the glm() estimation path. Fixed effects
  # are not supported alongside a glm family (se_compare() strips them before
  # calling), so the glm branch is always non-FE.
  is_glm <- !is.null(fam_obj)

  # Create a list of NAs to return for cases when bootstrapping fails
  fallback_list <- as.list(rep(NA, n_x + 1))

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
    # Estimate the model on a bootstrap resample of `sample_size` rows drawn
    # WITH replacement. slice_sample() is the dplyr replacement for the
    # superseded sample_n(). Resampling with replacement (rather than drawing a
    # subset without replacement) is what makes this a bootstrap: the standard
    # deviation of the resampled coefficients estimates the sampling
    # variability of the estimator at sample size `sample_size`, which is
    # rescaled to the full-sample standard error below.
    model <- tryCatch(
      {
        if(FE){
          if(is.null(weights)){
            suppressMessages(feols(as.formula(formula),
                                   slice_sample(data, n=sample_size, replace=TRUE)))
          }
          else{
            suppressMessages(feols(as.formula(formula),
                                   slice_sample(data, n=sample_size, replace=TRUE),
                                   weights=data[[weights]]))
          }
        }
        else{
          if(is.null(weights)){
            if(is_glm){
              suppressMessages(glm(as.formula(formula),
                                   slice_sample(data, n=sample_size, replace=TRUE),
                                   family=fam_obj))
            }
            else{
              suppressMessages(lm(as.formula(formula),
                                  slice_sample(data, n=sample_size, replace=TRUE)))
            }
          }
          else{
            fmla <- as.formula(formula)
            environment(fmla) <- environment()
            if(is_glm){
              suppressMessages(glm(fmla,
                                   slice_sample(data, n=sample_size, replace=TRUE),
                                   family=fam_obj,
                                   weights=get(weights)))
            }
            else{
              suppressMessages(lm(fmla,
                                  slice_sample(data, n=sample_size, replace=TRUE),
                                  weights=get(weights)))
            }
          }
        }
      },
      error=function(cond){
        if(FE){
          message(paste0("Estimation failed during bootstrap for fixed ",
                         "effects model with n_samples=",
                         n_samples, " and sample_size=", sample_size,
                         ".\nConsider respecifying bootstrap parameters or ",
                         "model.\n"))
        }
        else{
        message(paste0("Estimation failed during bootstrap for non-fixed ",
                       "effects model with n_samples=",
                       n_samples, " and sample_size=", sample_size,
                       ".\nConsider respecifying bootstrap parameters or ",
                       "model.\n"))
        }

        return(fallback_list)
      }
    )
    # If the model can't be estimated then return NULL
    if(is.null(model)) return(fallback_list)
    # Catch instances where coefficients aren't estimated due to
    # collinearity. We can't in good conscience estimate bootstrapped SEs
    # from different numbers of coefficients, it could bias the SEs.
    else if(FE & length(model$coefficients) != n_x){

      message(paste0("Estimation failed due to collinearity for ",
                     paste(model$collin.var, collapse=", "),
                     " during bootstrap of fixed effects model with n_samples=",
                     n_samples, " and sample_size=", sample_size,
                     ".\nConsider respecifying bootstrap parameters.\n"))
      return(fallback_list)
    }
    else if(!FE & length(model$coefficients) != n_x+1){

      message(paste0("Estimation failed due to collinearity for ",
                     paste(names(model$coefficients[is.na(model$coefficients)]),
                           collapse=", "),
                     " during bootstrap of non-fixed effects model
                     with n_samples=",
                     n_samples, " and sample_size=", sample_size,
                     ".\nConsider respecifying bootstrap parameters.\n"))
      return(fallback_list)
    }

    # Store the coefficient estimates in row i of our matrix
    coefs[i,] <- model$coefficients
  }

  # The standard deviation of the resampled coefficients estimates the
  # sampling variability of the estimator at sample size `sample_size`. For an
  # m-out-of-n bootstrap (sample_size = m < n), that overstates the
  # full-sample standard error by a factor of roughly sqrt(n / m), so we
  # rescale by sqrt(m / n) to recover the standard error of the estimator fit
  # to the full data. When sample_size == n this factor is 1 and the procedure
  # reduces to the ordinary nonparametric bootstrap.
  scale <- sqrt(sample_size / nrow(data))
  retVal <- apply(coefs, FUN=sd, MARGIN=2) * scale

  names(retVal) <- names(model$coefficients)

  if(FE) retVal <- c("(Intercept)"=NA, unlist(retVal))

  return(retVal)
}

# Internal: the colour mapping for the significance bins used across the
# package's plots. Distinct hues -- deep purple for the strongest evidence,
# through teal and amber, to a muted warm grey for non-significant -- so the
# bins are easy to tell apart on thin error bars, and the strength of evidence
# reads from the colour. The three significant hues are colour-blind
# distinguishable.
sca_sig_colors <- function(){
  c("p < .005" = "#4A1486",
    "p < .05"  = "#1D9E75",
    "p < .1"   = "#E69F00",
    "p >= .1"  = "#B4B2A9")
}

#' A clean, consistent ggplot2 theme for speccurvieR plots.
#'
#' @description
#' `theme_sca()` is the shared theme applied by the package's plotting
#' functions. It is exported so the same look can be reused or tweaked when
#' customising the plots they return.
#'
#' @param base_size Base font size, passed to [ggplot2::theme_minimal()].
#'                  Defaults to `11`.
#' @param base_family Base font family, passed to [ggplot2::theme_minimal()].
#'                    Defaults to the value of the `speccurvieR.base_family`
#'                    option, or `""` (the graphics device's default font) if
#'                    that option is unset, so plots stay portable across
#'                    machines. Set the option (e.g.
#'                    `options(speccurvieR.base_family = "Roboto")`) or pass a
#'                    family directly to use a specific font.
#'
#' @return A ggplot2 theme object.
#'
#' @export
#'
#' @examples
#' library(ggplot2)
#' ggplot(bottles, aes(T_degC, Salnty)) + geom_point() + theme_sca();
theme_sca <- function(base_size = 11,
                      base_family = getOption("speccurvieR.base_family", "")){
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
      plot.title       = element_text(face = "bold"),
      legend.position  = "top",
      legend.title     = element_blank(),
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold")
    )
}
