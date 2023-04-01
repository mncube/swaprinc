#' Swap in Principal Components
#'
#'
#'
#' @param data A dataframe
#' @param formula A quoted model formula
#' @param engine The engine for fitting the model.  Options are "stats", "lme4",
#' or "nlme".
#' @param pca_vars Variables to include in the principal component analysis.
#' These variables will be swapped out for principal components
#' @param n_pca_components The number of principal components to include in the
#' model
#' @param ... Pass additional arguments arguments
#'
#' @return A list with fitted models
#' @export
#'
#' @examples
#' data(iris)
#' res <- swaprinc(iris,
#' "Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width",
#' pca_vars = c("Sepal.Width", "Petal.Length", "Petal.Width"),
#' n_pca_components = 2)
swaprinc <- function(data, formula, engine = "stats", pca_vars,
                     n_pca_components, ...) {
  # Helper function for model fitting
  fit_model <- function(data, formula, engine, ...) {
    if (engine == "stats") {
      lm_model <- try(stats::lm(formula, data), silent = TRUE)
      if (inherits(lm_model, "lm")) {
        return(lm_model)
      } else {
        glm_model <- try(stats::glm(formula, data, ...), silent = TRUE)
        if (inherits(glm_model, "glm")) {
          return(glm_model)
        } else {
          stop("Neither lm nor glm from stats could fit the model.")
        }
      }
    } else if (engine == "lme4") {
      lmer_model <- try(lme4::lmer(formula, data, ...), silent = TRUE)
      if (inherits(lmer_model, "merMod")) {
        return(lmer_model)
      } else {
        glmer_model <- try(lme4::glmer(formula, data, ...), silent = TRUE)
        if (inherits(glmer_model, "merMod")) {
          return(glmer_model)
        } else {
          stop("Neither lmer nor glmer from lme4 could fit the model.")
        }
      }
    } else if (engine == "nlme") {
      lme_model <- try(nlme::lme(formula, data, ...), silent = TRUE)
      if (inherits(lme_model, "lme")) {
        return(lme_model)
      } else {
        nlme_model <- try(nlme::nlme(formula, data, ...), silent = TRUE)
        if (inherits(nlme_model, "nlme")) {
          return(nlme_model)
        } else {
          stop("Neither lme nor nlme from nlme could fit the model.")
        }
      }
    } else {
      stop("Invalid engine specified.")
    }
  }

  # Fit the regular model
  model_raw <- fit_model(data, formula, engine, ...)

  # Perform PCA
  pca_data <- data[, pca_vars]
  pca_result <- stats::princomp(pca_data, cor = TRUE)
  pca_scores <- pca_result$scores[, 1:n_pca_components]
  colnames(pca_scores) <- paste0("PC", 1:n_pca_components)

  # Replace the original variables with the principal components
  data_pca <- data %>%
    dplyr::select(-tidyselect::one_of(pca_vars)) %>%
    cbind(pca_scores)

  # # Generate the string with the principal components
  # pca_terms <- paste("PC", 1:n_pca_components, sep = "", collapse = " + ")
  #
  # # Create the new formula string
  # formula_pca_str <- paste(as.character(formula), "+", pca_terms)

  replace_pca_vars_in_formula <- function(formula, pca_vars, pca_terms) {
    formula_terms <- all.vars(stats::as.formula(formula))
    new_formula_terms <- base::setdiff(formula_terms, pca_vars)
    new_formula_terms <- c(new_formula_terms, unlist(strsplit(pca_terms, " \\+ ")))
    new_formula <- stats::as.formula(paste(formula_terms[1], "~", paste(new_formula_terms[-1], collapse = " + ")))
    return(new_formula)
  }

  # # Generate the string with the principal components
  # pca_terms <- paste("PC", 1:n_pca_components, sep = "", collapse = " + ")

  # Generate the string with the principal components
  pca_terms <- paste("PC", 1:n_pca_components, sep = "", collapse = " + ")

  # Replace the pca_vars in the formula with the pca_terms
  formula_pca <- replace_pca_vars_in_formula(formula, pca_vars, pca_terms)

  # # Create the new formula string
  # formula_pca_str <- paste(as.character(formula), "+", pca_terms, collapse = " ")

  # # Convert the string to a formula
  # formula_pca <- stats::formula(paste(formula_pca_str, collapse = " "))

  # Fit the PCA model
  model_pca <- fit_model(data_pca, formula_pca, engine, ...)

  # Compare models
  # ... (implement comparison metrics here)

  return(list(model_raw = model_raw, model_pca = model_pca))
}

