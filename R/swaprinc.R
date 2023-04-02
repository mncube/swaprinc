#' Swap in Principal Components
#'
#'
#'
#' @param data A dataframe
#' @param formula A quoted model formula
#' @param engine The engine for fitting the model.  Options are "stats" or"lme4".
#' @param pca_vars Variables to include in the principal component analysis.
#' These variables will be swapped out for principal components
#' @param n_pca_components The number of principal components to include in the
#' model
#' @param ... Pass additional arguments
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
  # # Helper function for model fitting
  # fit_model <- function(data, formula, engine, ...) {
  #   if (engine == "stats") {
  #     lm_model <- try(stats::lm(formula, data), silent = TRUE)
  #     if (inherits(lm_model, "lm")) {
  #       return(lm_model)
  #     } else {
  #       glm_model <- try(stats::glm(formula, data, ...), silent = TRUE)
  #       if (inherits(glm_model, "glm")) {
  #         return(glm_model)
  #       } else {
  #         stop("Neither lm nor glm from stats could fit the model.")
  #       }
  #     }
  #   } else if (engine == "lme4") {
  #     lmer_model <- try(lme4::lmer(formula, data, ...), silent = TRUE)
  #     if (inherits(lmer_model, "merMod")) {
  #       return(lmer_model)
  #     } else {
  #       glmer_model <- try(lme4::glmer(formula, data, ...), silent = TRUE)
  #       if (inherits(glmer_model, "merMod")) {
  #         return(glmer_model)
  #       } else {
  #         stop("Neither lmer nor glmer from lme4 could fit the model.")
  #       }
  #     }
  #   } else {
  #     stop("Invalid engine specified.")
  #   }
  # }

  # Helper function for model fitting
  fit_model <- function(data, formula, engine, ...) {
    if (engine == "stats") {
      glm_model <- try(stats::glm(formula, data, ...), silent = TRUE)
      if (inherits(glm_model, "glm")) {
        return(glm_model)
      } else {
        lm_model <- try(stats::lm(formula, data), silent = TRUE)
        if (inherits(lm_model, "lm")) {
          return(lm_model)
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
  # data_pca <- data %>%
  #   dplyr::select(-tidyselect::one_of(pca_vars), dplyr::everything()) %>%
  #   cbind(pca_scores)

  replace_pca_vars_in_formula <- function(formula, pca_vars, pca_terms) {
    original_formula <- stats::as.formula(formula)
    response_var <- all.vars(original_formula)[1]

    # Separate fixed and random effects using gsub and strsplit
    #https://stackoverflow.com/questions/62966793/how-to-extract-just-the-random-effects-part-of-the-formula-from-lme4
    inp <- deparse(original_formula)
    formula_terms <- gsub(" ", "", unlist(strsplit(inp, "+", fixed = T)), fixed = T)

    # Identify and separate fixed and random effects
    fixed_effects <- formula_terms[!grepl("\\|", formula_terms)]
    random_effects_terms <- formula_terms[grepl("\\|", formula_terms)]

    # Remove the response variable from the fixed_effects
    fixed_effects[1] <- gsub(paste0(response_var, "~"), "", fixed_effects[1])

    # Remove pca_vars from the fixed_effects
    fixed_effects <- base::setdiff(fixed_effects, pca_vars)

    # Combine fixed_effects with the PCA terms
    fixed_effects <- c(fixed_effects, unlist(strsplit(pca_terms, " \\+ ")))

    if (length(random_effects_terms) > 0) {

      # Check that random effect variables are not included as principal components
      term_matrix <- attr(stats::terms(original_formula), "factors")
      random_effects_vars <- rownames(term_matrix)[apply(term_matrix, 1, function(x) any(x == 1))]
      random_effects_vars <- grep("\\|", random_effects_vars, value = TRUE)
      random_effects_vars <- unlist(strsplit(gsub("\\||\\(|\\)", "", random_effects_vars), split = " "))

      if (any(random_effects_vars %in% pca_vars)) {
        rlang::abort("Using a random effect variable as one of the pca_vars is not allowed.")
      }


      new_formula <- paste(response_var, "~", paste(fixed_effects, collapse = " + "), "+", paste(random_effects_terms, collapse = " + "))
    } else {
      new_formula <- paste(response_var, "~", paste(fixed_effects, collapse = " + "))
    }

    return(stats::as.formula(new_formula))
  }

  # Generate the string with the principal components
  pca_terms <- paste("PC", 1:n_pca_components, sep = "", collapse = " + ")

  # Replace the pca_vars in the formula with the pca_terms
  formula_pca <- replace_pca_vars_in_formula(formula, pca_vars, pca_terms)

  # Fit the PCA model
  model_pca <- fit_model(data_pca, formula_pca, engine, ...)

  # Compare models
  # ... (implement comparison metrics here)

  return(list(model_raw = model_raw, model_pca = model_pca))
}

