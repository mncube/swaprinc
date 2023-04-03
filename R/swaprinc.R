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
#' @param norun_raw Include regression on raw variables if TRUE, exclude if FALSE.
#' @param center Set center parameter for prcomp
#' @param scale. Set scale. parameter for prcomp
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
                     n_pca_components, norun_raw = FALSE, center = TRUE,
                     scale. = FALSE, ...) {

  # Helper function for model fitting
  fit_model <- function(data, formula, engine, ...) {
    if (engine == "stats") {
      glm_model <- try(stats::glm(formula, data, ...), silent = TRUE)
      if (inherits(glm_model, "glm")) {
        return(glm_model)
      } else {
        lm_model <- try(stats::lm(formula, data, ...), silent = TRUE)
        if (inherits(lm_model, "lm")) {
          return(lm_model)
        } else {
          rlang::abort("Neither lm nor glm from stats could fit the model.")
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
          rlang::abort("Neither lmer nor glmer from lme4 could fit the model.")
        }
      }
    } else {
      rlang::abort("Invalid engine specified.")
    }
  }

  # Fit the regular model conditionally
  if (!norun_raw) {
    model_raw <- fit_model(data, formula, engine, ...)
  } else {
    model_raw <- NULL
  }

  # Perform PCA
  pca_data <- data[, pca_vars]
  pca_result <- stats::prcomp(pca_data, center = center, scale. = scale., ...)
  pca_scores <- pca_result$x[, 1:n_pca_components]
  colnames(pca_scores) <- paste0("PC", 1:n_pca_components)

  # Replace the original variables with the principal components
  data_pca <- data %>%
    dplyr::select(-tidyselect::one_of(pca_vars)) %>%
    cbind(pca_scores)

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

  #Compare Models
  compare_models <- function(model_raw, model_pca) {
    # Tidy model output
    if (is.null(model_raw)){
      raw_summary <- NULL
    } else if (inherits(model_raw, "merMod")) {
      raw_summary <- broom.mixed::glance(model_raw)
    } else {
      raw_summary <- broom::glance(model_raw)
    }

    if (inherits(model_pca, "merMod")) {
      pca_summary <- broom.mixed::glance(model_pca)
    } else {
      pca_summary <- broom::glance(model_pca)
    }

    # Create comparison metrics data frame
    if(is.null(model_raw) & inherits(model_pca, c("glm", "glmerMod"))) {
      # For glm and glmer models
      pca_mcfadden_pseudo_r_squared <- 1 - (pca_summary$logLik / pca_summary$null.deviance)

      comparison <- data.frame(
        model = c("PCA"),
        pseudo_r_squared = c(pca_mcfadden_pseudo_r_squared),
        AIC = c(pca_summary$AIC),
        BIC = c(pca_summary$BIC)
      )
    } else if (inherits(model_raw, c("glm", "glmerMod")) & inherits(model_pca, c("glm", "glmerMod"))) {
      # For glm and glmer models
      raw_mcfadden_pseudo_r_squared <- 1 - (raw_summary$logLik / raw_summary$null.deviance)
      pca_mcfadden_pseudo_r_squared <- 1 - (pca_summary$logLik / pca_summary$null.deviance)

      comparison <- data.frame(
        model = c("Raw", "PCA"),
        pseudo_r_squared = c(raw_mcfadden_pseudo_r_squared, pca_mcfadden_pseudo_r_squared),
        AIC = c(raw_summary$AIC, pca_summary$AIC),
        BIC = c(raw_summary$BIC, pca_summary$BIC)
      )
    } else if (is.null(model_raw) & inherits(model_pca, "lm")) {
      # For lm models only
      comparison <- data.frame(
        model = c("PCA"),
        r_squared = c(pca_summary$r.squared),
        adj_r_squared = c(pca_summary$adj.r.squared),
        AIC = c(pca_summary$AIC),
        BIC = c(pca_summary$BIC)
      )
    } else if (inherits(model_raw, "lm") & inherits(model_pca, "lm")) {
      # For lm models only
      comparison <- data.frame(
        model = c("Raw", "PCA"),
        r_squared = c(raw_summary$r.squared, pca_summary$r.squared),
        adj_r_squared = c(raw_summary$adj.r.squared, pca_summary$adj.r.squared),
        AIC = c(raw_summary$AIC, pca_summary$AIC),
        BIC = c(raw_summary$BIC, pca_summary$BIC)
      )
    } else if (is.null(model_raw) & inherits(model_pca, "lmerMod")){
      # For lmer models only
      comparison <- data.frame(
        model = c("PCA"),
        logLik = c(pca_summary$logLik),
        AIC = c(pca_summary$AIC),
        BIC = c(pca_summary$BIC)
      )
    } else if (inherits(model_raw, "lmerMod") & inherits(model_pca, "lmerMod")) {
      # For lmer models only
      comparison <- data.frame(
        model = c("Raw", "PCA"),
        logLik = c(raw_summary$logLik, pca_summary$logLik),
        AIC = c(raw_summary$AIC, pca_summary$AIC),
        BIC = c(raw_summary$BIC, pca_summary$BIC)
      )
    } else {
      rlang::abort("The two models must be of the same type, either linear (lm or lmer) or generalized linear (glm or glmer).")
    }

    return(comparison)
  }


  # Get comparison
  model_comparison <- compare_models(model_raw, model_pca)

  return(list(model_raw = model_raw, model_pca = model_pca, comparison = model_comparison))

}
