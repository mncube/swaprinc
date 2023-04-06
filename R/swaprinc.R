#' Swap in Principal Components
#'
#'
#'
#' @param data A dataframe
#' @param formula A quoted model formula
#' @param engine The engine for fitting the model.  Options are 'stats' or 'lme4'.
#' @param prc_eng Then engine or extracting principal components.  Options are
#' 'stats', 'Gifi', and 'stats_Gifi'.  The stats_Gifi engine uses
#' tidyselect::where(is.numeric) to select the pca_vars for stats::prcomp and
#' -tidyselect::where(is.numeric) to select the pca_vars for Gifi::princals.
#' @param pca_vars Variables to include in the principal component analysis.
#' These variables will be swapped out for principal components
#' @param n_pca_components The number of principal components to include in the
#' model. If using a complex prc_eng (i.e., stats_Gifi) then provide a named
#' vector (i.e., n_pca_components = c("stats" = 2, "Gifi" = 3)).
#' @param norun_raw Include regression on raw variables if TRUE, exclude if FALSE.
#' @param center Set center parameter for prcomp
#' @param scale. Set scale. parameter for prcomp
#' @param lpca_center Center data as in the LearnPCA Step-by-Step PCA vignette.  Only
#' numeric variables will be included in the centering.  Parameter takes values
#' 'all' to center raw and pca variables, 'raw' to only center variables for the
#' raw variable model fitting, 'pca' to only center pca_vars before pca regression
#' model fitting, and 'none' to skip lpca centering.
#' @param lpca_scale Scale data as in the LearnPCA Step-by-Step PCA vignette.  Only
#' numeric variables will be included in the scaling.  Parameter takes values
#' 'all' to scale raw and pca variables, 'raw' to only scale variables for the
#' raw variable model fitting, 'pca' to only scale pca_vars before pca regression
#' model fitting, and 'none' to skip lpca scaling.
#' @param lpca_undo Undo centering and scaling of pca_vars as in the LearnPCA
#' Step-by-Step PCA vignette.
#' @param gifi_transform Use Gifi optimal scaling to transform a set of variables.
#' Parameter takes values 'none', 'all', 'raw', and 'pca'
#' @param gifi_trans_vars A vector of variables to include in the Gifi optimal
#' scaling transformation
#' @param gifi_trans_dims Number of dimensions to extract in the Gifi optimal
#' scaling transformation algorithm
#' @param no_tresp When set to TRUE, no_tresp (No transform response) will exclude
#' the response variable from from pre-modeling and pre-pca transformations.
#' Specifically, setting no_tresp to TRUE will exclude the response variable from
#' the transformation specified in lpca_center and lpca_scale.
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
swaprinc <- function(data, formula, engine = "stats", prc_eng = "stats", pca_vars,
                     n_pca_components, norun_raw = FALSE, center = TRUE,
                     scale. = FALSE, lpca_center = "none", lpca_scale = "none",
                     lpca_undo = FALSE, gifi_transform = "none", gifi_trans_vars,
                     gifi_trans_dims, no_tresp = FALSE, ...) {
  # Test function parameters
  if (!(lpca_center == "none" | lpca_center == "all" | lpca_center == "raw" |
        lpca_center == "pca")){
    rlang::abort("lpca_center must be set to: 'none', 'all', 'raw', or 'pca'")
  }

  if (!(lpca_scale == "none" | lpca_scale == "all" | lpca_scale == "raw" |
        lpca_scale == "pca")){
    rlang::abort("lpca_center must be set to: 'none', 'all', 'raw', or 'pca'")
  }

  if(lpca_undo == TRUE & (lpca_scale == "none" | lpca_scale == "raw" |
                                lpca_center == "none" | lpca_center == "raw" )) {
    rlang::abort("To use lpca_undo, lpca_scale and lpca_center must be set
                 to 'all' or 'pca'")
  }

  # Helper function to get numerics
  get_nums <- function(df){
    dplyr::select(df, tidyselect::where(is.numeric))
  }

  # Helper function to bind processed numeric variables to non numeric variables
  bind_nums <- function(df_nums, df){
    dplyr::bind_cols(df_nums, dplyr::select(df, -tidyselect::where(is.numeric)))
  }

  # Helper function to get response for no_tresp
  get_resp <- function(df){
    df_resp <- dplyr::select(df, all.vars(stats::update(stats::as.formula(formula), . ~ 1)))
    df_pred <- dplyr::select(df, -all.vars(stats::update(stats::as.formula(formula), . ~ 1)))
    out <- list("df_resp" = df_resp,
                "df_pred" = df_pred)
  }

  # Helper function to bind response variable to processed variables
  bind_resp <- function(df_resp, df){
    dplyr::bind_cols(df_resp, df)
  }

  # Create helper function for lpca center and scale
  lpca_cs <- function(df, scl, cnt){
  if (no_tresp == TRUE){
    df <- get_resp(df)
    df_resp <- df$df_resp
    df <- df$df_pred
  }
    df <- scale(get_nums(df), scale = scl, center = cnt) %>%
      as.data.frame() %>%
      bind_nums(df)

    if (no_tresp == TRUE){
      df <- bind_resp(df_resp, df)
    }

    return(df)
  }

  # Create helper function to get os_trans vars
  gifi_trans <- function(df, gifi_trans_dims, ...){
    # Split data frame
    dftr <- dplyr::select(df, tidyselect::all_of(gifi_trans_vars))
    df_notr <- dplyr::select(df, -tidyselect::all_of(gifi_trans_vars))

    # Get transformed data
    gifi_trans <- Gifi::princals(dftr, ndim=gifi_trans_dims, ...)

    # Combine transformed and non-transformed data
    df_trans <- gifi_trans$transform %>%
      as.data.frame() %>%
      dplyr::mutate_all(as.numeric) %>%
      cbind(df_notr)

    return(df_trans)
  }

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

  # Transform All Data using Gifi::princals
  if(gifi_transform == "all"){
    data <- gifi_trans(data, gifi_trans_dims, ...)
  }

  # Scale All Data and Raw Data According to LearnPCA
  if(lpca_center == "all"){
    data <- lpca_cs(data, scl = FALSE, cnt = TRUE)
  }

  if(lpca_scale == "all"){
    data <- lpca_cs(data, scl = TRUE, cnt = FALSE)
  }

  # Fit the regular model conditionally
  if (!norun_raw) {
    # Copy data
    df_raw <- data

    # Gifi trans data
    if(gifi_transform == "raw"){
      df_raw <- gifi_trans(df_raw, gifi_trans_dims, ...)
    }

    # Scale raw data only according to LearnPCA
    if(lpca_center == "raw"){
      df_raw <- lpca_cs(df_raw, scl = FALSE, cnt = TRUE)
    }

    if(lpca_scale == "raw"){
      df_raw <- lpca_cs(df_raw, scl = TRUE, cnt = FALSE)
    }

    model_raw <- fit_model(df_raw, formula, engine, ...)
  } else {
    model_raw <- NULL
  }

  # Perform PCA

  # Gifi trans data
  if(gifi_transform == "pca"){
    data <- gifi_trans(data, gifi_trans_dims, ...)
  }

  if(lpca_center == "pca"){
    data <- lpca_cs(data, scl = FALSE, cnt = TRUE)
  }

  if(lpca_scale == "pca"){
    data <- lpca_cs(data, scl = TRUE, cnt = FALSE)
  }

  # Get PCA data

  #Split data based on prc_eng
  split_stats_Gifi <- FALSE
  if (sum(names(pca_vars) == c("stats", "Gifi")) == 2 |
      sum(names(pca_vars) == c("Gifi", "stats")) == 2){
    pca_stats_Gifi_vars <- pca_vars
    pca_vars <- unlist(pca_vars, use.names = FALSE)
    split_stats_Gifi <- TRUE
  }

  pca_data <- data[, pca_vars]

  # Extraction helper functions
  extract_stats <- function(df = pca_data, comps = n_pca_components, ...){
    pca_result <- stats::prcomp(df, center = center, scale. = scale., ...)

    if (lpca_undo == TRUE) {
      Xhat <- pca_result$x[, 1:comps] %*% t(pca_result$rotation[, 1:comps])
      Xhat <- scale(Xhat, center = FALSE, scale = 1/pca_result$scale)
      pca_scores <- scale(Xhat, center = -pca_result$center, scale = FALSE)
    } else {
      pca_scores <- pca_result$x[, 1:comps]
    }
  }

  extract_Gifi <- function(df = pca_data, comps = n_pca_components, ...){
    gifi_results <- Gifi::princals(df, ndim=comps, ...)
    pca_scores <- gifi_results$objectscores
  }

  # Run prc_eng
  if (prc_eng == "stats"){
    pca_scores <- extract_stats()
  } else if (prc_eng == "Gifi") {
    pca_scores <- extract_Gifi()
  } else if (prc_eng == "stats_Gifi"){

    # Split pca_data by
    if (split_stats_Gifi){
      pca_data_stats <- pca_data %>% dplyr::select(pca_stats_Gifi_vars[["stats"]])
      pca_data_Gifi <- pca_data %>% dplyr::select(pca_stats_Gifi_vars[["Gifi"]])
    } else {
      pca_data_stats <- pca_data %>% dplyr::select(tidyselect::where(is.numeric))
      pca_data_Gifi <- pca_data %>% dplyr::select(-tidyselect::where(is.numeric))
    }

    #stats
    pca_scores_stats <- extract_stats(df = pca_data_stats,
                                      comps = n_pca_components[["stats"]])
    #Gifi
    pca_scores_Gifi <- extract_Gifi(df = pca_data_Gifi,
                                    comps = n_pca_components[["Gifi"]])
    #Collapse
    pca_scores <- cbind(pca_scores_stats, pca_scores_Gifi)
    n_pca_components <- sum(n_pca_components)
  }else {
    rlang::abort("Must specify a valid per_engine.  Use 'stats' to call prcomp,
    or 'Gifi to call princals")
  }

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
