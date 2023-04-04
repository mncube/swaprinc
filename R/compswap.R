#' Compare swaprinc Models
#'
#' @param data A dataframe
#' @param formula A quoted model formula
#' @param engine The engine for fitting the model.  Options are "stats" or"lme4".
#' @param .pca_varlist A list of pca_vars (see swaprinc documentation)
#' @param .n_pca_list A list of n_pca_components (see swaprinc documentation)
#' @param .center_list A list of center values (see swaprinc documentation)
#' @param .scale._list A list of scale values (see swaprinc documentation)
#' @param .lpca_center_list A list of lpca_center values (see swaprinc documentation)
#' @param .lpca_scale_list A list of lpca_scale values (see swaprinc documentation)
#' @param .lpca_undo_list A list of lpca_undo values (see swaprinc documentation)
#' @param ... Pass additional arguments to the swaprinc function
#'
#' @return A list containing a list of fitted models and a comparison metrics
#' data frame.
#' @export
#'
#' @examples
#'# Load the iris dataset
#'data(iris)
#'
#'# Define the formula
#'formula <- "Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width"
#'
#'# Define the pca_varlist
#'pca_varlist <- list(c("Sepal.Width", "Petal.Length"),
#'                    c("Sepal.Width", "Petal.Width"))
#'
#'# Define the n_pca_list
#'n_pca_list <- list(2, 2)
#'
#'# Set scaling values
#'center_list <- list(TRUE, TRUE)
#'scale._list <- list(FALSE, FALSE)
#'lpca_center_list <- list("none", "none")
#'lpca_scale_list <- list("none", "none")
#'lpca_undo_list <- list(FALSE, FALSE)
#'
#'# Run compswap
#'compswap_results <- compswap(data = iris,
#'                             formula = formula,
#'                             engine = "stats",
#'                             .pca_varlist = pca_varlist,
#'                             .n_pca_list = n_pca_list,
#'                             .center_list = center_list,
#'                             .scale._list = scale._list,
#'                             .lpca_center_list = lpca_center_list,
#'                             .lpca_scale_list = lpca_scale_list,
#'                             .lpca_undo_list = lpca_undo_list)
compswap <- function(data, formula,
                     engine = "stats",
                     .pca_varlist,
                     .n_pca_list,
                     .center_list,
                     .scale._list,
                     .lpca_center_list,
                     .lpca_scale_list,
                     .lpca_undo_list,...) {

  if (length(.pca_varlist) != length(.n_pca_list)) {
    rlang::abort("Length of .pca_varlist and .n_pca_list must be the same.")
  }

  all_models <- list()
  all_comparisons <- data.frame()

  norun_raw <- FALSE

  for (i in seq_along(.pca_varlist)) {
    pca_vars <- .pca_varlist[[i]]
    n_pca_components <- .n_pca_list[[i]]
    center <- .center_list[[i]]
    scale. <- .scale._list[[i]]
    lpca_center <- .lpca_center_list[[i]]
    lpca_scale <- .lpca_scale_list[[i]]
    lpca_undo <- .lpca_undo_list[[i]]

    #swaprinc_result <- swaprinc(data, formula, engine, pca_vars, n_pca_components, ..., norun_raw = norun_raw)
    swaprinc_result <- swaprinc(data, formula, engine, pca_vars,
                                n_pca_components, norun_raw = norun_raw, center,
                                scale., lpca_center, lpca_scale, lpca_undo, ...)

    if (!norun_raw) {
      all_models$model_raw <- swaprinc_result$model_raw
      norun_raw <- TRUE
    }

    all_models[[paste0("model_pca_", i)]] <- swaprinc_result$model_pca

    comparison_df <- swaprinc_result$comparison
    comparison_df$model_set <- paste0("model_pca_", i)
    all_comparisons <- rbind(all_comparisons, comparison_df)

    # Update the model_set column for the Raw model row
    if (i == 1) {
      all_comparisons[all_comparisons$model == "Raw", "model_set"] <- "model_raw"
    }
  }

  return(list(all_models = all_models, all_comparisons = all_comparisons))
}


