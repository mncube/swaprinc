# Test case for compswap using lm engine
test_that("compswap works correctly with lm engine", {
  # Load the iris dataset
  data(iris)

  # Define the formula
  formula <- "Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width"

  # Define the pca_varlist
  pca_varlist <- list(c("Sepal.Width", "Petal.Length"),
                      c("Sepal.Width", "Petal.Width"))

  # Define the n_pca_list
  n_pca_list <- list(2, 2)

  # Run compswap
  compswap_results <- compswap(data = iris,
                               formula = formula,
                               engine = "stats",
                               .pca_varlist = pca_varlist,
                               .n_pca_list = n_pca_list)

  # Check if the output is a list with the expected length
  expect_type(compswap_results, "list")
  expect_equal(length(compswap_results), 2)

  # Check if the model_raw is included only once
  expect_s3_class(compswap_results$all_models$model_raw, "lm")
  expect_equal(sum(names(compswap_results$all_models) == "model_raw"), 1)

  # Check if the correct number of model_pca objects are included
  expect_equal(sum(grepl("model_pca", names(compswap_results$all_models))), 2)

  # Check if the correct number of comparisons dataframes are included
  expect_equal(sum(grepl("comparison", names(compswap_results))), 1)
})
