# Run tests for stats engine
test_that("swaprinc returns correct class for basic lm", {
  #Get iris data
  data(iris)

  #Run basic lm model using stats
  res <- swaprinc(iris,
                  "Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width",
                  pca_vars = c("Sepal.Width", "Petal.Length"),
                  n_pca_components = 2)

  #See if swaprinc returned a list
  expect_type(res, "list")


})

test_that("swaprinc works with logistic regression using stats engine", {
  suppressWarnings({
    # Get iris data
    data(iris)
    iris$Species_binary <- ifelse(iris$Species == "setosa", 1, 0)

    # Run logistic regression model using stats engine
    res_logistic <- swaprinc(data = iris,
                             formula = "Species_binary ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width",
                             engine = "stats",
                             pca_vars = c("Sepal.Length", "Petal.Width"),
                             n_pca_components = 2,
                             family = binomial(link = "logit"))

    # Check if swaprinc returned a list
    expect_type(res_logistic, "list")

    # Check if both models in the list are of class "glm"
    expect_s3_class(res_logistic$model_raw, "glm")
    expect_s3_class(res_logistic$model_pca, "glm")
  })
})

# Run tests for lme4 engine
test_that("A basic lme4 random intercept model works.", {

  # Set seed for reproducibility
  set.seed(42)

  # Define the number of subjects and observations per subject
  n_subjects <- 30
  n_obs_per_subject <- 10

  # Generate subject IDs
  subject_ids <- 1:n_subjects

  # Simulate the random effects for each subject
  random_intercepts <- rnorm(n_subjects, mean = 0, sd = 2)

  # Create an empty data frame to store the simulated data
  simulated_data <- data.frame()

  # Generate the data for each subject
  for (subject_id in subject_ids) {
    subject_data <- data.frame(
      Subject = rep(subject_id, n_obs_per_subject),
      Time = 1:n_obs_per_subject,
      X1 = rnorm(n_obs_per_subject, mean = 0, sd = 1),
      X2 = rnorm(n_obs_per_subject, mean = 0, sd = 1),
      Random_Intercept = rep(random_intercepts[subject_id], n_obs_per_subject)
    )

    # Simulate the response variable (e.g., Y) with fixed effects and the random intercept
    fixed_effects <- c(1.5, 2, -1)
    subject_data$Y <- fixed_effects[1] * subject_data$Time +
      fixed_effects[2] * subject_data$X1 +
      fixed_effects[3] * subject_data$X2 +
      subject_data$Random_Intercept +
      rnorm(n_obs_per_subject, mean = 0, sd = 1)

    # Combine the data for each subject into a single data frame
    simulated_data <- rbind(simulated_data, subject_data)
  }

  # Fit a random intercept model using lme4
  res_ri <- swaprinc(data = simulated_data,
                     formula = "Y ~ Time + X1 + X2 + (1 | Subject)",
                     pca_vars = c("Time", "X1", "X2"),
                     engine = "lme4",
                     n_pca_components = 2)

  #See if swaprinc returned a list
  expect_type(res_ri, "list")


})


test_that("Test that swaprinc throws error if a variable is included in both
          the random effects and pca_vars", {
  # Set seed for reproducibility
  set.seed(42)

  # Define the number of subjects and observations per subject
  n_subjects <- 30
  n_obs_per_subject <- 10

  # Generate subject IDs
  subject_ids <- 1:n_subjects

  # Simulate the random effects for each subject
  random_intercepts <- rnorm(n_subjects, mean = 0, sd = 2)
  random_slopes <- rnorm(n_subjects, mean = 0, sd = 1)

  # Create an empty data frame to store the simulated data
  simulated_data <- data.frame()

  # Generate the data for each subject
  for (subject_id in subject_ids) {
    subject_data <- data.frame(
      Subject = rep(subject_id, n_obs_per_subject),
      Time = 1:n_obs_per_subject,
      X1 = rnorm(n_obs_per_subject, mean = 0, sd = 1),
      X2 = rnorm(n_obs_per_subject, mean = 0, sd = 1),
      Random_Intercept = rep(random_intercepts[subject_id], n_obs_per_subject),
      Random_Slope = rep(random_slopes[subject_id], n_obs_per_subject)
    )

    # Simulate the response variable (e.g., Y) with fixed effects, random intercept, and random slope
    fixed_effects <- c(1.5, 2, -1)
    subject_data$Y <- fixed_effects[1] * subject_data$Time +
      fixed_effects[2] * subject_data$X1 +
      fixed_effects[3] * subject_data$X2 +
      subject_data$Random_Intercept +
      subject_data$Random_Slope * subject_data$Time +
      rnorm(n_obs_per_subject, mean = 0, sd = 1)

    # Combine the data for each subject into a single data frame
    simulated_data <- rbind(simulated_data, subject_data)
  }

  # Fit a random slope model using lme4
  expect_error(swaprinc(data = simulated_data,
                        formula = "Y ~ Time + X1 + X2 + (Time | Subject)",
                        pca_vars = c("Time", "X1", "X2"),
                        engine = "lme4",
                        n_pca_components = 2))

})

test_that("Test random slopes with variables not included in pca_vars", {
  # Set seed for reproducibility
  set.seed(42)

  # Define the number of subjects and observations per subject
  n_subjects <- 30
  n_obs_per_subject <- 10

  # Generate subject IDs
  subject_ids <- 1:n_subjects

  # Simulate the random effects for each subject
  random_intercepts <- rnorm(n_subjects, mean = 0, sd = 2)
  random_slopes <- rnorm(n_subjects, mean = 0, sd = 1)

  # Create an empty data frame to store the simulated data
  simulated_data <- data.frame()

  # Generate the data for each subject
  for (subject_id in subject_ids) {
    subject_data <- data.frame(
      Subject = rep(subject_id, n_obs_per_subject),
      Time = 1:n_obs_per_subject,
      X1 = rnorm(n_obs_per_subject, mean = 0, sd = 1),
      X2 = rnorm(n_obs_per_subject, mean = 0, sd = 1),
      Random_Intercept = rep(random_intercepts[subject_id], n_obs_per_subject),
      Random_Slope = rep(random_slopes[subject_id], n_obs_per_subject)
    )

    # Simulate the response variable (e.g., Y) with fixed effects, random intercept, and random slope
    fixed_effects <- c(1.5, 2, -1)
    subject_data$Y <- fixed_effects[1] * subject_data$Time +
      fixed_effects[2] * subject_data$X1 +
      fixed_effects[3] * subject_data$X2 +
      subject_data$Random_Intercept +
      subject_data$Random_Slope * subject_data$Time +
      rnorm(n_obs_per_subject, mean = 0, sd = 1)

    # Combine the data for each subject into a single data frame
    simulated_data <- rbind(simulated_data, subject_data)
  }

  # Fit a random slope model using lme4
  res_rs <- swaprinc(data = simulated_data,
                     formula = "Y ~ Time + X1 + X2 + (Time | Subject)",
                     pca_vars = c("X1", "X2"),
                     engine = "lme4",
                     n_pca_components = 2)

  # See if swaprinc returned a list
  expect_type(res_rs, "list")

  # Check if both models in the list are of class "merMod"
  expect_s4_class(res_rs$model_raw, "merMod")
  expect_s4_class(res_rs$model_pca, "merMod")
})

#Test norun_raw option
# Run tests for stats engine
test_that("swaprinc returns correct class for basic lm for norun_raw = TRUE", {
  #Get iris data
  data(iris)

  #Run basic lm model using stats
  res <- swaprinc(iris,
                  "Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width",
                  pca_vars = c("Sepal.Width", "Petal.Length"),
                  n_pca_components = 2,
                  norun_raw = TRUE)

  #See if swaprinc returned a list
  expect_type(res, "list")


})
