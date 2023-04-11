
<!-- README.md is generated from README.Rmd. Please edit that file -->

# swaprinc

<!-- badges: start -->
<!-- badges: end -->

The goal of swaprinc is to help compare a regression model with raw
variables to a regression model where some of the raw variables have
been swapped out for principal components.

## Installation

You can install the development version of swaprinc from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mncube/swaprinc")
```

## A Simple Example

In the example below, a regression model is used to estimate the
relation between x1 and y after controlling for x2 through x10.

Using the default engine of “stats” specifies that stats::lm will be
used to fit the statistical model and using the default prc_eng of
“stats” specifies that stats::prcomp will be used to extract principal
components.

The formula parameter specifies the “raw model” which will be passed to
stats::lm. The pca_vars and n_pca_components parameters specify that raw
variables x2 through x10 will be used to extract 3 principal components,
and then the following “pca model” will be passed to stats::lm: y \~
x1 + PC1 + PC2 + PC3

Setting the parameters lpca_center and lpca_scale to ‘pca’ specifies
that before passing the data in pca_vars to stats::prcomp, centering and
scaling will be executed as outlined in the LearnPCA Step-by-Step PCA
vignette:
(<https://cran.r-project.org/web/packages/LearnPCA/vignettes/Vig_03_Step_By_Step_PCA.pdf>).
Setting the miss_handler parameter to ‘omit’ specifies that the rows of
the data frame passed to swaprinc will be subset with
stats::complete.cases.

``` r
library(swaprinc)

  # Create a small simulated dataset
  set.seed(40)
  n <- 50
  x1 <- rnorm(n)
  x2 <- rnorm(n, 5, 15)
  x3 <- rnorm(n, -5.5, 20)
  x4 <- rnorm(n, 3, 3) + x3*1.5
  x5 <- rnorm(n, -2, 4) + x3*.25
  x6 <- rnorm(n, -5, 5) + x4
  x7 <- rnorm(n, -2, 6)
  x8 <- rnorm(n, 2, 7)
  x9 <- rnorm(n, -2, 3) +x2*.4
  x10 <- rnorm(n, 5, 4)
  y <- 1 + 2 * x1 + 3 * x2 + 2.5*x4 - 3.5*x5 + 2*x6 + 1.5*x7 + x8 + 2*x9 + x10 + rnorm(n)
  data <- data.frame(y, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)

  # Run swaprinc with
  swaprinc_result <- swaprinc(data,
                              formula = "y ~ x1 + x2 + x3 + x4 + x5 + x5 + x6 + x7 + x8 + x9 + x10",
                              pca_vars = c("x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10"),
                              n_pca_components = 3,
                              lpca_center = "pca", 
                              lpca_scale = "pca",
                              miss_handler = "omit")
  
  # Summarize raw model
  summary(swaprinc_result$model_raw)
#> 
#> Call:
#> stats::lm(formula = formula, data = data)
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -2.22717 -0.55113  0.00573  0.45296  3.03938 
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  1.75250    0.38418   4.562 4.94e-05 ***
#> x1           2.06070    0.17472  11.794 1.96e-14 ***
#> x2           2.93211    0.02976  98.525  < 2e-16 ***
#> x3           0.06272    0.07596   0.826    0.414    
#> x4           2.40907    0.06368  37.830  < 2e-16 ***
#> x5          -3.48453    0.04467 -78.007  < 2e-16 ***
#> x6           2.04578    0.03497  58.501  < 2e-16 ***
#> x7           1.49816    0.02621  57.158  < 2e-16 ***
#> x8           0.97008    0.02530  38.345  < 2e-16 ***
#> x9           2.12426    0.06502  32.670  < 2e-16 ***
#> x10          0.99118    0.03787  26.174  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 1.064 on 39 degrees of freedom
#> Multiple R-squared:  0.9999, Adjusted R-squared:  0.9999 
#> F-statistic: 7.237e+04 on 10 and 39 DF,  p-value: < 2.2e-16
  
  # Summarize pca model
  summary(swaprinc_result$model_pca)
#> 
#> Call:
#> stats::lm(formula = formula, data = data)
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -0.47301 -0.16369 -0.03991  0.15923  0.48937 
#> 
#> Coefficients:
#>               Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  4.130e-17  3.701e-02   0.000    1.000    
#> x1           6.056e-02  4.050e-02   1.495    0.142    
#> PC1          4.227e-01  1.955e-02  21.621  < 2e-16 ***
#> PC2          3.692e-01  2.781e-02  13.275  < 2e-16 ***
#> PC3         -1.513e-01  3.422e-02  -4.423 6.11e-05 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 0.2617 on 45 degrees of freedom
#> Multiple R-squared:  0.9371, Adjusted R-squared:  0.9315 
#> F-statistic: 167.6 on 4 and 45 DF,  p-value: < 2.2e-16
  
  # Get model comparisons
  print(swaprinc_result$comparison)
#>   model r_squared adj_r_squared      AIC       BIC
#> 1   Raw 0.9999461     0.9999323 159.7190 182.66323
#> 2   PCA 0.9370900     0.9314980  14.5812  26.05334
```

## The Motivating Example

A chronic difficulty in applied statistics and data science is logistic
regression with a set of categorical independent variables. In this
motivating example, swaprinc will be used to compare a ‘raw’ logistic
regression model with five categorical independent variables to a ‘pca’
logistic regression model where Gifi::princals will be used to swap out
six of the independent variables for the first three principal
components.

``` r
 # Create a small simulated dataset
  set.seed(42)
  n <- 50
  x1 <- rnorm(n, 0.5, 4)
  x2 <- rnorm(n, 3, 15)
  x3 <- rnorm(n, -2.5, 5)
  x4 <- -2.5*x2 + 3*x3 + rnorm(n, 0, 4)
  x5 <- x2*x3 + rnorm(n, -5, 5)*rnorm(n, 5, 10)
  x6 <- rnorm(n, -2, 4)*rnorm(n, 3, 5)
  x7 <- x4 + x6 + rnorm(n, 0, 3)
  y <- 1 + 2*x1 + 3*x2 + -2*x3 + .5*x4 + x5 + 1.5*x6 + x7 + rnorm(n)
  data <- data.frame(y, x1, x2, x3, x4, x5, x6, x7)

  # Categorize the variables
  yq <- stats::quantile(data$y,c(0,1/2, 1))
  x1q <- stats::quantile(data$x1,c(0,1/2, 1))
  x2q <- stats::quantile(data$x2,c(0,1/4,3/4,1))
  x3q <- stats::quantile(data$x3,c(0,2/5,3/5,1))
  x4q <- stats::quantile(data$x4,c(0,1/5,4/5,1))
  x5q <- stats::quantile(data$x5,c(0,2/5,3/5,1))
  x6q <- stats::quantile(data$x6,c(0,2/5,4/5,1))
  x7q <- stats::quantile(data$x7,c(0,2/5,3/5,1))


  data <- data %>% dplyr::mutate(
    y = cut(y, breaks=yq, labels=c("0", "1"),include.lowest = TRUE),
    x1 = cut(x1, breaks=x1q, labels=c("control", "treatment"),include.lowest = TRUE),
    x2 = cut(x2, breaks=x2q, labels=c("small","medium","large"),include.lowest = TRUE),
    x3 = cut(x3, breaks=x3q, labels=c("short","average","tall"),include.lowest = TRUE),
    x4 = cut(x4, breaks=x4q, labels=c("lowbit","most","highbit"),include.lowest = TRUE),
    x5 = cut(x5, breaks=x5q, labels=c("under","healthy","over"),include.lowest = TRUE),
    x6 = cut(x6, breaks=x6q, labels=c("small","medium","large"),include.lowest = TRUE),
    x7 = cut(x7, breaks=x7q, labels=c("small","medium","large"),include.lowest = TRUE)) %>%
    dplyr::mutate(y = as.numeric(ifelse(y == "0", 0, 1)))

  # Run swaprinc with prc_eng set to Gifi
  swaprinc_result <- swaprinc(data,
                              formula = "y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7",
                              pca_vars = c("x2", "x3", "x4", "x5", "x6", "x7"),
                              n_pca_components = 3,
                              prc_eng = "Gifi",
                              model_options = list(family = binomial(link = "logit")))
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
  
  # Summarize raw model
  summary(swaprinc_result$model_raw)
#> 
#> Call:
#> (function (formula, family = gaussian, data, weights, subset, 
#>     na.action, start = NULL, etastart, mustart, offset, control = list(...), 
#>     model = TRUE, method = "glm.fit", x = FALSE, y = TRUE, singular.ok = TRUE, 
#>     contrasts = NULL, ...) 
#> {
#>     cal <- match.call()
#>     if (is.character(family)) 
#>         family <- get(family, mode = "function", envir = parent.frame())
#>     if (is.function(family)) 
#>         family <- family()
#>     if (is.null(family$family)) {
#>         print(family)
#>         stop("'family' not recognized")
#>     }
#>     if (missing(data)) 
#>         data <- environment(formula)
#>     mf <- match.call(expand.dots = FALSE)
#>     m <- match(c("formula", "data", "subset", "weights", "na.action", 
#>         "etastart", "mustart", "offset"), names(mf), 0L)
#>     mf <- mf[c(1L, m)]
#>     mf$drop.unused.levels <- TRUE
#>     mf[[1L]] <- quote(stats::model.frame)
#>     mf <- eval(mf, parent.frame())
#>     if (identical(method, "model.frame")) 
#>         return(mf)
#>     if (!is.character(method) && !is.function(method)) 
#>         stop("invalid 'method' argument")
#>     if (identical(method, "glm.fit")) 
#>         control <- do.call("glm.control", control)
#>     mt <- attr(mf, "terms")
#>     Y <- model.response(mf, "any")
#>     if (length(dim(Y)) == 1L) {
#>         nm <- rownames(Y)
#>         dim(Y) <- NULL
#>         if (!is.null(nm)) 
#>             names(Y) <- nm
#>     }
#>     X <- if (!is.empty.model(mt)) 
#>         model.matrix(mt, mf, contrasts)
#>     else matrix(, NROW(Y), 0L)
#>     weights <- as.vector(model.weights(mf))
#>     if (!is.null(weights) && !is.numeric(weights)) 
#>         stop("'weights' must be a numeric vector")
#>     if (!is.null(weights) && any(weights < 0)) 
#>         stop("negative weights not allowed")
#>     offset <- as.vector(model.offset(mf))
#>     if (!is.null(offset)) {
#>         if (length(offset) != NROW(Y)) 
#>             stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
#>                 length(offset), NROW(Y)), domain = NA)
#>     }
#>     mustart <- model.extract(mf, "mustart")
#>     etastart <- model.extract(mf, "etastart")
#>     fit <- eval(call(if (is.function(method)) "method" else method, 
#>         x = X, y = Y, weights = weights, start = start, etastart = etastart, 
#>         mustart = mustart, offset = offset, family = family, 
#>         control = control, intercept = attr(mt, "intercept") > 
#>             0L, singular.ok = singular.ok))
#>     if (length(offset) && attr(mt, "intercept") > 0L) {
#>         fit2 <- eval(call(if (is.function(method)) "method" else method, 
#>             x = X[, "(Intercept)", drop = FALSE], y = Y, mustart = fit$fitted.values, 
#>             weights = weights, offset = offset, family = family, 
#>             control = control, intercept = TRUE))
#>         if (!fit2$converged) 
#>             warning("fitting to calculate the null deviance did not converge -- increase 'maxit'?")
#>         fit$null.deviance <- fit2$deviance
#>     }
#>     if (model) 
#>         fit$model <- mf
#>     fit$na.action <- attr(mf, "na.action")
#>     if (x) 
#>         fit$x <- X
#>     if (!y) 
#>         fit$y <- NULL
#>     structure(c(fit, list(call = cal, formula = formula, terms = mt, 
#>         data = data, offset = offset, control = control, method = method, 
#>         contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, 
#>             mf))), class = c(fit$class, c("glm", "lm")))
#> })(formula = "y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7", family = structure(list(
#>     family = "binomial", link = "logit", linkfun = function (mu) 
#>     .Call(C_logit_link, mu), linkinv = function (eta) 
#>     .Call(C_logit_linkinv, eta), variance = function (mu) 
#>     mu * (1 - mu), dev.resids = function (y, mu, wt) 
#>     .Call(C_binomial_dev_resids, y, mu, wt), aic = function (y, 
#>         n, mu, wt, dev) 
#>     {
#>         m <- if (any(n > 1)) 
#>             n
#>         else wt
#>         -2 * sum(ifelse(m > 0, (wt/m), 0) * dbinom(round(m * 
#>             y), round(m), mu, log = TRUE))
#>     }, mu.eta = function (eta) 
#>     .Call(C_logit_mu_eta, eta), initialize = {
#>         if (NCOL(y) == 1) {
#>             if (is.factor(y)) 
#>                 y <- y != levels(y)[1L]
#>             n <- rep.int(1, nobs)
#>             y[weights == 0] <- 0
#>             if (any(y < 0 | y > 1)) 
#>                 stop("y values must be 0 <= y <= 1")
#>             mustart <- (weights * y + 0.5)/(weights + 1)
#>             m <- weights * y
#>             if ("binomial" == "binomial" && any(abs(m - round(m)) > 
#>                 0.001)) 
#>                 warning(gettextf("non-integer #successes in a %s glm!", 
#>                   "binomial"), domain = NA)
#>         }
#>         else if (NCOL(y) == 2) {
#>             if ("binomial" == "binomial" && any(abs(y - round(y)) > 
#>                 0.001)) 
#>                 warning(gettextf("non-integer counts in a %s glm!", 
#>                   "binomial"), domain = NA)
#>             n <- (y1 <- y[, 1L]) + y[, 2L]
#>             y <- y1/n
#>             if (any(n0 <- n == 0)) 
#>                 y[n0] <- 0
#>             weights <- weights * n
#>             mustart <- (n * y + 0.5)/(n + 1)
#>         }
#>         else stop(gettextf("for the '%s' family, y must be a vector of 0 and 1's\nor a 2 column matrix where col 1 is no. successes and col 2 is no. failures", 
#>             "binomial"), domain = NA)
#>     }, validmu = function (mu) 
#>     all(is.finite(mu)) && all(mu > 0 & mu < 1), valideta = function (eta) 
#>     TRUE, simulate = function (object, nsim) 
#>     {
#>         ftd <- fitted(object)
#>         n <- length(ftd)
#>         ntot <- n * nsim
#>         wts <- object$prior.weights
#>         if (any(wts%%1 != 0)) 
#>             stop("cannot simulate from non-integer prior.weights")
#>         if (!is.null(m <- object$model)) {
#>             y <- model.response(m)
#>             if (is.factor(y)) {
#>                 yy <- factor(1 + rbinom(ntot, size = 1, prob = ftd), 
#>                   labels = levels(y))
#>                 split(yy, rep(seq_len(nsim), each = n))
#>             }
#>             else if (is.matrix(y) && ncol(y) == 2) {
#>                 yy <- vector("list", nsim)
#>                 for (i in seq_len(nsim)) {
#>                   Y <- rbinom(n, size = wts, prob = ftd)
#>                   YY <- cbind(Y, wts - Y)
#>                   colnames(YY) <- colnames(y)
#>                   yy[[i]] <- YY
#>                 }
#>                 yy
#>             }
#>             else rbinom(ntot, size = wts, prob = ftd)/wts
#>         }
#>         else rbinom(ntot, size = wts, prob = ftd)/wts
#>     }), class = "family"), data = structure(list(y = c(1, 1, 
#> 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 
#> 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 
#> 1, 1, 1, 1, 1, 1), x1 = structure(c(2L, 1L, 2L, 2L, 2L, 1L, 2L, 
#> 2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 
#> 2L, 2L, 1L, 1L, 1L, 2L, 1L, 2L, 2L, 2L, 1L, 2L, 1L, 1L, 1L, 1L, 
#> 2L, 2L, 1L, 2L, 1L, 1L, 2L, 1L, 2L, 1L, 2L), levels = c("control", 
#> "treatment"), class = "factor"), x2 = structure(c(2L, 1L, 3L, 
#> 2L, 2L, 2L, 3L, 2L, 1L, 2L, 2L, 2L, 2L, 3L, 1L, 3L, 2L, 3L, 3L, 
#> 3L, 1L, 2L, 2L, 1L, 1L, 2L, 3L, 2L, 1L, 1L, 3L, 2L, 2L, 2L, 1L, 
#> 2L, 2L, 2L, 3L, 3L, 3L, 2L, 2L, 3L, 1L, 1L, 1L, 1L, 2L, 2L), levels = c("small", 
#> "medium", "large"), class = "factor"), x3 = structure(c(3L, 3L, 
#> 1L, 3L, 1L, 3L, 2L, 2L, 3L, 3L, 3L, 3L, 1L, 1L, 1L, 2L, 1L, 3L, 
#> 1L, 3L, 1L, 1L, 3L, 1L, 3L, 2L, 1L, 1L, 1L, 3L, 3L, 1L, 3L, 3L, 
#> 3L, 1L, 2L, 3L, 1L, 2L, 2L, 1L, 1L, 2L, 2L, 3L, 1L, 2L, 3L, 1L
#> ), levels = c("short", "average", "tall"), class = "factor"), 
#>     x4 = structure(c(2L, 3L, 1L, 2L, 2L, 2L, 2L, 2L, 3L, 2L, 
#>     2L, 2L, 2L, 1L, 2L, 1L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 3L, 2L, 
#>     2L, 2L, 1L, 2L, 3L, 1L, 2L, 2L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 
#>     1L, 2L, 2L, 1L, 3L, 3L, 3L, 3L, 2L, 2L), levels = c("lowbit", 
#>     "most", "highbit"), class = "factor"), x5 = structure(c(3L, 
#>     2L, 1L, 3L, 1L, 3L, 2L, 1L, 3L, 2L, 1L, 3L, 1L, 2L, 3L, 1L, 
#>     1L, 3L, 2L, 1L, 3L, 2L, 1L, 3L, 3L, 1L, 1L, 1L, 1L, 3L, 1L, 
#>     2L, 3L, 1L, 1L, 1L, 3L, 3L, 2L, 1L, 1L, 3L, 1L, 2L, 3L, 2L, 
#>     3L, 3L, 3L, 3L), levels = c("under", "healthy", "over"), class = "factor"), 
#>     x6 = structure(c(1L, 2L, 1L, 2L, 2L, 1L, 2L, 3L, 1L, 3L, 
#>     2L, 3L, 2L, 2L, 2L, 1L, 1L, 2L, 2L, 3L, 1L, 2L, 1L, 3L, 3L, 
#>     2L, 2L, 2L, 1L, 1L, 3L, 3L, 1L, 1L, 2L, 2L, 1L, 1L, 1L, 1L, 
#>     1L, 2L, 1L, 2L, 2L, 1L, 3L, 1L, 2L, 3L), levels = c("small", 
#>     "medium", "large"), class = "factor"), x7 = structure(c(2L, 
#>     3L, 1L, 3L, 2L, 1L, 1L, 3L, 3L, 3L, 3L, 3L, 1L, 1L, 3L, 1L, 
#>     1L, 3L, 1L, 2L, 2L, 2L, 1L, 3L, 3L, 1L, 1L, 1L, 2L, 3L, 1L, 
#>     2L, 1L, 3L, 3L, 1L, 2L, 2L, 1L, 1L, 1L, 3L, 1L, 1L, 3L, 3L, 
#>     3L, 3L, 3L, 2L), levels = c("small", "medium", "large"), class = "factor")), class = "data.frame", row.names = c(NA, 
#> -50L)))
#> 
#> Deviance Residuals: 
#>     Min       1Q   Median       3Q      Max  
#> -1.7567   0.0000   0.0000   0.0000   0.9994  
#> 
#> Coefficients:
#>               Estimate Std. Error z value Pr(>|z|)
#> (Intercept) -2.897e+01  3.924e+04  -0.001    0.999
#> x1treatment  8.684e-01  2.203e+00   0.394    0.693
#> x2medium    -2.759e+01  3.604e+04  -0.001    0.999
#> x2large     -3.264e+01  3.809e+04  -0.001    0.999
#> x3average    2.511e+01  1.841e+04   0.001    0.999
#> x3tall       4.673e+00  1.769e+04   0.000    1.000
#> x4most       9.150e-12  2.507e+00   0.000    1.000
#> x4highbit   -2.219e+01  3.954e+04  -0.001    1.000
#> x5healthy    4.025e+01  1.611e+04   0.002    0.998
#> x5over       5.233e+01  2.104e+04   0.002    0.998
#> x6medium    -4.181e+00  1.968e+04   0.000    1.000
#> x6large      5.650e+01  2.313e+04   0.002    0.998
#> x7medium    -8.684e-01  2.203e+00  -0.394    0.693
#> x7large      2.995e+01  2.079e+04   0.001    0.999
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 69.3147  on 49  degrees of freedom
#> Residual deviance:  9.0793  on 36  degrees of freedom
#> AIC: 37.079
#> 
#> Number of Fisher Scoring iterations: 22
  
  # Summarize pca model
  summary(swaprinc_result$model_pca)
#> 
#> Call:
#> (function (formula, family = gaussian, data, weights, subset, 
#>     na.action, start = NULL, etastart, mustart, offset, control = list(...), 
#>     model = TRUE, method = "glm.fit", x = FALSE, y = TRUE, singular.ok = TRUE, 
#>     contrasts = NULL, ...) 
#> {
#>     cal <- match.call()
#>     if (is.character(family)) 
#>         family <- get(family, mode = "function", envir = parent.frame())
#>     if (is.function(family)) 
#>         family <- family()
#>     if (is.null(family$family)) {
#>         print(family)
#>         stop("'family' not recognized")
#>     }
#>     if (missing(data)) 
#>         data <- environment(formula)
#>     mf <- match.call(expand.dots = FALSE)
#>     m <- match(c("formula", "data", "subset", "weights", "na.action", 
#>         "etastart", "mustart", "offset"), names(mf), 0L)
#>     mf <- mf[c(1L, m)]
#>     mf$drop.unused.levels <- TRUE
#>     mf[[1L]] <- quote(stats::model.frame)
#>     mf <- eval(mf, parent.frame())
#>     if (identical(method, "model.frame")) 
#>         return(mf)
#>     if (!is.character(method) && !is.function(method)) 
#>         stop("invalid 'method' argument")
#>     if (identical(method, "glm.fit")) 
#>         control <- do.call("glm.control", control)
#>     mt <- attr(mf, "terms")
#>     Y <- model.response(mf, "any")
#>     if (length(dim(Y)) == 1L) {
#>         nm <- rownames(Y)
#>         dim(Y) <- NULL
#>         if (!is.null(nm)) 
#>             names(Y) <- nm
#>     }
#>     X <- if (!is.empty.model(mt)) 
#>         model.matrix(mt, mf, contrasts)
#>     else matrix(, NROW(Y), 0L)
#>     weights <- as.vector(model.weights(mf))
#>     if (!is.null(weights) && !is.numeric(weights)) 
#>         stop("'weights' must be a numeric vector")
#>     if (!is.null(weights) && any(weights < 0)) 
#>         stop("negative weights not allowed")
#>     offset <- as.vector(model.offset(mf))
#>     if (!is.null(offset)) {
#>         if (length(offset) != NROW(Y)) 
#>             stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
#>                 length(offset), NROW(Y)), domain = NA)
#>     }
#>     mustart <- model.extract(mf, "mustart")
#>     etastart <- model.extract(mf, "etastart")
#>     fit <- eval(call(if (is.function(method)) "method" else method, 
#>         x = X, y = Y, weights = weights, start = start, etastart = etastart, 
#>         mustart = mustart, offset = offset, family = family, 
#>         control = control, intercept = attr(mt, "intercept") > 
#>             0L, singular.ok = singular.ok))
#>     if (length(offset) && attr(mt, "intercept") > 0L) {
#>         fit2 <- eval(call(if (is.function(method)) "method" else method, 
#>             x = X[, "(Intercept)", drop = FALSE], y = Y, mustart = fit$fitted.values, 
#>             weights = weights, offset = offset, family = family, 
#>             control = control, intercept = TRUE))
#>         if (!fit2$converged) 
#>             warning("fitting to calculate the null deviance did not converge -- increase 'maxit'?")
#>         fit$null.deviance <- fit2$deviance
#>     }
#>     if (model) 
#>         fit$model <- mf
#>     fit$na.action <- attr(mf, "na.action")
#>     if (x) 
#>         fit$x <- X
#>     if (!y) 
#>         fit$y <- NULL
#>     structure(c(fit, list(call = cal, formula = formula, terms = mt, 
#>         data = data, offset = offset, control = control, method = method, 
#>         contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, 
#>             mf))), class = c(fit$class, c("glm", "lm")))
#> })(formula = y ~ x1 + PC1 + PC2 + PC3, family = structure(list(
#>     family = "binomial", link = "logit", linkfun = function (mu) 
#>     .Call(C_logit_link, mu), linkinv = function (eta) 
#>     .Call(C_logit_linkinv, eta), variance = function (mu) 
#>     mu * (1 - mu), dev.resids = function (y, mu, wt) 
#>     .Call(C_binomial_dev_resids, y, mu, wt), aic = function (y, 
#>         n, mu, wt, dev) 
#>     {
#>         m <- if (any(n > 1)) 
#>             n
#>         else wt
#>         -2 * sum(ifelse(m > 0, (wt/m), 0) * dbinom(round(m * 
#>             y), round(m), mu, log = TRUE))
#>     }, mu.eta = function (eta) 
#>     .Call(C_logit_mu_eta, eta), initialize = {
#>         if (NCOL(y) == 1) {
#>             if (is.factor(y)) 
#>                 y <- y != levels(y)[1L]
#>             n <- rep.int(1, nobs)
#>             y[weights == 0] <- 0
#>             if (any(y < 0 | y > 1)) 
#>                 stop("y values must be 0 <= y <= 1")
#>             mustart <- (weights * y + 0.5)/(weights + 1)
#>             m <- weights * y
#>             if ("binomial" == "binomial" && any(abs(m - round(m)) > 
#>                 0.001)) 
#>                 warning(gettextf("non-integer #successes in a %s glm!", 
#>                   "binomial"), domain = NA)
#>         }
#>         else if (NCOL(y) == 2) {
#>             if ("binomial" == "binomial" && any(abs(y - round(y)) > 
#>                 0.001)) 
#>                 warning(gettextf("non-integer counts in a %s glm!", 
#>                   "binomial"), domain = NA)
#>             n <- (y1 <- y[, 1L]) + y[, 2L]
#>             y <- y1/n
#>             if (any(n0 <- n == 0)) 
#>                 y[n0] <- 0
#>             weights <- weights * n
#>             mustart <- (n * y + 0.5)/(n + 1)
#>         }
#>         else stop(gettextf("for the '%s' family, y must be a vector of 0 and 1's\nor a 2 column matrix where col 1 is no. successes and col 2 is no. failures", 
#>             "binomial"), domain = NA)
#>     }, validmu = function (mu) 
#>     all(is.finite(mu)) && all(mu > 0 & mu < 1), valideta = function (eta) 
#>     TRUE, simulate = function (object, nsim) 
#>     {
#>         ftd <- fitted(object)
#>         n <- length(ftd)
#>         ntot <- n * nsim
#>         wts <- object$prior.weights
#>         if (any(wts%%1 != 0)) 
#>             stop("cannot simulate from non-integer prior.weights")
#>         if (!is.null(m <- object$model)) {
#>             y <- model.response(m)
#>             if (is.factor(y)) {
#>                 yy <- factor(1 + rbinom(ntot, size = 1, prob = ftd), 
#>                   labels = levels(y))
#>                 split(yy, rep(seq_len(nsim), each = n))
#>             }
#>             else if (is.matrix(y) && ncol(y) == 2) {
#>                 yy <- vector("list", nsim)
#>                 for (i in seq_len(nsim)) {
#>                   Y <- rbinom(n, size = wts, prob = ftd)
#>                   YY <- cbind(Y, wts - Y)
#>                   colnames(YY) <- colnames(y)
#>                   yy[[i]] <- YY
#>                 }
#>                 yy
#>             }
#>             else rbinom(ntot, size = wts, prob = ftd)/wts
#>         }
#>         else rbinom(ntot, size = wts, prob = ftd)/wts
#>     }), class = "family"), data = structure(list(y = c(1, 1, 
#> 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 
#> 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 
#> 1, 1, 1, 1, 1, 1), x1 = structure(c(2L, 1L, 2L, 2L, 2L, 1L, 2L, 
#> 2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 
#> 2L, 2L, 1L, 1L, 1L, 2L, 1L, 2L, 2L, 2L, 1L, 2L, 1L, 1L, 1L, 1L, 
#> 2L, 2L, 1L, 2L, 1L, 1L, 2L, 1L, 2L, 1L, 2L), levels = c("control", 
#> "treatment"), class = "factor"), PC1 = c(-0.388804811590457, 
#> -1.20744930399809, 1.58158319556763, -0.761302056137922, 0.38510884321309, 
#> -0.0568226410822124, 0.861603010772897, -0.184781684797962, -1.67306292427343, 
#> -0.278198610312821, -0.267664142410702, -0.76798457448232, 0.717091013721334, 
#> 1.55355890211575, -0.900102677374199, 1.36670008655428, 0.741263357615495, 
#> -0.398055000515289, 1.55355890211575, 0.437225295520837, -0.527605432826734, 
#> 0.381256893655369, 0.436815272645008, -1.3994697015415, -1.21123328068909, 
#> 0.502207904707985, 1.08033806934397, 1.19416379605084, -0.0339675190995146, 
#> -1.67306292427343, 1.24628024835859, 0.374574375310971, -0.0568226410822124, 
#> -0.243491798516541, -1.20359735444037, 1.19416379605084, -0.299239835633319, 
#> -0.881489317413358, 1.57773124600991, 0.889627304224779, 1.36670008655428, 
#> -0.456853971167434, 0.741263357615495, 1.3386757931024, -1.60767029221045, 
#> -1.18327696010393, -1.3994697015415, -1.58349794831629, -0.761302056137922, 
#> -0.115211588858529), PC2 = c(1.50301434407398, -0.404723238842585, 
#> 0.627435292253674, -0.200142788510268, -1.07389539194851, 1.73617204076097, 
#> -0.0922686881765645, -1.05949845825969, 1.2325067905731, -0.782063963990278, 
#> -0.380256669733227, -0.603355545285134, -0.840737695261517, -0.829680855939054, 
#> -1.24191578211839, 1.28968050459106, 0.617783915449027, -0.115324456280514, 
#> -0.829680855939054, -0.454015572404638, 0.461241350465851, -1.07248992943069, 
#> 1.55605815953801, -1.56750182100129, -0.70685429480428, -0.178492482924134, 
#> -0.755919363031762, -0.915904650686624, 0.281127469242892, 1.2325067905731, 
#> -0.296024831142756, -1.47570268620556, 1.73617204076097, 1.07826494097732, 
#> -0.4061287013604, -0.915904650686624, 1.22698531232238, 1.58064106196595, 
#> 0.628840754771489, 1.36484746001616, 1.28968050459106, -1.13841703259925, 
#> 0.617783915449027, -0.167435643601671, -0.502043851889039, 1.05379837186796, 
#> -1.56750182100129, 0.956477758821505, -0.200142788510268, -1.29699426750041
#> ), PC3 = c(0.266733485727984, 0.572345027148548, -0.87141409130504, 
#> 1.24013165551087, -0.625049492675555, 0.0448154964135076, 0.951660089843622, 
#> 1.08592992602929, -0.338705680172926, 1.61535557786255, 1.41196757020723, 
#> 1.44486053115305, -0.846967481990031, -0.132199298679922, -1.20657314051107, 
#> 0.401993942586714, -1.58752314260201, 1.76669206143963, -0.132199298679922, 
#> 1.91041434260725, -2.17997131029395, -0.626390360662411, 0.216651411109864, 
#> -1.19759569763064, 0.802330288842979, 0.426440551901723, -0.320407076061276, 
#> -0.65741883662182, -2.0081353955976, -0.338705680172926, 1.87804499866098, 
#> -0.421661485020234, 0.0448154964135076, 0.671411909595253, 0.573685895135404, 
#> -0.65741883662182, -0.264033034092131, 0.0709820529662292, -0.872754959291896, 
#> 0.212445297218504, 0.401993942586714, -0.564042898200998, -1.58752314260201, 
#> 1.14120873521183, -0.128916539381067, -0.168210633463426, -1.19759569763064, 
#> -0.869472199993041, 1.24013165551087, -0.592156531729734)), class = "data.frame", row.names = c("1", 
#> "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", 
#> "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", 
#> "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", 
#> "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", 
#> "47", "48", "49", "50")))
#> 
#> Deviance Residuals: 
#>      Min        1Q    Median        3Q       Max  
#> -2.45029  -0.48290   0.03456   0.53738   1.87631  
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  -0.1621     0.5669  -0.286 0.774968    
#> x1treatment   0.3030     0.8103   0.374 0.708422    
#> PC1          -2.0040     0.5391  -3.718 0.000201 ***
#> PC2          -0.5010     0.3869  -1.295 0.195358    
#> PC3           0.3391     0.3806   0.891 0.372944    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 69.315  on 49  degrees of freedom
#> Residual deviance: 40.765  on 45  degrees of freedom
#> AIC: 50.765
#> 
#> Number of Fisher Scoring iterations: 5
  
  # Get model comparisons
  print(swaprinc_result$comparison)
#>   model pseudo_r_squared      AIC      BIC
#> 1   Raw        0.8690138 37.07927 63.84759
#> 2   PCA        0.4118838 50.76511 60.32522
```

## Compare Multiple Models

Using the same data set for the logistic regression model above, it
would be useful to compare results for different ways of swapping
variables. In the example presented below, the compswap helper function
is used to compare results with 2, 3, 4, and 5 principal components
swapped in for the six raw independent variables.

``` r
  # Run swaprinc with prc_eng set to Gifi
  compswap_results <- compswap(data,
                              formula = "y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7",
                              .pca_varlist = list(c("x2", "x3", "x4", "x5", "x6", "x7")),
                              .n_pca_list = list(2, 3, 4, 5),
                              .prc_eng_list = list("Gifi"),
                              .model_options_list = list(list(family = binomial(link = "logit"))))
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

  # Show available models
  summary(compswap_results$all_models)
#>             Length Class Mode
#> model_raw   30     glm   list
#> model_pca_1 30     glm   list
#> model_pca_2 30     glm   list
#> model_pca_3 30     glm   list
#> model_pca_4 30     glm   list

  # Summarize raw model
  summary(compswap_results$all_models$model_raw)
#> 
#> Call:
#> (function (formula, family = gaussian, data, weights, subset, 
#>     na.action, start = NULL, etastart, mustart, offset, control = list(...), 
#>     model = TRUE, method = "glm.fit", x = FALSE, y = TRUE, singular.ok = TRUE, 
#>     contrasts = NULL, ...) 
#> {
#>     cal <- match.call()
#>     if (is.character(family)) 
#>         family <- get(family, mode = "function", envir = parent.frame())
#>     if (is.function(family)) 
#>         family <- family()
#>     if (is.null(family$family)) {
#>         print(family)
#>         stop("'family' not recognized")
#>     }
#>     if (missing(data)) 
#>         data <- environment(formula)
#>     mf <- match.call(expand.dots = FALSE)
#>     m <- match(c("formula", "data", "subset", "weights", "na.action", 
#>         "etastart", "mustart", "offset"), names(mf), 0L)
#>     mf <- mf[c(1L, m)]
#>     mf$drop.unused.levels <- TRUE
#>     mf[[1L]] <- quote(stats::model.frame)
#>     mf <- eval(mf, parent.frame())
#>     if (identical(method, "model.frame")) 
#>         return(mf)
#>     if (!is.character(method) && !is.function(method)) 
#>         stop("invalid 'method' argument")
#>     if (identical(method, "glm.fit")) 
#>         control <- do.call("glm.control", control)
#>     mt <- attr(mf, "terms")
#>     Y <- model.response(mf, "any")
#>     if (length(dim(Y)) == 1L) {
#>         nm <- rownames(Y)
#>         dim(Y) <- NULL
#>         if (!is.null(nm)) 
#>             names(Y) <- nm
#>     }
#>     X <- if (!is.empty.model(mt)) 
#>         model.matrix(mt, mf, contrasts)
#>     else matrix(, NROW(Y), 0L)
#>     weights <- as.vector(model.weights(mf))
#>     if (!is.null(weights) && !is.numeric(weights)) 
#>         stop("'weights' must be a numeric vector")
#>     if (!is.null(weights) && any(weights < 0)) 
#>         stop("negative weights not allowed")
#>     offset <- as.vector(model.offset(mf))
#>     if (!is.null(offset)) {
#>         if (length(offset) != NROW(Y)) 
#>             stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
#>                 length(offset), NROW(Y)), domain = NA)
#>     }
#>     mustart <- model.extract(mf, "mustart")
#>     etastart <- model.extract(mf, "etastart")
#>     fit <- eval(call(if (is.function(method)) "method" else method, 
#>         x = X, y = Y, weights = weights, start = start, etastart = etastart, 
#>         mustart = mustart, offset = offset, family = family, 
#>         control = control, intercept = attr(mt, "intercept") > 
#>             0L, singular.ok = singular.ok))
#>     if (length(offset) && attr(mt, "intercept") > 0L) {
#>         fit2 <- eval(call(if (is.function(method)) "method" else method, 
#>             x = X[, "(Intercept)", drop = FALSE], y = Y, mustart = fit$fitted.values, 
#>             weights = weights, offset = offset, family = family, 
#>             control = control, intercept = TRUE))
#>         if (!fit2$converged) 
#>             warning("fitting to calculate the null deviance did not converge -- increase 'maxit'?")
#>         fit$null.deviance <- fit2$deviance
#>     }
#>     if (model) 
#>         fit$model <- mf
#>     fit$na.action <- attr(mf, "na.action")
#>     if (x) 
#>         fit$x <- X
#>     if (!y) 
#>         fit$y <- NULL
#>     structure(c(fit, list(call = cal, formula = formula, terms = mt, 
#>         data = data, offset = offset, control = control, method = method, 
#>         contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, 
#>             mf))), class = c(fit$class, c("glm", "lm")))
#> })(formula = "y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7", family = structure(list(
#>     family = "binomial", link = "logit", linkfun = function (mu) 
#>     .Call(C_logit_link, mu), linkinv = function (eta) 
#>     .Call(C_logit_linkinv, eta), variance = function (mu) 
#>     mu * (1 - mu), dev.resids = function (y, mu, wt) 
#>     .Call(C_binomial_dev_resids, y, mu, wt), aic = function (y, 
#>         n, mu, wt, dev) 
#>     {
#>         m <- if (any(n > 1)) 
#>             n
#>         else wt
#>         -2 * sum(ifelse(m > 0, (wt/m), 0) * dbinom(round(m * 
#>             y), round(m), mu, log = TRUE))
#>     }, mu.eta = function (eta) 
#>     .Call(C_logit_mu_eta, eta), initialize = {
#>         if (NCOL(y) == 1) {
#>             if (is.factor(y)) 
#>                 y <- y != levels(y)[1L]
#>             n <- rep.int(1, nobs)
#>             y[weights == 0] <- 0
#>             if (any(y < 0 | y > 1)) 
#>                 stop("y values must be 0 <= y <= 1")
#>             mustart <- (weights * y + 0.5)/(weights + 1)
#>             m <- weights * y
#>             if ("binomial" == "binomial" && any(abs(m - round(m)) > 
#>                 0.001)) 
#>                 warning(gettextf("non-integer #successes in a %s glm!", 
#>                   "binomial"), domain = NA)
#>         }
#>         else if (NCOL(y) == 2) {
#>             if ("binomial" == "binomial" && any(abs(y - round(y)) > 
#>                 0.001)) 
#>                 warning(gettextf("non-integer counts in a %s glm!", 
#>                   "binomial"), domain = NA)
#>             n <- (y1 <- y[, 1L]) + y[, 2L]
#>             y <- y1/n
#>             if (any(n0 <- n == 0)) 
#>                 y[n0] <- 0
#>             weights <- weights * n
#>             mustart <- (n * y + 0.5)/(n + 1)
#>         }
#>         else stop(gettextf("for the '%s' family, y must be a vector of 0 and 1's\nor a 2 column matrix where col 1 is no. successes and col 2 is no. failures", 
#>             "binomial"), domain = NA)
#>     }, validmu = function (mu) 
#>     all(is.finite(mu)) && all(mu > 0 & mu < 1), valideta = function (eta) 
#>     TRUE, simulate = function (object, nsim) 
#>     {
#>         ftd <- fitted(object)
#>         n <- length(ftd)
#>         ntot <- n * nsim
#>         wts <- object$prior.weights
#>         if (any(wts%%1 != 0)) 
#>             stop("cannot simulate from non-integer prior.weights")
#>         if (!is.null(m <- object$model)) {
#>             y <- model.response(m)
#>             if (is.factor(y)) {
#>                 yy <- factor(1 + rbinom(ntot, size = 1, prob = ftd), 
#>                   labels = levels(y))
#>                 split(yy, rep(seq_len(nsim), each = n))
#>             }
#>             else if (is.matrix(y) && ncol(y) == 2) {
#>                 yy <- vector("list", nsim)
#>                 for (i in seq_len(nsim)) {
#>                   Y <- rbinom(n, size = wts, prob = ftd)
#>                   YY <- cbind(Y, wts - Y)
#>                   colnames(YY) <- colnames(y)
#>                   yy[[i]] <- YY
#>                 }
#>                 yy
#>             }
#>             else rbinom(ntot, size = wts, prob = ftd)/wts
#>         }
#>         else rbinom(ntot, size = wts, prob = ftd)/wts
#>     }), class = "family"), data = structure(list(y = c(1, 1, 
#> 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 
#> 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 
#> 1, 1, 1, 1, 1, 1), x1 = structure(c(2L, 1L, 2L, 2L, 2L, 1L, 2L, 
#> 2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 
#> 2L, 2L, 1L, 1L, 1L, 2L, 1L, 2L, 2L, 2L, 1L, 2L, 1L, 1L, 1L, 1L, 
#> 2L, 2L, 1L, 2L, 1L, 1L, 2L, 1L, 2L, 1L, 2L), levels = c("control", 
#> "treatment"), class = "factor"), x2 = structure(c(2L, 1L, 3L, 
#> 2L, 2L, 2L, 3L, 2L, 1L, 2L, 2L, 2L, 2L, 3L, 1L, 3L, 2L, 3L, 3L, 
#> 3L, 1L, 2L, 2L, 1L, 1L, 2L, 3L, 2L, 1L, 1L, 3L, 2L, 2L, 2L, 1L, 
#> 2L, 2L, 2L, 3L, 3L, 3L, 2L, 2L, 3L, 1L, 1L, 1L, 1L, 2L, 2L), levels = c("small", 
#> "medium", "large"), class = "factor"), x3 = structure(c(3L, 3L, 
#> 1L, 3L, 1L, 3L, 2L, 2L, 3L, 3L, 3L, 3L, 1L, 1L, 1L, 2L, 1L, 3L, 
#> 1L, 3L, 1L, 1L, 3L, 1L, 3L, 2L, 1L, 1L, 1L, 3L, 3L, 1L, 3L, 3L, 
#> 3L, 1L, 2L, 3L, 1L, 2L, 2L, 1L, 1L, 2L, 2L, 3L, 1L, 2L, 3L, 1L
#> ), levels = c("short", "average", "tall"), class = "factor"), 
#>     x4 = structure(c(2L, 3L, 1L, 2L, 2L, 2L, 2L, 2L, 3L, 2L, 
#>     2L, 2L, 2L, 1L, 2L, 1L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 3L, 2L, 
#>     2L, 2L, 1L, 2L, 3L, 1L, 2L, 2L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 
#>     1L, 2L, 2L, 1L, 3L, 3L, 3L, 3L, 2L, 2L), levels = c("lowbit", 
#>     "most", "highbit"), class = "factor"), x5 = structure(c(3L, 
#>     2L, 1L, 3L, 1L, 3L, 2L, 1L, 3L, 2L, 1L, 3L, 1L, 2L, 3L, 1L, 
#>     1L, 3L, 2L, 1L, 3L, 2L, 1L, 3L, 3L, 1L, 1L, 1L, 1L, 3L, 1L, 
#>     2L, 3L, 1L, 1L, 1L, 3L, 3L, 2L, 1L, 1L, 3L, 1L, 2L, 3L, 2L, 
#>     3L, 3L, 3L, 3L), levels = c("under", "healthy", "over"), class = "factor"), 
#>     x6 = structure(c(1L, 2L, 1L, 2L, 2L, 1L, 2L, 3L, 1L, 3L, 
#>     2L, 3L, 2L, 2L, 2L, 1L, 1L, 2L, 2L, 3L, 1L, 2L, 1L, 3L, 3L, 
#>     2L, 2L, 2L, 1L, 1L, 3L, 3L, 1L, 1L, 2L, 2L, 1L, 1L, 1L, 1L, 
#>     1L, 2L, 1L, 2L, 2L, 1L, 3L, 1L, 2L, 3L), levels = c("small", 
#>     "medium", "large"), class = "factor"), x7 = structure(c(2L, 
#>     3L, 1L, 3L, 2L, 1L, 1L, 3L, 3L, 3L, 3L, 3L, 1L, 1L, 3L, 1L, 
#>     1L, 3L, 1L, 2L, 2L, 2L, 1L, 3L, 3L, 1L, 1L, 1L, 2L, 3L, 1L, 
#>     2L, 1L, 3L, 3L, 1L, 2L, 2L, 1L, 1L, 1L, 3L, 1L, 1L, 3L, 3L, 
#>     3L, 3L, 3L, 2L), levels = c("small", "medium", "large"), class = "factor")), class = "data.frame", row.names = c(NA, 
#> -50L)))
#> 
#> Deviance Residuals: 
#>     Min       1Q   Median       3Q      Max  
#> -1.7567   0.0000   0.0000   0.0000   0.9994  
#> 
#> Coefficients:
#>               Estimate Std. Error z value Pr(>|z|)
#> (Intercept) -2.897e+01  3.924e+04  -0.001    0.999
#> x1treatment  8.684e-01  2.203e+00   0.394    0.693
#> x2medium    -2.759e+01  3.604e+04  -0.001    0.999
#> x2large     -3.264e+01  3.809e+04  -0.001    0.999
#> x3average    2.511e+01  1.841e+04   0.001    0.999
#> x3tall       4.673e+00  1.769e+04   0.000    1.000
#> x4most       9.150e-12  2.507e+00   0.000    1.000
#> x4highbit   -2.219e+01  3.954e+04  -0.001    1.000
#> x5healthy    4.025e+01  1.611e+04   0.002    0.998
#> x5over       5.233e+01  2.104e+04   0.002    0.998
#> x6medium    -4.181e+00  1.968e+04   0.000    1.000
#> x6large      5.650e+01  2.313e+04   0.002    0.998
#> x7medium    -8.684e-01  2.203e+00  -0.394    0.693
#> x7large      2.995e+01  2.079e+04   0.001    0.999
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 69.3147  on 49  degrees of freedom
#> Residual deviance:  9.0793  on 36  degrees of freedom
#> AIC: 37.079
#> 
#> Number of Fisher Scoring iterations: 22
  
  # Summarize pca model with 5 principal components
  summary(compswap_results$all_models$model_pca_4)
#> 
#> Call:
#> (function (formula, family = gaussian, data, weights, subset, 
#>     na.action, start = NULL, etastart, mustart, offset, control = list(...), 
#>     model = TRUE, method = "glm.fit", x = FALSE, y = TRUE, singular.ok = TRUE, 
#>     contrasts = NULL, ...) 
#> {
#>     cal <- match.call()
#>     if (is.character(family)) 
#>         family <- get(family, mode = "function", envir = parent.frame())
#>     if (is.function(family)) 
#>         family <- family()
#>     if (is.null(family$family)) {
#>         print(family)
#>         stop("'family' not recognized")
#>     }
#>     if (missing(data)) 
#>         data <- environment(formula)
#>     mf <- match.call(expand.dots = FALSE)
#>     m <- match(c("formula", "data", "subset", "weights", "na.action", 
#>         "etastart", "mustart", "offset"), names(mf), 0L)
#>     mf <- mf[c(1L, m)]
#>     mf$drop.unused.levels <- TRUE
#>     mf[[1L]] <- quote(stats::model.frame)
#>     mf <- eval(mf, parent.frame())
#>     if (identical(method, "model.frame")) 
#>         return(mf)
#>     if (!is.character(method) && !is.function(method)) 
#>         stop("invalid 'method' argument")
#>     if (identical(method, "glm.fit")) 
#>         control <- do.call("glm.control", control)
#>     mt <- attr(mf, "terms")
#>     Y <- model.response(mf, "any")
#>     if (length(dim(Y)) == 1L) {
#>         nm <- rownames(Y)
#>         dim(Y) <- NULL
#>         if (!is.null(nm)) 
#>             names(Y) <- nm
#>     }
#>     X <- if (!is.empty.model(mt)) 
#>         model.matrix(mt, mf, contrasts)
#>     else matrix(, NROW(Y), 0L)
#>     weights <- as.vector(model.weights(mf))
#>     if (!is.null(weights) && !is.numeric(weights)) 
#>         stop("'weights' must be a numeric vector")
#>     if (!is.null(weights) && any(weights < 0)) 
#>         stop("negative weights not allowed")
#>     offset <- as.vector(model.offset(mf))
#>     if (!is.null(offset)) {
#>         if (length(offset) != NROW(Y)) 
#>             stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
#>                 length(offset), NROW(Y)), domain = NA)
#>     }
#>     mustart <- model.extract(mf, "mustart")
#>     etastart <- model.extract(mf, "etastart")
#>     fit <- eval(call(if (is.function(method)) "method" else method, 
#>         x = X, y = Y, weights = weights, start = start, etastart = etastart, 
#>         mustart = mustart, offset = offset, family = family, 
#>         control = control, intercept = attr(mt, "intercept") > 
#>             0L, singular.ok = singular.ok))
#>     if (length(offset) && attr(mt, "intercept") > 0L) {
#>         fit2 <- eval(call(if (is.function(method)) "method" else method, 
#>             x = X[, "(Intercept)", drop = FALSE], y = Y, mustart = fit$fitted.values, 
#>             weights = weights, offset = offset, family = family, 
#>             control = control, intercept = TRUE))
#>         if (!fit2$converged) 
#>             warning("fitting to calculate the null deviance did not converge -- increase 'maxit'?")
#>         fit$null.deviance <- fit2$deviance
#>     }
#>     if (model) 
#>         fit$model <- mf
#>     fit$na.action <- attr(mf, "na.action")
#>     if (x) 
#>         fit$x <- X
#>     if (!y) 
#>         fit$y <- NULL
#>     structure(c(fit, list(call = cal, formula = formula, terms = mt, 
#>         data = data, offset = offset, control = control, method = method, 
#>         contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, 
#>             mf))), class = c(fit$class, c("glm", "lm")))
#> })(formula = y ~ x1 + PC1 + PC2 + PC3 + PC4 + PC5, family = structure(list(
#>     family = "binomial", link = "logit", linkfun = function (mu) 
#>     .Call(C_logit_link, mu), linkinv = function (eta) 
#>     .Call(C_logit_linkinv, eta), variance = function (mu) 
#>     mu * (1 - mu), dev.resids = function (y, mu, wt) 
#>     .Call(C_binomial_dev_resids, y, mu, wt), aic = function (y, 
#>         n, mu, wt, dev) 
#>     {
#>         m <- if (any(n > 1)) 
#>             n
#>         else wt
#>         -2 * sum(ifelse(m > 0, (wt/m), 0) * dbinom(round(m * 
#>             y), round(m), mu, log = TRUE))
#>     }, mu.eta = function (eta) 
#>     .Call(C_logit_mu_eta, eta), initialize = {
#>         if (NCOL(y) == 1) {
#>             if (is.factor(y)) 
#>                 y <- y != levels(y)[1L]
#>             n <- rep.int(1, nobs)
#>             y[weights == 0] <- 0
#>             if (any(y < 0 | y > 1)) 
#>                 stop("y values must be 0 <= y <= 1")
#>             mustart <- (weights * y + 0.5)/(weights + 1)
#>             m <- weights * y
#>             if ("binomial" == "binomial" && any(abs(m - round(m)) > 
#>                 0.001)) 
#>                 warning(gettextf("non-integer #successes in a %s glm!", 
#>                   "binomial"), domain = NA)
#>         }
#>         else if (NCOL(y) == 2) {
#>             if ("binomial" == "binomial" && any(abs(y - round(y)) > 
#>                 0.001)) 
#>                 warning(gettextf("non-integer counts in a %s glm!", 
#>                   "binomial"), domain = NA)
#>             n <- (y1 <- y[, 1L]) + y[, 2L]
#>             y <- y1/n
#>             if (any(n0 <- n == 0)) 
#>                 y[n0] <- 0
#>             weights <- weights * n
#>             mustart <- (n * y + 0.5)/(n + 1)
#>         }
#>         else stop(gettextf("for the '%s' family, y must be a vector of 0 and 1's\nor a 2 column matrix where col 1 is no. successes and col 2 is no. failures", 
#>             "binomial"), domain = NA)
#>     }, validmu = function (mu) 
#>     all(is.finite(mu)) && all(mu > 0 & mu < 1), valideta = function (eta) 
#>     TRUE, simulate = function (object, nsim) 
#>     {
#>         ftd <- fitted(object)
#>         n <- length(ftd)
#>         ntot <- n * nsim
#>         wts <- object$prior.weights
#>         if (any(wts%%1 != 0)) 
#>             stop("cannot simulate from non-integer prior.weights")
#>         if (!is.null(m <- object$model)) {
#>             y <- model.response(m)
#>             if (is.factor(y)) {
#>                 yy <- factor(1 + rbinom(ntot, size = 1, prob = ftd), 
#>                   labels = levels(y))
#>                 split(yy, rep(seq_len(nsim), each = n))
#>             }
#>             else if (is.matrix(y) && ncol(y) == 2) {
#>                 yy <- vector("list", nsim)
#>                 for (i in seq_len(nsim)) {
#>                   Y <- rbinom(n, size = wts, prob = ftd)
#>                   YY <- cbind(Y, wts - Y)
#>                   colnames(YY) <- colnames(y)
#>                   yy[[i]] <- YY
#>                 }
#>                 yy
#>             }
#>             else rbinom(ntot, size = wts, prob = ftd)/wts
#>         }
#>         else rbinom(ntot, size = wts, prob = ftd)/wts
#>     }), class = "family"), data = structure(list(y = c(1, 1, 
#> 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 
#> 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 
#> 1, 1, 1, 1, 1, 1), x1 = structure(c(2L, 1L, 2L, 2L, 2L, 1L, 2L, 
#> 2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 
#> 2L, 2L, 1L, 1L, 1L, 2L, 1L, 2L, 2L, 2L, 1L, 2L, 1L, 1L, 1L, 1L, 
#> 2L, 2L, 1L, 2L, 1L, 1L, 2L, 1L, 2L, 1L, 2L), levels = c("control", 
#> "treatment"), class = "factor"), PC1 = c(-0.174985567831333, 
#> -1.51738388060647, 1.47214289475174, -0.653351839178653, 0.643321490538164, 
#> 0.0509214014196, 0.697388287504449, -0.03899627224904, -1.80553625569436, 
#> -0.346075112611793, -0.12281376553672, -0.657034811896078, 0.869228459789096, 
#> 1.22975719619769, -0.878481442620089, 1.24906705183288, 0.892035783985493, 
#> -0.382537794186797, 1.22975719619769, 0.599876253888636, -0.400115171272768, 
#> 0.423743116180515, 0.581459475061534, -1.52145024368422, -1.19274072426147, 
#> 0.646152616870241, 1.14004250478095, 1.17852152556349, 0.130422902369166, 
#> -1.80553625569436, 1.13507628891396, 0.42006014346309, 0.0509214014196, 
#> -0.100006441340323, -1.29780550624882, 1.17852152556349, -0.0874851018262282, 
#> -0.814271396178037, 1.25256452039409, 0.939773986058494, 1.24906705183288, 
#> -0.342775530254694, 0.892035783985493, 1.00668135327884, -1.74084311388565, 
#> -1.49457655641007, -1.52145024368422, -1.71803578968925, -0.653351839178653, 
#> 0.109100444178805), PC2 = c(1.43220443157834, -0.291379536117545, 
#> 0.648795427435331, -0.432109433894531, -0.998624157524348, 1.58585626892584, 
#> -0.17990716202936, -0.914465515204791, 1.2503817945849, -0.670425582891932, 
#> -0.410438713236568, -0.683127255239056, -0.844972320176844, -0.914636623925082, 
#> -1.38752517716496, 1.29382275352271, 0.709490682872729, -0.403102553523582, 
#> -0.914636623925082, -0.322598791786851, 0.476788688307905, -1.00759320583519, 
#> 1.6075269895838, -1.45313533652027, -0.740506691798413, -0.199944994089469, 
#> -0.815965439805896, -0.934674455985191, 0.498459408965869, 1.2503817945849, 
#> -0.258649090247694, -1.25861102717971, 1.58585626892584, 1.14402428981301, 
#> -0.282410487806705, -0.934674455985191, 1.17919545095464, 1.61761209356755, 
#> 0.639826379124491, 1.38352488933105, 1.29382275352271, -1.3301457406056, 
#> 0.709490682872729, -0.269609297837707, -0.557090189088367, 1.26308346693203, 
#> -1.45313533652027, 0.997372813961206, -0.432109433894531, -1.27131269952684
#> ), PC3 = c(0.449785884690901, 0.308640006745756, -0.902416024113477, 
#> 1.45758826948625, -0.639743989196388, 0.284492917553269, 0.862910897615201, 
#> 1.08014368947902, -0.37450516649822, 1.57517295286018, 1.47237853932543, 
#> 1.56650405973069, -0.80503695633402, -0.234061120708687, -1.10538603669352, 
#> 0.37522975963663, -1.47951323644851, 1.85401171639506, -0.234061120708687, 
#> 1.64439167179783, -2.11318842148887, -0.645865365906085, 0.299283187392454, 
#> -1.36990951998025, 0.782326177391888, 0.472608827416086, -0.408613509425209, 
#> -0.624363190907802, -2.09839815164968, -0.37450516649822, 1.65977247008642, 
#> -0.53694957566164, 0.284492917553269, 0.797902259210947, 0.314761383455453, 
#> -0.624363190907802, -0.0513647553999538, 0.0763466111597259, 
#> -0.908537400823174, 0.194555994210411, 0.37522975963663, -0.321208154354713, 
#> -1.47951323644851, 1.04358466304142, -0.201179526474588, -0.365836273368731, 
#> -1.36990951998025, -0.875655806589075, 1.45758826948625, -0.545618468791128
#> ), PC4 = c(1.41043171825911, -1.66952763368889, 0.00806621520617041, 
#> 1.04739675792601, -0.883793949504214, 1.50645633422901, 0.251645238448398, 
#> -1.35005368172934, -0.0545824909985992, -0.42550858271601, -1.41890351166185, 
#> 1.02004269489766, -0.787769333534313, 0.859421342103733, 0.866776557209718, 
#> -0.23719835873413, -0.618375468457685, 1.31132691187483, 0.859421342103733, 
#> -0.988686325484912, 1.22981151754282, 0.136955042469977, -0.959843935358847, 
#> 0.0901380477975801, 0.4979540272802, -1.03303390747461, -0.523839179585493, 
#> -0.425257803819279, -1.23648875204504, -0.0545824909985992, -0.530150179799976, 
#> 0.109600979441622, 1.50645633422901, -1.24950964658522, -2.69027662566309, 
#> -0.425257803819279, 1.50663561121997, 0.661147271875326, 1.02881520718036, 
#> -0.599709888449165, -0.23719835873413, 1.38886522482717, -0.618375468457685, 
#> 0.614156768163433, -0.127772463114366, -1.50013376861227, 0.0901380477975801, 
#> 0.0416214019622631, 1.04739675792601, 1.55515225705529), PC5 = c(0.0348022564926567, 
#> -0.890986783794518, 1.05392622708464, 0.771152127648713, -0.0906445840156123, 
#> -0.880860990272764, -1.87737186652042, 1.27545855833962, -0.112258093969213, 
#> 0.923309354501443, 1.33661093166784, 0.591881869102698, -1.00630783078103, 
#> -0.2902612667598, 1.85485731891141, 0.752793581605774, 0.103848344443018, 
#> 0.435252056008574, -0.2902612667598, -1.02506544489842, 1.11850744775535, 
#> -0.324675902635999, -0.315402186253632, -0.982433997042629, 1.25633652966874, 
#> -1.3074404762599, -1.34220790242117, 0.279670123500725, 1.68396625177449, 
#> -0.112258093969213, -0.654750737382078, -0.503946161182014, -0.880860990272764, 
#> 2.44676710689189, -0.656955465174132, 0.279670123500725, 0.152920141710444, 
#> -2.62321880091537, 0.81989490846425, -0.533184372675984, 0.752793581605774, 
#> 1.19040265834536, 0.103848344443018, -0.591393912238663, -1.10429638397548, 
#> 0.219169391429533, -0.982433997042629, 0.00585979124857411, 0.771152127648713, 
#> -0.835373646580759)), class = "data.frame", row.names = c("1", 
#> "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", 
#> "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", 
#> "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", 
#> "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", 
#> "47", "48", "49", "50")))
#> 
#> Deviance Residuals: 
#>      Min        1Q    Median        3Q       Max  
#> -2.05463  -0.14828  -0.00598   0.13932   2.09125  
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)   
#> (Intercept)  -0.7371     0.8436  -0.874  0.38229   
#> x1treatment   1.9810     1.3851   1.430  0.15263   
#> PC1          -3.3874     1.2246  -2.766  0.00567 **
#> PC2          -1.6857     0.9528  -1.769  0.07685 . 
#> PC3           1.0024     0.7908   1.268  0.20493   
#> PC4           2.4352     0.9726   2.504  0.01229 * 
#> PC5           0.4204     0.5355   0.785  0.43238   
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 69.315  on 49  degrees of freedom
#> Residual deviance: 23.233  on 43  degrees of freedom
#> AIC: 37.233
#> 
#> Number of Fisher Scoring iterations: 7
  
  # Get model comparisons
  print(compswap_results$all_comparisons)
#>   model pseudo_r_squared      AIC      BIC   model_set
#> 1   Raw        0.8690138 37.07927 63.84759   model_raw
#> 2   PCA        0.3827987 50.78113 58.42923 model_pca_1
#> 3   PCA        0.4118838 50.76511 60.32522 model_pca_2
#> 4   PCA        0.6023670 39.56182 51.03395 model_pca_3
#> 5   PCA        0.6648237 37.23265 50.61681 model_pca_4
```
