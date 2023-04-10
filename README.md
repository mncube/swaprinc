
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
#>              y          x1          x2         x3          x4          x5
#> 1    39.475820  0.47773904  18.7383130  -8.773589 -10.6331796  -6.0540123
#> 2   -63.079847  0.49618282  13.6177720 -13.933367 -21.7238792  -2.0951518
#> 3   -12.182720 -0.85958430  -2.2810604  -1.773555  -0.7757329   0.3537918
#> 4   400.803255 -0.82905996   5.2361401  62.775202  97.0701178  16.0071664
#> 5   107.198006 -0.32157308  20.2132202  -7.899481  -6.4430556  -8.6888340
#> 6   -64.241618 -1.30377040   5.8722388 -17.541489 -24.4387685 -12.6625036
#> 7   129.510782 -1.42148660  39.0244323  -6.002883  -8.4235930  -5.8399478
#> 8   -53.124280  1.74491495  15.6297978 -23.468726 -33.5144449 -12.1440258
#> 9     0.761734 -0.28827936  28.5871818 -12.469999 -19.2738582  -4.6809061
#> 10   35.644342 -1.30886572  18.3080345 -12.286187 -11.8774819  -2.8944067
#> 11  -60.635240 -0.06945219 -16.1750371   1.317377   2.0170811  -3.2399018
#> 12  272.374932 -1.22492668  17.0566269  34.717773  52.1120135   0.5874081
#> 13 -199.534216  0.80899626  -3.2857059 -35.477599 -41.8802165  -6.0349577
#> 14  179.443269 -0.49215034  11.3860745  23.228872  38.8074789   8.3111588
#> 15  -75.935952  0.45269393  -6.5571196 -10.812509 -15.3091826  -2.2768390
#> 16  -88.888762  0.99963310  -4.4288023 -18.425761 -21.8841270 -11.8480187
#> 17  -83.432877  0.46702960   4.7241553 -14.012933 -22.7163856  -6.9538292
#> 18  -71.532618  0.37605199  -4.8996119 -11.818129 -12.5014386  -4.7731135
#> 19  -43.211607  1.70349095  -2.6564299  -6.134988 -10.7646334  -4.8909273
#> 20  -96.430833 -1.03546152  14.6955564 -24.506841 -32.9999212  -0.4274484
#> 21  -15.926773  1.32812210  -1.9895743  -8.528803  -9.8073460  -5.2829679
#> 22   55.760262 -0.59428711  25.6319046 -15.546777 -20.0432134 -11.4686783
#> 23  -61.244411  1.61128458  22.9639935 -31.064006 -42.9324710  -8.0258999
#> 24 -143.368296 -1.11267383  -4.2449467 -17.133744 -28.2176092  -4.0042549
#> 25 -122.840965 -1.46018298 -10.5841940 -13.499040 -19.4418277  -6.2930692
#> 26  100.127060  0.73215598  10.9932423   8.075351  14.4056833   2.9301229
#> 27    5.986796 -1.61033943   8.3087131  -8.205605  -3.3703453  -1.6914257
#> 28   49.229898  0.33206736  10.9593663   1.071770   3.3121805  -2.9690864
#> 29  160.573197  0.76085616   5.4054455  20.642732  31.2447633   1.7672148
#> 30  246.402685 -1.85366955  23.5110932  26.622061  41.3701420   9.4795387
#> 31    5.431842  0.79115157   5.4591328  -6.452399  -1.3654737  -3.9401572
#> 32  187.906892 -1.28174039 -15.4587114  32.115782  56.9063296   4.1195489
#> 33   31.941937 -0.77987734  18.5743065  -2.102877  -3.5158984   1.2766131
#> 34 -104.086055  1.43834228 -26.0414246   3.880662   6.7273872  -0.1211065
#> 35  -79.470175 -0.23482954   8.1918257 -19.407632 -25.5101523  -1.0977488
#> 36   70.398697  0.61219714   6.8603411   5.324253   6.1754848  -1.1545379
#> 37 -162.457201 -1.45847063   5.8029624 -35.248747 -46.2305324 -13.6369710
#> 38  324.965968  0.24276907   9.9507763  47.259884  65.9535423   6.9595572
#> 39  126.354990 -0.31739211  29.9739944  -4.100188  -1.5527377  -2.6602058
#> 40 -147.171578  0.85937333   0.4359594 -33.133698 -47.3157301 -11.6322443
#> 41  -60.073277  1.34415507 -16.1773607   1.224252   0.1601444  -3.5143958
#> 42   46.601524  0.31555720  26.3215622 -10.682936 -18.3218265  -7.9862870
#> 43   50.330363  2.04544486  -4.7911703   4.888043  10.2265600   0.1489307
#> 44  -24.756178 -0.18247053 -20.0750636   3.974270  10.1359559  -4.7950096
#> 45   36.303899 -1.48393823  -8.0869199   9.946346  13.5809210   4.6784509
#> 46  203.377610  0.39716563  17.5883375  14.899321  29.1934264  -5.0340896
#> 47  -23.554735  1.56965755 -12.5200583   7.718097  16.4851478   6.8018347
#> 48  -13.904639  1.27023296  21.3109947 -12.070845 -17.9835883  -2.7860730
#> 49   36.835912 -0.87144335  17.8876096 -15.503695 -19.6327828  -9.1540557
#> 50 -149.971367  0.36088326   3.4197420 -23.591043 -32.3348350  -6.6624585
#>             x6         x7          x8           x9         x10
#> 1  -13.3618582  -7.491364 -0.51567685  10.21341530  4.03424043
#> 2  -30.5287909  -8.952915  1.58818675   4.81066064  3.72181797
#> 3   -8.7356168   3.618402  4.53421788  -1.58406980  9.57477135
#> 4   98.6220412 -11.053882  8.61827429   0.87747664  7.69875397
#> 5   -4.8310215   7.137167 12.27702640   8.45320706  2.66501550
#> 6  -30.2004566  -7.038647  0.91699203  -0.87069957  7.62271151
#> 7  -12.0116446  -4.909563  0.02973866  17.71371724 10.95218853
#> 8  -39.8707328   3.019071  7.39229978  -2.79999008 12.24696558
#> 9  -34.2171646  -5.949003  1.97197720  12.68421841 -1.84775250
#> 10 -16.8766719   1.712549 17.35467749   5.71104019  7.49289279
#> 11  -9.2469139   3.768522  9.04674532 -12.08392586 -1.35198076
#> 12  51.5626099  -8.065459 -2.70265502   5.29816491 -3.63447029
#> 13 -44.2910590   4.976009 -6.33640714 -10.80645150  2.45948059
#> 14  34.7779677  -1.121521 -3.26591816   1.02825313 10.93384634
#> 15 -14.5782013  -1.848137  5.59900156  -3.20371056  4.77536227
#> 16 -33.5387803   8.164848 -4.99424662  -4.39788078  1.61423806
#> 17 -32.4237141  -9.129137 -1.23136895   0.82693074 11.16317481
#> 18 -17.4575787  -5.436237  4.00391445  -1.46050273 -3.51983384
#> 19 -19.5484265   1.536239  7.57952853  -3.47077049  5.90806539
#> 20 -42.0710284   1.218807  6.36636742   5.03655258  4.77814132
#> 21  -6.9711081   2.226684  1.86051776  -1.80614508  4.29690063
#> 22 -24.2616653   1.202681 10.75330394  10.01275418  4.90453857
#> 23 -40.8742366   3.187506 10.69822121   6.06804352 -1.12015838
#> 24 -28.6901478 -12.871297  3.48393798  -3.59837556  7.97608829
#> 25 -29.2993920 -10.298432 12.69196787  -3.50241381  5.81041086
#> 26   9.7596034   7.816931 -2.19063811  -0.08869435 10.03462948
#> 27  -9.2981005 -10.465826 12.44393779   0.72326807  6.99286252
#> 28   1.4769861  -6.729397 -3.11884555   1.64009093  2.51191164
#> 29  32.0322469  -3.194732 14.31880695  -3.54332938  4.48702090
#> 30  43.2950888   1.963718  0.73102413   7.33118281  2.51045583
#> 31   0.7504028  -6.350634 -4.17990698  -5.03847190 -1.98700431
#> 32  51.7432504   2.415718  2.69589432  -5.31335765  9.16496512
#> 33 -10.3333026   1.086745 -4.69005811   6.99718842 -0.26110600
#> 34 -11.5427398  -1.428252 -2.17551478 -11.70769366  5.57920639
#> 35 -25.4537456   5.773674 -8.47529357   2.34744411  1.31529391
#> 36   6.1305140   3.522552  9.64599638  -3.92488299  9.12538281
#> 37 -53.6524547  -6.092165  4.01330251  -2.70244700  7.62456392
#> 38  69.4738825  -3.797380  5.49656609   3.17886831  5.86885328
#> 39   1.3078906   1.408648  6.78984234   7.76375160  4.86361455
#> 40 -52.9899865  10.078518  6.38241956  -0.29448291 10.30999783
#> 41  -6.4241700  -5.085768  4.74006505  -8.38602834  3.78348587
#> 42 -18.1028712   3.206713 -2.34838068   9.55498522 -0.58694618
#> 43  12.1717436   4.158772  9.12651941  -4.62940907  4.15538266
#> 44   3.0478804   2.117315 -3.57159982  -9.76076879  6.17141166
#> 45  13.2442334  -4.654611 18.36363071  -2.03799963 11.22556066
#> 46  19.5758519   2.324029 -0.26589530   8.08271530  0.06877747
#> 47  10.2846536 -11.089204 -0.73572517 -10.67085672 10.11413548
#> 48 -17.3308509 -15.365838  1.49703122   7.72735495 -5.94848997
#> 49 -17.8869103   5.186855  3.84717412   8.57907723  8.17863365
#> 50 -38.5748828 -10.094115 -9.11893592  -1.75945401  1.16800493
#>              y          x1         PC1         PC2         PC3
#> 1   0.16868428  0.42678011 -0.72811937  1.57543900  0.47260744
#> 2  -0.62406801  0.44408399 -0.90504308  0.72517780  0.55688758
#> 3  -0.23063472 -0.82789086  0.37141700 -0.95942184 -0.92354152
#> 4   2.96173473 -0.79925305  6.70468267  0.62596833 -0.29692272
#> 5   0.69217479 -0.32313106 -0.83620073  1.39745632 -1.47828568
#> 6  -0.63304847 -1.24462433 -1.70541176 -0.47559207  0.11577096
#> 7   0.86465189 -1.35506516 -0.70669769  2.89273663 -1.02115819
#> 8  -0.54711177  1.61563907 -2.19682801 -0.65082942 -1.79932048
#> 9  -0.13057447 -0.29189503 -1.37894831  2.49950259  0.79241103
#> 10  0.13906706 -1.24940474 -0.55798501  0.73504151 -2.27226429
#> 11 -0.60517127 -0.08659233  0.26253403 -2.07295830  0.41908155
#> 12  1.96898757 -1.17065349  2.89473727  1.92480139  1.97033593
#> 13 -1.67885628  0.73756421 -2.34546966 -1.82411114  1.27062867
#> 14  1.25062854 -0.48316590  2.89059840  0.25308682 -0.14996077
#> 15 -0.72344533  0.40328291 -0.33614468 -1.09520312 -0.03365475
#> 16 -0.82357017  0.91641886 -1.95042961 -1.09486899  0.92350152
#> 17 -0.78139634  0.41673257 -1.12999725 -0.51104099 -0.01419982
#> 18 -0.68940769  0.33137777 -0.66312551 -0.35211540  1.39671425
#> 19 -0.47048711  1.57677526 -0.45895870 -1.03361045 -0.62502375
#> 20 -0.88187016 -0.99289810 -1.53048888  0.54042673 -0.65259949
#> 21 -0.25957614  1.22460582 -0.43087545 -0.69572211  0.11665052
#> 22  0.29456253 -0.57899018 -1.71391360  1.62513407 -1.40921412
#> 23 -0.60988015  1.49026764 -2.56891468  1.27586566 -0.56670018
#> 24 -1.24469537 -1.06533835 -0.87849133 -1.19372662  0.31825568
#> 25 -1.08601970 -1.39136994 -0.78470015 -1.39611921 -0.44988241
#> 26  0.63751658  0.66547299  1.16094495 -0.09545987 -0.63899601
#> 27 -0.09018490 -1.53224607  0.14996234 -0.02776078 -0.81494057
#> 28  0.24408301  0.29011158  0.29915405  0.54586576  1.19522771
#> 29  1.10476344  0.69239937  2.30504528 -0.22624697 -0.90950930
#> 30  1.76822286 -1.76053733  2.95973928  1.92855106  0.17096935
#> 31 -0.09447467  0.72082239 -0.04136065 -0.17322501  2.08520954
#> 32  1.31605210 -1.22395587  3.56666392 -1.49107610 -0.35810341
#> 33  0.11044758 -0.75311013 -0.03738083  1.48261135  1.15918175
#> 34 -0.94104480  1.32801392  0.76548384 -2.77706398  1.05027938
#> 35 -0.75076477 -0.24174865 -1.23280175  0.25998483  1.26309137
#> 36  0.40771721  0.55292812  0.78985204 -0.68034554 -1.44531107
#> 37 -1.39225204 -1.38976342 -2.94061144 -0.86197621 -0.23996766
#> 38  2.37551473  0.20633230  4.42360361  0.91255843 -0.26929213
#> 39  0.84025771 -0.31920849 -0.07927358  1.78652906 -0.92882246
#> 40 -1.27409463  0.78482776 -2.91919242 -1.16059355 -1.71643506
#> 41 -0.60082732  1.23964789  0.40092584 -1.88744440  0.56983793
#> 42  0.22376576  0.27462181 -1.41932139  1.99741763  0.64739476
#> 43  0.25258957  1.89759493  1.03012369 -0.99777932 -0.60909249
#> 44 -0.32782718 -0.19262565  0.59514110 -2.31409735  0.79583821
#> 45  0.14416542 -1.41365701  1.88593019 -1.21189739 -2.06999447
#> 46  1.43564029  0.35118649  0.96227329  1.70438797  0.53749591
#> 47 -0.31854006  1.45121334  2.03484293 -2.10245093  0.69010475
#> 48 -0.24394510  1.17029448 -0.76665357  2.02824867  2.09469407
#> 49  0.14827787 -0.83901696 -1.42969720  0.95220237 -1.22242617
#> 50 -1.29573691  0.31714653 -1.78061938 -0.30625690  2.30344911
  
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
