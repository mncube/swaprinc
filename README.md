
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
