Joint modeling of endpoints for RCT
================
Lars Mølgaard Saxhaug
2022-09-01

Recreation of analysis from [“Joint modeling of endpoints can be used to
answer various research questions in randomized clinical
trials”](https://www.sciencedirect.com/science/article/pii/S0895435622000701#bib33)
(Eijk et al. 2022) using {rstanarm}(Goodrich et al. 2020)

``` r
survival_data <- import(here("data/VPA_SURVIVAL_dataset.xlsx")) |> 
  janitor::clean_names()
str(survival_data)
```

    ## 'data.frame':    154 obs. of  7 variables:
    ##  $ avisit : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ total  : num  41 45 46 30 30 43 44 40 43 21 ...
    ##  $ maxtime: num  12 16 16 8 12 20 20 20 20 20 ...
    ##  $ trt    : num  0 0 1 1 0 1 0 1 1 0 ...
    ##  $ stime  : num  13.4 18.5 19.4 12.3 15.3 ...
    ##  $ status : num  1 1 1 1 1 0 0 0 0 0 ...
    ##  $ id     : num  1 2 3 4 5 6 7 8 9 10 ...

``` r
head(survival_data)
```

    ##   avisit total maxtime trt    stime status id
    ## 1      0    41      12   0 13.43737      1  1
    ## 2      0    45      16   0 18.49692      1  2
    ## 3      0    46      16   1 19.38398      1  3
    ## 4      0    30       8   1 12.25462      1  4
    ## 5      0    30      12   0 15.34292      1  5
    ## 6      0    43      20   1 20.63244      0  6

``` r
long_data <- import(here("data/VPA_ALSFRS_dataset.xlsx")) |> 
  janitor::clean_names()
str(long_data)
```

    ## 'data.frame':    701 obs. of  7 variables:
    ##  $ avisit : num  0 2 4 8 12 0 2 4 8 12 ...
    ##  $ total  : num  41 41 39 36 31 45 46 44 42 38 ...
    ##  $ maxtime: num  12 12 12 12 12 16 16 16 16 16 ...
    ##  $ trt    : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ stime  : num  13.4 13.4 13.4 13.4 13.4 ...
    ##  $ status : num  1 1 1 1 1 1 1 1 1 1 ...
    ##  $ id     : num  1 1 1 1 1 2 2 2 2 2 ...

``` r
head(long_data)
```

    ##   avisit total maxtime trt    stime status id
    ## 1      0    41      12   0 13.43737      1  1
    ## 2      2    41      12   0 13.43737      1  1
    ## 3      4    39      12   0 13.43737      1  1
    ## 4      8    36      12   0 13.43737      1  1
    ## 5     12    31      12   0 13.43737      1  1
    ## 6      0    45      16   0 18.49692      1  2

``` r
# check if file exists, only run model if it does not
if (file.exists(here("output/fits/model1.Rdata"))) { 
  model1 <- import(here("output/fits/model1.Rdata"))
  
} else {
  model1 <-
    stan_jm(
      formulaLong = total ~ avisit + avisit:trt + (avisit + I(avisit ^ 2) | id), # longitudinal mixed effects model
      dataLong = long_data,
      formulaEvent = Surv(stime, status) ~ trt, # survival model 
      dataEvent = survival_data,
      time_var = "avisit",
      priorEvent = rstanarm::normal(0, 1), # prior on treatment effect in survival model
      iter = 4000
    )
  save(model1, file = here("output/fits/model1.Rdata"))
  
}
```

``` r
prior_summary(model1) # 
```

    ## Priors for model 'model1' 
    ## ------
    ## Long1|Intercept (after predictors centered)
    ##   Specified prior:
    ##     ~ normal(location = 0, scale = 10)
    ##   Adjusted prior:
    ##     ~ normal(location = 0, scale = 89)
    ## 
    ## Long1|Coefficients
    ##   Specified prior:
    ##     ~ normal(location = [0,0], scale = [2.5,2.5])
    ##   Adjusted prior:
    ##     ~ normal(location = [0,0], scale = [3.94,4.61])
    ## 
    ## Long1|Auxiliary (sigma)
    ##   Specified prior:
    ##     ~ half-cauchy(location = 0, scale = 5)
    ##   Adjusted prior:
    ##     ~ half-cauchy(location = 0, scale = 44)
    ## 
    ## Event|Coefficients
    ##  ~ normal(location = 0, scale = 1)
    ## 
    ## Event|Auxiliary (B-spline-coefficients)
    ##  ~ cauchy(location = [0,0,0,...], scale = [20,20,20,...])
    ## 
    ## Association parameters
    ##   Specified prior:
    ##     ~ normal(location = 0, scale = 2.5)
    ##   Adjusted prior:
    ##     ~ normal(location = 0, scale = 0.29)
    ## 
    ## Covariance
    ##  ~ lkj(reg. = 1, df = [1,1,1], scale = [10,10,10])
    ##      **adjusted scale = [88.54,15.78, 0.90]
    ## ------
    ## See help('prior_summary.stanreg') for more details

``` r
summary(model1)
```

    ## 
    ## Model Info:
    ## 
    ##  function:         stan_jm
    ##  formula (Long1):  total ~ avisit + avisit:trt + (avisit + I(avisit^2) | id)
    ##  family  (Long1):  gaussian [identity]
    ##  formula (Event):  Surv(stime, status) ~ trt
    ##  baseline hazard:  bs
    ##  assoc:            etavalue (Long1)
    ##  algorithm:        sampling
    ##  sample:           8000 (posterior sample size)
    ##  priors:           see help('prior_summary')
    ##  observations:     701 (Long1)
    ##  subjects:         154
    ##  events:           33 (21.4%)
    ##  groups:           id (154)
    ##  runtime:          20.6 mins
    ## 
    ## Estimates:
    ##                                                 mean      sd        10%    
    ## Long1|(Intercept)                                40.071     0.461    39.490
    ## Long1|avisit                                     -0.966     0.086    -1.076
    ## Long1|avisit:trt                                 -0.117     0.120    -0.272
    ## Long1|sigma                                       2.093     0.088     1.982
    ## Long1|mean_PPD                                   34.660     0.112    34.516
    ## Event|(Intercept)                                -1.608     0.625    -2.403
    ## Event|trt                                         0.529     0.339     0.102
    ## Event|b-splines-coef1                           -14.578     8.014   -25.245
    ## Event|b-splines-coef2                             1.614     2.374    -1.339
    ## Event|b-splines-coef3                            -0.719     1.472    -2.582
    ## Event|b-splines-coef4                            -0.784     1.244    -2.394
    ## Event|b-splines-coef5                             1.750     1.545    -0.236
    ## Event|b-splines-coef6                            -3.691     2.851    -7.390
    ## Assoc|Long1|etavalue                             -0.082     0.018    -0.105
    ## Sigma[id:Long1|(Intercept),Long1|(Intercept)]    28.912     3.666    24.448
    ## Sigma[id:Long1|avisit,Long1|(Intercept)]          2.077     0.515     1.424
    ## Sigma[id:Long1|I(avisit^2),Long1|(Intercept)]    -0.095     0.024    -0.126
    ## Sigma[id:Long1|avisit,Long1|avisit]               0.857     0.172     0.647
    ## Sigma[id:Long1|I(avisit^2),Long1|avisit]         -0.023     0.008    -0.034
    ## Sigma[id:Long1|I(avisit^2),Long1|I(avisit^2)]     0.001     0.000     0.001
    ## log-posterior                                 -2350.503    27.843 -2386.214
    ##                                                 50%       90%    
    ## Long1|(Intercept)                                40.075    40.658
    ## Long1|avisit                                     -0.964    -0.858
    ## Long1|avisit:trt                                 -0.116     0.037
    ## Long1|sigma                                       2.091     2.206
    ## Long1|mean_PPD                                   34.660    34.801
    ## Event|(Intercept)                                -1.605    -0.815
    ## Event|trt                                         0.524     0.967
    ## Event|b-splines-coef1                           -13.263    -5.556
    ## Event|b-splines-coef2                             1.504     4.673
    ## Event|b-splines-coef3                            -0.693     1.126
    ## Event|b-splines-coef4                            -0.772     0.814
    ## Event|b-splines-coef5                             1.751     3.750
    ## Event|b-splines-coef6                            -3.348    -0.401
    ## Assoc|Long1|etavalue                             -0.082    -0.059
    ## Sigma[id:Long1|(Intercept),Long1|(Intercept)]    28.657    33.820
    ## Sigma[id:Long1|avisit,Long1|(Intercept)]          2.063     2.745
    ## Sigma[id:Long1|I(avisit^2),Long1|(Intercept)]    -0.094    -0.065
    ## Sigma[id:Long1|avisit,Long1|avisit]               0.844     1.080
    ## Sigma[id:Long1|I(avisit^2),Long1|avisit]         -0.022    -0.013
    ## Sigma[id:Long1|I(avisit^2),Long1|I(avisit^2)]     0.001     0.002
    ## log-posterior                                 -2350.161 -2314.674
    ## 
    ## Diagnostics:
    ##                                               mcse  Rhat  n_eff
    ## Long1|(Intercept)                             0.012 1.001  1478
    ## Long1|avisit                                  0.001 1.000  4062
    ## Long1|avisit:trt                              0.002 1.000  5123
    ## Long1|sigma                                   0.002 1.004  1729
    ## Long1|mean_PPD                                0.001 1.000  7737
    ## Event|(Intercept)                             0.005 1.000 15033
    ## Event|trt                                     0.003 1.000 14936
    ## Event|b-splines-coef1                         0.102 1.000  6135
    ## Event|b-splines-coef2                         0.029 1.000  6539
    ## Event|b-splines-coef3                         0.018 1.000  6452
    ## Event|b-splines-coef4                         0.015 1.000  6602
    ## Event|b-splines-coef5                         0.018 1.000  7206
    ## Event|b-splines-coef6                         0.032 1.000  8139
    ## Assoc|Long1|etavalue                          0.000 1.000 14378
    ## Sigma[id:Long1|(Intercept),Long1|(Intercept)] 0.069 1.001  2782
    ## Sigma[id:Long1|avisit,Long1|(Intercept)]      0.008 1.001  3952
    ## Sigma[id:Long1|I(avisit^2),Long1|(Intercept)] 0.000 1.001  5415
    ## Sigma[id:Long1|avisit,Long1|avisit]           0.004 1.003  1919
    ## Sigma[id:Long1|I(avisit^2),Long1|avisit]      0.000 1.002  1706
    ## Sigma[id:Long1|I(avisit^2),Long1|I(avisit^2)] 0.000 1.004  1325
    ## log-posterior                                 0.920 1.008   917
    ## 
    ## For each parameter, mcse is Monte Carlo standard error, n_eff is a crude measure of effective sample size, and Rhat is the potential scale reduction factor on split chains (at convergence Rhat=1).

``` r
ps_check(model1)
```

![](README_files/figure-gfm/ps_check-1.png)<!-- -->

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-vaneijk2022" class="csl-entry">

Eijk, Ruben P. A. van, Kit C. B. Roes, Leonard H. van den Berg, and Ying
Lu. 2022. “Joint Modeling of Endpoints Can Be Used to Answer Various
Research Questions in Randomized Clinical Trials.” *Journal of Clinical
Epidemiology* 147 (July): 32–39.
<https://doi.org/10.1016/j.jclinepi.2022.03.009>.

</div>

<div id="ref-rstanarm" class="csl-entry">

Goodrich, Ben, Jonah Gabry, Imad Ali, and Sam Brilleman. 2020.
“Rstanarm: Bayesian Applied Regression Modeling via Stan.”
<https://mc-stan.org/rstanarm>.

</div>

</div>
