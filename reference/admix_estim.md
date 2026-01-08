# Estimate the unknown weight in an admixture model

Estimate the unknown component weight (and possibly location shift
parameter in case of a symmetric unknown component density), using
different estimation techniques. We remind that the i-th admixture model
has probability density function (pdf) l_i such that: l_i = p_i \* f_i +
(1-p_i) \* g_i, where g_i is the known component density. The unknown
quantities p_i and f_i then have to be estimated.

## Usage

``` r
admix_estim(samples, admixMod, est_method = c("PS", "BVdk", "IBM"), ...)
```

## Arguments

- samples:

  A list of the K (K\>0) samples to be studied, all following admixture
  distributions.

- admixMod:

  A list of objects of class [admix_model](admix_model.md), containing
  useful information about distributions and parameters.

- est_method:

  The estimation method to be applied. Can be one of 'BVdk' (Bordes and
  Vandekerkhove estimator), 'PS' (Patra and Sen estimator), or 'IBM'
  (Inversion Best-Matching approach) in the continuous case (continuous
  random variable). Only 'IBM' for discrete random variables. The same
  estimation method is performed on each sample if several samples are
  provided.

- ...:

  Optional arguments to [estim_PS](estim_PS.md),
  [estim_BVdk](estim_BVdk.md) or [estim_IBM](estim_IBM.md) depending on
  the choice made by the user for the estimation method.

## Value

An object of class `estim_BVdk`, `estim_PS` or `estim_IBM` (that
inherits from class admix_estim), containing at least 5 attributes: 1)
the number of samples under study; 2) the information about the mixture
components (distributions and parameters); 3) the sizes of the samples;
4) the chosen estimation technique (one of 'BVdk', 'PS' or 'IBM'); 5)
the estimated mixing proportions (weights of the unknown component
distributions in the mixture model). In case of 'BVdk' estimation, one
additional attribute corresponding to the estimated location shift
parameter is included.

## Details

For further details on the different estimation techniques, see
references below on i) Patra and Sen estimator ; ii) Bordes and
Vandekerkhove estimator ; iii) Inversion Best-Matching approach.
Important note: estimation by 'IBM' requires at least two samples at
hand, and provides unbiased estimators only if the distributions of
unknown components are equal (meaning that it requires to perform
previously this test between the pairs of samples, see
[admix_test](admix_test.md).

## References

Patra RK, Sen B (2016). “Estimation of a two-component mixture model
with applications to multiple testing.” *Journal of the Royal
Statistical Society Series B*, **78**(4), 869-893. Bordes L, Delmas C,
Vandekerkhove P (2006). “Semiparametric Estimation of a Two-Component
Mixture Model Where One Component Is Known.” *Scandinavian Journal of
Statistics*, **33**(4), 733–752. ISSN 03036898, 14679469,
<http://www.jstor.org/stable/4616955>. Bordes L, Vandekerkhove P (2010).
“Semiparametric two-component mixture model with a known component: An
asymptotically normal estimator.” *Mathematical Methods of Statistics*,
**19**(1), 22–41.
[doi:10.3103/S1066530710010023](https://doi.org/10.3103/S1066530710010023)
. Milhaud X, Pommeret D, Salhi Y, Vandekerkhove P (2024). “Two-sample
contamination model test.” *Bernoulli*, **30**(1), 170–197.
[doi:10.3150/23-BEJ1593](https://doi.org/10.3150/23-BEJ1593) .

## Author

Xavier Milhaud <xavier.milhaud.research@gmail.com>

## Examples

``` r
## Simulate mixture data:
mixt1 <- twoComp_mixt(n = 300, weight = 0.7,
                      comp.dist = list("norm", "norm"),
                      comp.param = list(list("mean" = -2, "sd" = 0.5),
                                        list("mean" = 0, "sd" = 1)))
mixt2 <- twoComp_mixt(n = 250, weight = 0.85,
                      comp.dist = list("norm", "exp"),
                      comp.param = list(list("mean" = -2, "sd" = 0.5),
                                        list("rate" = 1)))
mixt3 <- twoComp_mixt(n = 500, weight = 0.5,
                      comp.dist = list("pois", "pois"),
                      comp.param = list(list("lambda" = 2),
                                        list("lambda" = 7)))
mixt4 <- twoComp_mixt(n = 1500, weight = 0.2, comp.dist = list("multinom", "multinom"),
                      comp.param = list(list("size"=1, "prob" = c(0.8,0.1,0.1)),
                                   list("size"=1, "prob" = c(0.1,0.2,0.7))))
data1 <- get_mixture_data(mixt1)
data2 <- get_mixture_data(mixt2)
data3 <- get_mixture_data(mixt3)
data4 <- get_mixture_data(mixt4)
## Define the admixture models:
admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
                         knownComp_param = mixt1$comp.param[[2]])
admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
                         knownComp_param = mixt2$comp.param[[2]])
admixMod3 <- admix_model(knownComp_dist = mixt3$comp.dist[[2]],
                         knownComp_param = mixt3$comp.param[[2]])
admixMod4 <- admix_model(knownComp_dist = mixt4$comp.dist[[2]],
                         knownComp_param = mixt4$comp.param[[2]])
# Estimation by different methods:
admix_estim(samples = list(data1), admixMod = list(admixMod1), est_method = "BVdk")
#> Mixing weight estimation using 'BVdk' assumes the unknown component
#> distribution to have a symmetric probability density function.
#> 
#> Call:admix_estim(samples = list(data1), admixMod = list(admixMod1), 
#>     est_method = "BVdk")
#> 
#> ******** Sample #1 ********
#> Estimated mixing weight: 0.685 / Estimated location shift: -2.052 
#> 
admix_estim(samples = list(data1, data2, data3, data4),
            admixMod = list(admixMod1, admixMod2, admixMod3, admixMod4), est_method = "PS")
#> 
#> Call:admix_estim(samples = list(data1, data2, data3, data4), admixMod = list(admixMod1, 
#>     admixMod2, admixMod3, admixMod4), est_method = "PS")
#> 
#> ******** Sample #1 ********
#>  Estimated mixing weight (of the unknown component): 0.655
#>  Selected c_n equals 0.174 in the penalization term. See ?estim_PS
#> 
#> ******** Sample #2 ********
#>  Estimated mixing weight (of the unknown component): 0.841
#>  Selected c_n equals 0.171 in the penalization term. See ?estim_PS
#> 
#> ******** Sample #3 ********
#>  Estimated mixing weight (of the unknown component): 0.458
#>  Selected c_n equals 0.183 in the penalization term. See ?estim_PS
#> 
#> ******** Sample #4 ********
#>  Estimated mixing weight (of the unknown component): 0.203
#>  Selected c_n equals 0.199 in the penalization term. See ?estim_PS
#> 
admix_estim(samples = list(data1,data2), admixMod = list(admixMod1,admixMod2), est_method = "IBM")
#>  IBM estimators of two unknown proportions are reliable only if the two corresponding
#>  unknown component distributions have previously been tested equal (see ?admix_test).
#> 
#> Call:admix_estim(samples = list(data1, data2), admixMod = list(admixMod1, 
#>     admixMod2), est_method = "IBM")
#> 
#> Pairwise estimation performed (IBM estimation method).
#> 
#> ******** Samples #1 with #2 ********
#> Estimated mixing weight of the unknown distribution in the 1st sample: 0.748 
#> Estimated mixing weight of the unknown distribution in the 2nd sample: 0.895 
#> 
```
