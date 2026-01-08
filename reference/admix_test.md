# Equality test for the unknown components in admixture models

Perform hypothesis test between unknown components of a list of
admixture models, where we remind that the i-th admixture model has
probability density function (pdf) l_i such that: l_i = p_i \* f_i +
(1-p_i) \* g_i, with g_i the known component density. The unknown
quantities p_i and f_i are thus estimated, leading to the test given by
the following null and alternative hypothesis: H0: f_i = f_j for all i
!= j against H1 : there exists at least i != j such that f_i differs
from f_j. The test can be performed using two methods, either the
comparison of coefficients obtained through polynomial basis expansions
of the component densities, or by the inner-convergence property
obtained using the IBM approach. See 'Details' below for further
information.

## Usage

``` r
admix_test(
  samples,
  admixMod,
  test_method = c("poly", "icv"),
  conf_level = 0.95,
  ...
)
```

## Arguments

- samples:

  A list of the K (K \> 0) samples to be studied, each one assumed to
  follow a mixture distribution.

- admixMod:

  A list of objects of class [admix_model](admix_model.md), containing
  useful information about distributions and parameters of the
  contamination / admixture models under study.

- test_method:

  The testing method to be applied. Can be either 'poly' (polynomial
  basis expansion) or 'icv' (inner convergence from IBM). The same
  testing method is performed between all samples. In the one-sample
  case, only 'poly' is available and the test is a gaussianity test. For
  further details, see section 'Details' below.

- conf_level:

  The confidence level of the K-sample test.

- ...:

  Depending on the choice made by the user for the test method ('poly'
  or 'icv'), optional arguments to
  [gaussianity_test](gaussianity_test.md),
  [orthobasis_test](orthobasis_test.md) (in case of 'poly'), or
  [IBM_k_samples_test](IBM_k_samples_test.md) in case of 'icv'.

## Value

An object of class `gaussianity_test`, `orthobasis_test`, or `IBM_test`
(that inherits from class `htest`), containing attributes specific to
the object class (in addition to classical attributes from `htest`).
Usually, the test decision (reject the null hypothesis or not); the
confidence level of the test (1-alpha, where alpha denotes the level of
the test or equivalently the type-I error); the number of samples under
study; the respective size of each sample; the information about known
mixture components.

## Details

For further details on implemented hypothesis tests, see the references
hereafter. .

## References

Milhaud X, Pommeret D, Salhi Y, Vandekerkhove P (2024).
“Contamination-source based K-sample clustering.” *Journal of Machine
Learning Research*, **25**(287), 1–32.
<https://jmlr.org/papers/v25/23-0914.html>. Milhaud X, Pommeret D, Salhi
Y, Vandekerkhove P (2022). “Semiparametric two-sample admixture
components comparison test: The symmetric case.” *Journal of Statistical
Planning and Inference*, **216**, 135-150. ISSN 0378-3758,
[doi:10.1016/j.jspi.2021.05.010](https://doi.org/10.1016/j.jspi.2021.05.010)
. Pommeret D, Vandekerkhove P (2019). “Semiparametric density testing in
the contamination model.” *Electronic Journal of Statistics*, 4743–4793.
[doi:10.1214/19-EJS1650](https://doi.org/10.1214/19-EJS1650) .

## Author

Xavier Milhaud <xavier.milhaud.research@gmail.com>

## Examples

``` r
####### Example with 2 samples
mixt1 <- twoComp_mixt(n = 380, weight = 0.7,
                      comp.dist = list("norm", "norm"),
                      comp.param = list(list("mean" = -2, "sd" = 0.5),
                                        list("mean" = 0, "sd" = 1)))
mixt2 <- twoComp_mixt(n = 350, weight = 0.85,
                      comp.dist = list("norm", "norm"),
                      comp.param = list(list("mean" = -2, "sd" = 0.5),
                                        list("mean" = -1, "sd" = 1)))
data1 <- get_mixture_data(mixt1)
data2 <- get_mixture_data(mixt2)
admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
                         knownComp_param = mixt1$comp.param[[2]])
admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
                         knownComp_param = mixt2$comp.param[[2]])
admix_test(samples = list(data1,data2), admixMod = list(admixMod1,admixMod2),
           conf_level = 0.95, test_method = "poly", ask_poly_param = FALSE, support = "Real")
#>   Default estimation method is 'BVdk' when testing with polynomial basis expansions (ensuring
#>   theoretical guarantees, but relying on symmetric unknown component densities). To consider
#>   other frameworks in the 2-sample case, use 'PS' estimator (setting argument 'est_method' to 'PS').
#> 
#>  Equality test of unknown distributions with polynomial expansions of
#>  pdfs
#> 
#> data:  samples
#> T = 1.4929, expansion order S = 1, p-value = 0.2218
#> alternative hypothesis: Distributions of unknown components involved 
#>                         in the contamination models are different
#> 
```
