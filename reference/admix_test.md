# Hypothesis test for the unknown component(s) in admixture model(s)

Perform hypothesis test on the unknown component(s) of a list of
admixture model(s), where we remind that the \\i\\-th admixture model
has probability density function (pdf) \\\ell_i\\ such that: \$\$ \ell_i
= p_i f_i + (1 - p_i) g_i, \$\$ with \\g_i\\ the known component
density, and where \\\ell_i\\ can be estimated consistently thanks to
the observations. The test is made on the \\f_i\\'s, can be performed
using two methods: either the comparison of coefficients obtained
through polynomial basis expansions of the component densities, or by
the inner-convergence property obtained using the IBM approach. See
'Details' below for further information.

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

  A list of objects of class [admix_model](admix_model.md), with
  information about known distributions and known parameters.

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
hereafter. When choosing the 'icv' testing method, it is recommended to
use parallel computing.

## References

Milhaud X, Pommeret D, Salhi Y, Vandekerkhove P (2024).
“Contamination-source based K-sample clustering.” *Journal of Machine
Learning Research*, **25**(287), 1–32.
<https://jmlr.org/papers/v25/23-0914.html>. Milhaud X, Pommeret D, Salhi
Y, Vandekerkhove P (2022). “Semiparametric two-sample admixture
components comparison test: The symmetric case.” *Journal of Statistical
Planning and Inference*, **216**, 135-150. ISSN 0378-3758.
[doi:10.1016/j.jspi.2021.05.010](https://doi.org/10.1016/j.jspi.2021.05.010)
. Pommeret D, Vandekerkhove P (2019). “Semiparametric density testing in
the contamination model.” *Electronic Journal of Statistics*, 4743–4793.
[doi:10.1214/19-EJS1650](https://doi.org/10.1214/19-EJS1650) .

## See also

[`gaussianity_test()`](gaussianity_test.md),
[`orthobasis_test()`](orthobasis_test.md),
[`IBM_k_samples_test()`](IBM_k_samples_test.md),
[`get_known_component()`](get_known_component.md),
[`get_mixing_weights()`](get_mixing_weights.md),
[`reject_nullHyp()`](reject_nullHyp.md), [`which_rank()`](which_rank.md)

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
#> T = 0.52976, expansion order S = 1, p-value = 0.4667
#> alternative hypothesis: Distributions of unknown components involved 
#>                         in the contamination models are different
#> 
```
