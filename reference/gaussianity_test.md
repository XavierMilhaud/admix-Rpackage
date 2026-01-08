# Gaussianity test in an admixture model

Performs an hypothesis test to check for the gaussianity of the unknown
mixture component. Recall that an admixture model has probability
density function (pdf) l = p\*f + (1-p)\*g, where g is the known pdf and
l is observed (others are unknown). This test requires optimization (to
estimate the unknown parameters) as defined by Bordes & Vandekerkhove
(2010), which means that the unknown mixture component must have a
symmetric density.

## Usage

``` r
gaussianity_test(
  sample,
  admixMod,
  conf_level = 0.95,
  ask_poly_param = FALSE,
  K = 3,
  s = 0.25,
  support = c("Real", "Integer", "Positive", "Bounded.continuous"),
  ...
)
```

## Arguments

- sample:

  (numeric) The sample under study.

- admixMod:

  An object of class [admix_model](admix_model.md), containing useful
  information about distributions and parameters.

- conf_level:

  (default to 0.95) The confidence level. Equals 1-alpha, where alpha is
  the level of the test (type-I error).

- ask_poly_param:

  (default to FALSE) If TRUE, ask the user to choose both the order 'K'
  of expansion coefficients in the orthonormal polynomial basis, and the
  penalization rate 's' involved on the penalization rule for the test.

- K:

  (default to 3) If not asked (see the previous argument), number of
  coefficients considered for the polynomial basis expansion.

- s:

  (in \]0,1/2\[, default to 0.25) If not asked (see the previous
  argument), normalization rate involved in the penalization rule for
  model selection. See the reference below.

- support:

  Support of the probability distributions, useful to choose the
  appropriate polynomial orthonormal basis. One of 'Real', 'Integer',
  'Positive', or 'Bounded.continuous'.

- ...:

  Optional arguments to [estim_BVdk](estim_BVdk.md).

## Value

An object of class gaussianity_test, inherited class from "htest".
Contains attributes like 1) the number of populations under study (1 in
this case); 2) the sample size; 3) the information about the known
component distribution; 4) the reject decision of the test; 5) the
confidence level of the test, 6) the p-value of the test; 7) the value
of the test statistic; 8) the variance of the test statistic at each
order in the polynomial orthobasis expansion; 9) the selected rank
(order) for the test statistic; 10) a list of estimates (mixing weight,
mean and standard deviation of the Gaussian unknown distribution).

## References

Pommeret D, Vandekerkhove P (2019). “Semiparametric density testing in
the contamination model.” *Electronic Journal of Statistics*, 4743–4793.
[doi:10.1214/19-EJS1650](https://doi.org/10.1214/19-EJS1650) .

## Author

Xavier Milhaud <xavier.milhaud.research@gmail.com>

## Examples

``` r
if (FALSE) { # \dontrun{
####### Under the null hypothesis H0.
## Simulate mixture data:
mixt1 <- twoComp_mixt(n = 250, weight = 0.4,
                      comp.dist = list("norm", "exp"),
                      comp.param = list(list("mean" = -2, "sd" = 0.5),
                                        list("rate" = 1)))
data1 <- get_mixture_data(mixt1)
## Define the admixture models:
admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
                         knownComp_param = mixt1$comp.param[[2]])
## Performs the test:
gaussianity_test(sample = data1, admixMod = admixMod1,
                 conf_level = 0.95, K = 3, s = 0.1, support = "Real")
} # }
```
