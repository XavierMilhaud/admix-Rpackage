# Estimation of the admixture parameters by Bordes & Vandekerkhove (2010)

Estimates parameters in an admixture model where the unknown component
is assumed to have a symmetric density. More precisely, estimates the
two parameters (mixture weight and location shift) in the admixture
model with pdf: l(x) = p\*f(x-mu) + (1-p)\*g(x), x in R, where g is the
known component, p is the proportion and f is the unknown component with
symmetric density. The localization shift parameter is denoted mu, and
the component weight p. See the reference below for further details.

## Usage

``` r
estim_BVdk(
  samples,
  admixMod,
  method = c("L-BFGS-B", "Nelder-Mead"),
  compute_var = FALSE
)
```

## Arguments

- samples:

  The observed sample under study.

- admixMod:

  An object of class [admix_model](admix_model.md), containing useful
  information about distributions and parameters.

- method:

  The method used throughout the optimization process, either 'L-BFGS-B'
  or 'Nelder-Mead' (see ?optim).

- compute_var:

  (default to FALSE) A boolean that indicates whether one computes the
  variance of the estimators of unknown mixing proportions and location
  shift parameter.

## Value

An object of class estim_BVdk, containing 8 attributes: 1) the number of
sample under study (set to 1 here); 2) the sample size; 3) the
information about mixture components (distributions and parameters); 4)
the estimation method (Bordes and Vandekerkhove here, see the given
reference); 5) the estimated mixing proportion (weight of the unknown
component distribution); 6) the estimated location parameter of the
unknown component distribution (with symetric density); 7) the variance
of the two estimators (respectively the mixing proportion and location
shift); 8) the optimization method that was used.

## References

Bordes L, Vandekerkhove P (2010). “Semiparametric two-component mixture
model with a known component: An asymptotically normal estimator.”
*Mathematical Methods of Statistics*, **19**(1), 22–41.
[doi:10.3103/S1066530710010023](https://doi.org/10.3103/S1066530710010023)
.

## See also

[`print.estim_BVdk()`](print.estim_BVdk.md) for printing a short version
of the results from this estimation method, and
[`summary.estim_BVdk()`](summary.estim_BVdk.md) for more comprehensive
results.

## Author

Xavier Milhaud <xavier.milhaud.research@gmail.com>

## Examples

``` r
if (FALSE) { # \dontrun{
## Simulate mixture data:
mixt1 <- twoComp_mixt(n = 200, weight = 0.4,
                      comp.dist = list("norm", "norm"),
                      comp.param = list(list("mean" = -2, "sd" = 0.5),
                                        list("mean" = 0, "sd" = 1)))
## Retrieves the mixture data:
data1 <- get_mixture_data(mixt1)
## Define the admixture model:
admixMod <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
                        knownComp_param = mixt1$comp.param[[2]])
## Perform the estimation of parameters in real-life:
ex <- estim_BVdk(samples = data1, admixMod = admixMod, method = 'L-BFGS-B')
print.estim_BVdk(ex)

## Second example:
mixt2 <- twoComp_mixt(n = 200, weight = 0.65,
                      comp.dist = list("norm", "exp"),
                      comp.param = list(list("mean" = -1, "sd" = 0.5),
                                        list("rate" = 1)))
data2 <- get_mixture_data(mixt2)
admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
                        knownComp_param = mixt2$comp.param[[2]])
## Perform the estimation of parameters in real-life:
ex <- estim_BVdk(samples = data2, admixMod = admixMod2, method = 'L-BFGS-B',
                 compute_var = TRUE)
print.estim_BVdk(ex)
} # }
```
