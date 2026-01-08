# Extractor for estimated mixing weights

Extracts the estimated mixing weights from fitted objects of class
`admix_estim`, `gaussianity_test` and `orthobasis_test`.

## Usage

``` r
get_mixing_weights(x)
```

## Arguments

- x:

  An object of class `admix_estim`, `gaussianity_test` or
  `orthobasis_test`.

## Value

A numeric vector of estimated mixing weight(s).

## Details

This is a generic extractor. The exact behavior depends on the class of
the input object:

- `admix_estim`: returns the estimated mixture proportions.

- `gaussianity_test`, `orthobasis_test`: returns weights derived from
  hypothesis testing results.

## Author

Xavier Milhaud <xavier.milhaud.research@gmail.com>

## Examples

``` r
## Simulate a two-component Gaussian mixture:
mixt1 <- twoComp_mixt(n = 380, weight = 0.7,
                      comp.dist = list("norm", "norm"),
                      comp.param = list(list("mean" = -2, "sd" = 0.5),
                                        list("mean" = 0, "sd" = 1)))
data1 <- get_mixture_data(mixt1)
admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
                         knownComp_param = mixt1$comp.param[[2]])
## Estimate the unknown quantities:
x <- admix_estim(samples = list(data1), admixMod = list(admixMod1), est_method = "BVdk")
#> Mixing weight estimation using 'BVdk' assumes the unknown component
#> distribution to have a symmetric probability density function.
## Extract the information about the known component:
get_mixing_weights(x)
#> [1] 0.7146939
```
