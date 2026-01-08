# Extractor for known component(s) in admixture model(s)

Get the known component of the admixture model considered for
estimation, test, or clustering.

## Usage

``` r
get_known_component(x)
```

## Arguments

- x:

  An object of class `admix_estim`, `gaussianity_test`,
  `orthobasis_test`, `IBM_test`, or `admix_cluster`.

## Value

A list providing information on the known component (distribution,
parameters).

## Details

This is a generic extractor, providing with the same information
whatever the object class.

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
get_known_component(x)
#> [[1]]
#> Call:admix_model(knownComp_dist = mixt1$comp.dist[[2]], knownComp_param = mixt1$comp.param[[2]])
#> 
#> Known component distribution:  norm 
#> Known component parameters: mean=0 sd=1
#> 
#> 
```
