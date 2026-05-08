# Extractor for simulated data from two-component mixture

Get the mixture data generated from method
[`twoComp_mixt()`](twoComp_mixt.md).

## Usage

``` r
get_mixture_data(x)
```

## Arguments

- x:

  An object of class `twoComp_mixt`.

## Value

A numeric vector of the simulated data.

## Author

Xavier Milhaud <xavier.milhaud.research@gmail.com>

## Examples

``` r
sim.X <- twoComp_mixt(n = 20, weight = 0.5,
                      comp.dist = list("norm", "norm"),
                      comp.param = list(list("mean"=3, "sd"=0.5),
                                        list("mean"=0, "sd"=1)))
get_mixture_data(sim.X)
#>  [1]  2.52054402  3.39532943  2.66234220  3.39846378 -0.31723667  2.90550389
#>  [7]  3.46344163  3.77090799  0.13322040  0.57573498  1.61661603  0.27934333
#> [13]  3.31501573  4.29501914 -0.01179492  0.83966119  1.33855061  3.30550843
#> [19]  2.29059740  2.96786824
```
