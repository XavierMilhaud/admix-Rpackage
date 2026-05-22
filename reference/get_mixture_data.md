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
#>  [1]  3.8141905  0.6450361  2.9824421  3.2862745  2.7484593  2.6916632
#>  [7]  2.9545718 -0.2731024  0.9818921  2.1035109 -0.5616319  0.4850168
#> [13]  1.1777993  1.0441347  2.3441173  0.9594856  2.6542894  3.1942147
#> [19] -0.5619008  2.8197637
```
