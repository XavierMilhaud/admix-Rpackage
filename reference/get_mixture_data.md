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
#>  [1]  0.14952991 -2.14804235 -1.05260879 -0.48597053  3.30978955  3.60907338
#>  [7] -0.24812569  3.28337064  3.34080056 -0.07020902 -0.17769259  2.79419545
#> [13]  3.90009150 -1.29015095  2.91253514  0.14652847  0.67151584  2.53934891
#> [19] -0.23538262  2.45098311
```
