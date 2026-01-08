# Extractor for the selected rank in the test statistic

Provide the selected rank of the test statistic (connected to the
expansion order of the densities in the orthonormal polynomial basis if
method 'poly' was chosen; or to the number of terms, i.e. discrepancies
between couples of samples, included in the test statistic with method
'icv').

## Usage

``` r
which_rank(x)
```

## Arguments

- x:

  An object of class `gaussianity_test`, `orthobasis_test` or
  `IBM_test`.

## Value

An integer corresponding to the selected rank in the test statistics,
i.e. how many terms were kept in the test statistic.

## Author

Xavier Milhaud <xavier.milhaud.research@gmail.com>

## Examples

``` r
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
x <- admix_test(samples = list(data1,data2), admixMod = list(admixMod1,admixMod2),
                conf_level = 0.95, test_method = "poly", ask_poly_param = FALSE, support = "Real")
#>   Default estimation method is 'BVdk' when testing with polynomial basis expansions (ensuring
#>   theoretical guarantees, but relying on symmetric unknown component densities). To consider
#>   other frameworks in the 2-sample case, use 'PS' estimator (setting argument 'est_method' to 'PS').
which_rank(x)
#> [1] 1
```
