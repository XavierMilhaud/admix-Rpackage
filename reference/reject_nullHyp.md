# Extractor for the test decision

Provide the decision of the statistical test: reject or do not reject
the null hypothesis.

## Usage

``` r
reject_nullHyp(x)
```

## Arguments

- x:

  An object of class `gaussianity_test`, `orthobasis_test` or
  `IBM_test`.

## Value

A boolean giving the result of the test, TRUE if the null hypothesis is
rejected, otherwise FALSE.

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
reject_nullHyp(x)
#> [1] FALSE
```
