# Extractor for tabulated distribution in the k-sample test

Provide (the list of) tabulated distribution(s) that allow to define the
extreme quantile(s) against which the test statistic(s) is compared.

## Usage

``` r
get_tabulated_dist(x)
```

## Arguments

- x:

  An object of class `IBM_test` or `admix_cluster`.

## Value

A numeric vector containing the simulated values of the tabulated
distribution, sorted in increasing order.

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
                conf_level = 0.95, test_method = "icv", n_sim_tab = 10)
get_tabulated_dist(x)
#> logical(0)
```
