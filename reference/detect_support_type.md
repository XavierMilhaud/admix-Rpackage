# Detect the type of support of some random variables

Given one or two sets of observations (samples), the function provides
with the most plausible type of support for the underlying random
variables to be studied. If less than 3 percents of the observations
have different values, we consider that the support is discrete.
Otherwise, we consider it as a continuous support.

## Usage

``` r
detect_support_type(sample1, sample2 = NULL)
```

## Arguments

- sample1:

  The first sample of observations under study.

- sample2:

  The second sample of observations under study.

## Value

The type of support, either discrete or continuous.

## Author

Xavier Milhaud <xavier.milhaud.research@gmail.com>

## Examples

``` r
## Simulate mixture data:
mixt1 <- twoComp_mixt(n = 1500, weight = 0.5,
                      comp.dist = list("norm", "norm"),
                      comp.param = list(list("mean" = 3, "sd" = 0.5),
                                        list("mean" = 0, "sd" = 1)))
data1 <- get_mixture_data(mixt1)
mixt2 <- twoComp_mixt(n = 2000, weight = 0.7,
                      comp.dist = list("norm", "norm"),
                      comp.param = list(list("mean" = 3, "sd" = 0.5),
                                        list("mean" = 5, "sd" = 2)))
data2 <- get_mixture_data(mixt2)
## Test the type of support:
detect_support_type(data1, data2)
#> [1] "Continuous"
```
