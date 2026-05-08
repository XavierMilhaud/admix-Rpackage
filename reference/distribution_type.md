# Determine the type of distribution under consideration

Determine the type of distribution under consideration

## Usage

``` r
distribution_type(dist)
```

## Arguments

- dist:

  A single string naming the distribution under consideration.

## Value

The type of distribution under study (continuous, discrete, or
multivariate).

## Author

Xavier Milhaud <xavier.milhaud.research@gmail.com>

## Examples

``` r
distribution_type("norm")
#> [1] "Continuous"
distribution_type("pois")
#> [1] "Discrete"
distribution_type("multinom")
#> [1] "Multivariate"
distribution_type("weibull")
#> [1] "Continuous"
```
