# Check the validity of the specified distributions

Check the validity of the specified distributions

## Usage

``` r
get_distribution_parameters(dist)
```

## Arguments

- dist:

  A single string naming the distribution under consideration.

## Value

A vector composed of the names of the parameters of the distribution.

## Author

Xavier Milhaud <xavier.milhaud.research@gmail.com>

## Examples

``` r
get_distribution_parameters("norm")
#> [1] "mean" "sd"  
get_distribution_parameters("gamma")
#> [1] "shape" "rate"  "scale"
get_distribution_parameters("weibull")
#> [1] "shape" "scale"
```
