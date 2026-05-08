# Check the validity of the specified distribution and parameter(s)

Check the validity of the specified distribution and parameter(s)

## Usage

``` r
validate_distribution(dist, params)
```

## Arguments

- dist:

  A single string naming the distribution under consideration.

- params:

  A character vector composed of the names of the parameters for a given
  distribution.

## Value

A vector composed of the names of the parameters of the distribution.

## Author

Xavier Milhaud <xavier.milhaud.research@gmail.com>

## Examples

``` r
validate_distribution(dist = "norm", params = c("mean" = 0, "sd" = 1))
```
