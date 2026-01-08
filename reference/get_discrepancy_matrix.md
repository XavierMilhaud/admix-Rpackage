# Extractor for discrepancies b/w unknown components

Provide the matrix storing pairwise discrepancies b/w unknown components
in admixture models, using Inversion-Best Matching approach.

## Usage

``` r
get_discrepancy_matrix(x)
```

## Arguments

- x:

  An object of class `IBM_test` or `admix_cluster`.

## Value

A matrix of pairwise discrepancies among the K (K\>2) samples under
study.

## Author

Xavier Milhaud <xavier.milhaud.research@gmail.com>
