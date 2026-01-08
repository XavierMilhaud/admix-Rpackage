# Extractor for pairwise discrepancy rankings

Provide the matrix storing the ranks of discrepancies using
Inversion-Best Matching approach between all couples among the K (K\>2)
samples under study.

## Usage

``` r
get_discrepancy_rank(x)
```

## Arguments

- x:

  An object of class `IBM_test`.

## Value

A matrix of ranks, from the closest couple (rank 1) in terms of
discrepancy measure to the most different one.

## Author

Xavier Milhaud <xavier.milhaud.research@gmail.com>
