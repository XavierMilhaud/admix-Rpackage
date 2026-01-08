# Print method for objects of class `admix_model`

Print an object of class 'admix_mod'. An admixture model has probability
density function (pdf) l_i such that: l_i = p_i \* f_i + (1-p_i) \* g_i,
with g_i the known component density. The unknown quantities are
therefore p_i and f_i.

## Usage

``` r
# S3 method for class 'admix_model'
print(x, ...)
```

## Arguments

- x:

  An object of class `admix_model`.

- ...:

  A list of additional parameters belonging to the default method.

## Author

Xavier Milhaud <xavier.milhaud.research@gmail.com>
