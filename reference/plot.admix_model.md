# Plot method for objects of class `admix_model`

Plots the probability density function of the known component of the
admixture model, where we recall that an admixture model has probability
density function (pdf) l_i such that: l_i = p_i \* f_i + (1-p_i) \* g_i,
with g_i the known component density. The unknown quantities are
therefore p_i and f_i.

## Usage

``` r
# S3 method for class 'admix_model'
plot(x, ...)
```

## Arguments

- x:

  An object of class `admix_model`.

- ...:

  A list of additional parameters belonging to the default method.

## Author

Xavier Milhaud <xavier.milhaud.research@gmail.com>

## Examples

``` r
plot(admix_model(knownComp_dist = "norm", knownComp_param = list("mean"=0, "sd"=1)))

plot(admix_model(knownComp_dist = "pois", knownComp_param = list("lambda"=1.5)))

```
