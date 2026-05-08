# Plot method for objects of class `admix_model`

Plots the probability density function of the known component of the
admixture model.

## Usage

``` r
# S3 method for class 'admix_model'
plot(x, n = 1000, main = "Known component distribution", ...)
```

## Arguments

- x:

  An object of class `admix_model`.

- n:

  The number of 'x' values to consider for plotting the pdf in the
  continuous case.

- main:

  The title of the plot.

- ...:

  A list of additional parameters belonging to the default method.

## Author

Xavier Milhaud <xavier.milhaud.research@gmail.com>

## Examples

``` r
plot(admix_model(knownComp_dist = "norm", knownComp_param = list("mean"=0, "sd"=1)))

plot(admix_model(knownComp_dist = "pois", knownComp_param = list("lambda"=1.5)))

```
