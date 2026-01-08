# Define the distribution/parameter(s) of the known component

Create an object of class `admix_model`, containing the information
about the known component distribution in the admixture model. An
admixture (aka contamination) model is a two-component mixture model
with one known component. Both the second component distribution and the
mixing weight are unknown.

## Usage

``` r
admix_model(knownComp_dist, knownComp_param)
```

## Arguments

- knownComp_dist:

  (Character) The name of the distribution (specified as in R glossary)
  of the known component of the admixture model

- knownComp_param:

  (Character) A vector of the names of the parameters (specified as in R
  glossary) involved in the chosen known distribution, with their
  values.

## Value

An object of class admix_model, containing 2 attributes: 1) a list that
gives the information about the distributions involved in the
two-component mixture model (the unknown and the known ones); 2) a list
that gives the information about the corresponding parameters of those
distributions.

## Author

Xavier Milhaud <xavier.milhaud.research@gmail.com>

## Examples

``` r
admix_model(knownComp_dist = "norm", knownComp_param = list("mean"=0, "sd"=1))
#> Call:admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, 
#>     sd = 1))
#> 
#> Known component distribution:  norm 
#> Known component parameters: mean=0 sd=1
#> 
admix_model(knownComp_dist = "exp", knownComp_param = list("rate"=2))
#> Call:admix_model(knownComp_dist = "exp", knownComp_param = list(rate = 2))
#> 
#> Known component distribution:  exp 
#> Known component parameters: rate=2
#> 
admix_model(knownComp_dist = "pois", knownComp_param = list("lambda"=5))
#> Call:admix_model(knownComp_dist = "pois", knownComp_param = list(lambda = 5))
#> 
#> Known component distribution:  pois 
#> Known component parameters: lambda=5
#> 
admix_model(knownComp_dist = "multinom", knownComp_param = list("size"=1, "prob"=c(0.2,0.8,0.1)))
#> Call:admix_model(knownComp_dist = "multinom", knownComp_param = list(size = 1, 
#>     prob = c(0.2, 0.8, 0.1)))
#> 
#> Known component distribution:  multinom 
#> Known component parameters: size=1 prob=c(0.2, 0.8, 0.1)
#> 
```
