# Equality of known components in two admixture models

Test if the known component distributions coming from two admixture
models are identical.

## Usage

``` r
is_equal_knownComp(admixMod1, admixMod2)
```

## Arguments

- admixMod1:

  An object of class `admix_model` related to the first admixture model.

- admixMod2:

  An object of class `admix_model` related to the second admixture
  model.

## Value

A boolean (TRUE if the known components are the same, otherwise FALSE).

## Author

Xavier Milhaud <xavier.milhaud.research@gmail.com>

## Examples

``` r
admixMod1 <- admix_model(knownComp_dist = "norm",
                         knownComp_param = list("mean"=0, "sd"=1))
admixMod2 <- admix_model(knownComp_dist = "norm",
                         knownComp_param = list("mean"=0, "sd"=1))
is_equal_knownComp(admixMod1, admixMod2)
#> [1] TRUE

admixMod1 <- admix_model(knownComp_dist = "multinom",
                         knownComp_param = list("size"=1, "prob"=c(0.2,0.5,0.3)))
admixMod2 <- admix_model(knownComp_dist = "multinom",
                         knownComp_param = list("size"=1, "prob"=c(0.2,0.4,0.4)))
is_equal_knownComp(admixMod1, admixMod2)
#> [1] FALSE
```
