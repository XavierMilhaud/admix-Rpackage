# Extractor for members of clusters

Extract the clusters that were discovered among K samples, where
belonging to one given cluster means having equal unknown component
distributions.

## Usage

``` r
get_cluster_members(x)
```

## Arguments

- x:

  An object of class `admix_cluster`.

## Value

The samples included in each detected cluster.

## Author

Xavier Milhaud <xavier.milhaud.research@gmail.com>

## Examples

``` r
# \donttest{
## Simulate mixture data:
mixt1 <- twoComp_mixt(n = 1600, weight = 0.8,
                      comp.dist = list("gamma", "exp"),
                      comp.param = list(list("shape" = 16, "scale" = 1/4),
                                        list("rate" = 1/3.5)))
mixt2 <- twoComp_mixt(n = 2000, weight = 0.7,
                      comp.dist = list("gamma", "exp"),
                      comp.param = list(list("shape" = 14, "scale" = 1/2),
                                        list("rate" = 1/5)))
mixt3 <- twoComp_mixt(n = 2500, weight = 0.6,
                      comp.dist = list("gamma", "gamma"),
                      comp.param = list(list("shape" = 16, "scale" = 1/4),
                                        list("shape" = 12, "scale" = 1/2)))
mixt4 <- twoComp_mixt(n = 3800, weight = 0.5,
                      comp.dist = list("gamma", "exp"),
                      comp.param = list(list("shape" = 14, "scale" = 1/2),
                                        list("rate" = 1/7)))
data1 <- get_mixture_data(mixt1) ; data2 <- get_mixture_data(mixt2)
data3 <- get_mixture_data(mixt3) ; data4 <- get_mixture_data(mixt4)
## Define the admixture models:
admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
                         knownComp_param = mixt1$comp.param[[2]])
admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
                         knownComp_param = mixt2$comp.param[[2]])
admixMod3 <- admix_model(knownComp_dist = mixt3$comp.dist[[2]],
                         knownComp_param = mixt3$comp.param[[2]])
admixMod4 <- admix_model(knownComp_dist = mixt4$comp.dist[[2]],
                         knownComp_param = mixt4$comp.param[[2]])
## Clustering procedure:
x <- admix_cluster(samples = list(data1, data2, data3, data4),
              admixMod = list(admixMod1, admixMod2, admixMod3, admixMod4),
              conf_level = 0.95, tune_penalty = TRUE, n_sim_tab = 10)
#>   |                                                          |                                                  |   0%  |                                                          |======                                            |  12%  |                                                          |======================================            |  75%  |                                                          |==================================================| 100%
get_cluster_members(x)
#>            [,1] [,2] [,3] [,4]
#> Id_sample     1    2    3    4
#> Id_cluster    1    2    1    2
get_cluster_sizes(x)
#> [1] 2 2
# }
```
