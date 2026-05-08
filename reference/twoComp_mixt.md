# Simulation of a two-component mixture model

Simulate a two-component mixture model following the probability density
function (pdf) \\\ell\\ such that \$\$ \ell = p f + (1 - p) g, \$\$ with
\\f\\ and \\g\\ the mixture component distributions, and \\p\\ the
mixing weight.

## Usage

``` r
twoComp_mixt(
  n = 1000,
  weight = 0.5,
  comp.dist = list("norm", "norm"),
  comp.param = list(list(mean = 0, sd = 1), list(mean = 2, sd = 1))
)
```

## Arguments

- n:

  Number of observations to be simulated.

- weight:

  Weight of the first component distribution (distribution f) in the
  mixture.

- comp.dist:

  A list of two elements corresponding to the component distributions
  (with available names listed in object 'Distribution.df' in package
  EnvStats) involved in the mixture model. These elements respectively
  refer to the two component distributions f and g. By convention, in
  the framework of admixture models where one of the two components is
  unknown, the first element of the list corresponds to the 'unknown'
  component distribution, whereas the second one refers to the known
  one.

- comp.param:

  A list of two elements corresponding to the parameters of the
  component distributions, each element being a list itself. The names
  used in each list must correspond to the available parameters listed
  in object 'Distribution.df' in package EnvStats. These elements
  respectively refer to the parameters of f and g distributions of the
  mixture model. By convention, in the framework of admixture models
  where one of the two components is unknown, the first element of the
  list corresponds to the 'unknown' component parameters, whereas the
  second one refers to the known ones.

## Value

An object of class twoComp_mixt, containing eight attributes: 1) the
number of simulated observations, 2) the simulated mixture data, 3) the
support of the distributions, 4) the name of the component
distributions, 5) the name of the parameters of the component
distributions and their values, 6) the mixing proportion, 7) the
observations coming from the first component, 8) the observations coming
from the second component.

## See also

[`get_mixture_data()`](get_mixture_data.md) to access the simulated
mixture data.

## Author

Xavier Milhaud <xavier.milhaud.research@gmail.com>

## Examples

``` r
## Mixture of continuous random variables:
sim.X <- twoComp_mixt(n = 1200, weight = 0.7,
                      comp.dist = list("norm", "exp"),
                      comp.param = list(list("mean"=-3, "sd"=0.5),
                                        list("rate"=1)))
print(sim.X)
#> 
#> Call:twoComp_mixt(n = 1200, weight = 0.7, comp.dist = list("norm", 
#>     "exp"), comp.param = list(list(mean = -3, sd = 0.5), list(rate = 1)))
#> 
#> Number of observations: 1200 
#> 
#> Simulated data (first 5 obs.): 
#>  0.5504287 -2.349743 -3.420734 -3.350772 -4.005017
#> Simulated observations coming from the 1st component (first 5 obs.): 
#>  -2.349743 -3.420734 -3.350772 -4.005017 -3.701633
#> Simulated observations coming from the 2nd component (first 5 obs.): 
#>  0.5504287 1.698026 0.265485 2.098534 1.058672
#> 
data.X <- get_mixture_data(sim.X)
plot(density(data.X))


## Mixture of discrete random variables:
sim.Y <- twoComp_mixt(n = 1800, weight = 0.7,
                      comp.dist = list("multinom", "multinom"),
                      comp.param = list(list("size"=1, "prob"=c(0.3,0.4,0.3)),
                                        list("size"=1, "prob"=c(0.6,0.2,0.2))))
print(sim.Y)
#> 
#> Call:twoComp_mixt(n = 1800, weight = 0.7, comp.dist = list("multinom", 
#>     "multinom"), comp.param = list(list(size = 1, prob = c(0.3, 
#>     0.4, 0.3)), list(size = 1, prob = c(0.6, 0.2, 0.2))))
#> 
#> Number of observations: 1800 
#> 
#> Obtained multinomial mixture distribution: 
#>  663 612 525 
#> 
#> 
```
