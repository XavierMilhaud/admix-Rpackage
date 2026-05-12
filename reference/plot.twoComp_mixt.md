# Plot the empirical mixture pdf

Plots the empirical densities of the samples provided, with optional
arguments to improve the visualization.

## Usage

``` r
# S3 method for class 'twoComp_mixt'
plot(
  x,
  add_plot = FALSE,
  offset = 0,
  bar_width = 0.2,
  main = "Mixture distribution (density or mass function)",
  ...
)
```

## Arguments

- x:

  Object of class `twoComp_mixt` from which the density will be plotted.

- add_plot:

  (default to FALSE) Option to plot another mixture distribution on the
  same graph.

- offset:

  Numeric. Position of the bars relative to the labels on the x-axis.

- bar_width:

  Width of bars to be plotted.

- main:

  The title of the plot.

- ...:

  further classical arguments and graphical parameters for methods plot
  and hist.

## Value

A plot with the densities of the samples provided as inputs.

## Author

Xavier Milhaud <xavier.milhaud.research@gmail.com>

## Examples

``` r
## Mixture of continuous random variables:
sim.X <- twoComp_mixt(n = 2000, weight = 0.5,
                      comp.dist = list("norm", "norm"),
                      comp.param = list(list("mean"=3, "sd"=0.5),
                                        list("mean"=0, "sd"=1)))
sim.Y <- twoComp_mixt(n = 1200, weight = 0.7,
                      comp.dist = list("norm", "exp"),
                      comp.param = list(list("mean"=-3, "sd"=0.5),
                                        list("rate"=1)))
plot(sim.X, xlim=c(-5,5), ylim=c(0,0.5))
plot(sim.Y, add_plot = TRUE, xlim=c(-5,5), ylim=c(0,0.5), col = "red")
legend("topright", legend = c("sim.X","sim.Y"), col = c("black","red"),
       lty = rep(1,2), bty = "n")


## Mixture of discrete random variables:
sim.X <- twoComp_mixt(n = 2000, weight = 0.5,
                      comp.dist = list("multinom", "multinom"),
                      comp.param = list(list("size"=1, "prob"=c(0.3,0.4,0.3)),
                                        list("size"=1, "prob"=c(0.1,0.2,0.7))))
sim.Y <- twoComp_mixt(n = 1800, weight = 0.7,
                      comp.dist = list("multinom", "multinom"),
                      comp.param = list(list("size"=1, "prob"=c(0.3,0.4,0.3)),
                                        list("size"=1, "prob"=c(0.6,0.2,0.2))))
sim.Z <- twoComp_mixt(n = 1800, weight = 0.3,
                      comp.dist = list("multinom", "multinom"),
                      comp.param = list(list("size"=1, "prob"=c(0.2,0.1,0.7)),
                                        list("size"=1, "prob"=c(1/3,1/3,1/3))))
plot(sim.X, offset = -0.05, bar_width = 0.05, col = "steelblue")
plot(sim.Y, add_plot = TRUE, offset = 0, bar_width = 0.05, col = "orange")
plot(sim.Z, add_plot = TRUE, offset = +0.05, bar_width = 0.05, col = "red")
legend("topleft", legend = c("sim.X","sim.Y","sim.Z"), col = c("steelblue","orange","red"),
       lty = rep(1,3), bty = "n")

```
