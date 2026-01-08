# Plot the empirical mixture pdf

Plots the empirical densities of the samples provided, with optional
arguments to improve the visualization.

## Usage

``` r
# S3 method for class 'twoComp_mixt'
plot(x, add_plot = FALSE, offset = 0, bar_width = 0.2, ...)
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

- ...:

  further classical arguments and graphical parameters for methods plot
  and hist.

## Value

A plot with the densities of the samples provided as inputs.

## Author

Xavier Milhaud <xavier.milhaud.research@gmail.com>
