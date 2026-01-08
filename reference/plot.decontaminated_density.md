# Plot method for object of class `decontaminated_density`

Plot the decontaminated density of the unknown component from some
admixture model, after inversion of the admixture cumulative
distribution functions.

## Usage

``` r
# S3 method for class 'decontaminated_density'
plot(x, x_val = NULL, add_plot = FALSE, offset = 0, bar_width = 0.3, ...)
```

## Arguments

- x:

  An object of class `decontaminated_density` (see
  ?decontaminated_density).

- x_val:

  Values at which to evaluate the decontaminated density.

- add_plot:

  Boolean, TRUE when a new plot is added to the existing one.

- offset:

  Numeric. Position of the bars relative to the labels on the x-axis.

- bar_width:

  Width of bars to be plotted.

- ...:

  Arguments to be passed to generic method `plot`, such as graphical
  parameters (see ?par).

## Value

The plot of the decontaminated density if one sample is provided, or the
comparison of decontaminated densities plotted on the same graph in the
case of multiple samples.

## Details

The decontaminated density is obtained by inverting the admixture
density, given by l = p\*f + (1-p)\*g, to isolate the unknown component
f after having estimated p and l.

## Author

Xavier Milhaud <xavier.milhaud.research@gmail.com>
