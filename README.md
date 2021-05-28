
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Welcome to R package admix

<!-- badges: start -->

[![R-CMD-check](https://github.com/XavierMilhaud/admix/workflows/R-CMD-check/badge.svg)](https://github.com/XavierMilhaud/admix/actions)

The goal of admix is to provide code for estimation, hypothesis testing
and clustering methods in admixture models.

We remind that an admixture model has the following cumulative
distribution function (cdf)
*L*(*x*) = *p**F*(*x*) + (1−*p*)*G*(*x*),   *x* ∈ ℝ,

where *G* is a perfectly known cdf, and *p* and *F* are unknown.

The cdf *F* relates to the contamination phenomenon that is added to the
well-known signal *G*, with proportion *p*.

The proportion of the unknown component in the two-component mixture
model can be easily estimated under weak nonparametric assumptions on
the related distribution. The decontaminated version of this unknown
component distribution can then be tested against some other specified
distribution (included another decontaminated unknown component).
Finally, clustering of *K* populations is made possible, based on
hypothesis tests that compare unknown component distributions. The
package is suited to one-sample as well as multi-samples analysis.

# Installation

<!-- You can install the released version of admix from [CRAN](https://CRAN.R-project.org) with: -->

You can install the released version of admix from
[Github](https://github.com/XavierMilhaud/admix) with:

``` r
#once on CRAN with : install.package("admix")
# from now on:
remotes::install_github(repo = "XavierMilhaud/admix@main", build_manual = TRUE, build_vignettes = FALSE)
```

The optional argument build\_vignettes can be set to TRUE to get
vignettes that help to understand the functionnalities of the package.

To get some help about the functionalities of the package, do once
installed:

``` r
help(package = 'admix')
```

More details can also be found through the vignettes, available in admix
github-pages (see <https://xaviermilhaud.github.io/admix/>, in Menu
Articles).

## Example

This is a basic example which shows you how to estimate the unknown
component proportion and the localization shift parameters in an
admixture model where the unknown component density is assumed to be
symmetric. In practice, the cdf *L* is given by
*L*(*x*) = *p**F*(*x*−*μ*) + (1−*p*)*G*(*x*),   *x* ∈ ℝ,
where *p* is the unknown component weight, and *μ* is the localization
shift parameter of the unknown cdf *F* with symmetric density.

The estimation would be made through the following commands:

``` r
library(admix)
## Simulate data:
list.comp <- list(f = 'norm', g = 'norm')
list.param <- list(f = list(mean = 3, sd = 0.5),
                   g = list(mean = 0, sd = 1))
data1 <- rsimmix(n = 1000, unknownComp_weight = 0.8, list.comp, list.param)[['mixt.data']]
## Perform the estimation of parameters in real-life:
list.comp <- list(f = NULL, g = 'norm')
list.param <- list(f = NULL, g = list(mean = 0, sd = 1))
BVdk_estimParam(data1, method = 'L-BFGS-B', list.comp, list.param)
#> [1] 0.7977239 3.0114174
```

<!-- badges: end -->
