
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Welcome to R package admix

<!-- badges: start -->

[![R-CMD-check](https://github.com/XavierMilhaud/admix-Rpackage/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/XavierMilhaud/admix-Rpackage/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/XavierMilhaud/admix-Rpackage/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/XavierMilhaud/admix-Rpackage/actions/workflows/test-coverage.yaml)
[![pkgdown](https://github.com/XavierMilhaud/admix-Rpackage/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/XavierMilhaud/admix-Rpackage/actions/workflows/pkgdown.yaml)
[![pages-build-deployment](https://github.com/XavierMilhaud/admix-Rpackage/actions/workflows/pages/pages-build-deployment/badge.svg)](https://github.com/XavierMilhaud/admix-Rpackage/actions/workflows/pages/pages-build-deployment)

The goal of admix is to provide code for estimation, hypothesis testing
and clustering methods in admixture models.

We remind that an admixture model has the following cumulative
distribution function (cdf) $$
  L(x) = pF(x) + (1-p)G(x), \qquad x \in \mathbb{R},
$$

where $G$ is a perfectly known cdf, and $p$ and $F$ are unknown.

The cdf $F$ relates to the contamination phenomenon that is added to the
well-known signal $G$, with proportion $p$.

The proportion of the unknown component in the two-component mixture
model can be easily estimated under weak nonparametric assumptions on
the related distribution. The decontaminated version of this unknown
component distribution can then be tested against some other specified
distribution (included another decontaminated unknown component).
Finally, clustering of $K$ populations is made possible, based on
hypothesis tests that compare unknown component distributions. The
package is suited to one-sample as well as multi-samples analysis.

# Installation

<!-- You can install the released version of admix from [CRAN](https://CRAN.R-project.org) with: -->

You can install the released version of admix from
[Github](https://github.com/XavierMilhaud/admix-Rpackage) with:

``` r
#once on CRAN with : install.package("admix")
# from now on:
remotes::install_git("git@github.com:XavierMilhaud/admix.git", build_manual = TRUE, build_vignettes = TRUE)
```

The optional argument build_vignettes can be set to TRUE to get
vignettes that help to understand the functionalities of the package.

To get some help about the functionalities of the package, do once
installed:

``` r
help(package = 'admix')
```

More details can also be found through the vignettes, available in admix
github-pages (see <https://xaviermilhaud.github.io/admix-Rpackage/>, in Menu
Articles).

## Example

This is a basic example which shows you how to estimate the unknown
component proportion and the localization shift parameters in an
admixture model where the unknown component density is assumed to be
symmetric. In practice, the cdf $L$ is given by $$
L(x) = p F(x-\mu) + (1-p) G(x), \qquad x \in \mathbb{R},
$$ where $p$ is the unknown component weight, and $\mu$ is the
localization shift parameter of the unknown cdf $F$ with symmetric
density.

The estimation would be made through the following commands:

``` r
library(admix)
#> Package 'admix' version 2.3.1
#> -------------------------------
#> Type 'citation("admix")' for citing this R package in publications.
#> -------------------------------
#> This work was partly conducted within the Research Chair DIALog under the aegis of the Risk Foundation, an initiative by CNP Assurances.
```

``` r
## Simulate mixture data:
mixt1 <- twoComp_mixt(n = 450, weight = 0.4,
                      comp.dist = list("norm", "norm"),
                      comp.param = list(list("mean" = -2, "sd" = 0.5),
                                        list("mean" = 0, "sd" = 1)))
data1 <- getmixtData(mixt1)
## Define the admixture models:
admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
                         knownComp_param = mixt1$comp.param[[2]])
## Estimation step:
admix_estim(samples = list(data1),
            admixMod = list(admixMod1),
            est.method = 'BVdk', sym.f = TRUE)
#> Call:
#> admix_estim(samples = list(data1), admixMod = list(admixMod1), 
#>     est.method = "BVdk", sym.f = TRUE)
#> 
#> Estimated mixing weight of the unknown component distribution in Sample 1: 0.4
#> 
#> Estimated location parameters of the unknown component distribution in Sample 1: -2.09
```

<!-- badges: end -->
