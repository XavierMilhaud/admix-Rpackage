# Hypothesis test in admixture models

``` r

library(admix)
```

We remind that a random variable $`X`$ following an admixture
distribution has cumulative distribution function (cdf) $`L`$ given by
``` math
L(x) = pF(x) + (1-p)G(x), \qquad x \in \mathbb{R},
```
where $`G`$ is a mixture component whose distribution is perfectly
known, whereas $`p`$ and $`F`$ are unknown. In this setting, if no
parametric assumption is made on the unknown component distribution
$`F`$, the mixture is considered as a semiparametric mixture. For an
overview on semiparametric extensions of finite mixture models, see
(Xiang and Yang 2018).

The goal of this vignette is to introduce the functionalities that
enable to perform hypothesis tests on the unknown component distribution
$`F`$. We aim to test whether $`F`$ belongs to certain parametric family
(e.g. the Gaussian one) in a 1-sample case, or if two different
decontaminated versions of $`F_1`$ and $`F_2`$ (obtained from two
observed samples $`X_1`$ and $`X_2`$) are similar in the K-sample case
($`K \geq 2`$). All the specific tests presented hereafter can be
performed using one single generic function for testing with appropriate
arguments, the so-called $`admix\_test`$ function.

## The one-sample case only available to symmetric unknown density

In this setting, the test to be performed is a parametric family
testing, i.e.
``` math
H_0: \, F\in \mathcal{F} \qquad \mbox{against} \qquad  H_1: \, F\notin \mathcal{F},
```
where $`\mathcal{F}=\left\{F_\theta:~\theta\in \Theta \right\}`$.

The support of the known component density has to be in line with the
one of the unknown component density. Such tests have been introduced in
(Pommeret and Vandekerkhove 2019), and the idea underlying this
hypothesis test follows these steps:

1.  decompose the observed density and the known density in an
    orthonormal polynomial basis,
2.  get the expansion coefficients of such densities,
3.  reformulate the null hypothesis of the test using these
    coefficients,
4.  adopt a $`\chi^2`$ test strategy that relies on Central Limit
    Theorem (CLT) results on estimators of the (unknown) weight related
    to the unknown component distribution $`F`$.

Because of the use of asymptotically normal estimators, it is not
possible to use the estimator provided in (Patra and Sen 2016) to
perform hypothesis testing. On the contrary, Bordes and Vandekerkhove
(2010) propose an estimator that can be used if the unknown component
density is assumed to be symmetric. More generally, Pommeret and
Vandekerkhove (2019) give more details about the distribution of the
test statistic under the null (hypothesis $`H_0`$), and under the
alternative $`H_1`$.

Here, the implemented function allows to perform the so-called
Gaussianity test, meaning that the parametric family against which the
unknown component is tested belongs to Gaussian distributions. Below is
an example of hypothesis testing in this 1-sample case:

``` r

####### Under the null hypothesis H0.
## Simulate mixture data:
mixt1 <- twoComp_mixt(n = 300, weight = 0.6,
                      comp.dist = list("norm", "norm"),
                      comp.param = list(c("mean" = 2, "sd" = 0.5),
                                        c("mean" = 0, "sd" = 1)))
data1 <- get_mixture_data(mixt1)
## Define the admixture model:
admixMod <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
                        knownComp_param = mixt1$comp.param[[2]])
admix_test(samples = list(data1), admixMod = list(admixMod), conf_level = 0.95,
           test_method = "poly", ask_poly_param = FALSE, support = "Real")
#> Call:admix_test(samples = list(data1), admixMod = list(admixMod), 
#>     test_method = "poly", conf_level = 0.95, ask_poly_param = FALSE, 
#>     support = "Real")
#> 
#> Is the null hypothesis H0 rejected? No
#> p-value of the test: 0.682
```

The result of the test is that we cannot reject the null hypothesis
$`H_0`$, which is in line with the specified distribution for the
unknown component. Indeed, simulated data is a Gaussian mixture with two
components, i.e. $`F \sim \mathcal{N}(\mu,\sigma)`$ where $`\mu=2`$ and
$`\sigma=0.5`$.

## The two-sample case

Let us introduce two random samples $`X_1`$ and $`X_2`$ following
admixture models, such that
``` math
\begin{align*}
  \left\{
    \begin{array}{l}
        L_1(x) = (1-p_1)G_1(x) + p_1F_1(x) \\
        L_2(x) = (1-p_2)G_2(x) + p_2F_2(x),
    \end{array}
    \right.
\end{align*}
```

The goal here is to perform the following hypothesis test:
``` math
H_0: ~ F_1=F_2 \qquad \mbox{against} \qquad H_1: F_1\neq F_2.
```

### Case of symmetric unknown densities

In this framework, we assume that $`F_1`$ and $`F_2`$ both have a
symmetric density. This way the normally-distributed estimator of
$`p_1`$ and $`p_2`$, proposed in (Bordes and Vandekerkhove 2010), can be
used together with the testing strategy suggested in (Milhaud et al.
2022). This testing strategy is closely connected to (Pommeret and
Vandekerkhove 2019), where the computation of the expansion coefficients
is duplicated on each of the two samples under study.

In what follows, we simulate two samples under the null and check
whether the test provides satisfactory results.

``` r

mixt1 <- twoComp_mixt(n = 600, weight = 0.8,
                      comp.dist = list("norm", "norm"),
                      comp.param = list(list("mean" = 3, "sd" = 0.5),
                                        list("mean" = 0, "sd" = 1)))
mixt2 <- twoComp_mixt(n = 800, weight = 0.7,
                      comp.dist = list("norm", "norm"),
                      comp.param = list(list("mean" = 3, "sd" = 0.5),
                                        list("mean" = 0, "sd" = 1)))
data1 <- get_mixture_data(mixt1)
data2 <- get_mixture_data(mixt2)
## Define the admixture models:
admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
                         knownComp_param = mixt1$comp.param[[2]])
admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
                         knownComp_param = mixt2$comp.param[[2]])
## Using expansion coefficients in orthonormal polynomial basis:
admix_test(samples = list(data1,data2), admixMod = list(admixMod1,admixMod2), conf_level = 0.95,
           test_method = "poly", ask_poly_param = FALSE, support = "Real")
#> Call:admix_test(samples = list(data1, data2), admixMod = list(admixMod1, 
#>     admixMod2), test_method = "poly", conf_level = 0.95, ask_poly_param = FALSE, 
#>     support = "Real")
#> 
#> Is the null hypothesis H0 rejected? No
#> p-value of the test: 0.096
```

The hypothesis test concludes that the null hypothesis cannot be
rejected, once again in line with what was expected given the specified
parameters when simulating the data.

Note that the following arguments, involved in subroutines, were set to
default values (but that the user can choose them setting parameter
‘ask_poly_param’ to TRUE):

- ‘est_method’ is set to ‘BVdk’ to tell the program to estimate the
  unknown proportions $`p_1`$ and $`p_2`$ using the estimator proposed
  in (Bordes and Vandekerkhove 2010),
- ‘K’ equals 3 to mention that such expansions are computed up to the
  third order of the decomposition in the polynomial basis,
- ‘s’ equals 0.25 as the penalization rate involved in the penalization
  rule used in (Milhaud et al. 2022).

When the unknown component distributions are not supposed to have
symmetric densities, another solution to perform the test is to set
‘est_method’ to ‘PS’ following ‘ask_poly_param’ set to TRUE, but keeping
in mind that plugging in such estimators for the test procedure should
not be allowed in theory. However, that works quite well in practice.
Another solution to perform this test in full generality is to use the
IBM method (see below).

### Case of fully unknown densities

Estimation of the unknown quantities is made by the Inversion - Best
Matching approach, see (Milhaud et al. 2024b). In this case, one can
still use the function with same first arguments except ‘method’, and
‘method’ should be set to ‘icv’. The user also has to define the number
of simulated Gaussian processes used to tabulate the test statistic
distribution (‘n_sim_tab’), and can accelerate computations using
parallel computations and choosing an adequate number of cpus. Other
arguments such as $`support`$ are useless.

## The K-sample case

We introduce hereafter a natural extension of the two-sample case to the
K-sample one, see (Milhaud et al. 2024a). In what follows, the K-sample
test is illustrated within the framework of the IBM approach, i.e. using
the associated inner convergence property. Of course, in the case when
all the unknown component densities are assumed to be symmetric, one
could use a pairwise version of the two sample test using the comparison
of expansion coefficients in a polynomial orthonormal basis, associated
to the estimation method provided by (Bordes and Vandekerkhove 2010).

Consider $`K`$ samples. For $`i=1,...,K`$, sample
$`X^{(i)} = (X_1^{(i)}, ..., X_{n_i}^{(i)})`$ follows
``` math
L_i(x) = p_i F_i(x) + (1-p_i) G_i, \qquad x \in \mathbb{R}.
```
The test to perform is given by
``` math
H_0 : \; F_1 = ... = F_K \qquad \mbox{against} \qquad H_1: \; F_i \neq F_j \quad \mbox{for some} \quad i \neq j.
```
We use the IBM approach to do so, where assumptions are
(straightforwardly) adapted to deal with the $`K`$ samples.

Basically, we apply the theoretical results of IBM for each pair of
populations $`(i,j)`$, and then build a series of embedded statistics.

Consider the set of pair indices: \${\cal S}(K) = \\(i,j)\in
\mathbb{N}^2 ; \\ 1\leq i\<j \leq K\\\$.\\ Order \${\cal S}(K)\$
lexicographically, and denote $`r_K[(i,j)]`$ the rank of $`(i,j)`$ in
the set $`S(K)`$.

Then, $`\forall i\neq j \in \{1,...,K\}`$,

1.  Estimate
    $`\hat{\theta}_{n}(i,j) = \arg\min_{\theta\in \Theta_{i,j}}d_n[i,j](\theta)`$,
2.  Compute the statistic
    $`T_{i,j} = n \, d_n[i,j](\hat{\theta}_n(i,j))`$.

We then obtain $`d(K) = K(K-1)/2`$ comparisons that we embed in a series
of statistics:
``` math
\begin{eqnarray*}
    U_1 &= &T_{1,2} \\
    U_2 & = & T_{1,2}+T_{1,3} \\
        &\vdots& \\
    U_{d(K)} & = &T_{1,2} + \cdots + T_{K-1,K},
\end{eqnarray*}
```

To choose automatically the right order $`k`$ for testing, consider the
penalization rule (mimicking Schwarz criteria procedure, see (Schwarz
1978)):
``` math
S(n) = \min \left\{\arg\max_{1 \leq k \leq d(K)} \left ({U}_{k} - k \sum_{(i,j)\in S(K)} l_n(i,j) \; 1_{\left\{r_K(i,j)=k\right\}} \right )  \right\}.
```

Our data-driven test statistic is given by
``` math
\tilde{U}_n =  U_{S(n)}.
```

It can be shown that under $`H_0`$ and appropriate assumptions, $`S(n)`$
converges in probablity towards 1 as $`n \rightarrow +\infty`$; meaning
that we asymptotically choose the first element of \${\cal S}(K)\$.\\
Moreover, under $`H_0`$, $`U_{S(n)}`$ converges in law towards
$`U^{0}(1,2)`$, which is exactly the null limit distribution studied in
the two-sample case. Finally, we thus consider the $`H_0`$-rejection
rule:
``` math
\tilde{U}_n \quad \geq \quad \hat q_{1-\alpha} \qquad \Rightarrow \qquad H_0 \; \mbox{is rejected}.
```

We now provide the way to perform this test with the package $`admix`$
with Gaussian mixtures. First, let us study the case where we are under
the null hypothesis $`H_0`$, considering $`K=3`$ different populations.

``` r

mixt1 <- twoComp_mixt(n = 450, weight = 0.4,
                      comp.dist = list("norm", "norm"),
                      comp.param = list(c("mean" = -2, "sd" = 0.5),
                                        c("mean" = 0, "sd" = 1)))
mixt2 <- twoComp_mixt(n = 600, weight = 0.24,
                      comp.dist = list("norm", "norm"),
                      comp.param = list(c("mean" = -2, "sd" = 0.5),
                                        c("mean" = -1, "sd" = 1)))
mixt3 <- twoComp_mixt(n = 400, weight = 0.53,
                      comp.dist = list("norm", "norm"),
                      comp.param = list(c("mean" = -2, "sd" = 0.5),
                                        c("mean" = 2, "sd" = 1)))
data1 <- get_mixture_data(mixt1)
data2 <- get_mixture_data(mixt2)
data3 <- get_mixture_data(mixt3)
admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
                         knownComp_param = mixt1$comp.param[[2]])
admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
                         knownComp_param = mixt2$comp.param[[2]])
admixMod3 <- admix_model(knownComp_dist = mixt3$comp.dist[[2]],
                         knownComp_param = mixt3$comp.param[[2]])
admix_test(samples = list(data1, data2, data3), 
           admixMod = list(admixMod1, admixMod2, admixMod3), 
           conf_level = 0.95, test_method = "icv", n_sim_tab = 4, 
           tune_penalty = FALSE, parallel = FALSE, n_cpu = 2)
#> Call:admix_test(samples = list(data1, data2, data3), admixMod = list(admixMod1, 
#>     admixMod2, admixMod3), test_method = "icv", conf_level = 0.95, 
#>     n_sim_tab = 4, tune_penalty = FALSE, parallel = FALSE, n_cpu = 2)
#> 
#> Is the null hypothesis H0 rejected? No
#> p-value of the test: 0.667
```

## References

Bordes, L., and P. Vandekerkhove. 2010. “Semiparametric Two-Component
Mixture Model with a Known Component: An Asymptotically Normal
Estimator.” *Mathematical Methods of Statistics* 19 (1): 22–41.
https://doi.org/<https://doi.org/10.3103/S1066530710010023>.

Milhaud, Xavier, Denys Pommeret, Yahia Salhi, and Pierre Vandekerkhove.
2022. “Semiparametric Two-Sample Admixture Components Comparison Test:
The Symmetric Case.” *Journal of Statistical Planning and Inference*
216: 135–50.
https://doi.org/<https://doi.org/10.1016/j.jspi.2021.05.010>.

Milhaud, Xavier, Denys Pommeret, Yahia Salhi, and Pierre Vandekerkhove.
2024a. “Contamination-Source Based k-Sample Clustering.” *Journal of
Machine Learning Research* 25 (287): 1–32.
<https://jmlr.org/papers/v25/23-0914.html>.

Milhaud, Xavier, Denys Pommeret, Yahia Salhi, and Pierre Vandekerkhove.
2024b. “Two-sample contamination model test.” *Bernoulli* 30 (1):
170–97. <https://doi.org/10.3150/23-BEJ1593>.

Patra, Rohit Kumar, and Bodhisattva Sen. 2016. “Estimation of a
two-component mixture model with applications to multiple testing.”
*Journal of the Royal Statistical Society Series B* 78 (4): 869–93.

Pommeret, Denys, and Pierre Vandekerkhove. 2019. “Semiparametric Density
Testing in the Contamination Model.” *Electronic Journal of Statistics*,
no. 13: 4743–93. https://doi.org/<https://doi.org/10.1214/19-EJS1650>.

Schwarz, G. 1978. “Estimating the Dimension of a Model.” *The Annals of
Statistics* 6 (2): 461–64.

Xiang, Yao, S., and G. Yang. 2018. “An Overview of Semiparametric
Extensions of Finite Mixture Models.” *Statistica Scinica* 34: 391–404.
