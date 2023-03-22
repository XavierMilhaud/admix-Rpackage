#' Two-samples hypothesis test on the unknown component in admixture models
#'
#' Test hypothesis on the unknown component of admixture models using different estimation techniques, and different
#' testing strategies.
#'
#' @param samples A list of the two observed samples, where each sample follows the mixture distribution given by l = p*f + (1-p)*g,
#'                with f and p unknown and g known.
#' @param known.p (default to NULL) The true component weights p1 and p2 if known, only useful in simulation studies.
#' @param comp.dist A list with four elements corresponding to the component distributions (specified with R native names for these distributions)
#'                  involved in the two admixture models. The two first elements refer to the unknown and known components of the 1st admixture model,
#'                  and the last two ones to those of the second admixture model. If there are unknown elements, they must be specified as 'NULL' objects.
#'                  For instance, 'comp.dist' could be specified as follows: list(f1=NULL, g1='norm', f2=NULL, g2='norm').
#' @param comp.param A list with four elements corresponding to the parameters of the component distributions, each element being a list
#'                   itself. The names used in this list must correspond to the native R argument names for these distributions.
#'                   The two first elements refer to the parameters of unknown and known components of the 1st admixture model, and the last
#'                   two ones to those of the second admixture model. If there are unknown elements, they must be specified as 'NULL' objects.
#'                   For instance, 'comp.param' could be specified as follows: : list(f1=NULL, g1=list(mean=0,sd=1), f2=NULL, g2=list(mean=3,sd=1.1)).
#' @param method Method used for testing. Choose either 'Poly' or 'ICV'. 'Poly' refers to comparison of expansion coefficients in polynomial
#'               orthonormal basis, whereas 'ICV' refers to the Inner Convergence property obtained when using the IBM approach.
#'               More details are provided below in 'Details'.
#' @param n_sim_tab (Only with 'ICV' method) Number of simulated gaussian processes used for the tabulation of the Inner Convergence distribution in IBM approach.
#' @param K (Only for 'Poly' method) Number of coefficients considered for the polynomial basis expansion.
#' @param support (Only for 'Poly' method) Support of the densities under consideration, useful to choose the polynomial orthonormal basis.
#'                 One of 'Real', 'Integer', 'Positive', or 'Bounded.continuous'.
#' @param est.method (Only for 'Poly' method) Either 'BVdk' (Bordes and Valdekerkhove estimation technique) or 'PS' (Patra and Sen estimation
#'                    technique). The latter should not be used since the estimators plugged into the test statistic are not square-root n
#'                    consistent. More details are given in Section 'Details' below.
#' @param s (Only for 'Poly' method) Rate at which the normalization factor is set in the penalization rule for model selection (in ]0,1/2[).
#' @param nb.ssEch (Only with 'Poly' method) Number of subsamples created from original data to decorrelate the estimation of the parameters.
#' @param var.explicit (Only with 'Poly' method) Boolean that enables to choose between explicit evaluation of the variance of the test
#'                     statistic or not (FALSE=bootstrap). FIXME: it seems that bootstrap procedure does not work in the context of admixtures.
#' @param nb.echBoot (Only with 'Poly' method) Number of bootstrap samples if 'var.explicit' is set to FALSE.
#' @param bounds.supp (Only with 'Poly' method) default to NULL. Useful if support = 'bounded.continuous', a list of minimum and maximum bounds,
#'                    specified as follows: list( list(min.f1,min.g1,min.f2,min.g2) , list(max.f1,max.g1,max.f2,max.g2) )
#' @param parallel Boolean to indicate whether parallel computations are performed (speed-up the tabulation).
#' @param n_cpu Number of cores used when parallelizing.
#'
#' @details Here as some details concerning the different methods that can be choosen: i) 'Poly' relies on two-sample testing strategy
#'          where each unknown component density is decomposed in an orthonormal polynomial basis, and the estimation of the component
#'          weights related to the two two-component admixture models can be performed either using Patra and Sen estimator (despite
#'          the latter is not square-root n consistent and thus should not be used in such hypothesis tests), or by Bordes and Vandekerkhove
#'          estimation technique (if the unknown component density is symmetric); ii) 'ICV' refers to Inversion - Best Matching strategy
#'          which has no constraints except that we need to handle two samples.
#'
#' @return The decision of the test with further information such as p-value and others, depending on the method used.
#'
#' @examples
#' \donttest{
#' ##### Under the null hypothesis H0 :
#' ## Simulate data:
#' list.comp <- list(f1 = "norm", g1 = "norm",
#'                   f2 = "norm", g2 = "norm")
#' list.param <- list(f1 = list(mean = 3, sd = 0.5), g1 = list(mean = 0, sd = 1),
#'                    f2 = list(mean = 3, sd = 0.5), g2 = list(mean = 6, sd = 1.2))
#' sample1 <- rsimmix(n=250, unknownComp_weight=0.85, comp.dist = list(list.comp$f1,list.comp$g1),
#'                    comp.param = list(list.param$f1,list.param$g1))[['mixt.data']]
#' sample2 <- rsimmix(n=300, unknownComp_weight=0.8, comp.dist = list(list.comp$f2,list.comp$g2),
#'                    comp.param = list(list.param$f2,list.param$g2))[['mixt.data']]
#' plot_mixt_density(samples = list(sample1,sample2), user.bounds=NULL, support='continuous')
#' ##### Performs the test by the different methods :
#' list.comp <- list(f1 = NULL, g1 = "norm",
#'                   f2 = NULL, g2 = "norm")
#' list.param <- list(f1 = NULL, g1 = list(mean = 0, sd = 1),
#'                    f2 = NULL, g2 = list(mean = 6, sd = 1.2))
#' ## Using expansion coefficients in orthonormal polynomial basis:
#' two_samples_test(samples = list(sample1, sample2), comp.dist=list.comp, comp.param=list.param,
#'                  method = 'Poly', K = 3, support = 'Real', est.method = 'BVdk', s = 0.4,
#'                  nb.ssEch = 2, var.explicit = TRUE)
#' }
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

two_samples_test <- function(samples, known.p = NULL, comp.dist = NULL, comp.param = NULL, method = c("ICV","Poly"),
                             n_sim_tab = NULL, K = 3, support = c('Real','Positive','Integer','Bounded.continuous'),
                             est.method = c("BVdk","PS"), s = 0.49, nb.ssEch = 2, var.explicit = F, nb.echBoot = NULL,
                             bounds.supp = NULL, parallel = FALSE, n_cpu = 2)
{
  if (length(samples) != 2) stop("Please consider TWO samples in the list of samples.")

  meth <- match.arg(method)

  if (meth == "ICV") {
    U <- IBM_tabul_stochasticInteg(n.sim = n_sim_tab, n.varCovMat = 100, sample1 = samples[[1]], sample2 = samples[[2]], min_size = NULL,
                                   comp.dist = comp.dist, comp.param = comp.param, parallel = parallel, n_cpu = n_cpu)
    test_res <- IBM_test_H0(samples = samples, known.p = known.p, comp.dist = comp.dist, comp.param = comp.param,
                            sim_U = U[["U_sim"]], min_size = NULL, parallel = parallel, n_cpu = n_cpu)

  } else if (meth == "Poly") {
    test_res <- orthoBasis_test_H0(samples = samples, known.p = known.p, comp.dist = comp.dist, comp.param = comp.param, known.coef = NULL,
                                   K=K, nb.ssEch = nb.ssEch, s=s, var.explicit = var.explicit, nb.echBoot = nb.echBoot,
                                   support = support, bounds.supp = bounds.supp, est.method = est.method)

  } else stop("Please specify an appropriate method to perform the test hypothesis.")

  return(test_res)
}

