#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix estimVarCov_empProcess_Rcpp(NumericVector t, NumericVector obs_data){
  int nobs = obs_data.size();
  std::sort(obs_data.begin(), obs_data.end());
  NumericMatrix Sigma(t.size(),t.size());
  for (int i = 0; i < Sigma.nrow(); i++) {
    for (int j = 0; j < Sigma.ncol(); j++) {
      double a = (std::upper_bound(obs_data.begin(), obs_data.end(),
                                   std::min(t[i], t[j])) - obs_data.begin()) / ((double) nobs);
      double b = (std::upper_bound(obs_data.begin(), obs_data.end(),
                                   std::max(t[i], t[j])) - obs_data.begin()) / ((double) nobs);
      Sigma(i,j) = a * (1 - b);
    }
  }
  return Sigma;
}
