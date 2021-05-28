#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double Donsker_correl_cpp(double u, double v, NumericVector sorted_obs_data) {
  int nobs = sorted_obs_data.size();
  double Sigma;
  double a = (std::upper_bound(sorted_obs_data.begin(), sorted_obs_data.end(), std::min(u,v)) - sorted_obs_data.begin()) / ((double) nobs);
  double b = (std::upper_bound(sorted_obs_data.begin(), sorted_obs_data.end(), std::max(u,v)) - sorted_obs_data.begin()) / ((double) nobs);
  Sigma = a * (1 - b);
  return Sigma;
}
