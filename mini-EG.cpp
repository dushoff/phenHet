#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector Rrunif_eg(int n) {
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    double r = R::runif(0, 1); 
    out[i] = r;
  }
  return out;
}

// [[Rcpp::export]]
int CppSample_eg(IntegerVector x) {
  IntegerVector samp = Rcpp::sample(x, 1, false);
  int out = samp[0];
  return out;
}



