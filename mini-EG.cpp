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
IntegerVector CppSample_eg(IntegerVector x, int n) {
  IntegerVector samp = Rcpp::sample(x, n, false);
  return samp;
}

// [[Rcpp::export]]
NumericVector test_conj(IntegerVector x, int n, int a) {
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    double r = R::runif(0, 1); 
    out[i] = r;
    if (i==a) {
      IntegerVector temp = Rcpp::sample(x, 1, false);
    }
  }
  return out;
}

