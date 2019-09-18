#include <Rcpp.h>
using namespace Rcpp;

//' This function multiply a numeric vector by 5
//' @param x a numeric vector
//' @export
// [[Rcpp::export]]
NumericVector timesFive(NumericVector x) {
  return x * 5;
}


//' This function uses mvrnorm from MASS
//' @param m a numeric vector
//' @param s the covariance matrix
//' @export
// [[Rcpp::export]]
NumericVector my_mvrnorm(NumericVector m, NumericMatrix s) {
    Environment mvrnorm_env("package:MASS");
    Function mvrnorm = mvrnorm_env["mvrnorm"];
    return (mvrnorm(1, m, s));
}
