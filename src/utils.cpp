#include <Rcpp.h>
using namespace Rcpp;

// from Rcpp gallery https://gallery.rcpp.org/articles/stl-random-shuffle/
// wrapper around R's RNG such that we get a uniform distribution over
// [0,n) as required by the STL algorithm
inline int randWrapper(const int n) { return floor(unif_rand()*n); }

// [[Rcpp::export]]
NumericMatrix row_mean_grouped_dgcmatrix(S4 matrix, IntegerVector group, 
                                         bool shuffle) {
  NumericVector x = matrix.slot("x");
  IntegerVector i = matrix.slot("i");
  IntegerVector p = matrix.slot("p");
  IntegerVector dim = matrix.slot("Dim");
  int rows = dim[0];
  int cols = dim[1];
  CharacterVector levs = group.attr("levels");
  int groups = levs.length();
  NumericMatrix ret(rows, groups);
  IntegerVector groupsize(groups, 0);
  int x_length = x.length();
  
  if (shuffle) {
    group = clone(group);
    std::random_shuffle(group.begin(), group.end(), randWrapper);
  }
  
  int col = 0;
  for (int k=0; k<x_length; ++k) {
    while (k>=p[col]) {
      ++col;
      ++groupsize[group[col-1]-1];
    }
    ret(i[k], group[col-1]-1) += x[k];
  }
  while (col < cols) {
    ++col;
    ++groupsize[group[col-1]-1];
  }
  
  for (int j=0; j<groups; ++j) {
    if (groupsize[j] == 0) {
      ret(_, j) = rep(NumericVector::get_na(), rows);
    } else{
      ret(_, j) = ret(_, j) / groupsize[j];
    }
  }
  colnames(ret) = levs;
  List dn = matrix.slot("Dimnames");
  if (dn[0] != R_NilValue) {
    rownames(ret) = as<CharacterVector>(dn[0]);
  }
  return ret;
}

// [[Rcpp::export]]
NumericMatrix row_gmean_grouped_dgcmatrix(S4 matrix, IntegerVector group, 
                                          double eps, bool shuffle) {
  NumericVector x = matrix.slot("x");
  IntegerVector i = matrix.slot("i");
  IntegerVector p = matrix.slot("p");
  IntegerVector dim = matrix.slot("Dim");
  int rows = dim[0];
  int cols = dim[1];
  CharacterVector levs = group.attr("levels");
  int groups = levs.length();
  NumericMatrix ret(rows, groups);
  IntegerVector groupsize(groups, 0);
  int x_length = x.length();
  IntegerMatrix nonzero(rows, groups);
  double log_eps = log(eps);
  
  if (shuffle) {
    group = clone(group);
    std::random_shuffle(group.begin(), group.end(), randWrapper);
  }
  
  int col = 0;
  for (int k=0; k<x_length; ++k) {
    while (k>=p[col]) {
      ++col;
      ++groupsize[group[col-1]-1];
    }
    ret(i[k], group[col-1]-1) += log(x[k] + eps);
    ++nonzero(i[k], group[col-1]-1);
  }
  while (col < cols) {
    ++col;
    ++groupsize[group[col-1]-1];
  }
  
  for (int j=0; j<groups; ++j) {
    for (int k=0; k<rows; ++k) {
      ret(k, j) = exp((ret(k, j) + log_eps * (groupsize[j] - nonzero(k, j))) / groupsize[j]) - eps;
    }
  }
  colnames(ret) = levs;
  List dn = matrix.slot("Dimnames");
  if (dn[0] != R_NilValue) {
    rownames(ret) = as<CharacterVector>(dn[0]);
  }
  return ret;
}

// [[Rcpp::export]]
IntegerMatrix row_nonzero_count_grouped_dgcmatrix(S4 matrix, IntegerVector group) {
  IntegerVector p = matrix.slot("p");
  IntegerVector i = matrix.slot("i");
  int i_length = i.length();
  IntegerVector dim = matrix.slot("Dim");
  int rows = dim[0];
  CharacterVector levs = group.attr("levels");
  int groups = levs.length();
  IntegerMatrix ret(rows, groups);
  
  int col = 0;
  for (int k=0; k<i_length; ++k) {
    while (k>=p[col]) {
      ++col;
    }
    ret(i[k], group[col-1]-1)++;
  }
  
  colnames(ret) = levs;
  List dn = matrix.slot("Dimnames");
  if (dn[0] != R_NilValue) {
    rownames(ret) = as<CharacterVector>(dn[0]);
  }
  return ret;
}
