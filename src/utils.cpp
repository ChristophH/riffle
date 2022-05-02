#include <Rcpp.h>
using namespace Rcpp;

// from Rcpp gallery https://gallery.rcpp.org/articles/stl-random-shuffle/
// wrapper around R's RNG such that we get a uniform distribution over
// [0,n) as required by the STL algorithm
inline int randWrapper(const int n) { return floor(unif_rand()*n); }

// [[Rcpp::export]]
NumericMatrix row_mean_grouped_dgcmatrix(S4 matrix, IntegerVector group_labels, 
                                         bool shuffle) {
  NumericVector x = matrix.slot("x");
  IntegerVector i = matrix.slot("i");
  IntegerVector p = matrix.slot("p");
  IntegerVector dim = matrix.slot("Dim");
  int rows = dim[0];
  int cols = dim[1];
  CharacterVector levs = group_labels.attr("levels");
  int groups = levs.length();
  NumericMatrix ret(rows, groups);
  IntegerVector groupsize(groups, 0);
  int x_length = x.length();
  
  if (shuffle) {
    group_labels = clone(group_labels);
    std::random_shuffle(group_labels.begin(), group_labels.end(), randWrapper);
  }
  
  int col = 0;
  for (int k=0; k<x_length; ++k) {
    while (k>=p[col]) {
      ++col;
      ++groupsize[group_labels[col-1]-1];
    }
    ret(i[k], group_labels[col-1]-1) += x[k];
  }
  while (col < cols) {
    ++col;
    ++groupsize[group_labels[col-1]-1];
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
NumericMatrix row_mean_grouped_dense(NumericMatrix matrix, IntegerVector group_labels, 
                                     bool shuffle) {
  int nrows = matrix.nrow();
  int ncols = matrix.ncol();
  List dn = matrix.attr("dimnames");
  CharacterVector levs = group_labels.attr("levels");
  int groups = levs.length();
  IntegerVector groupsize(groups, 0);
  NumericMatrix ret(nrows, groups);
  
  if (shuffle) {
    group_labels = clone(group_labels);
    std::random_shuffle(group_labels.begin(), group_labels.end(), randWrapper);
  }
  
  for (int j = 0; j < ncols; j++) {
    ++groupsize[group_labels[j] - 1];
    for (int i = 0; i < nrows; i++) {
      ret(i, group_labels[j] - 1) += matrix(i, j);
    }
  }
  
  for (int j = 0; j < groups; j++) {
    if (groupsize[j] == 0) {
      ret(_, j) = rep(NumericVector::get_na(), nrows);
    } else{
      ret(_, j) = ret(_, j) / groupsize[j];
    }
  }
  colnames(ret) = levs;
  if (dn[0] != R_NilValue) {
    rownames(ret) = as<CharacterVector>(dn[0]);
  }
  
  return ret;
}

// [[Rcpp::export]]
IntegerMatrix row_nonzero_count_grouped_dgcmatrix(S4 matrix, IntegerVector group_labels) {
  IntegerVector p = matrix.slot("p");
  IntegerVector i = matrix.slot("i");
  int i_length = i.length();
  IntegerVector dim = matrix.slot("Dim");
  int rows = dim[0];
  CharacterVector levs = group_labels.attr("levels");
  int groups = levs.length();
  IntegerMatrix ret(rows, groups);
  
  int col = 0;
  for (int k=0; k<i_length; ++k) {
    while (k>=p[col]) {
      ++col;
    }
    ret(i[k], group_labels[col-1]-1)++;
  }
  
  colnames(ret) = levs;
  List dn = matrix.slot("Dimnames");
  if (dn[0] != R_NilValue) {
    rownames(ret) = as<CharacterVector>(dn[0]);
  }
  return ret;
}
