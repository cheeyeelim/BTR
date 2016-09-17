#include <Rcpp.h>
#include <vector>
#include <exception>
#include <stdexcept> //for std::invalid_argument

using namespace std;

//external C++ functions. exported.

//' @title Calculating validation scores between two adjacency matrices
//' 
//' @description
//' This function calculates the validation scores between two adjacency matrices.
//' 
//' @param inf_mat matrix. It should be adjacency matrix of inferred network.
//' @param true_mat matrix. It should be adjacency matrix of true network.
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_validate(Rcpp::NumericMatrix inf_mat, Rcpp::NumericMatrix true_mat) {
  if(inf_mat.ncol() != true_mat.ncol()) {
    throw std::invalid_argument( "Two input matrices should have the same number of columns." );
  }
  if(inf_mat.nrow() != true_mat.nrow()) {
    throw std::invalid_argument( "Two input matrices should have the same number of rows." );
  }
  
  int tp=0;
  int tn=0;
  int fp=0;
  int fn=0;
  for(signed int i=0; i<inf_mat.nrow(); i++) {
    //Convert R objects into C++ objects.
    Rcpp::NumericVector xr = inf_mat.row(i);
    Rcpp::NumericVector yr = true_mat.row(i);
    std::vector<int> x = Rcpp::as<std::vector <int> >(xr);
    std::vector<int> y = Rcpp::as<std::vector <int> >(yr);
    std::vector<int> z;
    
    //Calculate the frequency of numbers.
    //tp=true positive [1,1], tn=true negative [0,0], fp=false positive [1,0], fn=false negative [0,1].
    for(unsigned int k=0; k<x.size(); k++) {
      z.push_back(x[k] + y[k]); //Calculate the summation of x and y between each element.
                
      if(z[k] == 2) {
        tp += 1;
      } else if(z[k] == 0) {
        tn += 1;
      } else if(z[k] == 1) {
        if(x[k] == 0) {
          fp += 1;
        } else {
          fn += 1;
        }
      } else {
        throw std::invalid_argument("Error in calculating the contigency table.");
      }
    }
  }
  
  //std::vector<int> output{tp, tn, fp, fn}; c++11 only
  int tmp_arr[4] = {tp, tn, fp, fn};
  std::vector<int> output(&tmp_arr[0], &tmp_arr[0]+4);

  return Rcpp::wrap(output);
}
