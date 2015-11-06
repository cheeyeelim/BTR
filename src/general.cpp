#include <Rcpp.h>
#include <string>
#include <iostream>
#include <sstream>
#include <regex>
#include <vector>
#include <iterator>
#include <map>
#include <exception>
#include <stdexcept>

using namespace std;

//internal C++ functions. not exported.
bool all_cpp(Rcpp::LogicalVector x) {
  /*#Function to test for equality between two vectors.*/
 
  return is_true(Rcpp::all(x==TRUE)); //#is_true is used to return a bool type. all() is an Rcpp sugar.
}

//external C++ functions. exported.
//' @title Find a match between two data frames.
//' 
//' @description
//' (&&&Not for public use&&&)This function finds a match between two df of states. Used in match_state(). Return an row index vector indicating which row of mstate matches the rows in xstate.
//' 
//' @param mstate data frame. It should be a state(row) x gene(column) df.
//' @param xstate data frame. It should be a state(row) x gene(column) df.
// [[Rcpp::export]] //#must be called in front of each C++ function to be exported.
Rcpp::NumericVector match_state_loop(Rcpp::NumericMatrix mstate, Rcpp::NumericMatrix xstate) {
  int nrow_ms = mstate.nrow();
  int nrow_xs = xstate.nrow();
  Rcpp::NumericVector ind(nrow_ms); //#must specify length if using vector.

  int turn = 0;
  for(auto i=0; i<nrow_ms; i++) { //#Pick each row of mstate.
    Rcpp::NumericVector row_ms = mstate.row(i);

    //#if there is intention to speed up the computation, do as below:
    //#2 steps check to speed up computation.
    //#1st step is to sum up all rows, and look for matching sums.
    //#2nd step is to only compare the rows with the same sums.

    for(auto j=0; j<nrow_xs; j++) { //#Pick each row of xstate.
      Rcpp::NumericVector row_xs = xstate.row(j);

  	  bool check = all_cpp(row_ms == row_xs);
      if(check) { //#Make note of matching state.
        ind[turn] = j+1;
		    turn += 1;
      }
    }
  }
  
  return ind;
}
