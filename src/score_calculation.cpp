#include <Rcpp.h>
#include <string>
#include <iostream>
#include <sstream>
#include <regex>
#include <vector>
#include <iterator>
#include <map>
#include <algorithm>
#include <exception>
#include <stdexcept> //for std::invalid_argument
#include <cmath> //for std::abs(float)
#include <numeric> //for std::accumulate

using namespace std;

//external C++ functions. exported.

//Rcpp::NumericVector rcpp_test(Rcpp::NumericVector x) {
//  std::vector<double> y = Rcpp::as<std::vector <double> >(x);
//  return(Rcpp::wrap(y));
//}

//' @title Calculating pairwise scores between model and data states.
//' 
//' @description
//' This function calculates the pairwise scores between each row of model and data states. The score is calculated using a custom binary distance measure.
//' 
//' @param x_df matrix. It should be numerical matrix of model states.
//' @param y_df matrix. It should be numerical matrix of data states.
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_man_dist(Rcpp::NumericMatrix x_df, Rcpp::NumericMatrix y_df) {
  int xdf_row = x_df.nrow();
  int ydf_row = y_df.nrow();
  
  if(x_df.ncol() != y_df.ncol()) {
    throw std::invalid_argument( "Two input matrices should have the same number of columns." );
  }
  
  std::vector<double> tmpout(xdf_row * ydf_row);
  int ind = 0;
  for(auto i=0; i<xdf_row; i++) {
    for(auto j=0; j<ydf_row; j++) {
      //Convert R objects into C++ objects.
      Rcpp::NumericVector xr = x_df.row(i);
      Rcpp::NumericVector yr = y_df.row(j);
      std::vector<double> x = Rcpp::as<std::vector <double> >(xr);
      std::vector<double> y = Rcpp::as<std::vector <double> >(yr);
      std::vector<double> z;
      
      //Calculate the distance between the two vectors. (Manhattan distance)
      //Calculates the absolute distance between two of each variables.
      std::transform(x.begin(), x.end(), y.begin(), std::back_inserter(z),
      [](double x, double y) { return std::abs(x-y); });
      
      //Get the sum of all differences.
      double o = std::accumulate(z.begin(), z.end(), 0.0);
      
      //Add the value to the output vector.
      tmpout[ind] = o;
      ind += 1; 
      
      //std::copy(x.begin(), x.end(), std::ostream_iterator<double>(std::cout, ","));
      //Rcpp::Rcout << std::endl;
    }
  }

  //Convert Cpp vector into R matrix.
  Rcpp::NumericVector output = Rcpp::wrap(tmpout);
  output.attr("dim") = Rcpp::Dimension(ydf_row, xdf_row);
  
  return output;
}

//' @title Calculating Hamming pairwise scores between model and data states.
//' 
//' @description
//' This function calculates the pairwise scores between each row of model and data states. The score is calculated using a custom binary distance measure.
//' 
//' @param x_df matrix. It should be logical matrix of model states.
//' @param y_df matrix. It should be logical matrix of data states.
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_ham_dist(Rcpp::LogicalMatrix x_df, Rcpp::LogicalMatrix y_df) {
  int xdf_row = x_df.nrow();
  int ydf_row = y_df.nrow();
  
  std::vector<double> m_vec(xdf_row * ydf_row);

  int ind = 0;
  for(auto i=0; i<xdf_row; i++) { //note that this is an R object. So start counting from 1.
    for(auto j=0; j<ydf_row; j++) {
      //Convert R objects into C++ objects.
      Rcpp::LogicalVector xr = x_df.row(i);
      Rcpp::LogicalVector yr = y_df.row(j);
      std::vector<bool> x = Rcpp::as<std::vector <bool> >(xr);
      std::vector<bool> y = Rcpp::as<std::vector <bool> >(yr);
      
      //std::copy(x.begin(), x.end(), std::ostream_iterator<bool>(std::cout, ","));
      //Rcpp::Rcout << "After conversion!" << std::endl;
      
      //Calculate the summation of x and y between each element.
      std::vector<int> z;
      if(x.size() == y.size()) {
        for(unsigned k=0; k<x.size(); k++) {
          z.push_back(x[k] + y[k]); //conversion from bool to int is implicit.
        }
      } else {
        throw std::invalid_argument("row x and row y have different size.");
      }
      
      //Calculate the frequency of numbers.
      //The expected numbers should be 0,1 or 2. a=positive match where both 1s, d=positive match where both 0s, bc=negative match, either 01 or 10.
      double a = 0;
      double d = 0;
      double bc = 0;
      for(unsigned k=0; k<z.size(); k++) {
        if(z[k] == 2) {
          a += 1;
        } else if(z[k] == 0) {
          d += 1;
        } else if(z[k] == 1) {
          bc += 1;
        } else {
          throw std::invalid_argument("Error in calculating the contigency table.");
        }
      }
      
      //Calculate the distance between the two vectors.
      //custom distance measure. it is also known as the mean manhattan distance.
      //bc is included in denominator to restrict the possible range of scores. the score will always be equal to or less than 1.
      //m=0 means exact match, while m=1 means complete difference.
      
      double m = bc/(a+d+bc);
      m_vec[ind] = m;
      ind += 1;
    }
  }
  
  //Convert Cpp vector into R matrix.
  Rcpp::NumericVector output = Rcpp::wrap(m_vec);
  output.attr("dim") = Rcpp::Dimension(ydf_row, xdf_row);
  
  return output;
}

//' @title Inner core for m_score()
//' 
//' @description
//' This function takes in a df with columns ranked wrt each row, and try to assign each row to a unique column without repetition.
//' 
//' @param x_df matrix. Matrix with columns ranked wrt each row.
// [[Rcpp::export]]
Rcpp::IntegerVector rcpp_m_score(Rcpp::IntegerMatrix x_df) {
  int xdf_col = x_df.ncol();
  int xdf_row = x_df.nrow();
  
  Rcpp::IntegerVector yr(xdf_row);
  for(auto i=0; i<xdf_row; i++) { //note that this is an R object. So start counting from 1.
    //Convert R objects into C++ objects.
    Rcpp::IntegerVector xr = x_df.row(i);
    std::vector<int> x = Rcpp::as<std::vector <int> >(xr);
    
    //Obtain the position for the minimum value.
    std::vector<int>::iterator min_val = std::min_element(std::begin(x), std::end(x)); //return the first minimum element. Note that the input result should not contain multiple minimum.
    yr[i] = std::distance(std::begin(x), min_val);
    
    //Check if the value is already present, then pick the next min.
    Rcpp::IntegerVector unique_yr = Rcpp::unique(yr);
    Rcpp::IntegerVector sub_i = Rcpp::seq_len(i);
    Rcpp::IntegerVector test_yr = yr[sub_i];
    bool dup_check = Rcpp::as<bool >(Rcpp::any(Rcpp::duplicated(test_yr)));
  	int ind = 2; //start with finding 2nd min value.
  	// if the model state space is smaller than data state space, it is not possible to assign a unique model state to each data state.
    while (dup_check && unique_yr.size() < xdf_col && ind < xdf_row) {
      //Rcpp::Rcout << i << " : " << ind << std::endl;
      //Rcpp::Rcout << dup_check << std::endl;
      //std::copy(yr.begin(), yr.end(), std::ostream_iterator<int>(Rcpp::Rcout, " "));
      //Rcpp::Rcout << std::endl;
      
  		std::vector<int>::iterator alt_val = std::find(std::begin(x), std::end(x), ind); //find 2nd minimum.
  		yr[i] = std::distance(std::begin(x), alt_val);
  		
      sub_i = Rcpp::seq_len(i);
      test_yr = yr[sub_i];
      dup_check = Rcpp::as<bool >(Rcpp::any(Rcpp::duplicated(test_yr)));
  		unique_yr = Rcpp::unique(yr);
  		ind += 1; //for next iteration, look for the 3rd (and subsequent) values.
    }
  	
  	//Checking results.
    sub_i = Rcpp::seq_len(i);
    test_yr = yr[sub_i];
    dup_check = Rcpp::as<bool >(Rcpp::any(Rcpp::duplicated(test_yr)));
  	if(xdf_col >= xdf_row && dup_check) { //if there are more columns than rows, there should not be any duplicated values.
  		throw std::invalid_argument( "Duplicated values are still present." );
  	}
  }

  return yr + 1; //change C++ positions to R positions.
}

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
  for(auto i=0; i<inf_mat.nrow(); i++) {
    //Convert R objects into C++ objects.
    Rcpp::NumericVector xr = inf_mat.row(i);
    Rcpp::NumericVector yr = true_mat.row(i);
    std::vector<int> x = Rcpp::as<std::vector <int> >(xr);
    std::vector<int> y = Rcpp::as<std::vector <int> >(yr);
    std::vector<int> z;
    
    //Calculate the frequency of numbers.
    //tp=true positive [1,1], tn=true negative [0,0], fp=false positive [1,0], fn=false negative [0,1].
    for(unsigned k=0; k<x.size(); k++) {
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
  
  std::vector<int> output{tp, tn, fp, fn};

  return Rcpp::wrap(output);
}
