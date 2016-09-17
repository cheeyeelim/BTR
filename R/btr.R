#' @title BTR: A package for studying asynchronous Boolean models
#' 
#' @description
#' This package contains tools for Boolean model manipulation, as well as the search for the best Boolean model.
#' 
#' @docType package
#' @name BTR
NULL

## All the Roxygen codes below are for generating the correct NAMESPACE file.
#' @import methods
#' @import parallel
#' @import foreach
#' @import doParallel
NULL

#' @useDynLib BTR
#' @importFrom Rcpp sourceCpp evalCpp
NULL

#' @importFrom grDevices rgb
#' @importFrom stats aggregate kmeans median rbinom runif
#' @importFrom utils combn write.csv
NULL