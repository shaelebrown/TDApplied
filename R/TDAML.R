#' Permutation testing for groups of persistence barcodes
#'
#' This package implements a permutation test which can be used to detect differences
#' in multiple groups of persistence barcodes - the most common output of the most
#' common tool in Topological Data Analysis (TDA), persistent homology. Additional
#' functionality is included to account for between-group dependencies if they are
#' present. Also, a fast version of wasserstein and bottleneck distances is provided.
#'
#' @docType package
#' @name TDAML
#' @import TDA
#' @importFrom clue solve_LSAP
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel clusterEvalQ clusterExport detectCores makeCluster stopCluster
#' @importFrom rdist cdist
#' @importFrom stats complete.cases
#' @importFrom utils combn
NULL
