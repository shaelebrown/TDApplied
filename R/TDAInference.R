#' Statistical Inference for Persistence Diagrams
#'
#' This package provides the functionality to infer when
#' there are topological differences between an arbitrary
#' number of groups of persistence diagrams, and to quickly
#' compute distances between persistence diagrams.
#'
#' @useDynLib TDAInference
#' @importFrom clue solve_LSAP
#' @importFrom doParallel, registerDoParallel
#' @importFrom foreach "%dopar%" foreach
#' @importFrom parallel clusterEvalQ clusterExport detectCores makeCluster stopCluster
#' @importFrom rdist cdist
#' @importFrom stats complete.cases
#' @importFrom utils combn
#' @name TDAInference
#' @docType package
NULL
