#' Machine learning and inference for persistence diagrams
#'
#' This package aims to bridge topological data analysis (TDA) with data, statistical
#' and machine learning practitioners so that more analyses may benefit from the
#' power of TDA. The main tool of TDA is persistent homology, which computes a 
#' shape descriptor of a dataset, called a persistence diagram. There are four
#' goals of this package: (1) convert the output from the persistent homology
#' calculations in the popular TDA package for R, 'TDA', into a data frame which can
#' easily be used in personalized analyses, (2) provide a fast method for computing
#' distances between persistence diagrams, (3) implement kernel machine learning
#' method for persistence diagrams (currently kernel multidimensional scaling and
#' kernel SVM), and (4) provide methods for inference on groups of persistence diagrams.
#'
#' @docType package
#' @name TDAML
#' @import TDA
#' @importFrom clue solve_LSAP
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel clusterEvalQ clusterExport detectCores makeCluster stopCluster
#' @importFrom rdist cdist
#' @importFrom stats complete.cases pgamma
#' @importFrom utils combn
#' @importFrom iterators iter
NULL
