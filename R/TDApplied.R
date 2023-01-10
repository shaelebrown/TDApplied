#' Machine learning and inference for persistence diagrams
#'
#' Topological data analysis is a powerful tool for finding non-linear global structure
#' in whole datasets. 'TDApplied' aims to bridge topological data analysis with data, statistical
#' and machine learning practitioners so that more analyses may benefit from the
#' power of topological data analysis. The main tool of topological data analysis is
#' persistent homology, which computes a shape descriptor of a dataset, called
#' a persistence diagram. There are five goals of this package: (1) deliver a fast implementation
#' of persistent homology via a python interface, (2) convert persistence diagrams
#' computed using the two main R packages for topological data analysis into a data frame, 
#' (3) implement fast versions of both distance and kernel calculations
#' for pairs of persistence diagrams, (4) contribute tools for the interpretation of
#' persistence diagrams, and (5) provide parallelized methods for machine learning
#' and inference for persistence diagrams.
#'
#' @docType package
#' @name TDApplied
#' @importFrom clue solve_LSAP
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom iterators iter
#' @importFrom kernlab as.kernelMatrix kkmeans kpca ksvm predict
#' @importFrom parallel clusterEvalQ clusterExport detectCores makeCluster stopCluster
#' @importFrom parallelly availableCores
#' @importFrom rdist cdist
#' @importFrom stats cmdscale complete.cases pgamma lm quantile
#' @importFrom utils combn
#' @importFrom iterators iter
#' @importFrom graphics abline legend
NULL
