#' Machine learning and inference for persistence diagrams
#'
#' Topological data analysis is a powerful tool for finding non-linear global structure
#' in whole datasets. The main tool of topological data analysis is persistent homology, which computes
#' a topological shape descriptor of a dataset called a persistence diagram. 'TDApplied' provides 
#' useful and scalable methods for analyzing groups of persistence diagrams with machine learning and statistical inference,
#' and these functions can also interface with other data science packages to form flexible and integrated
#' topological data analysis pipelines.
#'
#' @useDynLib TDApplied, .registration = TRUE
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
#' @keywords internal
"_PACKAGE"
NULL
