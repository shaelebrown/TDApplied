#### PERMUTATION TEST FOR DIFFERENCES ####
#' Permutation test for finding group differences between persistence diagrams.
#' 
#' A non-parametric ANOVA-like test for persistence diagrams 
#' (see \url{https://link.springer.com/article/10.1007/s41468-017-0008-7} for details). In each
#' desired dimension a test statistic (loss) is calculated, then the group labels are shuffled
#' for some number of iterations and the loss is recomputed each time thereby generating a null
#' distribution for the test statistic. This test generates a p-value in each desired dimension.
#' 
#' The test is carried out in parallel and optimized in order to not recompute already-calculated distances. As such, memory issues
#' may occur when the number of persistence diagrams is very large. 
#' Like in (\url{https://github.com/hassan-abdallah/Statistical_Inference_PH_fMRI/blob/main/Abdallah_et_al_Statistical_Inference_PH_fMRI.pdf})
#' an option is provided for pairing diagrams between groups to reduce variance (in order to boost statistical power), and
#' like it was suggested in the original paper functionality is provided for an arbitrary number of groups (not just 2).
#' A small p-value in a dimension suggests that the groups are different (separated) in that dimension.
#' If `distance` is "fisher" then `sigma` must not be NULL. TDAstats also has a `permutation_test` function
#' so care should be taken to use the desired function when using TDApplied with TDAstats.
#'
#' @param ... lists of persistence diagrams which are either the output of a TDA/TDAstats calculations like \code{\link[TDA]{ripsDiag}}/\code{\link[TDAstats]{calculate_homology}}, or \code{\link{diagram_to_df}}. Each list must contain at least 2 diagrams.
#' @param iterations the number of iterations for permuting group labels, default 20.
#' @param p a positive number representing the wasserstein power parameter, a number at least 1 (and Inf if using the bottleneck distance) and default 2.
#' @param q  a finite number at least 1 for exponentiation in the Turner loss function, default 2.
#' @param dims a non-negative integer vector of the homological dimensions in which the test is to be carried out, default c(0,1).
#' @param paired a boolean flag for if there is a second-order pairing between diagrams at the same index in different groups, default FALSE
#' @param distance a string which determines which type of distance calculation to carry out, either "wasserstein" (default) or "fisher".
#' @param sigma the positive bandwidth for the Fisher information metric, default NULL.
#' @param num_workers the number of cores used for parallel computation, default is one less than the number of cores on the machine.
#' @param verbose a boolean flag for if the time duration of the function call should be printed, default FALSE
#'
#' @return a list with the following elements:
#' \describe{
#' 
#'  \item{dimensions}{the input `dims` argument.}
#' 
#'  \item{permvals}{a numeric vector of length `iterations` with the permuted loss value for each iteration (permutation)}
#'  
#'  \item{test_statisics}{a numeric vector of the test statistic value in each dimension.}
#'  
#'  \item{p_values}{a numeric vector of the p-values in each dimension.}
#'  
#'  \item{run_time}{the run time of the function call, containing time units.}
#' 
#' }
#' 
#' @importFrom parallelly availableCores
#' @export
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @references
#' Robinson T, Turner K (2017). "Hypothesis testing for topological data analysis." \url{https://link.springer.com/article/10.1007/s41468-017-0008-7}.
#' 
#' Abdallah H et al. (2021). "Statistical Inference for Persistent Homology applied to fMRI." \url{https://github.com/hassan-abdallah/Statistical_Inference_PH_fMRI/blob/main/Abdallah_et_al_Statistical_Inference_PH_fMRI.pdf}.
#' @examples
#'
#' # create two groups of diagrams
#' D1 <- TDAstats::calculate_homology(TDA::circleUnif(n = 10,r = 1),
#'                                    dim = 0,threshold = 2)
#' D2 <- TDAstats::calculate_homology(TDA::circleUnif(n = 10,r = 1),
#'                                    dim = 0,threshold = 2)
#' g1 <- list(D1,D2)
#' g2 <- list(D1,D2)
#'
#' # run test in dimension 0 with 1 iteration
#' perm_test <- permutation_test(g1,g2,iterations = 1,
#'                               num_workers = 2,
#'                               dims = c(0))

permutation_test <- function(...,iterations = 20,p = 2,q = 2,dims = c(0,1),paired = FALSE,distance = "wasserstein",sigma = NULL,num_workers = parallelly::availableCores(omit = 1),verbose = FALSE){

  # function to test whether or not multiple groups of persistence diagrams come from the same geometric process
  # ... are the groups of diagrams, either stored as lists or vectors
  # iterations is the number of permutations we will calculate for group labels
  # p is the wasserstein distance parameter, p >= 1
  # q is the finite distance exponential, q >= 1
  # dims is a vector of desired homological dimensions
  # paired is a boolean which determines if dependencies exist between diagrams of the same indices in different groups
  # distance is either "wasserstein" or "fisher" and determines how distances will be computed between diagrams
  # sigma is the positive bandwidth for the Fisher information metric, NULL by default
  # num_workers is the number of cores used for parallelization, default available number of cores minus 1.
  # verbose is either TRUE or FALSE (default), printing runtime of function call

  # retrieve diagram groups
  diagram_groups <- list(...)

  # make sure there are at least two groups
  check_param("diagram_groups",diagram_groups,min_length = 2)
  
  # make sure each group has at least two elements
  lapply(X = diagram_groups,FUN = function(X){
    
    if(length(X) < 2)
    {
      stop("Each group of diagrams must have at least 2 elements.")
    }
    
  })

  # check each diagram, converting each to a data frame and storing their indices in all the diagrams
  all_diagrams(diagram_groups,inference = "difference")
  diagram_groups <- all_diagrams(diagram_groups,inference = "difference")

  # error check function parameters
  check_param("iterations",iterations,whole_numbers = T,at_least_one = T,numeric = T,finite = T)
  check_param("p",p,finite = F,at_least_one = T,numeric = T,multiple = F)
  check_param("q",q,at_least_one = T,numeric = T,multiple = F)
  check_param("dims",dims,multiple = T,whole_numbers = T,non_negative = T,positive = F,numeric = T,finite = T)
  check_param("paired",paired,numeric = F,multiple = F)
  check_param("distance",distance)
  if(distance == "fisher")
  {
    check_param("sigma",sigma,non_negative = T,positive = F,numeric = T,finite = T,multiple = F)
  }
  check_param("num_workers",num_workers,whole_numbers = T,at_least_one = T,numeric = T,multiple = F)
  if(num_workers > parallelly::availableCores())
  {
    warning("num_workers is greater than the number of available cores - setting to maximum value.")
    num_workers <- parallelly::availableCores()
  }

  # make sure that if paired == T then all groups have the same number of diagrams
  if(paired)
  {
    if(length(unique(unlist(lapply(diagram_groups,FUN = length)))) != 1)
    {
      stop("If paired is true then all groups of diagrams must have the same number of elements.")
    }
  }

  # time computations
  start_time <- Sys.time()

  # make distance matrix for all diagrams for each dimension, store total number of diagrams across
  # all groups in variable n
  n <- sum(unlist(lapply(diagram_groups,FUN = length)))
  dist_mats <- lapply(X = dims,FUN = function(X){return(matrix(data = -1,nrow = n,ncol = n))})

  # compute loss function on observed data and update dist_mats
  test_loss <- loss(diagram_groups = diagram_groups,dist_mats = dist_mats,dims = dims,p = p,q = q,distance = distance,sigma = sigma,num_workers = num_workers)
  dist_mats <- test_loss$dist_mats
  test_statistics <- test_loss$statistics

  # compute cumulative sum of diagram group sizes for inverting diagram indices
  csum_group_sizes <- c(0,cumsum(unlist(lapply(X = diagram_groups,FUN = length))))

  # create empty list in each dimension to store loss function statistic for each permutation
  perm_values <- lapply(X = dims,FUN = function(X){return(list())})

  # permute group labels and recalculate the loss function each time
  for(perm in 1:iterations)
  {

    if(paired == F)
    {
      # sample groups from their union, maintaining group sizes
      perm <- split(x = 1:n,f = sample(unlist(lapply(X = 1:length(diagram_groups),FUN = function(X){return(rep(X,length = length(diagram_groups[[X]])))})),size = n,replace = F))

      permuted_groups <- lapply(X = perm,FUN = function(X){

        res <- list()

        # for each index in that permuted group
        for(ind in 1:length(X))
        {
          # invert diagram indices
          g <- min(which(csum_group_sizes >= X[ind])) - 1

          # append to permuted group
          res[[length(res) + 1]] <- diagram_groups[[g]][[X[[ind]] - csum_group_sizes[g]]]
        }
        return(res)})
    }else
    {
      # permute all the 1st diagrams between the groups, 2nd diagrams between the groups etc.
      perm <- lapply(X = 1:length(diagram_groups[[1]]),FUN = function(X){return(sample(1:length(diagram_groups),size = length(diagram_groups),replace = F))})

      permuted_groups <- lapply(X = 1:length(diagram_groups),FUN = function(X){

        # get the Xth element of each list
        groups <- sapply(perm,"[[",X)

        res <- list()

        # for each group g in the permuted group labels
        for(g in 1:length(groups))
        {
          # append correct diagram to the end of the list
          res[[length(res) + 1]] <- diagram_groups[[groups[[g]]]][[X]]
        }
        return(res)})
    }

    # compute loss function, add to permutation values and updated distance matrices
    permuted_loss <- loss(permuted_groups,dist_mats = dist_mats,dims = dims,p = p,q = q,distance = distance,sigma = sigma,num_workers = num_workers)
    dist_mats <- permuted_loss$dist_mats
    for(d in 1:length(dims))
    {
      perm_values[[d]][[(length(perm_values[[d]]) + 1)]] <- permuted_loss$statistics[[d]]
    }

  }

  # set up return list
  names(perm_values) <- as.character(dims)
  for(i in 1:length(perm_values))
  {
    perm_values[[i]] <- unlist(perm_values[[i]])
  }

  names(test_statistics) <- as.character(dims)

  pval <- lapply(X = 1:length(dims),FUN = function(X){

    return((sum(perm_values[[X]] <= test_statistics[[X]]) + 1)/(iterations + 1))

  })
  names(pval) <- as.character(dims)
  pval <- unlist(pval)

  runtime = Sys.time() - start_time

  results <- list(dimensions = dims,
                  permvals = perm_values,
                  test_statistics = test_statistics,
                  p_values = pval,
                  run_time = runtime)

  if(verbose == T)
  {
    # print time duration of function call
    print(paste0("Time taken: ",runtime))
  }

  return(results)

}

#### INDEPENDENCE TEST FOR PERSISTENCE DIAGRAMS ####
#' Independence test for two groups of persistence diagrams.
#'
#' Carries out inference to determine if two groups of persistence diagrams are independent or not
#' based on kernel calculations (see 
#' (\url{https://proceedings.neurips.cc/paper/2007/file/d5cfead94f5350c12c322b5b664544c1-Paper.pdf}) for details).
#' A small p-value in a certain dimension suggests that the groups are not independent in that dimension.
#' 
#' The test is carried out with a parametric null distribution, making it much faster than non-parametric
#' approaches. If all of the diagrams in either g1 or g2 are the same in some dimension, then some p-values may be NaN.
#'
#' @param g1 the first group of persistence diagrams, where each diagram was either the output from a TDA/TDAstats calculation like \code{\link[TDA]{ripsDiag}}/\code{\link[TDAstats]{calculate_homology}}, or \code{\link{diagram_to_df}}.
#' @param g2 the second group of persistence diagrams, where each diagram was either the output from a TDA/TDAstats calculation like \code{\link[TDA]{ripsDiag}}/\code{\link[TDAstats]{calculate_homology}}, or \code{\link{diagram_to_df}}.
#' @param dims a non-negative integer vector of the homological dimensions in which the test is to be carried out, default c(0,1).
#' @param sigma a positive number representing the bandwidth for the Fisher information metric, default 1.
#' @param t a positive number representing the scale for the persistence Fisher kernel, default 1.
#' @param num_workers the number of cores used for parallel computation, default is one less than the number of cores on the machine.
#' @param verbose a boolean flag for if the time duration of the function call should be printed, default FALSE
#'
#' @return a list with the following elements:
#' \describe{
#' 
#'  \item{dimensions}{the input `dims` argument.}
#'  
#'  \item{test_statisics}{a numeric vector of the test statistic value in each dimension.}
#'  
#'  \item{p_values}{a numeric vector of the p-values in each dimension.}
#'  
#'  \item{run_time}{the run time of the function call, containing time units.}
#' 
#' }
#' 
#' @importFrom stats pgamma
#' @importFrom parallelly availableCores
#' @export
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @references
#' Gretton A et al. (2007). "A Kernel Statistical Test of Independence." \url{https://proceedings.neurips.cc/paper/2007/file/d5cfead94f5350c12c322b5b664544c1-Paper.pdf}.
#' @examples
#'
#' # create two independent groups of diagrams of length 6, which
#' # is the minimum length
#' D1 <- TDAstats::calculate_homology(TDA::circleUnif(n = 10,r = 1),
#'                                    dim = 0,threshold = 2)
#' D2 <- TDAstats::calculate_homology(TDA::circleUnif(n = 10,r = 1),
#'                                    dim = 0,threshold = 2)
#' g1 <- list(D1,D2,D2,D2,D2,D2)
#' g2 <- list(D2,D1,D1,D1,D1,D1)
#' 
#' # do independence test with sigma = t = 1 in dimension 0
#' indep_test <- independence_test(g1,g2,dims = c(0),num_workers = 2)

independence_test <- function(g1,g2,dims = c(0,1),sigma = 1,t = 1,num_workers = parallelly::availableCores(omit = 1),verbose = FALSE){
  
  # function to test whether or not two groups of persistence diagrams are independent
  # g1 and g2 are the groups of diagrams, either stored as lists or vectors
  # dims is a vector of desired homological dimensions
  # sigma is the positive bandwidth for the Fisher information metric, 1 by default
  # t is the positive scale for the persistence Fisher kernel, 1 by default
  # num_workers is the number of cores used in parallelization, maximum minus 1 by default.
  # verbose is either TRUE or FALSE (default), printing runtime of function call
  
  # set internal variables to NULL to avoid build issues
  r <- NULL
  
  # retrieve diagram groups
  diagram_groups <- list(g1,g2)
  
  # make sure there are at least two groups
  check_param("diagram_groups",diagram_groups,min_length = 2)
  
  # check each diagram, converting each to a data frame
  all_diagrams(diagram_groups,inference = "independence")
  diagram_groups <- all_diagrams(diagram_groups,inference = "independence")
  
  # error check function parameters
  check_param("dims",dims,multiple = T,whole_numbers = T,non_negative = T,positive = F,numeric = T,finite = T)
  check_param("sigma",sigma,non_negative = F,positive = T,multiple = F,finite = T,numeric = T)
  check_param("t",t,non_negative = F,positive = T,numeric = T,multiple = F,finite = T)
  
  # make sure that the two groups have the same number of diagrams
  if(length(g1) != length(g2))
  {
    stop("g1 and g2 must be the same length.")
  }
  
  # make sure that the two groups each have at least 6 elements
  m <- length(g1)
  if(m < 6)
  {
    stop("g1 and g2 must have at least 6 elements each.")
  }
  
  # time computations
  start_time <- Sys.time()
  
  # conduct test in each dimension
  test_statistics <- c()
  p_value <- c()
  
  for(dim in dims)
  {
    
    # compute test statistics in parallel
    K <- gram_matrix(diagrams = diagram_groups[[1]],dim = dim,t = t,sigma = sigma,num_workers = num_workers)
    L <- gram_matrix(diagrams = diagram_groups[[2]],dim = dim,t = t,sigma = sigma,num_workers = num_workers)

    H <- matrix(data = -1/m,nrow = m,ncol = m)
    diag(H) <- rep((m-1)/m,nrow(H))
    
    HSIC <- sum(diag(K %*% H %*% L %*% H))/(m^2) # normalized trace
    test_statistics <- c(test_statistics,HSIC)
    
    # compute null distribution parameters
    mu_x_sq <- mean(K[upper.tri(K)])
    mu_y_sq <- mean(L[upper.tri(L)])
    mu <- (1 + mu_x_sq*mu_y_sq - mu_x_sq - mu_y_sq)/m
    B <- (H %*% K %*% H) * (H %*% L %*% H)
    B <- B * B
    diag(B) <- rep(0,m)
    var <- 2*(m-4)*(m-5)*factorial(m - 4)*(sum(colSums(B)))/factorial(m)
    if(var == 0)
    {
      stop("A zero variance was calculated, please make sure that both g1 and g2 contain at least 2 distinct diagrams.")
    }
    
    # compute p-value
    p_val <- stats::pgamma(q = HSIC,rate = mu/var,shape = mu^2/var,lower.tail = F)
    p_value <- c(p_value,p_val)
    
  }
  
  # set up return lists
  names(test_statistics) <- as.character(dims)
  names(p_value) <- as.character(dims)
  runtime = Sys.time() - start_time
  
  results <- list(dimensions = dims,
                  test_statistics = test_statistics,
                  p_values = p_value,
                  run_time = runtime)
  
  if(verbose == T)
  {
    # print time duration of function call
    print(runtime)
  }
  
  return(results)
  
}
