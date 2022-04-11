#### PERMUTATION TEST FOR DIFFERENCES ####
#' Permutation test for persistence diagrams
#' 
#' Permutation test for finding differences between groups of persistence diagrams,
#' based on the paper \url{https://link.springer.com/article/10.1007/s41468-017-0008-7}. The test is
#' carried out in parallel and optimized in order to not recompute already-calculated distances.
#' Like in \url{https://github.com/hassan-abdallah/Statistical_Inference_PH_fMRI/blob/main/Abdallah_et_al_Statistical_Inference_PH_fMRI.pdf}
#' an option is provided for pairing diagrams between groups to reduce variance (boost statistical power), and
#' like it was suggested in the original paper functionality is provided for an arbitrary number of groups (not just 2).
#'
#' The `...` parameter should be a number of lists of persistence diagrams, outputted from a
#' TDA calculation like \code{\link[TDA]{ripsDiag}} or a \code{\link{diagram_to_df}} function call. The `iterations` parameter should be the number of permutations
#' desired for generating the null distribution. The `p` parameter is the wasserstein power, and `q`
#' is the exponent for distances. `dims` is a numeric vector of the homological dimensions in which
#' to carry out the test. The `paired` parameter is a boolean flag for whether there are correspondences
#' between diagrams at the same location across groups, as this affects which permutations are permissible
#' when generating the null distribution. The `distance` parameter determines which distance metric
#' should be used between persistence diagrams. The `sigma` parameter is the positive bandwith for the
#' Fisher information metric `verbose` determines if the time duration of the function call should be printed.
#'
#' @param ... groups of persistence diagrams, outputted from a homology calculation in TDA.
#' @param iterations the number of iterations for permuting group labels, default 100.
#' @param p the wasserstein parameter, number at least 1 (and Inf if using the bottleneck distance), default 2.
#' @param q  a finite number at least 1 for exponentiation in the Turner loss function, default 2.
#' @param dims a numeric vector of the homological dimensions in which the test is to be carried out, default c(0,1).
#' @param paired a boolean flag for if there is a second-order pairing between diagrams at the same index in different groups. Default value is False.
#' @param distance a string which determines which type of distance calculation to carry out, either "wasserstein" (default) or "fisher".
#' @param sigma the positive bandwith for the Fisher information metric, default NULL.
#' @param verbose a boolean flag for if the time duration of the function call should be printed, default False.
#'
#' @return list with dimensions used (named vector), permutation loss values in each dimension (named list), test statistics in each dimension (named vector)
#'                   a p-value for each dimension (named vector) and the time duration of the function call.
#' @export
#' @examples
#'
#' # create three groups of persistence diagrams on 2D Gaussians using TDA
#' g1 <- lapply(X = 1:3,FUN = function(X){
#'
#' diag <- TDA::ripsDiag(data.frame(x = rnorm(100,mean = 0,sd = 1),
#' y = rnorm(100,mean = 0,sd = 1)),
#' maxscale = 1,
#' maxdimension = 1)
#' df <- diagram_to_df(d = diag)
#' return(df)
#'
#' })
#'
#' g2 <- lapply(X = 1:3,FUN = function(X){
#'
#' diag <- TDA::ripsDiag(data.frame(x = rnorm(100,mean = 0,sd = 1),
#' y = rnorm(100,mean = 0,sd = 1)),
#' maxscale = 1,
#' maxdimension = 1)
#' df <- diagram_to_df(d = diag)
#' return(df)
#'
#' })
#'
#' g3 <- lapply(X = 1:3,FUN = function(X){
#'
#' diag <- TDA::ripsDiag(data.frame(x = rnorm(100,mean = 0,sd = 1),
#' y = rnorm(100,mean = 0,sd = 1)),
#' maxscale = 1,
#' maxdimension = 1)
#' df <- diagram_to_df(d = diag)
#' return(df)
#'
#' })
#'
#' # do permutation test with 20 iterations, p,q = 2, in dimensions 0 and 1, with
#' # no pairing using persistence Fisher distance, sigma = 1, and printing the time duration
#' perm_test = permutation_test(g1,g2,g3,
#' iterations = 20,
#' distance = "fisher",
#' sigma = 1, 
#' verbose = TRUE)

permutation_test <- function(...,iterations = 100,p = 2,q = 2,dims = c(0,1),paired = F,distance = "wasserstein",sigma = NULL,verbose = FALSE){

  # function to test whether or not multiple groups of persistence diagrams come from the same geometric process
  # ... are the groups of diagrams, either stored as lists or vectors
  # iterations is the number of permutations we will calculate for group labels
  # p is the wasserstein distance parameter, p >= 1
  # q is the finite distance exponential, q >= 1
  # dims is a vector of desired homological dimensions
  # paired is a boolean which determines if dependencies exist between diagrams of the same indices in different groups
  # distance is either "wasserstein" or "fisher" and determines how distances will be computed between diagrams
  # sigma is the positive bandwith for the Fisher information metric, NULL by default
  # verbose is either TRUE or FALSE (default), printing runtime of function call

  # retrieve diagram groups
  diagram_groups <- list(...)

  # make sure there are at least two groups
  if(length(diagram_groups) < 2)
  {
    stop("At least two groups of persistence diagrams must be supplied.")
  }

  # check each diagram, converting each to a data frame and storing their indices in all the diagrams
  diagram_groups <- all_diagrams(diagram_groups,inference = "difference")

  # error check function parameters
  check_params(iterations,p,q,dims,paired,distance,sigma)

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
  test_loss <- loss(diagram_groups = diagram_groups,dist_mats = dist_mats,dims = dims,p = p,q = q,distance = distance,sigma = sigma)
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
    permuted_loss <- loss(permuted_groups,dist_mats = dist_mats,dims = dims,p = p,q = q,distance = distance,sigma = sigma)
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
                  test_statistic = test_statistics,
                  p_value = pval,
                  run_time = runtime)

  if(verbose == T)
  {
    # print time duration of function call
    print(paste0("Time taken: ",runtime))
  }

  return(results)

}

#### INDEPENDENCE TEST FOR PERSISTENCE DIAGRAMS ####
#' Independence test for persistence diagrams
#'
#' Calculates (an estimate of) the Hilbert-Schmidt independence criteria for 
#' two groups of paired persistence diagrams, the approximate null distribution
#' and a p-value for each desired homological dimension. See 
#' \url{https://doi.org/10.1007/s41468-017-0008-7} for details.
#'
#' The `g1` and `g2` parameters should be lists of persistence diagrams, outputted from a
#' TDA calculation like \code{\link[TDA]{ripsDiag}} or a \code{\link{diagram_to_df}} function call. `dims` is a numeric vector of the homological dimensions in which
#' to carry out the test. The `sigma` parameter is the positive bandwith for the
#' Fisher information metric, `t` is the scale parameter for the persistence Fisher kernel. 
#' `verbose` determines if the time duration of the function call should be printed.
#'
#' @param g1 the first group of persistence diagrams, outputted from a TDA calculation or \code{\link{diagram_to_df}}.
#' @param g2 the second group of persistence diagrams, outputted from a TDA calculation or \code{\link{diagram_to_df}}.
#' @param dims a numeric vector of the homological dimensions in which the test is to be carried out, default c(0,1).
#' @param sigma the positive bandwith for the Fisher information metric, default 1.
#' @param t the positive scale for the persistence Fisher kernel, default 1.
#' @param verbose a boolean flag for if the time duration of the function call should be printed, default False.
#'
#' @return a list with dimensions used (named vector), test statistics in each dimension (named vector)
#'                   a p-value for each dimension (named vector) and the time duration of the function call.
#' @importFrom stats pgamma
#' @export
#' @examples
#'
#' # create two groups of persistence diagrams on 2D Gaussians using TDA
#' g1 <- lapply(X = 1:6,FUN = function(X){
#'
#' diag <- TDA::ripsDiag(data.frame(x = rnorm(100,mean = 0,sd = 1),
#' y = rnorm(100,mean = 0,sd = 1)),
#' maxscale = 1,
#' maxdimension = 1)
#' df <- diagram_to_df(d = diag)
#' return(df)
#'
#' })
#'
#' g2 <- lapply(X = 1:6,FUN = function(X){
#'
#' diag <- TDA::ripsDiag(data.frame(x = rnorm(100,mean = 0,sd = 1),
#' y = rnorm(100,mean = 0,sd = 1)),
#' maxscale = 1,
#' maxdimension = 1)
#' df <- diagram_to_df(d = diag)
#' return(df)
#'
#' })
#'
#' # do independence test with sigma = 1, t = 1, in dimensions 0 and 1, printing the time duration
#' ind_test = independence_test(g1,g2,verbose = TRUE)

independence_test <- function(g1,g2,dims = c(0,1),sigma = 1,t = 1,verbose = FALSE){
  
  # function to test whether or not two groups of persistence diagrams are independent
  # g1 and g2 are the groups of diagrams, either stored as lists or vectors
  # dims is a vector of desired homological dimensions
  # sigma is the positive bandwith for the Fisher information metric, 1 by default
  # t is the positive scale for the persistence Fisher kernel, 1 by default
  # verbose is either TRUE or FALSE (default), printing runtime of function call
  
  # set internal variables to NULL to avoid build issues
  r <- NULL
  
  # retrieve diagram groups
  diagram_groups <- list(g1,g2)
  
  # make sure there are at least two groups
  if(length(diagram_groups) < 2)
  {
    stop("At least two groups of persistence diagrams must be supplied.")
  }
  
  # check each diagram, converting each to a data frame and storing their indices in all the diagrams
  diagram_groups <- all_diagrams(diagram_groups,inference = "independence")
  
  # error check function parameters
  check_params(iterations = 10,p = 2,q = 2,dims,paired = T,distance = "fisher",sigma)
  if(is.null(t))
  {
    stop("t must not be NULL.")
  }
  if(!is.numeric(t) | length(t) > 1 | is.na(t) | is.nan(t) | t <= 0)
  {
    stop("t must be a positive number.")
  }
  
  # make sure that the two groups have the same number of diagrams
  if(length(g1) != length(g2))
  {
    stop("g1 and g2 must be the same length.")
  }
  
  # make sure that the two groups each have at least 5 elements
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
    K <- gram_matrix(diagrams = diagram_groups[[1]],dim = dim,t = t,sigma = sigma)
    L <- gram_matrix(diagrams = diagram_groups[[2]],dim = dim,t = t,sigma = sigma)

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
    
    # compute p-value
    p_val <- stats::pgamma(q = HSIC,rate = mu/var,shape = mu^2/var,lower.tail = F)
    p_value <- c(p_value,p_val)
    
  }
  
  # set up return lists
  names(test_statistics) <- as.character(dims)
  names(p_value) <- as.character(dims)
  runtime = Sys.time() - start_time
  
  results <- list(dimensions = dims,
                  test_statistic = test_statistics,
                  p_value = p_value,
                  run_time = runtime)
  
  if(verbose == T)
  {
    # print time duration of function call
    print(runtime)
  }
  
  return(results)
  
}
