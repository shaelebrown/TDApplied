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
#' so care should be taken to use the desired function when using TDApplied with TDAstats. If `dist_mats` is supplied
#' then the sum of the elements of `group_sizes` must equal the number of rows and columns of each of its elements.
#'
#' @param ... lists of persistence diagrams which are either the output of persistent homology calculations like ripsDiag/\code{\link[TDAstats]{calculate_homology}}/\code{\link{PyH}}, or \code{\link{diagram_to_df}}. Each list must contain at least 2 diagrams.
#' @param iterations the number of iterations for permuting group labels, default 20.
#' @param p a positive number representing the wasserstein power parameter, a number at least 1 (and Inf if using the bottleneck distance) and default 2.
#' @param q  a finite number at least 1 for exponentiation in the Turner loss function, default 2.
#' @param dims a non-negative integer vector of the homological dimensions in which the test is to be carried out, default c(0,1).
#' @param dist_mats an optional list of precomputed distances matrices, one for each dimension, where the rows and columns would correspond to the unlisted groups of diagrams (in order), default NULL. If not NULL then no lists of diagrams need to be supplied.
#' @param group_sizes a vector of group sizes, one for each group, when `dist_mats` is not NULL.
#' @param paired a boolean flag for if there is a second-order pairing between diagrams at the same index in different groups, default FALSE
#' @param distance a string which determines which type of distance calculation to carry out, either "wasserstein" (default) or "fisher".
#' @param sigma the positive bandwidth for the Fisher information metric, default NULL.
#' @param rho an optional positive number representing the heuristic for Fisher information metric approximation, see \code{\link{diagram_distance}}. Default NULL. If supplied, code execution is sequential.
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
#' @seealso \code{\link{independence_test}} for an inferential test of independence for two groups of persistence diagrams.
#' @references
#' Robinson T, Turner K (2017). "Hypothesis testing for topological data analysis." \url{https://link.springer.com/article/10.1007/s41468-017-0008-7}.
#' 
#' Abdallah H et al. (2021). "Statistical Inference for Persistent Homology applied to fMRI." \url{https://github.com/hassan-abdallah/Statistical_Inference_PH_fMRI/blob/main/Abdallah_et_al_Statistical_Inference_PH_fMRI.pdf}.
#' @examples
#'
#' if(require("TDAstats"))
#' {
#'   # create two groups of diagrams
#'   D1 <- TDAstats::calculate_homology(TDAstats::circle2d[sample(1:100,10),],
#'                                      dim = 0,threshold = 2)
#'   D2 <- TDAstats::calculate_homology(TDAstats::circle2d[sample(1:100,10),],
#'                                      dim = 0,threshold = 2)
#'   g1 <- list(D1,D2)
#'   g2 <- list(D1,D2)
#'
#'   # run test in dimension 0 with 1 iteration, note that the TDA package function
#'   # "permutation_test" can mask TDApplied's function, so we will specify explicitly
#'   # which function we are using
#'   perm_test <- TDApplied::permutation_test(g1,g2,iterations = 1,
#'                                            num_workers = 2,
#'                                            dims = c(0))
#'                                  
#'   # repeat with precomputed distance matrix, gives similar results
#'   # (same but the randomness of the permutations can give small differences)
#'   # just much faster
#'   D <- distance_matrix(diagrams = list(D1,D2,D1,D2),dim = 0,
#'                        num_workers = 2)
#'   perm_test <- TDApplied::permutation_test(dist_mats = list(D),group_sizes = c(2,2),
#'                                            dims = c(0))
#' }

permutation_test <- function(...,iterations = 20,p = 2,q = 2,dims = c(0,1),dist_mats = NULL,group_sizes = NULL,paired = FALSE,distance = "wasserstein",sigma = NULL,rho = NULL,num_workers = parallelly::availableCores(omit = 1),verbose = FALSE){

  # function to test whether or not multiple groups of persistence diagrams come from the same geometric process

  if(is.null(dist_mats))
  {
    # dist_mats not precomputed, retrieve diagram groups
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
    diagram_groups <- all_diagrams(diagram_groups,inference = "difference") 
    
    # make sure that if paired == T then all groups have the same number of diagrams
    if(paired)
    {
      if(length(unique(unlist(lapply(diagram_groups,FUN = length)))) != 1)
      {
        stop("If paired is true then all groups of diagrams must have the same number of elements.")
      }
    }
    
    # make distance matrix for all diagrams for each dimension, store total number of diagrams across
    # all groups in variable n
    n <- sum(unlist(lapply(diagram_groups,FUN = length)))
    dist_mats <- lapply(X = dims,FUN = function(X){return(matrix(data = -1,nrow = n,ncol = n))})
    
  }else
  {
    # distance matrices supplied
    
    # error check parameters
    check_param(group_sizes,param_name = "group_sizes",numeric = T,whole_numbers = T,multiple = T,finite = T,at_least_one = T)
    
    if(!is.list(dist_mats))
    {
      stop("dist_mats must be a list.")
    }
    if(length(dist_mats)!=length(dims))
    {
      stop("dist_mats must have the same length as dims.")
    }
    dist_mat_dims <- lapply(X = dist_mats,FUN = function(X){
      
      check_matrix(M = X,name = "distance matrices",type = "matrix")
      return(nrow(X))
      
    })
    if(length(unique(unlist(dist_mat_dims))) > 1)
    {
      stop("dist_mats must all have the same size.")
    }
    
    if(unique(unlist(dist_mat_dims)) != sum(group_sizes))
    {
      stop("The sum of group_sizes must equal the dimenion (i.e. number of rows/columns) of all dist_mats.")
    }
    
    if(paired)
    {
      if(length(unique(group_sizes)) != 1)
      {
        stop("If paired is true then all groups must have the same number of elements.")
      }
    }
    
    # split the vector 1:n (i.e. diagram indices) into groups according to group_sizes
    n <- sum(group_sizes)
    diagram_groups <- split(x = 1:n,f = unlist(lapply(X = 1:length(group_sizes),FUN = function(X){return(rep(X,length = group_sizes[[X]]))})))
  }

  # error check general parameters
  check_param("iterations",iterations,whole_numbers = T,at_least_one = T,numeric = T,finite = T)
  check_param("p",p,finite = F,at_least_one = T,numeric = T,multiple = F)
  check_param("q",q,at_least_one = T,numeric = T,multiple = F)
  check_param("dims",dims,multiple = T,whole_numbers = T,non_negative = T,positive = F,numeric = T,finite = T)
  check_param("paired",paired,numeric = F,multiple = F)
  check_param("distance",distance)
  if(distance == "fisher")
  {
    check_param("sigma",sigma,non_negative = T,positive = F,numeric = T,finite = T,multiple = F)
    if(!is.null(rho))
    {
      check_param("rho",rho,positive = T,non_negative = T,numeric = T,finite = T,multiple = F)
    }
  }
  check_param("num_workers",num_workers,whole_numbers = T,at_least_one = T,numeric = T,multiple = F)
  if(num_workers > parallelly::availableCores())
  {
    warning("num_workers is greater than the number of available cores - setting to maximum value.")
    num_workers <- parallelly::availableCores()
  }

  # time computations
  start_time <- Sys.time()

  # compute loss function on observed data and update dist_mats
  test_loss <- loss(diagram_groups = diagram_groups,dist_mats = dist_mats,dims = dims,p = p,q = q,distance = distance,sigma = sigma,rho = rho,num_workers = num_workers,group_sizes = group_sizes)
  dist_mats <- test_loss$dist_mats
  test_statistics <- test_loss$statistics

  # compute cumulative sum of diagram group sizes for inverting diagram indices
  csum_group_sizes <- c(0,cumsum(unlist(lapply(X = diagram_groups,FUN = length))))

  # create empty list in each dimension to store loss function statistic for each permutation
  perm_values <- lapply(X = dims,FUN = function(X){return(list())})

  # permute group labels and recalculate the loss function each time
  for(perm in 1:iterations)
  {

    if(is.list(diagram_groups[[1]][[1]]))
    {
      # distance matrices not supplied, use actual diagrams when permuting
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
            res[[length(res) + 1]] <- diagram_groups[[groups[[g]]]][[g]] # CHANGED THIS!!
          }
          return(res)})
      }
    }else
    {
      # distance matrices supplied, use indices when permuting
      if(paired == F)
      {
        # sample groups from their union, maintaining group sizes
        permuted_groups <- split(x = 1:n,f = sample(unlist(lapply(X = 1:length(group_sizes),FUN = function(X){return(rep(X,length = group_sizes[[X]]))})),size = n,replace = F))

      }else
      {
        # permute all the 1st diagrams between the groups, 2nd diagrams between the groups etc.
        perm <- lapply(X = 1:group_sizes[[1]],FUN = function(X){return(sample(1:length(group_sizes),size = length(group_sizes),replace = F))})
        
        permuted_groups <- lapply(X = 1:length(group_sizes),FUN = function(X){
          
          # get the Xth element of each list
          groups <- sapply(perm,"[[",X)
          
          # for each group g in the permuted group labels
          for(g in 1:length(groups))
          {
            if(groups[[g]] == 1)
            {
              ind <- 0
            }else
            {
              ind <- sum(group_sizes[1:(groups[[g]] - 1)])
            }
            groups[[g]] <- ind + g
          }
          return(groups)})
      }
    }
    
    # compute loss function, add to permutation values and updated distance matrices
    permuted_loss <- loss(diagram_groups = permuted_groups,dist_mats = dist_mats,dims = dims,p = p,q = q,distance = distance,sigma = sigma,rho = rho,num_workers = num_workers,group_sizes = group_sizes)
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
#' @param g1 the first group of persistence diagrams, where each diagram was either the output from a persistent homology calculation like ripsDiag/\code{\link[TDAstats]{calculate_homology}}/\code{\link{PyH}}, or \code{\link{diagram_to_df}}.
#' @param g2 the second group of persistence diagrams, where each diagram was either the output from a persistent homology calculation like ripsDiag/\code{\link[TDAstats]{calculate_homology}}/\code{\link{PyH}}, or \code{\link{diagram_to_df}}.
#' @param dims a non-negative integer vector of the homological dimensions in which the test is to be carried out, default c(0,1).
#' @param sigma a positive number representing the bandwidth for the Fisher information metric, default 1.
#' @param rho an optional positive number representing the heuristic for Fisher information metric approximation, see \code{\link{diagram_distance}}. Default NULL. If supplied, calculation of Gram matrices is sequential.
#' @param t a positive number representing the scale for the persistence Fisher kernel, default 1.
#' @param num_workers the number of cores used for parallel computation, default is one less than the number of cores on the machine.
#' @param Ks an optional list of precomputed Gram matrices for the first group of diagrams, with one element for each dimension. If not NULL and `Ls` is not NULL then `g1` and `g2` do not need to be supplied.
#' @param Ls an optional list of precomputed Gram matrices for the second group of diagrams, with one element for each dimension. If not NULL and `Ks` is not NULL then `g1` and `g2` do not need to be supplied.
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
#' @seealso \code{\link{permutation_test}} for an inferential group difference test for groups of persistence diagrams.
#' @references
#' Gretton A et al. (2007). "A Kernel Statistical Test of Independence." \url{https://proceedings.neurips.cc/paper/2007/file/d5cfead94f5350c12c322b5b664544c1-Paper.pdf}.
#' @examples
#'
#' if(require("TDAstats"))
#' {
#'   # create two independent groups of diagrams of length 6, which
#'   # is the minimum length
#'   D1 <- TDAstats::calculate_homology(TDAstats::circle2d[sample(1:100,10),],
#'                                      dim = 0,threshold = 2)
#'   D2 <- TDAstats::calculate_homology(TDAstats::circle2d[sample(1:100,10),],
#'                                      dim = 0,threshold = 2)
#'   g1 <- list(D1,D2,D2,D2,D2,D2)
#'   g2 <- list(D2,D1,D1,D1,D1,D1)
#' 
#'   # do independence test with sigma = t = 1 in dimension 0, using
#'   # precomputed Gram matrices
#'   K = gram_matrix(diagrams = g1,dim = 0,t = 1,sigma = 1,num_workers = 2)
#'   L = gram_matrix(diagrams = g2,dim = 0,t = 1,sigma = 1,num_workers = 2)
#'   indep_test <- independence_test(Ks = list(K),Ls = list(L),dims = c(0))
#'   
#' }

independence_test <- function(g1,g2,dims = c(0,1),sigma = 1,rho = NULL,t = 1,num_workers = parallelly::availableCores(omit = 1),verbose = FALSE,Ks = NULL,Ls = NULL){
  
  # function to test whether or not two groups of persistence diagrams are independent
  
  # set internal variables to NULL to avoid build issues
  r <- NULL
  
  if(is.null(Ks) | is.null(Ls))
  {
    # Gram matrices are not supplied, retrieve diagram groups
    diagram_groups <- list(g1,g2)
    
    # make sure there are at least two groups
    check_param("diagram_groups",diagram_groups,min_length = 2)
    
    # check each diagram, converting each to a data frame
    all_diagrams(diagram_groups,inference = "independence")
    diagram_groups <- all_diagrams(diagram_groups,inference = "independence") 
    
    # check params
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
    
  }else
  {
    # Gram matrices are supplied
    
    # set g1 and g2 to NULL to avoid issues
    if(missing(g1))
    {
      g1 <- NULL
    }
    if(missing(g2))
    {
      g2 <- NULL
    }
    
    # error check Ks and Ls
    if(!is.list(Ks) | !is.list(Ls))
    {
      stop("Ks and Ls must be lists.")
    }
    
    if(length(Ks) != length(dims) | length(Ls) != length(dims))
    {
      stop("Ks and Ls must have the same length as dims.")
    }
    
    Ks_dims <- lapply(X = Ks,FUN = function(X){
      
      check_matrix(M = X,name = "Gram matrices")
      if(nrow(X) < 6)
      {
        stop("Gram matrices must have at least 6 rows and columns.")
      }
      return(nrow(X))
      
    })
    
    Ls_dims <- lapply(X = Ls,FUN = function(X){
      
      check_matrix(M = X,name = "Gram matrices")
      if(nrow(X) < 6)
      {
        stop("Gram matrices must have at least 6 rows and columns.")
      }
      return(nrow(X))
      
    })
    
    if(length(unique(c(unlist(Ks_dims),unlist(Ls_dims)))) > 1)
    {
      stop("Gram matrices in Ks and Ls must all have the same dimensions.")
    }
    
    m <- unique(c(unlist(Ks_dims),unlist(Ls_dims)))[[1]]
    
  }
  
  # error check dims
  check_param("dims",dims,multiple = T,whole_numbers = T,non_negative = T,positive = F,numeric = T,finite = T)
  
  # error check rho
  if(!is.null(rho))
  {
    check_param("rho",rho,positive = T,non_negative = T,numeric = T,finite = T,multiple = F)
  }
  
  # time computations
  start_time <- Sys.time()
  
  # conduct test in each dimension
  test_statistics <- c()
  p_value <- c()
  
  for(dim in dims)
  {
    
    # compute test statistics
    if(is.null(g1))
    {
      K <- Ks[[which(dims == dim)[[1]]]]
      L <- Ls[[which(dims == dim)[[1]]]]
    }else
    {
      K <- gram_matrix(diagrams = diagram_groups[[1]],dim = dim,t = t,sigma = sigma,rho = rho,num_workers = num_workers)
      L <- gram_matrix(diagrams = diagram_groups[[2]],dim = dim,t = t,sigma = sigma,rho = rho,num_workers = num_workers) 
    }

    H <- matrix(data = -1/m,nrow = m,ncol = m)
    diag(H) <- rep((m-1)/m,nrow(H))
    
    HSIC <- sum(diag(K %*% H %*% L %*% H))/(m^2) # normalized trace
    test_statistics <- c(test_statistics,HSIC)
    
    # compute null distribution parameters
    mu_x_sq <- mean(K[upper.tri(K)])
    mu_y_sq <- mean(L[upper.tri(L)])
    mu <- (1 + mu_x_sq*mu_y_sq - mu_x_sq - mu_y_sq)/m
    if(mu <= 0)
    {
      p_val <- 1
    }else
    {
      B <- (H %*% K %*% H) * (H %*% L %*% H)
      B <- B * B
      diag(B) <- rep(0,m)
      var <- 2*(m-4)*(m-5)*(sum(colSums(B)))/((m-3)*(m-2)*(m-1)*m)
      # error check var
      if(var == 0)
      {
        stop("A zero variance was calculated, please make sure that both g1 and g2 contain at least 2 distinct diagrams.")
      }
      
      # compute p-value
      p_val <- stats::pgamma(q = HSIC,rate = mu/var,shape = mu^2/var,lower.tail = F)
    }
    
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

#### MODEL INFERENCE ####
#' Model inference with permutation test.
#'
#' An inference procedure to determine if two datasets were unlikely to be generated by the same process (i.e. if
#' the persistence diagram of one dataset is a good model of the persistence diagram of the other dataset). 
#' 
#' Inference is carried out by generating bootstrap resampled persistence diagrams from the two datasets and carrying out a permutation test
#' on the resulting two groups. A small p-value in a certain dimension suggests that the datasets are not good models of each other. `samp` should
#' only be provided when `paired`is TRUE in order to generate the same row samplings of `D1` and `D2` for the bootstrapped persistence diagrams.
#' This makes a paired permutation test more appropriate, which has higher statistical power for detecting topological differences. See the examples
#' for how to properly supply `samp`.
#'
#' @param D1 the first dataset (a data frame).
#' @param D2 the second dataset (a data frame).
#' @param iterations the number of iterations for permuting group labels, default 20.
#' @param num_samples the number of bootstrap iterations, default 30.
#' @param dims a non-negative integer vector of the homological dimensions in which the test is to be carried out, default c(0,1).
#' @param samp an optional list of row-number samples of `D1`, default NULL. See details and examples for more information. Ignored when `paired` is FALSE.
#' @param paired a boolean flag for if there is a second-order pairing between diagrams at the same index in different groups, default FALSE.
#' @param num_workers the number of cores used for parallel computation, default is one less than the number of cores on the machine.
#' @param verbose a boolean flag for if the time duration of the function call should be printed, default FALSE
#' @param FUN_boot a string representing the persistent homology function to use for calculating the bootstrapped persistence diagrams, either
#' 'calculate_homology' (the default), 'PyH' or 'ripsDiag'.
#' @param thresh the positive numeric maximum radius of the Vietoris-Rips filtration.
#' @param distance_mat a boolean representing if `X` is a distance matrix (TRUE) or not (FALSE, default).
#' dimensions together (TRUE, the default) or if one threshold should be calculated for each dimension separately (FALSE).
#' @param ripser the imported ripser module when `FUN_boot` is `PyH`.
#' @param return_diagrams whether or not to return the two lists of bootstrapped persistence diagrams, default FALSE.
#'
#' @return a list which contains the output of the call to \code{\link{permutation_test}} and the two groups of bootstrapped
#' persistence diagrams if desired, in entries called `diagrams1` and `diagrams2`.
#' 
#' @importFrom parallelly availableCores
#' @export
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @seealso \code{\link{permutation_test}} for an inferential group difference test for groups of persistence diagrams and \code{\link{bootstrap_persistence_thresholds}} for computing confidence sets for persistence diagrams.
#' @references
#' Robinson T, Turner K (2017). "Hypothesis testing for topological data analysis." \url{https://link.springer.com/article/10.1007/s41468-017-0008-7}.
#' 
#' Chazal F et al (2017). "Robust Topological Inference: Distance to a Measure and Kernel Distance." \url{https://www.jmlr.org/papers/volume18/15-484/15-484.pdf}.
#' 
#' Abdallah H et al. (2021). "Statistical Inference for Persistent Homology applied to fMRI." \url{https://github.com/hassan-abdallah/Statistical_Inference_PH_fMRI/blob/main/Abdallah_et_al_Statistical_Inference_PH_fMRI.pdf}.
#' @examples
#'
#' if(require("TDAstats"))
#' {
#'   # create two datasets
#'   D1 <- TDAstats::calculate_homology(TDAstats::circle2d[sample(1:100,10),],
#'                                      dim = 0,threshold = 2)
#'   D2 <- TDAstats::calculate_homology(TDAstats::circle2d[sample(1:100,10),],
#'                                      dim = 0,threshold = 2)
#' 
#'   # do model inference test with 1 iteration (for speed, more
#'   # iterations should be used in practice)
#'   model_test <- permutation_model_inference(D1, D2, iterations = 1,
#'                                             thresh = 1.75,num_samples = 3,
#'                                             num_workers = 2L)
#'   
#'   # with more iterations, p-values show a difference in the 
#'   # clustering of points but not in the arrangement of loops
#'   model_test$p_values
#'   
#'   # to supply samp, when we believe there is a correspondence between
#'   # the rows in D1 and the rows in D2
#'   # note that the number of entries of samp (3 in this case) must
#'   # match the num_samples parameter to the function call
#'   samp <- lapply(X = 1:3,FUN = function(X){
#' 
#'            return(unique(sample(1:nrow(D1),size = nrow(D1),replace = TRUE)))
#' 
#'           })
#'   
#'   # model inference will theoretically have higher power now for a
#'   # paired test 
#'   model_test2 <- permutation_model_inference(D1, D2, iterations = 1,
#'                                              thresh = 1.75,num_samples = 3,
#'                                              paired = TRUE,samp = samp,
#'                                              num_workers = 2L)
#'   model_test2$p_values
#' }
permutation_model_inference <- function(D1, D2, iterations, num_samples, dims = c(0,1),samp = NULL,paired = F, num_workers = parallelly::availableCores(omit = 1), verbose = F,FUN_boot = "calculate_homology",thresh,distance_mat = FALSE,ripser = NULL,return_diagrams = FALSE){
  
  # do error checks for parameters
  if(!inherits(D1,"data.frame") & !inherits(D1,"matrix"))
  {
    stop("D1 must either be a dataframe or a matrix.")
  }
  if(nrow(D1) < 2 | ncol(D1) < 1)
  {
    stop("D1 must have at least two rows and one column.")
  }
  if(length(which(stats::complete.cases(D1) == F)) > 0)
  {
    stop("D1 must not contain any missing values.")
  }
  if(distance_mat == T & (ncol(D1) != nrow(D1) | !inherits(D1,"matrix")))
  {
    stop("if distance_mat is TRUE then D1 must be a square matrix.")
  }
  if((inherits(D1,"matrix") & !inherits(D1[1,1],"numeric")) | (inherits(D1,"data.frame") & length(which(unlist(lapply(D1,is.numeric)))) < ncol(D1)))
  {
    stop("D1 must have only numeric entries.")
  }
  
  if(!inherits(D2,"data.frame") & !inherits(D2,"matrix"))
  {
    stop("D2 must either be a dataframe or a matrix.")
  }
  if(nrow(D2) < 2 | ncol(D2) < 1)
  {
    stop("D2 must have at least two rows and one column.")
  }
  if(length(which(stats::complete.cases(D2) == F)) > 0)
  {
    stop("D2 must not contain any missing values.")
  }
  if(distance_mat == T & (ncol(D2) != nrow(D2) | !inherits(D2,"matrix")))
  {
    stop("if distance_mat is TRUE then D2 must be a square matrix.")
  }
  if((inherits(D2,"matrix") & !inherits(D2[1,1],"numeric")) | (inherits(D2,"data.frame") & length(which(unlist(lapply(D2,is.numeric)))) < ncol(D2)))
  {
    stop("D2 must have only numeric entries.")
  }
  
  check_param(param_name = 'num_samples',param = num_samples,numeric = T,whole_numbers = T,multiple = F,finite = T,positive = T)
  
  if(!is.null(samp))
  {
    if(!is.list(samp))
    {
      stop("samp must be a list.")
    }
    if(unique(unlist(lapply(samp, class))) != "integer")
    {
      stop("samp must be a list of integer vectors (i.e. row number samples).")
    }
    if(max(unlist(samp)) > nrow(D1))
    {
      stop("No entry of samp should be larger than the number of rows of D1.")
    }
  }
  
  if(is.null(paired))
  {
    stop("paired must not be NULL.")
  }
  if(length(paired) > 1 | !inherits(paired,"logical"))
  {
    stop("paired must be a single logical (i.e. T or F).")
  }
  if(is.na(paired) | is.nan(paired) )
  {
    stop("paired must not be NA/NAN.")
  }
  
  if(is.null(verbose))
  {
    stop("verbose must not be NULL.")
  }
  if(length(verbose) > 1 | !inherits(verbose,"logical"))
  {
    stop("verbose must be a single logical (i.e. T or F).")
  }
  if(is.na(verbose) | is.nan(verbose) )
  {
    stop("verbose must not be NA/NAN.")
  }
  
  if(is.null(return_diagrams))
  {
    stop("return_diagrams must not be NULL.")
  }
  if(length(return_diagrams) > 1 | !inherits(return_diagrams,"logical"))
  {
    stop("return_diagrams must be a single logical (i.e. T or F).")
  }
  if(is.na(return_diagrams) | is.nan(return_diagrams) )
  {
    stop("return_diagrams must not be NA/NAN.")
  }
  
  # generate bootstrap samples
  # if paired then generate shared samples, otherwise different samples
  if(is.null(samp))
  {
    samp <- lapply(X = 1:num_samples,FUN = function(X){
      
      return(unique(sample(1:nrow(D1),size = nrow(D1),replace = T)))
      
    }) 
  }
  if(paired)
  {
    samp2 <- samp
  }else
  {
    samp2 <- lapply(X = 1:num_samples,FUN = function(X){
      
      return(unique(sample(1:nrow(D2),size = nrow(D2),replace = T)))
      
    })
  }
  
  # use bootstrap samples to get bootstrapped diagrams
  if(verbose)
  {
    message(paste0("Starting to calculate bootstrapped diagrams at ",Sys.time()))
  }
  bootstrapped_diagrams1 <- bootstrap_persistence_thresholds(X = D1,maxdim = max(dims),distance_mat = distance_mat,thresh = thresh,num_workers = num_workers,ripser = ripser,FUN_boot = FUN_boot,bootstrap_samples = samp,num_samples = num_samples,return_diag = T)
  diag1 <- bootstrapped_diagrams1$diag
  bootstrapped_diagrams1 <- bootstrapped_diagrams1$diagrams
  bootstrapped_diagrams1 <- lapply(X = 1:num_samples,FUN = function(X){
    
    return(data.frame(dimension = bootstrapped_diagrams1[[(X-1)*3 + 1]],birth = bootstrapped_diagrams1[[(X-1)*3 + 2]],death = bootstrapped_diagrams1[[X*3]]))
    
  })
  bootstrapped_diagrams2 <- bootstrap_persistence_thresholds(X = D2,maxdim = max(dims),distance_mat = distance_mat,thresh = thresh,num_workers = num_workers,ripser = ripser,FUN_boot = FUN_boot,bootstrap_samples = samp2,num_samples = num_samples,return_diag = T)
  diag2 <- bootstrapped_diagrams2$diag
  bootstrapped_diagrams2 <- bootstrapped_diagrams2$diagrams
  bootstrapped_diagrams2 <- lapply(X = 1:num_samples,FUN = function(X){
    
    return(data.frame(dimension = bootstrapped_diagrams2[[(X-1)*3 + 1]],birth = bootstrapped_diagrams2[[(X-1)*3 + 2]],death = bootstrapped_diagrams2[[X*3]]))
    
  })
  
  # carry out permutation test on pairs of same subsetted diagrams
  if(verbose)
  {
    message(paste0("Starting permutation test at ",Sys.time()))
  }
  res <- permutation_test(bootstrapped_diagrams1,bootstrapped_diagrams2,iterations = iterations,dims = dims,paired = FALSE,num_workers = num_workers,p = Inf,q = 2)
  
  if(verbose)
  {
    message(paste0("Finished permutation test at ",Sys.time()))
  }
  if(return_diagrams)
  {
    res$diagrams1 <- bootstrapped_diagrams1
    res$diagrams2 <- bootstrapped_diagrams2
  }
  return(res)
  
}

#### UNIVERSAL NULL ####
#' Filtering topological features with the universal null distribution.
#'
#' An inference procedure to determine which topological features (if any) of a datasets are likely signal (i.e. significant)
#' vs noise (not). 
#' 
#' For each feature in a diagram we compute its persistence ratio \eqn{\pi = death/birth}, and a
#' test statistic \eqn{A log log \pi + B} (where \eqn{A} and \eqn{B} are constants). This statistic is compared to a left-skewed Gumbel distribution
#' to get a p-value. A Bonferroni correction is applied to all the p-values across all features, so when `return_pvals` is TRUE a list of 
#' p-value thresholds is also returned, one for each dimension, which is `alpha` divided by the number of features in that dimension.
#' If desired, infinite cycles (i.e. cycles whose death value is equal to the maximum distance threshold parameter for the persistent homology calculation) 
#' can be anaylzed for significance by determining their minimum distance thresholds where they might be significant (using the Gumbel distribution again),
#' calculating the persistence diagram up to those thresholds and seeing if they are still infinite (i.e. significant) or not.
#' This function is significantly faster than the \code{\link{bootstrap_persistence_thresholds}} function. Note that the `calculate_homology`
#' function does not seem to store infinite cycles (i.e. cycles that have death value equal to `thresh`).
#'
#' @param X the input dataset, must either be a matrix or data frame.
#' @param FUN_diag a string representing the persistent homology function to use for calculating the full persistence diagram, either
#' 'calculate_homology' (the default), 'PyH' or 'ripsDiag'.
#' @param maxdim the integer maximum homological dimension for persistent homology, default 0.
#' @param thresh the positive numeric maximum radius of the Vietoris-Rips filtration.
#' @param distance_mat a boolean representing if `X` is a distance matrix (TRUE) or not (FALSE, default).
#' dimensions together (TRUE, the default) or if one threshold should be calculated for each dimension separately (FALSE).
#' @param ripser the imported ripser module when `FUN_diag` is `PyH`.
#' @param ignore_infinite_cluster a boolean indicating whether or not to ignore the infinitely lived cluster when `FUN_diag` is `PyH`. If infinite cycle inference is to be performed,
#' this parameter should be set to FALSE.
#' @param calculate_representatives a boolean representing whether to calculate representative (co)cycles, default FALSE. Note that representatives cant be
#' calculated when using the 'calculate_homology' function. Note that representatives cannot be computed for (significant) infinite cycles.
#' @param alpha the type-1 error threshold, default 0.05.
#' @param return_pvals a boolean representing whether or not to return p-values for features in the subsetted diagram as well as a list of p-value thresholds, default FALSE.
#' Infinite cycles that are significant (see below) will have p-value NA in this list, as the true value is unknown but less than its dimension's p-value threshold.
#' @param infinite_cycle_inference a boolean representing whether or not to perform inference for features with infinite (i.e. `thresh`) death values, default FALSE. If `FUN_diag` is `calculate_homology` (the
#' default) then no infinite cycles will be returned by the persistent homology calculation at all.
#' @return a list containing the full persistence diagram, the subsetted diagram, representatives and/or subsetted representatives if desired, the p-values of subsetted features and the Bonferroni p-value thresholds in each dimension if desired. 
#' @export
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @references
#' Bobrowski O, Skraba P (2023). "A universal null-distribution for topological data analysis." \url{https://www.nature.com/articles/s41598-023-37842-2}.
#' @examples
#'
#' if(require("TDA"))
#' {
#'   # create dataset
#'   theta <- runif(n = 100,min = 0,max = 2*pi)
#'   x <- cos(theta)
#'   y <- sin(theta)
#'   circ <- data.frame(x = x,y = y)
#'
#'   # add noise
#'   x_noise <- -0.1 + 0.2*stats::runif(n = 100)
#'   y_noise <- -0.1 + 0.2*stats::runif(n = 100)
#'   circ$x <- circ$x + x_noise
#'   circ$y <- circ$y + y_noise
#'
#'   # determine significant topological features
#'   library(TDA)
#'   res <- universal_null(circ, thresh = 2,alpha = 0.1,return_pvals = TRUE,FUN_diag = "ripsDiag")
#'   res$subsetted_diag
#'   res$pvals
#'   res$alpha_thresh
#'
#'   # at a lower threshold we can check for 
#'   # infinite cycles
#'   res2 <- universal_null(circ, thresh = 1.1, 
#'                          infinite_cycle_inference = TRUE,
#'                          alpha = 0.1,
#'                          FUN_diag = "ripsDiag")
#'   res2$subsetted_diag
#' }
universal_null <- function(X,FUN_diag = "calculate_homology",maxdim = 1,thresh,distance_mat = FALSE,ripser = NULL,ignore_infinite_cluster = TRUE,calculate_representatives = FALSE,alpha = 0.05,return_pvals = FALSE,infinite_cycle_inference = FALSE){
  
  # function for thresholding the points computed in a diagram
  # based on the universal null distribution
  
  # error check parameters
  if(is.null(distance_mat))
  {
    stop("distance_mat must not be NULL.")
  }
  if(length(distance_mat) > 1 | !inherits(distance_mat,"logical"))
  {
    stop("distance_mat must be a single logical (i.e. T or F).")
  }
  if(is.na(distance_mat) | is.nan(distance_mat) )
  {
    stop("distance_mat must not be NA/NAN.")
  }
  
  if(is.null(return_pvals))
  {
    stop("return_pvals must not be NULL.")
  }
  if(length(return_pvals) > 1 | !inherits(return_pvals,"logical"))
  {
    stop("return_pvals must be a single logical (i.e. T or F).")
  }
  if(is.na(return_pvals) | is.nan(return_pvals) )
  {
    stop("return_pvals must not be NA/NAN.")
  }
  
  check_param(param = maxdim,param_name = "maxdim",numeric = T,whole_numbers = T,multiple = F,finite = T,non_negative = T)
  check_param(param = thresh,param_name = "thresh",numeric = T,whole_numbers = F,multiple = F,finite = T,non_negative = T,positive = T)
  
  if(!inherits(X,"data.frame") & !inherits(X,"matrix"))
  {
    stop("X must either be a dataframe or a matrix.")
  }
  if(nrow(X) < 2 | ncol(X) < 1)
  {
    stop("X must have at least two rows and one column.")
  }
  if(length(which(stats::complete.cases(X) == F)) > 0)
  {
    stop("X must not contain any missing values.")
  }
  if(distance_mat == T & (ncol(X) != nrow(X) | !inherits(X,"matrix")))
  {
    stop("if distance_mat is TRUE then X must be a square matrix.")
  }
  if((inherits(X,"matrix") & !inherits(X[1,1],"numeric")) | (inherits(X,"data.frame") & length(which(unlist(lapply(X,is.numeric)))) < ncol(X)))
  {
    stop("X must have only numeric entries.")
  }
  
  if(is.null(FUN_diag))
  {
    stop("FUN_diag must not be NULL.")
  }
  if(length(FUN_diag) > 1 | !inherits(FUN_diag,"character"))
  {
    stop("FUN_diag must be a single string.")
  }
  if(is.na(FUN_diag) | is.nan(FUN_diag) )
  {
    stop("FUN_diag must not be NA/NAN.")
  }
  if(FUN_diag %in% c("calculate_homology","ripsDiag","PyH") == F)
  {
    stop("FUN_diag must be either \'calculate_homology\', \'PyH\' or \'ripsDiag\'.")
  }
  if(FUN_diag == 'calculate_homology')
  {
    if(requireNamespace("TDAstats",quietly = T) == F)
    {
      stop("To use the \'calculate_homology\' function the package TDAstats must be installed.")
    }
  }
  if(FUN_diag == 'ripsDiag')
  {
    tryCatch(expr = {ripsDiag <- get("ripsDiag",envir = globalenv())},
             error = function(e){
               
               stop("To use the \'ripsDiag\' function the package TDA must be attached in the global environment.")
               
             })
  }
  
  if(FUN_diag == "PyH")
  {
    
    if(is.null(ignore_infinite_cluster))
    {
      stop("ignore_infinite_cluster must not be NULL.")
    }
    if(length(ignore_infinite_cluster) > 1 | !inherits(ignore_infinite_cluster,"logical"))
    {
      stop("ignore_infinite_cluster must be a single logical (i.e. T or F).")
    }
    if(is.na(ignore_infinite_cluster) | is.nan(ignore_infinite_cluster) )
    {
      stop("ignore_infinite_cluster must not be NA/NAN.")
    }
    
    check_ripser(ripser)
    
  }
  
  if(is.null(calculate_representatives))
  {
    stop("calculate_representatives must not be NULL.")
  }
  if(length(calculate_representatives) > 1 | !inherits(calculate_representatives,"logical"))
  {
    stop("calculate_representatives must be a single boolean value.")
  }
  if(is.na(calculate_representatives) | is.nan(calculate_representatives) )
  {
    stop("calculate_representatives must not be NA/NAN.")
  }
  check_param(param_name = "alpha",param = alpha,numeric = T,whole_numbers = F,multiple = F,finite = T,non_negative = T)
  if(alpha <= 0 | alpha >= 1)
  {
    stop("alpha must be between 0 and 1 (non-inclusive).")
  }
  if(is.null(infinite_cycle_inference))
  {
    stop("infinite_cycle_inference must not be NULL.")
  }
  if(length(infinite_cycle_inference) > 1 | !inherits(infinite_cycle_inference,"logical"))
  {
    stop("infinite_cycle_inference must be a single boolean value.")
  }
  if(is.na(infinite_cycle_inference) | is.nan(infinite_cycle_inference) )
  {
    stop("infinite_cycle_inference must not be NA/NAN.")
  }
  
  # first calculate the "real" persistence diagram, storing representative cycles if need be
  if(FUN_diag == "PyH")
  {
    diag <- PyH(X = X,maxdim = maxdim,thresh = thresh,distance_mat = distance_mat,ripser = ripser,ignore_infinite_cluster = ignore_infinite_cluster,calculate_representatives = calculate_representatives)
    if(calculate_representatives == T)
    {
      representatives = diag$representatives
      diag <- diag$diagram
    }
  }
  if(FUN_diag == "calculate_homology")
  {
    diag <- diagram_to_df(TDAstats::calculate_homology(mat = X,dim = maxdim,threshold = thresh,format = ifelse(test = distance_mat == T,yes = "distmat",no = "cloud")))
    representatives <- NULL
  }
  if(FUN_diag == "ripsDiag")
  {
    diag <- ripsDiag(X = X,maxdimension = maxdim,maxscale = thresh,dist = ifelse(test = distance_mat == F,yes = "euclidean",no = "arbitrary"),library = "dionysus",location = calculate_representatives,printProgress = F)
    if(calculate_representatives == T)
    {
      representatives <- diag$cycleLocation
    }
    diag <- diagram_to_df(diag)
  }
  
  # make return list
  ret_list <- list(diag = diag)
  if(calculate_representatives)
  {
    ret_list$representatives <- representatives
  }
  
  # make sure there is work to be done
  if(max(diag$dimension) > 0)
  {
    # subset diagram and representatives
    diag_highdim <- diag[which(diag$dimension > 0),]
    if(calculate_representatives)
    {
      if(FUN_diag == "calculate_homology")
      {
        representatives_highdim <- NULL
      }
      if(FUN_diag == "ripsDiag")
      {
        representatives_highdim <- representatives[which(diag$dimension > 0)]
      }
      if(FUN_diag == "PyH")
      {
        representatives_highdim <- representatives[2:(max(diag$dimension) + 1)]
      }
    }
    # compute test statistics
    dims <- c(1:maxdim)
    A <- 1 # for VR filtrations
    lambda <- -digamma(1)
    pi <- unlist(lapply(X = 1:nrow(diag_highdim),FUN = function(X){
      
      return(return(diag_highdim[X,3L]/diag_highdim[X,2L]))
      
    }))
    loglog_pi <- log(log(pi))
    L_bar <- unlist(lapply(X = dims, FUN = function(X){
      
      dim_inds <- which(diag_highdim$dimension == X)
      return(mean(loglog_pi[dim_inds]))
      
    }))
    B <- -1*lambda - A*L_bar
    test_stats <- unlist(lapply(X = 1:length(loglog_pi),FUN = function(X){
      
      return(A*loglog_pi[[X]] + B[[diag_highdim[X, 1L]]])
    
    }))
    
    # compute p-values
    pvals <- exp(-1*exp(test_stats))
    
    # get Bonferroni thresholds in each dimension
    alpha_thresh <- unlist(lapply(dims,FUN = function(X){
      
      L <- length(which(diag_highdim$dimension == X))
      if(L == 0)
      {
        return(Inf)
      }
      return(alpha/L)
      
    }))
    
    # subset in each dimension
    inds <- c()
    non_inf_inds <- c()
    for(d in 1:maxdim)
    {
      inds <- c(inds,which(diag_highdim$dimension == d & pvals < alpha_thresh[[d]])) 
      non_inf_inds <- c(non_inf_inds,which(diag_highdim$dimension == d & pvals < alpha_thresh[[d]] & diag_highdim$death < thresh))
    }
    ret_list$subsetted_diag <- diag_highdim[inds,]
    
    if(return_pvals)
    {
      if(length(non_inf_inds) > 0)
      {
        ret_list$pvals <- pvals[inds]
      }else
      {
        ret_list$pvals <- numeric()
      }
      ret_list$alpha_thresh <- alpha_thresh
    }
    
    # subset representatives if desired
    # offset inds
    if(length(non_inf_inds) > 0)
    {
      non_inf_inds <- non_inf_inds + length(which(diag$dimension == 0))
      if(calculate_representatives == T & FUN_diag == "ripsDiag")
      {
        ret_list$subsetted_representatives = representatives[non_inf_inds]
      }
      if(calculate_representatives == T & FUN_diag == "PyH")
      {
        ret_list$subsetted_representatives = representatives
        if(maxdim > 0)
        {
          for(d in 1:maxdim)
          {
            if(length(which(diag$dimension == d)) > 0)
            {
              min_ind <- min(which(diag$dimension == d))
              max_ind <- max(which(diag$dimension == d))
              ret_list$subsetted_representatives[[d + 1]] <- ret_list$subsetted_representatives[[d + 1]][non_inf_inds[which(non_inf_inds %in% c(min_ind:max_ind))] - min_ind + 1]
            }else
            {
              ret_list$subsetted_representatives[[d + 1]] <- list()
            }
          }
        }
      }
    }
    
    if(infinite_cycle_inference)
    {
      # perform inference for infinite cycles
      infinite_inds <- which(diag_highdim$death == thresh)
      if(length(infinite_inds) > 0)
      {
        diag_infinite <- diag_highdim[infinite_inds,]
        pvals_infinite <- pvals[infinite_inds]
        infinite_cycle_inference <- data.frame(dimension = integer(), birth = numeric(), death = numeric())
        unknown_death_cycles <- data.frame(dimension = integer(), birth = numeric(), death = numeric())
        for(i in 1:nrow(diag_infinite))
        {
          if(pvals_infinite[[i]] >= alpha_thresh[[diag_infinite[i,1L]]])
          {
            unknown_death_cycles <- rbind(unknown_death_cycles, data.frame(dimension = diag_infinite[i,1L], birth = diag_infinite[i,2L], death = diag_infinite[i,3L]))
          }
        }
        alpha_cutoffs <- unlist(lapply(1:length(alpha_thresh), FUN = function(X){
          
          quant <- log(log(1/alpha_thresh[[X]]))
          pi_min <- exp(exp((quant - B[[X]])/A))
          return(pi_min)
          
        }))
        # take lowest birth value point, determine its candidate death radius and test
        while(nrow(unknown_death_cycles) > 0)
        {
          min_birth_ind <- which(unknown_death_cycles$birth == min(unknown_death_cycles$birth))[[1]]
          min_birth_pt <- unknown_death_cycles[min_birth_ind,]
          unknown_death_cycles <- unknown_death_cycles[setdiff(1:nrow(unknown_death_cycles), min_birth_ind),]
          new_death <- alpha_cutoffs[[min_birth_pt[[1]]]]*min_birth_pt[[2]]
          if(FUN_diag == "PyH")
          {
            new_diag <- PyH(X = X,maxdim = maxdim,thresh = new_death,distance_mat = distance_mat,ripser = ripser,ignore_infinite_cluster = ignore_infinite_cluster,calculate_representatives = FALSE)
          }
          if(FUN_diag == "calculate_homology")
          {
            new_diag <- diagram_to_df(TDAstats::calculate_homology(mat = X,dim = maxdim,threshold = new_death,format = ifelse(test = distance_mat == T,yes = "distmat",no = "cloud")))
          }
          if(FUN_diag == "ripsDiag")
          {
            new_diag <- diagram_to_df(ripsDiag(X = X,maxdimension = maxdim,maxscale = new_death,dist = ifelse(test = distance_mat == F,yes = "euclidean",no = "arbitrary"),library = "dionysus",location = FALSE,printProgress = FALSE))
          }
          # check if feature with same birth value is still infinite or not
          # if so then add to subsetted_diagram
          if(length(which(new_diag$dimension == min_birth_pt[[1]] & new_diag$birth == min_birth_pt[[2]] & new_diag$death == new_death)) > 0)
          {
            infinite_cycle_inference <- rbind(infinite_cycle_inference, data.frame(dimension = min_birth_pt[[1]],birth = min_birth_pt[[2]],death = min_birth_pt[[3]]))
          }
        }
        # update subsetted diagram and pvalues
        ret_list$subsetted_diag <- rbind(ret_list$subsetted_diag,infinite_cycle_inference)
        ret_list$pvals <- c(ret_list$pvals,NA)
      }
      
    }
  }
  else
  {
    ret_list$subsetted_diag <- diagram_to_df(data.frame(dimension = integer(),birth = numeric(),death = numeric()))
    if(calculate_representatives)
    {
      ret_list$subsetted_representatives <- list()
    }
    if(return_pvals)
    {
      ret_list$pvals <- numeric()
      ret_list$alpha_thresh <- numeric()
    }
  }
  
  return(ret_list)
  
}


