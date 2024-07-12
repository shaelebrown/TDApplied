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

#### NEW STUFF ####
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
  
  # # now calculate distance metrics and confidence intervals if desired
  # if(paired & return_dist_CI)
  # {
  #   CI_data <- list(estimates = unlist(lapply(X = dims,FUN = function(X){
  #     
  #     return(diagram_distance(diag1, diag2, dim = X, p = Inf))
  #     
  #   })))
  #   sampling_distrs <- list()
  #   for(d in dims)
  #   {
  #     sampling_distrs[[length(sampling_distrs) + 1]] <- lapply(X = 1:num_samples,FUN = function(X){
  #       
  #       return(diagram_distance(bootstrapped_diagrams1[[X]],bootstrapped_diagrams2[[X]], dim = d, p = Inf))
  #       
  #     })
  #   }
  #   CI_data$sampling_distrs <- sampling_distrs
  #   res$CI_data <- CI_data
  # }
  
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

#### NEW STUFF ####
#' Model inference with bootstrap.
#'
#' Carries out inference to determine if two datasets were likely generated from the same process or not (i.e. if
#' the persistence diagram of one dataset is a good model for the persistence diagram of the other dataset). Inference is
#' carried out by generating bootstrap resampled persistence diagrams from the two datasets and carrying out a permutation test
#' on the resulting two groups.
#' A small p-value in a certain dimension suggests that the datasets are not good models of each other.
#'
#' @param D1 the first dataset (a data frame).
#' @param D2 the second dataset (a data frame).
#' @param num_samples the number of bootstrap iterations, default 30.
#' @param dims a non-negative integer vector of the homological dimensions in which the test is to be carried out, default c(0,1).
#' @param paired a boolean flag for if there is a second-order pairing between diagrams at the same index in different groups, default FALSE
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
#' 
#'   # do model inference test
#'   model_test <- bootstrap_model_inference(D1, D2, num_samples = 20,thresh = 1.75,num_samples = 20)
#'   
#'   # p-values show difference in the clustering of points but not in the
#'   # arrangement of loops:
#'   model_test$p_values
#'   
#' }

# bootstrap_model_inference <- function(D1, D2, num_samples, dims = c(0,1), paired = F, num_workers = parallelly::availableCores(omit = 1), verbose = F,FUN_diag = "calculate_homology",FUN_boot = "calculate_homology",thresh,distance_mat = FALSE,ripser = NULL,return_diagrams = FALSE){
#   
#   # do error checks for D1, D2 and paired
#   
#   # generate bootstrap samples
#   # if paired then generate shared samples, otherwise different samples
#   s1 <- lapply(X = 1:num_samples,FUN = function(X){
#     
#     return(unique(sample(1:nrow(D1),size = nrow(D1),replace = T)))
#     
#   })
#   if(paired)
#   {
#     s2 <- s1
#   }else
#   {
#     s2 <- lapply(X = 1:num_samples,FUN = function(X){
#       
#       return(unique(sample(1:nrow(D2),size = nrow(D2),replace = T)))
#       
#     })
#   }
#   
#   # store maximum dimension
#   m <- max(dims)
#   
#   # compute actual distance between the two full persistence diagrams
#   if(verbose)
#   {
#     message(paste0("Starting to calculate full persistence diagrams at ",Sys.time()))
#   }
#   if(FUN_diag == "PyH")
#   {
#     diag1 <- PyH(X = D1,maxdim = m,thresh = thresh,distance_mat = distance_mat,ripser = ripser,ignore_infinite_cluster = F)
#     diag2 <- PyH(X = D2,maxdim = m,thresh = thresh,distance_mat = distance_mat,ripser = ripser,ignore_infinite_cluster = F)
#   }
#   if(FUN_diag == "calculate_homology")
#   {
#     diag1 <- diagram_to_df(TDAstats::calculate_homology(mat = D1,dim = m,threshold = thresh,format = ifelse(test = distance_mat == T,yes = "distmat",no = "cloud")))
#     diag2 <- diagram_to_df(TDAstats::calculate_homology(mat = D2,dim = m,threshold = thresh,format = ifelse(test = distance_mat == T,yes = "distmat",no = "cloud")))
#   }
#   if(FUN_diag == "ripsDiag")
#   {
#     diag1 <- ripsDiag(X = D1,maxdimension = m,maxscale = thresh,dist = ifelse(test = distance_mat == F,yes = "euclidean",no = "arbitrary"),library = "dionysus",printProgress = F)
#     diag1 <- diagram_to_df(diag1)
#     diag2 <- ripsDiag(X = D2,maxdimension = m,maxscale = thresh,dist = ifelse(test = distance_mat == F,yes = "euclidean",no = "arbitrary"),library = "dionysus",printProgress = F)
#     diag2 <- diagram_to_df(diag2)
#   }
#   
#   # use bootstrap samples to get distances and perform inference
#   if(verbose)
#   {
#     message(paste0("Starting to calculate bootstrap distances at ",Sys.time()))
#   }
#   distances1 <- bootstrap_persistence_thresholds(X = D1,maxdim = m,distance_mat = distance_mat,thresh = thresh,num_workers = num_workers,ripser = ripser,FUN_boot = FUN_boot,bootstrap_samples = s1,return_distances = TRUE,num_samples = num_samples,FUN_diag = FUN_diag)
#   distances2 <- bootstrap_persistence_thresholds(X = D2,maxdim = m,distance_mat = distance_mat,thresh = thresh,num_workers = num_workers,ripser = ripser,FUN_boot = FUN_boot,bootstrap_samples = s2,return_distances = TRUE,num_samples = num_samples,FUN_diag = FUN_diag)
#   res <- lapply(X = 0:m,FUN = function(X){
#     
#     bottleneck_dist <- diagram_distance(D1 = diag1,D2 = diag2,dim = X,p = Inf)
#     boot_vals <- c(as.numeric(distances1[seq(X + 1,length(distances1),m + 1)]),as.numeric(distances2[seq(X + 1,length(distances2),m + 1)]))
#     Z <- length(which(boot_vals >= bottleneck_dist))
#     p <- (Z + 1)/(length(boot_vals) + 1)
#     return(list(boot_vals,bottleneck_dist,p))
#     
#   })
#   res <- list(dimensions = dims,permvals = lapply(res,"[[",1),test_statistics = lapply(res,"[[",2),p_values = lapply(res,"[[",3))
#   names(res$permvals) <- as.character(0:m)
#   res$test_statistics <- as.numeric(res$test_statistics)
#   names(res$test_statistics) <- as.character(0:m)
#   res$p_values <- as.numeric(res$p_values)
#   names(res$p_values) <- as.character(0:m)
#   
#   if(return_diagrams)
#   {
#     res$diag1 <- diag1
#     res$diag2 <- diag2
#   }
#   
#   if(verbose)
#   {
#     message(paste0("Finished inference test at ",Sys.time()))
#   }
#   
#   return(res)
#   
# }

#### NEW STUFF ####
#' Model inference with independence test.
#'
#' Carries out inference to determine if two datasets were likely generated from the same process or not (i.e. if
#' the persistence diagram of one dataset is a good model for the persistence diagram of the other dataset). Inference is
#' carried out by generating bootstrap resampled persistence diagrams from the two datasets and carrying out a permutation test
#' on the resulting two groups.
#' A small p-value in a certain dimension suggests that the datasets are not good models of each other.
#'
#' @param D1 the first dataset (a data frame).
#' @param D2 the second dataset (a data frame).
#' @param iterations the number of iterations for permuting group labels, default 20.
#' @param num_samples the number of bootstrap iterations, default 30.
#' @param dims a non-negative integer vector of the homological dimensions in which the test is to be carried out, default c(0,1).
#' @param paired a boolean flag for if there is a second-order pairing between diagrams at the same index in different groups, default FALSE
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
#' 
#'   # do model inference test
#'   model_test <- permutation_model_inference(D1, D2, iterations = 20,thresh = 1.75,num_samples = 20)
#'   
#'   # p-values show difference in the clustering of points but not in the
#'   # arrangement of loops:
#'   model_test$p_values
#'   
#' }
# 
# independence_model_inference <- function(D1, D2, num_samples, dims = c(0,1), sigma = 1, t = 1, rho = NULL, paired = F, num_workers = parallelly::availableCores(omit = 1), verbose = F,FUN_boot = "calculate_homology",thresh,distance_mat = FALSE,ripser = NULL,return_diagrams = FALSE){
#   
#   # do error checks for D1, D2 and paired
#   
#   # generate bootstrap samples
#   # if paired then generate shared samples, otherwise different samples
#   s1 <- lapply(X = 1:num_samples,FUN = function(X){
#     
#     return(unique(sample(1:nrow(D1),size = nrow(D1),replace = T)))
#     
#   })
#   if(paired)
#   {
#     s2 <- s1
#   }else
#   {
#     s2 <- lapply(X = 1:num_samples,FUN = function(X){
#       
#       return(unique(sample(1:nrow(D2),size = nrow(D2),replace = T)))
#       
#     })
#   }
#   
#   # use bootstrap samples to get bootstrapped diagrams
#   if(verbose)
#   {
#     message(paste0("Starting to calculate bootstrapped diagrams at ",Sys.time()))
#   }
#   bootstrapped_diagrams1 <- bootstrap_persistence_thresholds(X = D1,maxdim = max(dims),distance_mat = distance_mat,thresh = thresh,num_workers = num_workers,ripser = ripser,FUN_boot = FUN_boot,bootstrap_samples = s1,num_samples = num_samples)
#   bootstrapped_diagrams1 <- lapply(X = 1:num_samples,FUN = function(X){
#     
#     return(data.frame(dimension = bootstrapped_diagrams1[[(X-1)*3 + 1]],birth = bootstrapped_diagrams1[[(X-1)*3 + 2]],death = bootstrapped_diagrams1[[X*3]]))
#     
#   })
#   bootstrapped_diagrams2 <- bootstrap_persistence_thresholds(X = D2,maxdim = max(dims),distance_mat = distance_mat,thresh = thresh,num_workers = num_workers,ripser = ripser,FUN_boot = FUN_boot,bootstrap_samples = s2,num_samples = num_samples)
#   bootstrapped_diagrams2 <- lapply(X = 1:num_samples,FUN = function(X){
#     
#     return(data.frame(dimension = bootstrapped_diagrams2[[(X-1)*3 + 1]],birth = bootstrapped_diagrams2[[(X-1)*3 + 2]],death = bootstrapped_diagrams2[[X*3]]))
#     
#   })
#   
#   # carry out independence test
#   if(verbose)
#   {
#     message(paste0("Starting independence test at ",Sys.time()))
#   }
#   p_values <- c()
#   for(d in dims)
#   {
#     vals <- unlist(lapply(1:num_samples,FUN = function(X){
#       
#       return(cos(diagram_distance(bootstrapped_diagrams1[[X]],bootstrapped_diagrams2[[X]],dim = 0,distance = "fisher",sigma = sigma)))
#       
#     }))
#     test <- wilcox.test(vals,mu = 1,alternative = "less",conf.int = T)
#     p_values <- c(p_values,test$p.value)
#   }
#   res <- list(p_values = p_values)
#   
#   # res <- independence_test(bootstrapped_diagrams1,bootstrapped_diagrams2,dims = dims,sigma = sigma,t = t,rho = rho,num_workers = num_workers)
#   if(verbose)
#   {
#     message(paste0("Finished independence test at ",Sys.time()))
#   }
#   if(return_diagrams)
#   {
#     res$diagrams1 <- bootstrapped_diagrams1
#     res$diagrams2 <- bootstrapped_diagrams2
#   }
#   return(res)
#   
# }


