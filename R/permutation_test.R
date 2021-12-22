
#' Permutation test for persistence diagrams
#'
#' Calculates the turner loss function for a number of groups of persistence diagrams
#' and then permutes group labels several times and recomputes the loss function
#' to generate a null distribution and calculate a p value for each desired dimension.
#' This function keeps track of distance calculations as to not repeat calculations
#' which can individually take a long time.
#'
#' The `mat` parameter should be a numeric matrix with each row corresponding
#' to a single point, and each column corresponding to a single dimension. Thus,
#' if `mat` has 50 rows and 5 columns, it represents a point cloud with 50 points
#' in 5 dimensions. The `dim` parameter should be a positive integer.
#' Alternatively, the `mat` parameter could be a distance matrix (upper
#' triangular half is ignored); note: `format` should be specified as "ldm".
#'
#' @param ... groups of persistence diagrams, outputted from a homology calculation in TDA or TDAStats
#' @param iterations number of iterations for permuting group labels
#' @param p wasserstein parameter, number at least 1 (and bottleneck distance if == Inf)
#' @param q  finite number at least 1 for exponentiation in Turner loss function
#' @param dims homological dimensions in which the test is to be carried out
#' @param paired if there is a second-order pairing between diagrams at the same index in different groups
#' @param lib either "TDA" or "TDAStats" for consistency
#' @param distance either "wasserstein" or "Turner" for determining which distance of diagrams to use
#' @param verbose if the time duration of the function call should be printed
#'
#' @return list with dimensions used (named vector), permutation loss values in each dimension (named list), test statistics in each dimension (named vector)
#'                   a p-value for each dimension (all stored in a named vector) and the time duration of the function call.
#' @importFrom stats complete.cases
#' @export
#' @examples
#'
#' # create three groups of persistence diagrams on 2D Gaussians using TDA
#' g1 <- lapply(X = 1:3,FUN = function(X){
#'
#' diag <- ripsDiag(data.frame(x = rnorm(100,mean = 0,sd = 1),y = rnorm(100,mean = 0,sd = 1)),maxscale = 1,maxdimension = 1)
#' df <- TDA_diagram_to_df(d = diag)
#' return(list(diag = df,ind = X))
#'
#' })
#'
#' g2 <- lapply(X = 1:3,FUN = function(X){
#'
#' diag <- ripsDiag(data.frame(x = rnorm(100,mean = 0,sd = 1),y = rnorm(100,mean = 0,sd = 1)),maxscale = 1,maxdimension = 1)
#' df <- TDA_diagram_to_df(d = diag)
#' return(list(diag = df,ind = X + 3))
#'
#' })
#'
#' g3 <- lapply(X = 1:3,FUN = function(X){
#'
#' diag <- ripsDiag(data.frame(x = rnorm(100,mean = 0,sd = 1),y = rnorm(100,mean = 0,sd = 1)),maxscale = 1,maxdimension = 1)
#' df <- TDA_diagram_to_df(d = diag)
#' return(list(diag = df,ind = X + 6))
#'
#' })
#'
#' # do permutation test with 20 iterations, p,q = 2, in dimensions 0 and 1, with no pairing using Turner distance and printing the time duration
#' perm_test = permutation_test(g1,g2,g3,iterations = 20,distance = "Turner",verbose = TRUE)

permutation_test <- function(...,iterations = 100,p = 2,q = 2,dims = c(0,1),paired = F,lib = "TDA",distance = "wasserstein",verbose = FALSE){

  # function to test whether or not multiple groups of persistence diagrams come from the same geometric process
  # ... are the groups of diagrams, either stored as lists or vectors
  # iterations is the number of permutations we will calculate for group labels
  # p is the wasserstein distance parameter, p >= 1
  # q is the finite distance exponential, q >= 1
  # dims is a vector of desired homological dimensions
  # paired is a boolean which determines if dependencies exist between diagrams of the same indices in different groups
  # lib is the TDA library used for homological calculations, either "TDA" (default) or "TDAStats"
  # distance is either "wasserstein" or "Turner" and determines how distances will be computed between diagrams
  # verbose is either TRUE or FALSE (default), printing runtime of function call

  # retrieve diagram groups
  diagram_groups <- list(...)

  # make sure there are at least two groups
  if(length(diagram_groups) < 2)
  {
    stop("At least two groups of persistence diagrams must be supplied.")
  }

  # check each diagram, converting each to a data frame and storing their indices in all the diagrams
  diagram_groups <- all_diagrams(diagram_groups,lib)

  # error check function parameters
  check_params(iterations,p,q,dims,paired,lib,distance)

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
  test_loss <- loss(diagram_groups = diagram_groups,dist_mats = dist_mats,dims = dims,p = p,q = q,distance = distance)
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
    permuted_loss <- loss(permuted_groups,dist_mats = dist_mats,dims = dims,p = p,q = q,distance = distance)
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
