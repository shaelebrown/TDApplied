# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'


# Script for Statistical Inference using Persistent Hom
# not based on functions in the TDA package, because
# it uses the Dionysus wasserstein distance, which is by default
# an infinity-norm distance

# CHANGES TO MAKE:
# include option for infinity norm in Wasserstein metrics
# make sure test statistic is correct

# IMPORT LIBRARIES ----
library(TDA)
library(TDAStats)
library(parallel)
library(doParallel)
library(foreach)
library(clue)
library(rdist)

# FUNCTIONS ----
TDA_diagram_to_df <- function(d){
  # function to convert d to a dataframe
  # d is a diagram from library TDA

  d = d[[1]]
  class(d) = "matrix"
  return(as.data.frame(d))
}

TDAStats_diagram_to_df <- function(d){
  # function to convert d to a dataframe
  # d is a diagram from library TDAStats

  return(as.data.frame(d))
}

check_diagram <- function(d){
  # error checks for a diagram
  # d is a diagram, stored as a dataframe

  if(nrow(d) == 0)
  {
    stop("Every diagram must be non-empty.")
  }

  if(ncol(d) != 3)
  {
    stop("Every diagram must have three columns.")
  }

  if(class(d[,1]) != "numeric" | class(d[,2]) != "numeric" | class(d[,3]) != "numeric")
  {
    stop("Diagrams must have numeric columns.")
  }

  if(!all.equal(d[,1],as.integer(d[,1])))
  {
    stop("Homology dimensions must be whole numbers.")
  }

  if(length(which(d[,1]) < 0) > 0)
  {
    stop("Homology dimensions must be >= 0.")
  }

  if(length(which(d[,2]) < 0) > 0 | length(which(d[,3]) < 0) > 0)
  {
    stop("Birth and death radii must be >= 0.")
  }

  if(length(which(complete.cases(d))) != nrow(d))
  {
    stop("Diagrams can't have missing values.")
  }

}

all_diagrams <- function(groups,lib){
  # function to make sure all groups are lists or vectors of diagrams, to convert the diagrams to dataframes
  # and to error check each diagram
  # groups is a vector or list of vectors or lists of diagrams
  # lib is either "TDA" or "TDAStats

  # compute cumulative sums of groups lengths in order to correctly compute diagram indices
  csum_group_sizes = cumsum(unlist(lapply(diagram_groups,FUN = length)))
  csum_group_sizes = c(0,csum_group_sizes)

  for(g in 1:length(groups))
  {
    for(diag in 1:length(groups[g]))
    {
      if(class(groups[g][diag]) != ifelse(test = lib == "TDA",yes = "diagram",no = c("matrix","array")))
      {
        stop("Every diagram must be the output from a homology calculation from either TDA or TDAStats.")
      }else
      {
        # convert each diagram to a data frame according to which library they were computed with
        # and size index of each diagram in all diagrams
        if(lib == "TDA")
        {
          groups[g][diag] = list(diag = TDA_diagram_to_df(groups[g][diag]),ind = csum_group_sizes[g - 1] + d)
        }else
        {
          groups[g][diag] = list(diag = TDAStats_diagram_to_df(groups[g][diag]),ind = csum_group_sizes[g - 1] + d)
        }
      }

      check_diagram(groups[g][diag]$diag)
    }
  }

  return(groups)

}

check_params <- function(iterations,p,q,dims,paired,lib){
  # error checks on the parameters for the main function

  if(!is.numeric(iterations))
  {
    stop("Iterations must be numeric.")
  }

  if(iterations != as.integer(iterations))
  {
    stop("Iterations must be a whole number.")
  }

  if(iterations <= 1 | is.infinite(iterations))
  {
    stop("Iterations must be a finite number at least 1.")
  }

  if(!as.numeric(p))
  {
    stop("p must be numeric.")
  }

  if(p < 1)
  {
    stop("p must be at least 1.")
  }

  if(!is.numeric(q))
  {
    stop("q must be numeric.")
  }

  if(q < 1 | is.infinite(q))
  {
    stop("q must be a finite number at least 1.")
  }

  if(!is.numeric(dims))
  {
    stop("dims must be a numeric vector or value.")
  }

  if(!all.equal(dims,as.integer(dims)))
  {
    stop("dims must be whole numbers.")
  }

  if(length(which(dims < 0 | is.infinite(dims))) > 0)
  {
    stop("dims must be positive and finite.")
  }

  if(!is.logical(paired))
  {
    stop("paired must be T or F.")
  }

  if(!is.character(lib) | is.vector(lib) | lib %in% c("TDA","TDAStats") == F)
  {
    stop("lib must be a single character, either TDA or TDAStats.")
  }

}

d_wasserstein <- function(D1,D2,dim,p){
  # function to compute the wasserstein/bottleneck metric between two diagrams
  # D1 and D2 are diagrams stored as data frames
  # dim is the dimension to subset
  # p is the power of the wasserstein distance, p >= 1

  # subset both diagrams by dimension X
  D1_subset = D1[which(D1$dimension == dim),]
  D2_subset = D2[which(D2$dimension == dim),]
  D1_subset = D1_subset[,2:3]
  D2_subset = D2_subset[,2:3]

  # create empty diagonals for the persistence landscapes
  diag1 = D1_subset[0,]
  diag2 = D2_subset[0,]

  # if both subsets are empty then set their distance to 0
  if(nrow(D1_subset) == 0 & nrow(D2_subset) == 0)
  {
    return(0)
  }

  # remove diagonal entries from D1_subset and D2_subset
  D1_subset = D1_subset[which(D1_subset[,2] != D1_subset[,3]),]
  D2_subset = D2_subset[which(D2_subset[,2] != D2_subset[,3]),]

  if(nrow(D1_subset) > 0)
  {
    for(i in 1:nrow(D1_subset))
    {
      # for each non-trivial element in D1_subset we add its projection onto the diagonal in diag1
      proj_diag = mean(as.numeric(D1_subset[i,]))
      diag1 = rbind(diag1,data.frame(birth = proj_diag,death = proj_diag))
    }
  }

  if(nrow(D2_subset) > 0)
  {
    for(i in 1:nrow(D2_subset))
    {
      # for each non-trivial element in D2_subset we add its projection onto the diagonal in diag2
      proj_diag = mean(as.numeric(D2_subset[i,]))
      diag2 = rbind(diag2,data.frame(birth = proj_diag,death = proj_diag))
    }
  }

  # since an element b of D1_subset is either matched to an element of D2 or to the projection of b onto the diagonal
  # we form the two sets to be matched by row binding D1_subset with diag2 and D2_subset with diag1
  D1_subset = rbind(D1_subset,diag2)
  D2_subset = rbind(D2_subset,diag1)

  # compute the distance matrix between rows of D1_subset and D2_subset
  if(is.finite(p))
  {
    # compute the p-wasserstein distance
    dist_mat = as.matrix(cdist(D1_subset,D2_subset,metric = "minkowski",p = p))
  }else
  {
    # compute the bottleneck distance
    dist_mat = as.matrix(cdist(D1_subset,D2_subset,metric = "maximum"))
  }

  # use the Hungarian algorithm from the clue package to find the minimal weight matching
  best_match = solve_LSAP(x = dist_mat,maximum = F)

  # return the distance between D1 and D2
  if(is.finite(p))
  {
    return(sum(dist_mat[cbind(seq_along(best_match), best_match)]^(p))^(1/p))
  }
  return(max(dist_mat[cbind(seq_along(best_match), best_match)]))

}

loss <- function(diagram_groups,dist_mats,dims,p,q){
  # function to compute the F_{p,q} loss between groups of diagrams
  # diagram_groups are the (possibly permuted) groups of diagrams
  # dist mats stores the current distance matrices for the diagrams in each dimension
  # dims are the homological dimensions of diagrams to consider
  # p is a number >=1
  # q is a finite number >= 1

  # initialize a cluster for computing distances between diagrams in parallel
  cl = makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  clusterEvalQ(cl,c(library(clue),library(rdist)))
  clusterExport(cl,c("d_wasserstein"))
  clusterExport(cl,c("diagram_groups","dist_mats"),envir = environment())

  # initialize return vector of statistics, one for each dimension
  statistics = c()

  # compute loss function and update distance matrices in each dimension
  for(dim in dims)
  {

    # to store sum of test statistics within each group
    d_tots = unlist(lapply(X = diagram_groups,FUN = function(X){

      # compute distances in parallel for each unique pair of diagrams in each group
      d = foreach(i = combn(x = length(X),m = 2,simplify = F),.combine = c) %dopar%
        {
          # get two diagrams within the group to compute their distance
          i = unlist(i)

          # if the distance matrix has not already stored the distance calculation
          if(dist_mats[which.min(dims == dim)][X[i[[1]]]$ind,X[i[[2]]$ind]] == -1)
          {
            return(d_wasserstein(D1 = X[i[[1]]]$diag,D2 = X[i[[2]]]$diag,dim = dim,p = p)^q)
          }

          # else return the already stored distance value
          return(dist_mats[which.min(dims == dim)][X[i[[1]]]$ind,X[i[[2]]$ind]])
        }

      return(d)

    }))

    # update the upper triangle of dist_mat
    for(X in diagram_groups)
    {
      for(i in 1:length(combn(x = length(X),m = 2,simplify = F)))
      {
        j = comb1[[i]]
        dist_mats[which.min(dims == dim)][X[combn(x = length(X),m = 2,simplify = F)[i][[1]]]$ind,X[combn(x = length(X),m = 2,simplify = F)[i][[2]]$ind]] = d_tots[which.min(diagram_groups == X)][i]
      }
    }

    # append calculated loss statistic to statistics vector
    statistics = c(statistics,sum(unlist(lapply(X = 1:length(d_tots),FUN = function(X){

      return(sum(d_tots[X])/(length(diagram_groups[X])*(length(diagram_groups[X]) - 1)))

    }))))

  }

  stopCluster(cl)

  return(list(statistics = statistics,dist_mats = dist_mats))

}

permutation_test <- function(...,iterations = 100,p = 2,q = 2,dims = c(0,1),paired = F,lib = "TDA"){

  # function to test whether or not multiple groups of persistence diagrams come from the same geometric process
  # ... are the groups of diagrams, either stored as lists or vectors
  # iterations is the number of permutations we will calculate for group labels
  # p is the wasserstein distance parameter, p >= 1
  # q is the finite distance exponential, q >= 1
  # dims is a vector of desired homological dimensions
  # paired is a boolean which determines if dependencies exist between diagrams of the same indices in different groups

  # retrieve diagram groups
  diagram_groups <- list(...)

  # make sure there are at least two groups
  if(length(diagram_groups) < 2)
  {
    stop("At least two groups of persistence diagrams must be supplied.")
  }

  # check each diagram, converting each to a data frame and storing index in all diagrams
  diagram_groups = all_diagrams(diagram_groups,lib)

  # error check all parameters
  check_params(iterations,p,q,dims,paired,lib)

  # make sure that if paired == T then all groups have the same number of diagrams
  if(paired)
  {
    if(length(unique(unlist(lapply(diagrams_groups,FUN = length)))) != 1)
    {
      stop("If paired is true then all groups of diagrams must have the same number of elements.")
    }
  }

  # time computations
  s = Sys.time()

  # make distance matrix for all diagrams for each dimension
  n = sum(unlist(lapply(groups,FUN = length)))
  dist_mats = lapply(X = dims,FUN = function(X){return(matrix(data = -1,nrow = n,ncol = n))})

  # compute loss function on observed data
  test_loss = loss(diagram_groups = diagram_groups,dist_mat = dist_mat,dims = dims,p = p,q = q)
  dist_mats = test_loss$dist_mats
  test_statistics = test_loss$statistics

  # get permutation values
  perm_values = lapply(X = dims,FUN = function(X){return(list())})
  for(perm in 1:iterations)
  {
    if(paired == F)
    {
      # sample groups from their union, maintaining group sizes
      perm = split(x = 1:n,f = sample(unlist(lapply(X = 1:length(diagram_groups),FUN = function(X){return(rep(X,length = X))})),size = n,replace = F))

      permuted_groups = lapply(X = perm,FUN = function(X){

        res = list()
        for(i in 1:length(X))
        {
          # invert diagram indices
          g = min.which(csum_group_sizes > X[i])

          # append to permuted group
          res[[length(res) + 1]] <- diagram_groups[g][i - csum_group_sizes[g]]
        }
        return(res)})
    }else
    {
      # pairings are maintained
      perm = lapply(X = 1:length(diagram_groups[[1]]),FUN = function(X){return(sample(1:length(diagram_groups),size = length(diagram_groups),replace = F))})

      permuted_groups = lapply(X = length(diagram_groups),FUN = function(X){

        # get the Xth element of each list
        groups = do.call("[[",X,perm)
        res = list()
        for(i in 1:length(groups))
        {
          # append correct diagram to the end of the list
          res[[length(res) + 1]] <- diagram_groups[groups][i]
        }
        return(res)})
    }

    # compute loss function, add to permutation values and updated distance matrices
    permuted_loss = loss(permuted_groups,dist_mats = dist_mats,dims = dims,p = p,q = q)
    dist_mats = permuted_loss$dist_mats
    for(d in dims)
    {
      perm_values[d][(length(perm_values[d]) + 1):(length(perm_values[d]) + iterations)] = permuted_loss$statistics[d]
    }

  }

  # set up return list
  names(perm_values) = as.character(dims)
  names(test_statistics) = as.character(dims)
  pval = lapply(X = 1:length(dims),FUN = function(X){

    return((sum(perm_values[X] <= test_statistics[X]) + 1)/(iterations + 1))

  })
  names(pval) = as.character(dims)

  results = list(dimensions = dims,
                 permvals = perm_values,
                 test_statistic = test_statistics,
                 p_value = pval)

  # print time duration
  print(paste0("Time taken: ",Sys.time()-s))

  return(results)

}

