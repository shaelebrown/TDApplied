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
# include TDAStats and TDA functionalities
# remove dependencies on unnecessary libraries like data.table
# include option for infinity norm in Wasserstein metrics
# make sure test statistic is correct
# keep dependence functionality and cite paper
# allow for vector of desired dimensions
# allow for multiple groups of diagrams and cite papers

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
  d = d[[1]]
  class(d) = "matrix"
  return(as.data.frame(d))
}

TDAStats_diagram_to_df <- function(d){
  return(as.data.frame(d))
}

d_wasserstein <- function(D1,D2,dim,p){
  # function to compute the wasserstein metric between two barcodes
  # B1 and B2 are diagrams as outputted from the function ripsDiag
  # dim is the maximum dimension to consider
  # p is the finite power of the wasserstein distance, p >= 1

  # calculate the wasserstein metric in each dimension
  # subset both barcodes by dimension X

  B1_subset = B1[which(B1$dimension == dim),]
  B2_subset = B2[which(B2$dimension == dim),]
  B1_subset = B1_subset[,2:3]
  B2_subset = B2_subset[,2:3]

  # create empty diagonals for the persistence landscapes
  diag1 = B1_subset[0,]
  diag2 = B2_subset[0,]

  # if both subsets are empty then set their distance to 0
  if(nrow(B1_subset) == 0 & nrow(B2_subset) == 0)
  {
    return(0)
  }

  # remove diagonal entries from B1_subset and B2_subset
  B1_subset = B1_subset[which(B1_subset[,2] != B1_subset[,3]),]
  B2_subset = B2_subset[which(B2_subset[,2] != B2_subset[,3]),]

  if(nrow(B1_subset) > 0)
  {
    for(i in 1:nrow(B1_subset))
    {
      # for each non-trivial element in B1_subset we add its projection onto the diagonal in diag1
      proj_diag = mean(as.numeric(B1_subset[i,]))
      diag1 = rbind(diag1,data.frame(birth = proj_diag,death = proj_diag))
    }
  }

  if(nrow(B2_subset) > 0)
  {
    for(i in 1:nrow(B2_subset))
    {
      # for each non-trivial element in B2_subset we add its projection onto the diagonal in diag2
      proj_diag = mean(as.numeric(B2_subset[i,]))
      diag2 = rbind(diag2,data.frame(birth = proj_diag,death = proj_diag))
    }
  }

  # since an element b of B1_subset is either matched to an element of B2 or to the projection of b onto the diagonal
  # we form the two sets to be matched by row binding B1_subset with diag2 and B2_subset with diag1
  B1_subset = rbind(B1_subset,diag2)
  B2_subset = rbind(B2_subset,diag1)

  # use cdist function from the rdist package to find the wasserstein distance matrix between rows of the updated B1_subset and B2_subset
  dist_mat = as.matrix(cdist(B1_subset,B2_subset,metric = "minkowski",p = p))

  # use the Hungarian algorithm from the clue package to find the minimal weight matching
  best_match = solve_LSAP(x = dist_mat,maximum = F)

  # return the distance for dimension X
  return(sum(dist_mat[cbind(seq_along(best_match), best_match)]^(p))^(1/p))

}

loss <- function(diagrams_1,diagrams_2,dist_mat,perm,dim,p){
  # function to compute the F_{p,p} loss between two sets of diagrams
  # environment e contains diagrams_1, diagrams_2 and the current distmat
  # perm is the list of indices between 1 and n1+n2 which consistute the first set of diagrams
  # dim is the maximum dimension of diagrams to consider
  # p is a finite number >=1

  # get all possible pairs of distinct members of each list of diagrams
  comb1 = combn(perm,m = 2,simplify = F)
  comb2 = combn(setdiff(1:(length(diagrams_1) + length(diagrams_2)),perm),m = 2,simplify = F)

  # compute the pairwise barcode distances for both lists in parallel
  cl = makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  clusterEvalQ(cl,c(library(clue),library(rdist)))
  clusterExport(cl,c("d_wasserstein"))
  clusterExport(cl,c("diagrams_1","diagrams_2","dist_mat"),envir = environment())

  d1_tot = foreach(i = comb1,.combine = c) %dopar%
    {
      i = unlist(i)
      if(dist_mat[i[[1]],i[[2]]] == -1)
      {
        if(i[[1]] > length(diagrams_1))
        {
          D1 = diagrams_2[[i[[1]] - length(diagrams_1)]]
        }else
        {
          D1 = diagrams_1[[i[[1]]]]
        }

        if(i[[2]] > length(diagrams_1))
        {
          D2 = diagrams_2[[i[[2]] - length(diagrams_1)]]
        }else
        {
          D2 = diagrams_1[[i[[2]]]]
        }
        return(d_wasserstein(D1 = D1,D2 = D2,dim = dim,p = p))
      }
      return(dist_mat[i[[1]],i[[2]]])
    }

  d2_tot = foreach(i = comb2,.combine = c) %dopar%
    {
      if(dist_mat[i[[1]],i[[2]]] == -1)
      {
        if(i[[1]] > length(diagrams_1))
        {
          D1 = diagrams_2[[i[[1]] - length(diagrams_1)]]
        }else
        {
          D1 = diagrams_1[[i[[1]]]]
        }

        if(i[[2]] > length(diagrams_1))
        {
          D2 = diagrams_2[[i[[2]] - length(diagrams_1)]]
        }else
        {
          D2 = diagrams_1[[i[[2]]]]
        }
        return(d_wasserstein(D1 = D1,D2 = D2,dim = dim,p = p))
      }
      return(dist_mat[i[[1]],i[[2]]])
    }

  stopCluster(cl)

  # update the upper triang of dist_mat
  for(i in 1:length(comb1))
  {
    j = comb1[[i]]
    dist_mat[j[[1]],j[[2]]] = d1_tot[[i]]
  }
  for(i in 1:length(comb2))
  {
    j = comb2[[i]]
    dist_mat[j[[1]],j[[2]]] = d2_tot[[i]]
  }

  # return sum of within group distances
  return(list(res = (1/(length(diagrams_1)*(length(diagrams_1) - 1))) * sum(d1_tot) + (1/(length(diagrams_2)*(length(diagrams_2) - 1))) * sum(d2_tot),dist_mat = dist_mat))

}

permutation_test <- function(diagrams_1,diagrams_2,iterations,p,dim,subject_dependence){

  # function to test whether or not two sets of labelled persistence diagrams come from the same geometric process
  # diagrams_1 is a list of diagrams in the first group and likewise for diagrams_2 and the second group
  # iterations is the number of permutations we will calculate for class labels
  # dim is the maximum dimension in which to calculate homology for barcodes
  # p is the finite wasserstein distance parameter, p >= 1

  # time computations
  s = Sys.time()

  # make distance matrix for both sets of diagrams
  dist_mat = matrix(data = -1,nrow = length(diagrams_1) + length(diagrams_2),ncol = length(diagrams_1) + length(diagrams_2))

  # compute loss function on observed data
  test_loss = loss(diagrams_1 = diagrams_1,diagrams_2 = diagrams_2,dist_mat = dist_mat,perm = 1:length(diagrams_1),dim = dim,p = p)
  dist_mat = test_loss$dist_mat
  test_loss = test_loss$res

  # get permutation values
  perm_values = c()
  for(X in 1:iterations)
  {
    if(subject_dependence == F | length(diagrams_1) != length(diagrams_2))
    {
      # sample two groups from the union of the two sets of diagrams
      perm = sample(1:(length(diagrams_1) + length(diagrams_2)),size = n1,replace = F)
    }else
    {
      # group dependencies much be maintained
      perm = rbinom(n = length(diagrams_1),size = 1,prob = 0.5)
      perm = c(which(perm == 1),length(diagrams_1) + which(perm == 0))
    }

    # return loss function
    res = loss(diagrams_1 = diagrams_1,diagrams_2 = diagrams_2,dist_mat = dist_mat,perm = perm,dim = dim,p = p)
    dist_mat = res$dist_mat
    perm_values = c(perm_values,res$res)

  }

  # return results
  answer = list(dimension = dim,
                permvals = perm_values,
                test_statistic = test_loss,
                p_value = (sum(perm_values <= test_loss) + 1)/(length(perm_values) + 1))

  # print time duration
  print(paste0("Time taken: ",Sys.time()-s))

  return(answer)

}

