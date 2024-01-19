
# test_that("bootstrap_persistence_thresholds can do homology calculation with all three packages",{
# 
#   skip_if_not_installed("TDA")
#   skip_if_not_installed("TDAstats")
# 
#   D <- TDA::circleUnif(n = 50,r = 1)
#   # ripser = import_ripser()
#   expect_length(bootstrap_persistence_thresholds(X = D,maxdim = 1,thresh = 2.1,num_workers = 2,num_samples = 3),2)
#   expect_length(bootstrap_persistence_thresholds(X = D,FUN_diag = "ripsDiag",FUN_boot = "ripsDiag",maxdim = 1,thresh = 2.1,num_workers = 2,num_samples = 3),2)
#   # expect_length(bootstrap_persistence_thresholds(X = D,FUN_diag = "PyH",FUN_boot = "PyH",maxdim = 1,thresh = 2,ripser = ripser,num_workers = 2,num_samples = 3),2)
#   expect_length(bootstrap_persistence_thresholds(X = D,FUN_boot = "ripsDiag",maxdim = 1,thresh = 2,calculate_representatives = T,num_workers = 2,num_samples = 3),2)
#   expect_length(bootstrap_persistence_thresholds(X = D,FUN_diag = "ripsDiag",maxdim = 1,thresh = 2,calculate_representatives = T,num_workers = 2,num_samples = 3),3)
#   # expect_length(bootstrap_persistence_thresholds(X = D,FUN_boot = "PyH",maxdim = 1,thresh = 2,ripser = ripser,calculate_representatives = T,num_workers = 2,num_samples = 3),2)
#   expect_length(bootstrap_persistence_thresholds(X = as.matrix(dist(D)),maxdim = 1,thresh = 2,distance_mat = T,num_workers = 2,num_samples = 3),2)
#   expect_length(bootstrap_persistence_thresholds(X = as.matrix(dist(D)),FUN_boot = "ripsDiag",maxdim = 1,thresh = 2,distance_mat = T,num_workers = 2,num_samples = 3),2)
#   # expect_length(bootstrap_persistence_thresholds(X = as.matrix(dist(D)),FUN_boot = "PyH",maxdim = 1,thresh = 2,ripser = ripser,distance_mat = T,num_workers = 2,num_samples = 3),2)
#   expect_length(bootstrap_persistence_thresholds(X = as.matrix(dist(D)),FUN_boot = "ripsDiag",maxdim = 1,thresh = 2,calculate_representatives = T,distance_mat = T,num_workers = 2,num_samples = 3),2)
#   expect_length(bootstrap_persistence_thresholds(X = as.matrix(dist(D)),FUN_diag = "ripsDiag",maxdim = 1,thresh = 2,calculate_representatives = T,distance_mat = T,num_workers = 2,num_samples = 3),3)
#   # expect_length(bootstrap_persistence_thresholds(X = as.matrix(dist(D)),FUN_boot = "PyH",FUN_diag = "PyH",maxdim = 1,thresh = 2,ripser = ripser,calculate_representatives = T,distance_mat = T,num_workers = 2,num_samples = 3),3)
# 
# })


# test_that("bootstrap_persistence_thresholds can detect incorrect parameters correctly",{
# 
#   skip_if_not_installed("TDA")
#   skip_if_not_installed("TDAstats")
# 
#   # X, FUN, maxdim, thresh, distance_mat, ripser, ignore_infinite_cluster, calculate_representatives, num_samples, alpha, return_subsetted, return_diag
#   expect_error(bootstrap_persistence_thresholds(data.frame(),maxdim = 1,thresh = 2),"X")
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1),maxdim = 1,thresh = 2),"X")
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = c(1,2),y = c("1","2")),maxdim = 1,thresh = 2),"X")
#   expect_error(bootstrap_persistence_thresholds(X = NULL,maxdim = 1,thresh = 2),"X")
#   expect_error(bootstrap_persistence_thresholds(X = data.frame(x = c(1,NA)),maxdim = 1,thresh = 2),"missing")
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = c(1),y = c(2)),maxdim = 1,thresh = 2),"two")
# 
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN_diag = "calculatehomology",maxdim = 1,thresh = 2),"calculate_homology")
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN_diag = NULL,maxdim = 1,thresh = 2),"NULL")
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN_diag = 2,maxdim = 1,thresh = 2),"string")
# 
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN_boot = "calculatehomology",maxdim = 1,thresh = 2),"calculate_homology")
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN_boot = NULL,maxdim = 1,thresh = 2),"NULL")
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN_boot = 2,maxdim = 1,thresh = 2),"string")
# 
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),maxdim = 1.1,thresh = 2),"whole")
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),maxdim = NA,thresh = 2),"NA")
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),maxdim = -1,thresh = 2),"negative")
# 
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),maxdim = 1,thresh = 0),"positive")
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),maxdim = 1,thresh = NaN),"NaN")
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),maxdim = 1,thresh = "2"),"numeric")
# 
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN_diag = "calculate_homology",maxdim = 1,thresh = 2,distance_mat = "F"),"logical")
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN_diag = "ripsDiag",maxdim = 1,thresh = 2,distance_mat = c(F,T)),"single")
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN_diag = "ripsDiag",maxdim = 1,thresh = 2,distance_mat = NULL),"NULL")
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN_diag = "ripsDiag",maxdim = 1,thresh = 2,distance_mat = NA),"NA")
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN_diag = "ripsDiag",maxdim = 1,thresh = 2,distance_mat = T),"square")
# 
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN_diag = "calculate_homology",maxdim = 1,thresh = 2,p_less_than_alpha = "F"),"logical")
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN_diag = "ripsDiag",maxdim = 1,thresh = 2,p_less_than_alpha = c(F,T)),"single")
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN_diag = "ripsDiag",maxdim = 1,thresh = 2,p_less_than_alpha = NULL),"NULL")
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN_diag = "ripsDiag",maxdim = 1,thresh = 2,p_less_than_alpha = NA),"NA")
# 
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN_diag = "calculate_homology",maxdim = 1,thresh = 2,return_pvals = "F"),"logical")
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN_diag = "ripsDiag",maxdim = 1,thresh = 2,return_pvals = c(F,T)),"single")
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN_diag = "ripsDiag",maxdim = 1,thresh = 2,return_pvals = NULL),"NULL")
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN_diag = "ripsDiag",maxdim = 1,thresh = 2,return_pvals = NA),"NA")
# 
#   # PyH parameters (ripser,ignore_infinite_cluster,calculate_representatives) are all working, tested in test-python.R
# 
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN_diag = "ripsDiag",maxdim = 1,thresh = 2,calculate_representatives = c(T,F)),"boolean")
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN_diag = "calculate_homology",maxdim = 1,thresh = 2,calculate_representatives = NULL),"NULL")
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN_diag = "calculate_homology",maxdim = 1,thresh = 2,calculate_representatives = NA),"NA")
# 
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN_diag = "ripsDiag",maxdim = 1,thresh = 2,num_samples = 0),"one")
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN_diag = "calculate_homology",maxdim = 1,thresh = 2,num_samples = Inf),"finite")
# 
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN_diag = "ripsDiag",maxdim = 1,thresh = 2,alpha = NA),"NA")
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN_diag = "calculate_homology",maxdim = 1,thresh = 2,alpha = 2),"1")
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN_diag = "calculate_homology",maxdim = 1,thresh = 2,alpha = c(0.5,0.4)),"single")
# 
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN_diag = "calculate_homology",maxdim = 1,thresh = 2,return_subsetted = c(T,F)),"single")
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN_diag = "calculate_homology",maxdim = 1,thresh = 2,return_subsetted = NULL),"NULL")
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN_diag = "ripsDiag",maxdim = 1,thresh = 2,return_diag = Inf),"logical")
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN_diag = "ripsDiag",maxdim = 1,thresh = 2,return_diag = NULL),"NULL")
#   expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN_diag = "ripsDiag",maxdim = 1,thresh = 2,return_diag = NA),"NA")
# 
# })

# test_that("bootstrap_persistence_thresholds is computing properly",{
# 
#   skip_if_not_installed("TDA")
#   skip_if_not_installed("TDAstats")
# 
#   D <- TDA::circleUnif(n = 50,r = 1)
# 
#   # ripsDiag
#   bs <- bootstrap_persistence_thresholds(X = D,FUN_diag = "ripsDiag",maxdim = 1,thresh = 2,calculate_representatives = T,return_diag = T,num_workers = 2,num_samples = 3,return_subsetted = T)
#   expect_length(bs$representatives,nrow(bs$diag))
#   expect_lte(length(bs$thresholds),2)
#   expect_gt(bs$thresholds[[1]],0)
#   expect_gt(bs$thresholds[[2]],0)
#   expect_length(bs$subsetted_representatives,nrow(bs$subsetted_diag))
#   if(length(which(bs$subsetted_diag$dimension == 0)) > 0)
#   {
#     expect_true(min(bs$subsetted_diag[which(bs$subsetted_diag$dimension == 0),]$death - bs$subsetted_diag[which(bs$subsetted_diag$dimension == 0),]$birth) > bs$thresholds[[1]])
#   }
#   if(length(which(bs$subsetted_diag$dimension == 1)) > 0)
#   {
#     expect_true(min(bs$subsetted_diag[which(bs$subsetted_diag$dimension == 1),]$death - bs$subsetted_diag[which(bs$subsetted_diag$dimension == 1),]$birth) > bs$thresholds[[2]])
#   }
#   bs <- bootstrap_persistence_thresholds(X = D,FUN_diag = "ripsDiag",maxdim = 1,thresh = 2,calculate_representatives = T,return_diag = T,num_workers = 2,num_samples = 3,return_subsetted = T)
#   expect_length(bs$representatives,nrow(bs$diag))
#   expect_lte(length(bs$thresholds),2)
#   expect_gt(bs$thresholds[[1]],0)
#   expect_gt(bs$thresholds[[2]],0)
#   expect_length(bs$subsetted_representatives,nrow(bs$subsetted_diag))
#   if(length(which(bs$subsetted_diag$dimension == 0)) > 0)
#   {
#     expect_true(min(bs$subsetted_diag[which(bs$subsetted_diag$dimension == 0),]$death - bs$subsetted_diag[which(bs$subsetted_diag$dimension == 0),]$birth) > bs$thresholds[[1]])
#   }
#   if(length(which(bs$subsetted_diag$dimension == 1)) > 0)
#   {
#     expect_true(min(bs$subsetted_diag[which(bs$subsetted_diag$dimension == 1),]$death - bs$subsetted_diag[which(bs$subsetted_diag$dimension == 1),]$birth) > bs$thresholds[[2]])
#   }
#   bs <- bootstrap_persistence_thresholds(X = D,FUN_diag = "ripsDiag",FUN_boot = "ripsDiag",maxdim = 1,thresh = 2,calculate_representatives = T,return_diag = T,num_workers = 2,num_samples = 3,return_subsetted = T)
#   expect_length(bs$representatives,nrow(bs$diag))
#   expect_lte(length(bs$thresholds),2)
#   expect_gt(bs$thresholds[[1]],0)
#   expect_gt(bs$thresholds[[2]],0)
#   expect_length(bs$subsetted_representatives,nrow(bs$subsetted_diag))
#   expect_true(min(bs$subsetted_diag[which(bs$subsetted_diag$dimension == 0),]$death - bs$subsetted_diag[which(bs$subsetted_diag$dimension == 0),]$birth) > bs$thresholds[[1]])
#   expect_true(min(bs$subsetted_diag[which(bs$subsetted_diag$dimension == 1),]$death - bs$subsetted_diag[which(bs$subsetted_diag$dimension == 1),]$birth) > bs$thresholds[[2]])
# 
#   # calculate_homology
#   bs <- bootstrap_persistence_thresholds(X = D,maxdim = 1,thresh = 2,return_diag = T,num_workers = 2,num_samples = 3,return_subsetted = T)
#   expect_length(bs$thresholds,2)
#   expect_gt(bs$thresholds[[1]],0)
#   expect_gt(bs$thresholds[[2]],0)
#   if(length(which(bs$subsetted_diag$dimension == 1)) > 0)
#   {
#     expect_true(min(bs$subsetted_diag[which(bs$subsetted_diag$dimension == 1),]$death - bs$subsetted_diag[which(bs$subsetted_diag$dimension == 1),]$birth) > bs$thresholds[[2]])
#   }
# 
#   # check on circle:
#   bs <- bootstrap_persistence_thresholds(X = D,FUN_diag = "ripsDiag",maxdim = 1,thresh = 2,return_diag = T,num_workers = 2,num_samples = 3)
#   expect_lte(length(bs$subsetted_diag$dimension),2)
#   bs <- bootstrap_persistence_thresholds(X = D,maxdim = 1,thresh = 2,return_diag = T,num_workers = 2,num_samples = 3)
#   expect_lte(length(bs$subsetted_diag$dimension),1)
# 
#   # one example by hand to verify thresholds
#   D <- data.frame(x = c(0,0,0),y = c(1,2,2.5))
#   diag <- diagram_to_df(TDA::ripsDiag(D,maxdimension = 0,maxscale = 2.5,library = "dionysus"))
#   permuted_diag <- function(D,s)
#   {
#     return(TDA::ripsDiag(X = D[s,],maxdimension = 0,maxscale = 2.5,library = "dionysus"))
#   }
#   d1 <- diagram_distance(diag,permuted_diag(D = D,s = c(1,1,1)),dim = 0,p = Inf)
#   d2 <- diagram_distance(diag,permuted_diag(D = D,s = c(2,2,2)),dim = 0,p = Inf)
#   d3 <- diagram_distance(diag,permuted_diag(D = D,s = c(3,3,3)),dim = 0,p = Inf)
#   d4 <- diagram_distance(diag,permuted_diag(D = D,s = c(1,1,2)),dim = 0,p = Inf)
#   d5 <- diagram_distance(diag,permuted_diag(D = D,s = c(1,1,3)),dim = 0,p = Inf)
#   d6 <- diagram_distance(diag,permuted_diag(D = D,s = c(2,2,1)),dim = 0,p = Inf)
#   d7 <- diagram_distance(diag,permuted_diag(D = D,s = c(2,2,3)),dim = 0,p = Inf)
#   d8 <- diagram_distance(diag,permuted_diag(D = D,s = c(3,3,1)),dim = 0,p = Inf)
#   d9 <- diagram_distance(diag,permuted_diag(D = D,s = c(3,3,2)),dim = 0,p = Inf)
#   d10 <- diagram_distance(diag,diag,dim = 0,p = Inf)
#   unique_vals <- unique(c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10))
#   thresholds <- c()
#   for(i in 1:3)
#   {
#     for(j in 1:3)
#     {
#       for(k in 1:3)
#       {
#         thresholds <- c(thresholds,2*stats::quantile(unique_vals[c(i,j,k)],probs = 0.95)[[1]])
#       }
#     }
#   }
#   thresholds <- unique(thresholds)
#   expect_true(bootstrap_persistence_thresholds(X = D,FUN_diag = "ripsDiag",FUN_boot = "ripsDiag",maxdim = 0,thresh = 2.5,num_workers = 2,num_samples = 3)$thresholds %in% thresholds)
# 
#   # check p-values
#   D <- TDA::circleUnif(n = 100,r = 1)
#   bs <- bootstrap_persistence_thresholds(X = D,maxdim = 1,thresh = 2,return_diag = T,num_workers = 2,return_subsetted = T,return_pvals = T)
#   expect_lte(length(bs$pvals),2L)
#   expect_true(bs$pvals[[1]] < 0.1)
#   if(length(bs$pvals) == 2L)
#   {
#     expect_true(bs$pvals[[2]] < 0.1)
#   }
# 
#   bs <- bootstrap_persistence_thresholds(X = D,maxdim = 1,thresh = 2,return_diag = T,num_workers = 2,return_subsetted = T,return_pvals = T,p_less_than_alpha = T,alpha = 1/31)
#   expect_identical(length(bs$pvals),0L)
# 
#   bs <- bootstrap_persistence_thresholds(X = D,FUN_diag = "ripsDiag",FUN_boot = "ripsDiag",maxdim = 1,thresh = 2,return_diag = T,num_workers = 2,return_subsetted = T,return_pvals = T,calculate_representatives = T,alpha = 1/31)
#   expect_lte(length(bs$pvals),2L)
#   expect_true(bs$pvals[[1]] < 0.1)
#   if(length(bs$pvals) == 2L)
#   {
#     expect_true(bs$pvals[[2]] < 0.1)
#   }
#   bs <- bootstrap_persistence_thresholds(X = D,FUN_diag = "ripsDiag",FUN_boot = "ripsDiag",maxdim = 1,thresh = 2,return_diag = T,num_workers = 2,return_subsetted = T,return_pvals = T,calculate_representatives = T,alpha = 1/31,p_less_than_alpha = T)
#   expect_identical(length(bs$subsetted_representatives),0L)
# 
# })
