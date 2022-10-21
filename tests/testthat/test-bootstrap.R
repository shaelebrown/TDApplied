
test_that("bootstrap_persistence_thresholds can do homology calculation with all three packages",{
  
  D <- TDA::circleUnif(n = 50,r = 1)
  # ripser = import_ripser()
  expect_length(bootstrap_persistence_thresholds(X = D,FUN = "calculate_homology",maxdim = 1,thresh = 2.1),3)
  expect_length(bootstrap_persistence_thresholds(X = D,FUN = "ripsDiag",maxdim = 1,thresh = 2.1),3)
  # expect_length(bootstrap_persistence_thresholds(X = D,FUN = "PyH",maxdim = 1,thresh = 2,ripser = ripser),3)
  expect_length(bootstrap_persistence_thresholds(X = D,FUN = "calculate_homology",maxdim = 1,thresh = 2,calculate_representatives = T),3)
  expect_length(bootstrap_persistence_thresholds(X = D,FUN = "ripsDiag",maxdim = 1,thresh = 2,calculate_representatives = T),5)
  # expect_length(bootstrap_persistence_thresholds(X = D,FUN = "PyH",maxdim = 1,thresh = 2,ripser = ripser,calculate_representatives = T),5)
  expect_length(bootstrap_persistence_thresholds(X = as.matrix(dist(D)),FUN = "calculate_homology",maxdim = 1,thresh = 2,distance_mat = T),3)
  expect_length(bootstrap_persistence_thresholds(X = as.matrix(dist(D)),FUN = "ripsDiag",maxdim = 1,thresh = 2,distance_mat = T),3)
  # expect_length(bootstrap_persistence_thresholds(X = as.matrix(dist(D)),FUN = "PyH",maxdim = 1,thresh = 2,ripser = ripser,distance_mat = T),3)
  expect_length(bootstrap_persistence_thresholds(X = as.matrix(dist(D)),FUN = "calculate_homology",maxdim = 1,thresh = 2,calculate_representatives = T,distance_mat = T),3)
  expect_length(bootstrap_persistence_thresholds(X = as.matrix(dist(D)),FUN = "ripsDiag",maxdim = 1,thresh = 2,calculate_representatives = T,distance_mat = T),5)
  # expect_length(bootstrap_persistence_thresholds(X = as.matrix(dist(D)),FUN = "PyH",maxdim = 1,thresh = 2,ripser = ripser,calculate_representatives = T,distance_mat = T),5)
  
})


test_that("bootstrap_persistence_thresholds can detect incorrect parameters correctly",{
  
  # X, FUN, maxdim, thresh, distance_mat, global_threshold, ripser, ignore_infinite_cluster, calculate_representatives, num_samples, alpha, return_subsetted, return_diag
  expect_error(bootstrap_persistence_thresholds(data.frame(),FUN = "calculate_homology",maxdim = 1,thresh = 2),"X")
  expect_error(bootstrap_persistence_thresholds(data.frame(x = 1),FUN = "calculate_homology",maxdim = 1,thresh = 2),"X")
  expect_error(bootstrap_persistence_thresholds(data.frame(x = c(1,2),y = c("1","2")),FUN = "calculate_homology",maxdim = 1,thresh = 2),"X")
  expect_error(bootstrap_persistence_thresholds(X = NULL,FUN = "calculate_homology",maxdim = 1,thresh = 2),"X")
  expect_error(bootstrap_persistence_thresholds(X = data.frame(x = c(1,NA)),FUN = "calculate_homology",maxdim = 1,thresh = 2),"missing")
  expect_error(bootstrap_persistence_thresholds(data.frame(x = c(1),y = c(2)),FUN = "calculate_homology",maxdim = 1,thresh = 2),"two")
  
  expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN = "calculatehomology",maxdim = 1,thresh = 2),"calculate_homology")
  expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN = NULL,maxdim = 1,thresh = 2),"NULL")
  expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN = 2,maxdim = 1,thresh = 2),"string")
  
  expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN = "calculate_homology",maxdim = 1.1,thresh = 2),"whole")
  expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN = "ripsDiag",maxdim = NA,thresh = 2),"NA")
  expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN = "calculate_homology",maxdim = -1,thresh = 2),"negative")
  
  expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN = "calculate_homology",maxdim = 1,thresh = 0),"positive")
  expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN = "ripsDiag",maxdim = 1,thresh = NaN),"NaN")
  expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN = "calculate_homology",maxdim = 1,thresh = "2"),"numeric")
  
  expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN = "calculate_homology",maxdim = 1,thresh = 2,distance_mat = "F"),"logical")
  expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN = "ripsDiag",maxdim = 1,thresh = 2,distance_mat = c(F,T)),"single")
  
  expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN = "calculate_homology",maxdim = 1,thresh = 2,global_threshold = 2),"boolean")
  
  # PyH parameters (ripser,ignore_infinite_cluster,calculate_representatives) are all working, tested in test-python.R
  
  expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN = "ripsDiag",maxdim = 1,thresh = 2,calculate_representatives = c(T,F)),"boolean")
  expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN = "calculate_homology",maxdim = 1,thresh = 2,calculate_representatives = NULL),"NULL")
  
  expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN = "ripsDiag",maxdim = 1,thresh = 2,num_samples = 0),"one")
  expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN = "calculate_homology",maxdim = 1,thresh = 2,num_samples = Inf),"finite")
  
  expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN = "ripsDiag",maxdim = 1,thresh = 2,alpha = NA),"NA")
  expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN = "calculate_homology",maxdim = 1,thresh = 2,alpha = 2),"1")
  expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN = "calculate_homology",maxdim = 1,thresh = 2,alpha = c(0.5,0.4)),"single")
  
  expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN = "calculate_homology",maxdim = 1,thresh = 2,return_subsetted = c(T,F)),"single")
  expect_error(bootstrap_persistence_thresholds(data.frame(x = 1:10,y = 1:10),FUN = "ripsDiag",maxdim = 1,thresh = 2,return_diag = Inf),"logical")
  
})

test_that("bootstrap_persistence_thresholds is computing properly",{
  
  D <- TDA::circleUnif(n = 50,r = 1)
  
  # ripsDiag with global threshold
  bs <- bootstrap_persistence_thresholds(X = D,FUN = "ripsDiag",maxdim = 1,thresh = 2,calculate_representatives = T,return_diag = T)
  expect_length(bs$representatives,nrow(bs$diag))
  expect_length(bs$thresholds,1)
  expect_gt(bs$thresholds,0)
  expect_length(bs$subsetted_representatives,nrow(bs$subsetted_diag))
  expect_true(min(bs$subsetted_diag$death - bs$subsetted_diag$birth) > bs$thresholds)
  
  # ripsDiag with multiple thresholds
  bs <- bootstrap_persistence_thresholds(X = D,FUN = "ripsDiag",maxdim = 1,thresh = 2,calculate_representatives = T,return_diag = T,global_threshold = F)
  expect_length(bs$representatives,nrow(bs$diag))
  expect_length(bs$thresholds,2)
  expect_gt(bs$thresholds[[1]],0)
  expect_gt(bs$thresholds[[2]],0)
  expect_length(bs$subsetted_representatives,nrow(bs$subsetted_diag))
  expect_true(min(bs$subsetted_diag[which(bs$subsetted_diag$dimension == 0),]$death - bs$subsetted_diag[which(bs$subsetted_diag$dimension == 0),]$birth) > bs$thresholds[[1]])
  expect_true(min(bs$subsetted_diag[which(bs$subsetted_diag$dimension == 1),]$death - bs$subsetted_diag[which(bs$subsetted_diag$dimension == 1),]$birth) > bs$thresholds[[2]])
  
  # calculate_homology with global threshold
  bs <- bootstrap_persistence_thresholds(X = D,FUN = "calculate_homology",maxdim = 1,thresh = 2,return_diag = T)
  expect_length(bs$thresholds,1)
  expect_gt(bs$thresholds,0)
  expect_true(min(bs$subsetted_diag$death - bs$subsetted_diag$birth) > bs$thresholds)
  
  # calculate_homology with multiple thresholds
  bs <- bootstrap_persistence_thresholds(X = D,FUN = "calculate_homology",maxdim = 1,thresh = 2,return_diag = T,global_threshold = F)
  expect_length(bs$thresholds,2)
  expect_gt(bs$thresholds[[1]],0)
  expect_gt(bs$thresholds[[2]],0)
  expect_true(min(bs$subsetted_diag[which(bs$subsetted_diag$dimension == 1),]$death - bs$subsetted_diag[which(bs$subsetted_diag$dimension == 1),]$birth) > bs$thresholds[[2]])
  
  # ripser = import_ripser()
  # # PyH with global threshold
  # bs <- bootstrap_persistence_thresholds(X = D,FUN = "PyH",maxdim = 1,thresh = 2,calculate_representatives = T,return_diag = T,ripser = ripser)
  # expect_length(bs$representatives[[2]],length(which(bs$diag$dimension == 1)))
  # expect_length(bs$thresholds,1)
  # expect_gt(bs$thresholds,0)
  # expect_length(bs$subsetted_representatives[[2]],nrow(bs$subsetted_diag))
  # expect_true(min(bs$subsetted_diag$death - bs$subsetted_diag$birth) > bs$thresholds)
  
  # # PyH with multiple thresholds
  # bs <- bootstrap_persistence_thresholds(X = D,FUN = "ripsDiag",maxdim = 1,thresh = 2,calculate_representatives = T,return_diag = T,global_threshold = F,ripser = ripser)
  # expect_length(bs$representatives,nrow(bs$diag))
  # expect_length(bs$thresholds,2)
  # expect_gt(bs$thresholds[[1]],0)
  # expect_gt(bs$thresholds[[2]],0)
  # expect_length(bs$subsetted_representatives,nrow(bs$subsetted_diag))
  # expect_true(min(bs$subsetted_diag[which(bs$subsetted_diag$dimension == 0),]$death - bs$subsetted_diag[which(bs$subsetted_diag$dimension == 0),]$birth) > bs$thresholds[[1]])
  # expect_true(min(bs$subsetted_diag[which(bs$subsetted_diag$dimension == 1),]$death - bs$subsetted_diag[which(bs$subsetted_diag$dimension == 1),]$birth) > bs$thresholds[[2]])
  
  # check on circle:
  bs <- bootstrap_persistence_thresholds(X = D,FUN = "ripsDiag",maxdim = 1,thresh = 2,return_diag = T)
  expect_equal(bs$subsetted_diag$dimension,c(0,1))
  bs <- bootstrap_persistence_thresholds(X = D,FUN = "calculate_homology",maxdim = 1,thresh = 2,return_diag = T)
  expect_equal(bs$subsetted_diag$dimension,c(1))
  # bs <- bootstrap_persistence_thresholds(X = D,FUN = "PyH",maxdim = 1,thresh = 2,return_diag = T,ripser = ripser)
  # expect_equal(bs$subsetted_diag$dimension,c(1))
  # bs <- bootstrap_persistence_thresholds(X = D,FUN = "PyH",maxdim = 1,thresh = 2,return_diag = T,ripser = ripser,ignore_infinite_cluster = F)
  # expect_equal(bs$subsetted_diag$dimension,c(0,1))
  bs <- bootstrap_persistence_thresholds(X = D,FUN = "ripsDiag",maxdim = 1,thresh = 2,return_diag = T,global_threshold = F)
  expect_equal(bs$subsetted_diag$dimension,c(0,1))
  bs <- bootstrap_persistence_thresholds(X = D,FUN = "calculate_homology",maxdim = 1,thresh = 2,return_diag = T,global_threshold = F)
  expect_equal(bs$subsetted_diag$dimension,c(1))
  # bs <- bootstrap_persistence_thresholds(X = D,FUN = "PyH",maxdim = 1,thresh = 2,return_diag = T,ripser = ripser,global_threshold = F)
  # expect_equal(bs$subsetted_diag$dimension,c(1))
  # bs <- bootstrap_persistence_thresholds(X = D,FUN = "PyH",maxdim = 1,thresh = 2,return_diag = T,ripser = ripser,ignore_infinite_cluster = F,global_threshold = F)
  # expect_equal(bs$subsetted_diag$dimension,c(0,1))
  
  # one example by hand to verify thresholds
  D <- data.frame(x = c(0,0,0),y = c(1,2,2.5))
  diag <- diagram_to_df(TDA::ripsDiag(D,maxdimension = 0,maxscale = 1,library = "dionysus"))
  d1 <- diagram_distance(diag,rbind(diag[1,],diag[1,],diag[1,]),dim = 0,p = Inf)
  d2 <- diagram_distance(diag,rbind(diag[2,],diag[2,],diag[2,]),dim = 0,p = Inf)
  d3 <- diagram_distance(diag,rbind(diag[3,],diag[3,],diag[3,]),dim = 0,p = Inf)
  d4 <- diagram_distance(diag,rbind(diag[1,],diag[1,],diag[2,]),dim = 0,p = Inf)
  d5 <- diagram_distance(diag,rbind(diag[1,],diag[1,],diag[3,]),dim = 0,p = Inf)
  d6 <- diagram_distance(diag,rbind(diag[2,],diag[1,],diag[2,]),dim = 0,p = Inf)
  d7 <- diagram_distance(diag,rbind(diag[2,],diag[3,],diag[2,]),dim = 0,p = Inf)
  d8 <- diagram_distance(diag,rbind(diag[3,],diag[1,],diag[3,]),dim = 0,p = Inf)
  d9 <- diagram_distance(diag,rbind(diag[3,],diag[2,],diag[3,]),dim = 0,p = Inf)
  d10 <- diagram_distance(diag,diag,dim = 0,p = Inf)
  expect_equal(bootstrap_persistence_thresholds(X = D,FUN = "ripsDiag",maxdim = 0,thresh = 1)$thresholds,1)
  
})

