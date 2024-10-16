
# test_that("permutation_test detects incorrect parameters correctly",{
# 
#   skip_if_not_installed("TDA")
#   skip_if_not_installed("TDAstats")
#   circle = TDA::ripsDiag(X = TDA::circleUnif(n = 20,r = 1),maxdimension = 2,maxscale = 2)
#   sphere = TDA::ripsDiag(X = TDA::sphereUnif(n = 20,d = 2,r = 1),maxdimension = 2,maxscale = 2)
#   expect_error(permutation_test(list(circle,2,circle),list(sphere,sphere,sphere),iterations = 5,p = 2,q = 2,dims = c(0,1),paired = F,distance = "wasserstein",sigma = NULL,verbose = F,num_workers = 2),"Diagrams must")
#   expect_error(permutation_test(list(circle,circle,circle),list(sphere,"2",sphere),iterations = 5,p = 2,q = 2,dims = c(0,1),paired = F,distance = "wasserstein",sigma = NULL,verbose = F,num_workers = 2),"Diagrams must")
#   expect_error(permutation_test(list(circle,circle,circle),list(sphere,sphere,sphere),iterations = NA,p = 2,q = 2,dims = c(0,1),paired = F,distance = "wasserstein",sigma = NULL,verbose = F,num_workers = 2),"NA")
#   expect_error(permutation_test(list(circle,circle,circle),list(sphere,sphere,sphere),iterations = 5,p = NULL,q = 2,dims = c(0,1),paired = F,distance = "wasserstein",sigma = NULL,verbose = F,num_workers = 2),"NULL")
#   expect_error(permutation_test(list(circle,circle,circle),list(sphere,sphere,sphere),iterations = 5,p = 2,q = NA,dims = c(0,1),paired = F,distance = "wasserstein",sigma = NULL,verbose = F,num_workers = 2),"NA")
#   expect_error(permutation_test(list(circle,circle,circle),list(sphere,sphere,sphere),iterations = 5,p = 2,q = 2,dims = c(-1,1),paired = F,distance = "wasserstein",sigma = NULL,verbose = F,num_workers = 2),"non-negative")
#   expect_error(permutation_test(list(circle,circle),list(sphere,sphere,sphere),iterations = 5,p = 2,q = 2,dims = c(0,1),paired = T,distance = "wasserstein",sigma = NULL,verbose = F,num_workers = 2),"same number")
#   expect_error(permutation_test(list(circle,circle,circle),list(sphere,sphere,sphere),iterations = 5,p = 2,q = 2,dims = c(0,1),paired = F,distance = "Wasserstein",sigma = NULL,verbose = F,num_workers = 2),"wasserstein")
#   expect_error(permutation_test(list(circle,circle,circle),list(sphere,sphere,sphere),iterations = 5,p = 2,q = 2,dims = c(0,1),paired = F,distance = "fisher",sigma = NULL,verbose = F,num_workers = 2),"sigma")
#   expect_error(permutation_test(list(circle,circle,circle),list(sphere,sphere),iterations = 5,p = 2,q = 2,dims = c(0,1),paired = T,distance = "fisher",sigma = 1,verbose = F,num_workers = 2),"paired")
#   expect_error(permutation_test(list(circle,circle,circle),list(sphere,sphere),num_workers = NA),"num_workers")
#   expect_error(permutation_test(list(circle,circle,circle),list(sphere,sphere),num_workers = NULL),"num_workers")
#   expect_error(permutation_test(list(circle,circle,circle),list(sphere),num_workers = 2),"2")
#   expect_error(permutation_test(list(circle,circle,circle),list(sphere,sphere),num_workers = 2,dims = c(1),distance = "fisher",sigma = 1,rho = NaN),"rho")
# 
# })

# test_that("permutation_test can accept inputs from TDA, TDAstats and diagram_to_df",{
# 
#   skip_if_not_installed("TDA")
#   skip_if_not_installed("TDAstats")
#   circle = TDA::ripsDiag(X = TDA::circleUnif(n = 20,r = 1),maxdimension = 1,maxscale = 2)
#   sphere = TDAstats::calculate_homology(TDA::sphereUnif(n = 20,d = 2,r = 1),threshold = 2)
#   expect_length(permutation_test(list(circle,circle,diagram_to_df(circle)),list(sphere,sphere,sphere),iterations = 1,dims = c(1),num_workers = 2)$permvals[[1]],1)
# 
# })

# test_that("permutation_test can accept precomputed distance matrices",{
# 
#   skip_if_not_installed("TDA")
#   skip_if_not_installed("TDAstats")
# 
#   circle = TDA::ripsDiag(X = TDA::circleUnif(n = 20,r = 1),maxdimension = 1,maxscale = 2)
#   sphere = TDAstats::calculate_homology(TDA::sphereUnif(n = 20,d = 2,r = 1),threshold = 2)
#   D0 = distance_matrix(diagrams = list(circle,sphere),dim = 0,num_workers = 2)
#   D1 = distance_matrix(diagrams = list(circle,sphere),dim = 1,num_workers = 2)
#   D2 = distance_matrix(diagrams = list(circle,sphere,circle,sphere),dim = 0,num_workers = 2)
#   D3 = distance_matrix(diagrams = list(circle,sphere,circle,sphere),dim = 1,num_workers = 2)
#   expect_error(permutation_test(dist_mats = list(D0,D1),iterations = 1,dims = c(1),num_workers = 2,group_sizes = NULL),"group_sizes")
#   expect_error(permutation_test(dist_mats = list(D0,D1),iterations = 1,dims = c(1),num_workers = 2,group_sizes = c(1,1)),"dims")
#   expect_error(permutation_test(dist_mats = list(D0,D1),iterations = 1,dims = c(0,1),num_workers = 2,group_sizes = c(1,2)),"sum")
#   expect_error(permutation_test(dist_mats = NULL,iterations = 1,dims = c(1),num_workers = 2,group_sizes = c(1,1)),"list")
#   expect_error(permutation_test(dist_mats = list(D0,matrix(data = c(0),nrow = 1)),iterations = 1,dims = c(0,1),num_workers = 2,group_sizes = c(1,1)),"size")
#   expect_error(permutation_test(dist_mats = list(D2,D3),iterations = 1,dims = c(0,1),paired = T,num_workers = 2,group_sizes = c(3,1)),"paired")
# 
#   D0 = distance_matrix(diagrams = list(circle,circle,sphere,sphere),dim = 0,num_workers = 2)
#   D1 = distance_matrix(diagrams = list(circle,circle,sphere,sphere),dim = 1,num_workers = 2)
#   expect_length(permutation_test(dist_mats = list(D0,D1),group_sizes = c(2,2),iterations = 1,dims = c(0,1),paired = T,num_workers = 2)$test_statistics,2)
#   expect_length(permutation_test(dist_mats = list(D0,D1),group_sizes = c(2,2),iterations = 3,dims = c(0,1),paired = T,num_workers = 2)$permvals[[1]],3)
#   expect_length(permutation_test(dist_mats = list(D0,D1),group_sizes = c(2,2),iterations = 1,dims = c(0,1),paired = F,num_workers = 2)$test_statistics,2)
#   expect_length(permutation_test(dist_mats = list(D0,D1),group_sizes = c(2,2),iterations = 3,dims = c(0,1),paired = F,num_workers = 2)$permvals[[1]],3)
# 
#   circle2 <- TDA::ripsDiag(X = TDA::circleUnif(n = 50,r = 1),maxdimension = 2,maxscale = 2,library = "dionysus",location = T)
#   sphere2 <- TDA::ripsDiag(X = TDA::sphereUnif(n = 50,d = 2,r = 1),maxdimension = 2,maxscale = 2,library = "dionysus",location = T)
#   d <- diagram_distance(circle2,sphere2,dim = 1) # wasserstein distance in dimension 1
#   D1 <- distance_matrix(diagrams = list(circle2,circle2,circle2,circle2,circle2,sphere2),dim = 1,num_workers = 2)
#   D2 <- distance_matrix(diagrams = list(circle2,circle2,circle2,sphere2,sphere2,sphere2),dim = 2,num_workers = 2)
#   expect_length(permutation_test(dist_mats = list(D1),iterations = 1,dims = c(1),num_workers = 2,group_sizes = c(3,3))$permvals[[1]],1)
#   expect_length(permutation_test(dist_mats = list(D2),iterations = 1,dims = c(2),num_workers = 2,group_sizes = c(3,3))$permvals[[1]],1)
#   expect_length(permutation_test(dist_mats = list(D1,D2),iterations = 1,dims = c(1,2),num_workers = 2,group_sizes = c(3,3))$permvals,2)
#   expect_lte(abs(permutation_test(dist_mats = list(D1),iterations = 1,dims = c(1),num_workers = 2,group_sizes = c(3,3))$test_statistics[[1]]-d^2/3),0.003)
#   expect_equal(permutation_test(dist_mats = list(D2),iterations = 1,dims = c(2),num_workers = 2,group_sizes = c(3,3))$test_statistics[[1]],0)
# 
#   expect_lte(abs(permutation_test(dist_mats = list(D1),iterations = 1,dims = c(1),num_workers = 2,group_sizes = c(3,3))$permvals[[1]][[1]]-d^2/3),0.003)
#   v <- permutation_test(dist_mats = list(distance_matrix(diagrams = list(sphere2,sphere2,circle2,circle2,circle2,sphere2),dim = 1,num_workers = 2)),iterations = 1,dims = c(1),num_workers = 2,group_sizes = c(3,3))$permvals[[1]][[1]]
#   expect_lte(abs(v - 2*d^2/3)*v,0.002)
# 
# })

# test_that("permutation_test is working",{
# 
#   skip_if_not_installed("TDA")
#   skip_if_not_installed("TDAstats")
# 
#   circle <- TDA::ripsDiag(X = TDA::circleUnif(n = 50,r = 1),maxdimension = 2,maxscale = 2)
#   sphere <- TDA::ripsDiag(X = TDA::sphereUnif(n = 50,d = 2,r = 1),maxdimension = 2,maxscale = 2)
#   circle2 <- TDA::ripsDiag(X = TDA::circleUnif(n = 50,r = 1),maxdimension = 2,maxscale = 2,library = "dionysus",location = T)
#   sphere2 <- TDA::ripsDiag(X = TDA::sphereUnif(n = 50,d = 2,r = 1),maxdimension = 2,maxscale = 2,library = "dionysus",location = T)
#   d <- diagram_distance(circle,sphere,dim = 1) # wasserstein distance in dimension 1
#   # one check to make sure when cycle locations are calculated there are no errors
#   expect_length(permutation_test(list(circle2,circle2,circle2),list(sphere2,sphere2,sphere2),iterations = 1,dims = c(1),num_workers = 2)$permvals[[1]],1)
#   expect_length(permutation_test(list(circle,circle,circle),list(sphere,sphere,sphere),iterations = 3,dims = c(1),num_workers = 2)$permvals[[1]],3)
#   expect_length(permutation_test(list(circle,circle,circle),list(sphere,sphere,sphere),iterations = 5,dims = c(1),num_workers = 2)$permvals[[1]],5)
#   expect_equal(permutation_test(list(circle,circle,circle),list(sphere,sphere,sphere),iterations = 1,dims = c(1),num_workers = 2)$test_statistics[[1]],0)
#   expect_equal(permutation_test(list(circle,circle,circle),list(sphere,sphere,sphere),iterations = 1,dims = c(2),num_workers = 2)$test_statistics[[1]],0)
#   expect_lte(abs(permutation_test(list(circle,circle,circle),list(sphere,circle,circle),iterations = 1,dims = c(1),num_workers = 2)$permvals[[1]][[1]] - d^2/3),0.01)
#   v <- permutation_test(list(circle,sphere,circle),list(sphere,circle,sphere),iterations = 1,dims = c(1),num_workers = 2)$permvals[[1]][[1]]
#   expect_lte(abs(v - 2*d^2/3)*v,0.02)
#   expect_length(unique(permutation_test(list(circle,sphere,circle),list(circle,sphere,circle),paired = T,iterations = 3,dims = c(1),num_workers = 2)$permvals[[1]]),1)
#   expect_length(unique(permutation_test(list(sphere,sphere,circle),list(sphere,sphere,circle),paired = T,iterations = 3,dims = c(1),num_workers = 2)$permvals[[1]]),1)
#   # expect_length(permutation_test(list(sphere,sphere,circle),list(sphere,sphere,circle),paired = T,iterations = 3,dims = c(1),num_workers = 2,distance = "fisher",sigma = 1,rho = 0.001)$permvals[[1]],3)
# 
# })

# test_that("independence_test detects incorrect parameters correctly",{
# 
#   skip_if_not_installed("TDA")
#   skip_if_not_installed("TDAstats")
#   g1 <- lapply(X = 1:6,FUN = function(X){return(TDA::ripsDiag(X = TDA::circleUnif(n = 50,r = 1),maxdimension = 1,maxscale = 2))})
#   g2 <- lapply(X = 1:6,FUN = function(X){return(TDA::ripsDiag(X = TDA::sphereUnif(n = 50,d = 2,r = 1),maxdimension = 1,maxscale = 2))})
#   expect_error(independence_test(g1,g2,dims = c(0,1),sigma = 1,t = NA,num_workers = 2),"NA")
#   expect_error(independence_test(g1,g2,dims = c(0,1),sigma = NaN,t = 1,num_workers = 2),"NaN")
#   expect_error(independence_test(g1,g2,dims = c(0,1),sigma = c(1,2),t = 1,num_workers = 2),"single")
#   expect_error(independence_test(g1,g2,dims = c(0,1.1),sigma = 1,t = 1,num_workers = 2),"whole")
#   expect_error(independence_test(g1,g2[1:5],dims = c(0,1),sigma = 1,t = 1,num_workers = 2),"same length")
#   expect_error(independence_test(g1[1:5],g2[1:5],dims = c(0,1),sigma = 1,t = 1,num_workers = 2),"6")
#   # expect_error(independence_test(list(g1[[1]],g1[[2]],g1[[3]],g1[[4]],g1[[5]],data.frame(dimension = numeric(),birth = numeric(),death = numeric())),g2,dims = c(0,1),sigma = 1,t = 1,num_workers = 2),"empty")
#   expect_error(independence_test(g1,g2,dims = c(0,1),num_workers = 2.3),"num_workers")
#   expect_error(independence_test(g1,g2,dims = c(0,1),num_workers = -2),"num_workers")
# 
# })

# test_that("independence_test can accept inputs from TDA, TDAstats and diagram_to_df",{
# 
#   skip_if_not_installed("TDA")
#   skip_if_not_installed("TDAstats")
#   circle = TDA::ripsDiag(X = TDA::circleUnif(n = 20,r = 1),maxdimension = 1,maxscale = 2)
#   sphere = TDAstats::calculate_homology(TDA::sphereUnif(n = 20,d = 2,r = 1),threshold = 2)
#   expect_length(independence_test(list(circle,circle,diagram_to_df(circle),circle,circle,circle),list(sphere,sphere,sphere,sphere,sphere,circle),dims = c(1),num_workers = 2)$p_values,1)
# 
# })

# test_that("independence_test can accept precomputed Gram matrices",{
# 
#   skip_if_not_installed("TDA")
#   skip_if_not_installed("TDAstats")
#   circle = TDA::ripsDiag(X = TDA::circleUnif(n = 20,r = 1),maxdimension = 1,maxscale = 2)
#   sphere = TDAstats::calculate_homology(TDA::sphereUnif(n = 20,d = 2,r = 1),threshold = 2)
#   K0 = gram_matrix(list(circle,circle,circle,circle,circle,circle),dim = 0,num_workers = 2)
#   K1 = gram_matrix(list(circle,circle,circle,circle,circle,circle),dim = 1,num_workers = 2)
#   L0 = gram_matrix(list(sphere,sphere,sphere,sphere,sphere,sphere),dim = 0,num_workers = 2)
#   L1 = gram_matrix(list(sphere,sphere,sphere,sphere,sphere,sphere),dim = 1,num_workers = 2)
#   expect_length(independence_test(Ks = list(K0,K1),Ls = list(L0,L1),dims = c(0,1),num_workers = 2)$p_values,2)
#   expect_length(independence_test(g1 = list(circle,circle,circle,circle,circle,circle),g2 = list(sphere,sphere,sphere,sphere,sphere,sphere),dims = c(0,1),num_workers = 2)$p_values,2)
#   expect_identical(independence_test(Ks = list(K0),Ls = list(L0),dims = c(0),num_workers = 2)$p_values,independence_test(g1 = list(circle,circle,circle,circle,circle,circle),g2 = list(sphere,sphere,sphere,sphere,sphere,sphere),dims = c(0),num_workers = 2)$p_values)
# 
#   expect_error(independence_test(Ks = list(K0,K1),Ls = 2,dims = c(0,1),num_workers = 2),"list")
#   expect_error(independence_test(Ks = list(K0,K1),Ls = list(L0),dims = c(0,1),num_workers = 2),"length")
#   K2 = gram_matrix(list(circle,circle,circle,circle,circle,circle,circle),dim = 1,num_workers = 2)
#   K3 = gram_matrix(list(circle,circle,circle,circle,circle),dim = 1,num_workers = 2)
#   expect_error(independence_test(Ks = list(K0,K2),Ls = list(L0,L1),dims = c(0,1),num_workers = 2),"dim")
#   expect_error(independence_test(Ks = list(K0,K1),Ls = list(L0,K3),dims = c(0,1),num_workers = 2),"6")
# 
# })

test_that("independence_test is working",{
  
  D1 <- data.frame(dimension = 0,birth = 2,death = 3)
  D2 <- data.frame(dimension = 0,birth = 2,death = 3.1)
  D3 <- data.frame(dimension = 0,birth = c(2,5),death = c(3.1,6))
  k12 <- diagram_kernel(D1,D2,dim = 0) # sigma = 1,t = 1
  k13 <- diagram_kernel(D1,D3,dim = 0)
  k23 <- diagram_kernel(D3,D2,dim = 0)
  HSIC <- 100*(1 - k13)*(1 - k23)/36^2
  mu_x_sq <- (10 + 5*k13)/15
  mu_y_sq <- (10 + 5*k23)/15
  mu <- (1 + mu_x_sq*mu_y_sq - mu_x_sq - mu_y_sq)/6
  v <- 209*((k13 - 1)^2)*((k23 - 1)^2)/314928
  expect_equal(independence_test(g1 = list(D1,D1,D1,D1,D1,D3),g2 = list(D2,D2,D2,D2,D2,D3),dims = c(0),num_workers = 2)$test_statistics[[1]],HSIC)
  expect_equal(independence_test(g1 = list(D1,D1,D1,D1,D1,D3),g2 = list(D2,D2,D2,D2,D2,D3),dims = c(0),num_workers = 2)$p_values[[1]],stats::pgamma(q = HSIC,rate = mu/v,shape = mu^2/v,lower.tail = F))
  
})

# test_that("universal_null can detect incorrect parameters",{
# 
#     skip_if_not_installed("TDA")
#     skip_if_not_installed("TDAstats")
# 
#     library(TDA)
# 
#     # X, FUN_diag, maxdim, thresh, distance_mat, ripser, ignore_infinite_cluster, calculate_representatives, alpha, return_pvals, infinite_cycle_inference
#     expect_error(universal_null(data.frame(),maxdim = 1,thresh = 2),"X")
#     expect_error(universal_null(data.frame(x = 1),maxdim = 1,thresh = 2),"X")
#     expect_error(universal_null(data.frame(x = c(1,2),y = c("1","2")),maxdim = 1,thresh = 2),"X")
#     expect_error(universal_null(X = NULL,maxdim = 1,thresh = 2),"X")
#     expect_error(universal_null(X = data.frame(x = c(1,NA)),maxdim = 1,thresh = 2),"missing")
#     expect_error(universal_null(data.frame(x = c(1),y = c(2)),maxdim = 1,thresh = 2),"two")
# 
#     expect_error(universal_null(data.frame(x = 1:10,y = 1:10),FUN_diag = "calculatehomology",maxdim = 1,thresh = 2),"calculate_homology")
#     expect_error(universal_null(data.frame(x = 1:10,y = 1:10),FUN_diag = NULL,maxdim = 1,thresh = 2),"NULL")
#     expect_error(universal_null(data.frame(x = 1:10,y = 1:10),FUN_diag = 2,maxdim = 1,thresh = 2),"string")
# 
#     expect_error(universal_null(data.frame(x = 1:10,y = 1:10),maxdim = 1.1,thresh = 2),"whole")
#     expect_error(universal_null(data.frame(x = 1:10,y = 1:10),maxdim = NA,thresh = 2),"NA")
#     expect_error(universal_null(data.frame(x = 1:10,y = 1:10),maxdim = -1,thresh = 2),"negative")
# 
#     expect_error(universal_null(data.frame(x = 1:10,y = 1:10),maxdim = 1,thresh = 0),"positive")
#     expect_error(universal_null(data.frame(x = 1:10,y = 1:10),maxdim = 1,thresh = NaN),"NaN")
#     expect_error(universal_null(data.frame(x = 1:10,y = 1:10),maxdim = 1,thresh = "2"),"numeric")
# 
#     expect_error(universal_null(data.frame(x = 1:10,y = 1:10),FUN_diag = "calculate_homology",maxdim = 1,thresh = 2,distance_mat = "F"),"logical")
#     expect_error(universal_null(data.frame(x = 1:10,y = 1:10),FUN_diag = "ripsDiag",maxdim = 1,thresh = 2,distance_mat = c(F,T)),"single")
#     expect_error(universal_null(data.frame(x = 1:10,y = 1:10),FUN_diag = "ripsDiag",maxdim = 1,thresh = 2,distance_mat = NULL),"NULL")
#     expect_error(universal_null(data.frame(x = 1:10,y = 1:10),FUN_diag = "ripsDiag",maxdim = 1,thresh = 2,distance_mat = NA),"NA")
#     expect_error(universal_null(data.frame(x = 1:10,y = 1:10),FUN_diag = "ripsDiag",maxdim = 1,thresh = 2,distance_mat = T),"square")
# 
#     expect_error(universal_null(data.frame(x = 1:10,y = 1:10),FUN_diag = "calculate_homology",maxdim = 1,thresh = 2,return_pvals = "F"),"logical")
#     expect_error(universal_null(data.frame(x = 1:10,y = 1:10),FUN_diag = "ripsDiag",maxdim = 1,thresh = 2,return_pvals = c(F,T)),"single")
#     expect_error(universal_null(data.frame(x = 1:10,y = 1:10),FUN_diag = "ripsDiag",maxdim = 1,thresh = 2,return_pvals = NULL),"NULL")
#     expect_error(universal_null(data.frame(x = 1:10,y = 1:10),FUN_diag = "ripsDiag",maxdim = 1,thresh = 2,return_pvals = NA),"NA")
# 
#     expect_error(universal_null(data.frame(x = 1:10,y = 1:10),FUN_diag = "ripsDiag",maxdim = 1,thresh = 2,calculate_representatives = c(T,F)),"boolean")
#     expect_error(universal_null(data.frame(x = 1:10,y = 1:10),FUN_diag = "calculate_homology",maxdim = 1,thresh = 2,calculate_representatives = NULL),"NULL")
#     expect_error(universal_null(data.frame(x = 1:10,y = 1:10),FUN_diag = "calculate_homology",maxdim = 1,thresh = 2,calculate_representatives = NA),"NA")
# 
#     expect_error(universal_null(data.frame(x = 1:10,y = 1:10),FUN_diag = "ripsDiag",maxdim = 1,thresh = 2,alpha = NA),"NA")
#     expect_error(universal_null(data.frame(x = 1:10,y = 1:10),FUN_diag = "calculate_homology",maxdim = 1,thresh = 2,alpha = 2),"1")
#     expect_error(universal_null(data.frame(x = 1:10,y = 1:10),FUN_diag = "calculate_homology",maxdim = 1,thresh = 2,alpha = c(0.5,0.4)),"single")
# 
#     expect_error(universal_null(data.frame(x = 1:10,y = 1:10),FUN_diag = "ripsDiag",maxdim = 1,thresh = 2,infinite_cycle_inference = c(T,F)),"boolean")
#     expect_error(universal_null(data.frame(x = 1:10,y = 1:10),FUN_diag = "calculate_homology",maxdim = 1,thresh = 2,infinite_cycle_inference = NULL),"NULL")
#     expect_error(universal_null(data.frame(x = 1:10,y = 1:10),FUN_diag = "calculate_homology",maxdim = 1,thresh = 2,infinite_cycle_inference = NA),"NA")
# 
# })

# test_that("universal_null is working properly",{
# 
#   theta <- runif(n = 100, min = 0,max = 2*3.14)
#   x <- cos(theta)
#   y <- sin(theta)
#   circ <- data.frame(x = x,y = y)
#   res <- universal_null(circ, maxdim = 1, thresh = 2)
#   expect_equal(nrow(res$subsetted_diag), 0) # when there is only one feature, same for TDA::sphereUnif
# 
#   circ$x <- circ$x + rnorm(n = 100,sd = 0.1)
#   circ$y <- circ$y + rnorm(n = 100,sd = 0.1)
#   res <- universal_null(circ, maxdim = 1, thresh = 2,alpha = 0.1) # with enough leeway
#   expect_equal(nrow(res$subsetted_diag), 1L)
#   expect_equal(res$subsetted_diag[1,1L], 1) # dim 1
# 
#   # now trying with and without infinite cycle inference at a smaller radius
#   res <- universal_null(circ, maxdim = 1, thresh = 1.1,alpha = 0.05)
#   expect_equal(nrow(res$subsetted_diag), 0L)
#   library(TDA)
#   res <- universal_null(circ, FUN_diag = 'ripsDiag', maxdim = 1, thresh = 1.1,alpha = 0.05,infinite_cycle_inference = T)
#   expect_equal(nrow(res$subsetted_diag), 1L)
# 
# })

# test_that("universal_null subsets representatives properly",{
#   
#     # circle has only representative
#     theta <- runif(n = 100, min = 0,max = 2*3.14)
#     x <- cos(theta)
#     y <- sin(theta)
#     circ <- data.frame(x = x,y = y)
#     circ$x <- circ$x + rnorm(n = 100,sd = 0.1)
#     circ$y <- circ$y + rnorm(n = 100,sd = 0.1)
#     library(TDA)
#     res <- universal_null(circ, FUN_diag = "ripsDiag",maxdim = 1, thresh = 2,alpha = 0.1,calculate_representatives = TRUE,return_pvals = TRUE) # with enough leeway
#     expect_equal(length(res$subsetted_representatives),1L)
#     expect_equal(length(res$pvals),1L)
#     
#     # circle with infinite cycle still has one representative
#     res2 <- universal_null(circ, FUN_diag = "ripsDiag",maxdim = 1, thresh = 1.1,alpha = 0.1,calculate_representatives = TRUE,infinite_cycle_inference = TRUE,return_pvals = TRUE) # with enough leeway
#     expect_equal(length(res2$subsetted_representatives),0L)
#     expect_equal(length(res2$pvals),1L)
#     
# })

