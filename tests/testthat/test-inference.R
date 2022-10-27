
test_that("permutation_test detects incorrect parameters correctly",{
  
  circle = TDA::ripsDiag(X = TDA::circleUnif(n = 100,r = 1),maxdimension = 2,maxscale = 2)
  sphere = TDA::ripsDiag(X = TDA::sphereUnif(n = 100,d = 2,r = 1),maxdimension = 2,maxscale = 2)
  expect_error(permutation_test(list(circle,2,circle),list(sphere,sphere,sphere),iterations = 5,p = 2,q = 2,dims = c(0,1),paired = F,distance = "wasserstein",sigma = NULL,verbose = F,num_workers = 2),"Diagrams must")
  expect_error(permutation_test(list(circle,circle,circle),list(sphere,"2",sphere),iterations = 5,p = 2,q = 2,dims = c(0,1),paired = F,distance = "wasserstein",sigma = NULL,verbose = F,num_workers = 2),"Diagrams must")
  expect_error(permutation_test(list(circle,circle,circle),list(sphere,sphere,sphere),iterations = NA,p = 2,q = 2,dims = c(0,1),paired = F,distance = "wasserstein",sigma = NULL,verbose = F,num_workers = 2),"NA")
  expect_error(permutation_test(list(circle,circle,circle),list(sphere,sphere,sphere),iterations = 5,p = NULL,q = 2,dims = c(0,1),paired = F,distance = "wasserstein",sigma = NULL,verbose = F,num_workers = 2),"NULL")
  expect_error(permutation_test(list(circle,circle,circle),list(sphere,sphere,sphere),iterations = 5,p = 2,q = NA,dims = c(0,1),paired = F,distance = "wasserstein",sigma = NULL,verbose = F,num_workers = 2),"NA")
  expect_error(permutation_test(list(circle,circle,circle),list(sphere,sphere,sphere),iterations = 5,p = 2,q = 2,dims = c(-1,1),paired = F,distance = "wasserstein",sigma = NULL,verbose = F,num_workers = 2),"non-negative")
  expect_error(permutation_test(list(circle,circle),list(sphere,sphere,sphere),iterations = 5,p = 2,q = 2,dims = c(0,1),paired = T,distance = "wasserstein",sigma = NULL,verbose = F,num_workers = 2),"same number")
  expect_error(permutation_test(list(circle,circle,circle),list(sphere,sphere,sphere),iterations = 5,p = 2,q = 2,dims = c(0,1),paired = F,distance = "Wasserstein",sigma = NULL,verbose = F,num_workers = 2),"wasserstein")
  expect_error(permutation_test(list(circle,circle,circle),list(sphere,sphere,sphere),iterations = 5,p = 2,q = 2,dims = c(0,1),paired = F,distance = "fisher",sigma = NULL,verbose = F,num_workers = 2),"sigma")
  expect_error(permutation_test(list(circle,circle,circle),list(sphere,sphere),iterations = 5,p = 2,q = 2,dims = c(0,1),paired = T,distance = "fisher",sigma = 1,verbose = F,num_workers = 2),"paired")
  expect_error(permutation_test(list(circle,circle,circle),list(sphere,sphere),num_workers = NA),"num_workers")
  expect_error(permutation_test(list(circle,circle,circle),list(sphere,sphere),num_workers = NULL),"num_workers")
  expect_error(permutation_test(list(circle,circle,circle),list(sphere),num_workers = 2),"2")
  
})

test_that("permutation_test can accept inputs from TDA, TDAstats and diagram_to_df",{
  
  circle = TDA::ripsDiag(X = TDA::circleUnif(n = 20,r = 1),maxdimension = 1,maxscale = 2)
  sphere = TDAstats::calculate_homology(TDA::sphereUnif(n = 20,d = 2,r = 1),threshold = 2)
  expect_length(permutation_test(list(circle,circle,diagram_to_df(circle)),list(sphere,sphere,sphere),iterations = 1,dims = c(1),num_workers = 2)$permvals[[1]],1)
  
})

test_that("permutation_test is working",{
  
  circle <- TDA::ripsDiag(X = TDA::circleUnif(n = 50,r = 1),maxdimension = 2,maxscale = 2)
  sphere <- TDA::ripsDiag(X = TDA::sphereUnif(n = 50,d = 2,r = 1),maxdimension = 2,maxscale = 2)
  circle2 <- TDA::ripsDiag(X = TDA::circleUnif(n = 50,r = 1),maxdimension = 2,maxscale = 2,library = "dionysus",location = T)
  sphere2 <- TDA::ripsDiag(X = TDA::sphereUnif(n = 50,d = 2,r = 1),maxdimension = 2,maxscale = 2,library = "dionysus",location = T)
  d <- diagram_distance(circle,sphere,dim = 1) # wasserstein distance in dimension 1
  # one check to make sure when cycle locations are calculated there are no errors
  expect_length(permutation_test(list(circle2,circle2,circle2),list(sphere2,sphere2,sphere2),iterations = 1,dims = c(1),num_workers = 2)$permvals[[1]],1)
  expect_length(permutation_test(list(circle,circle,circle),list(sphere,sphere,sphere),iterations = 3,dims = c(1),num_workers = 2)$permvals[[1]],3)
  expect_length(permutation_test(list(circle,circle,circle),list(sphere,sphere,sphere),iterations = 5,dims = c(1),num_workers = 2)$permvals[[1]],5)
  expect_equal(permutation_test(list(circle,circle,circle),list(sphere,sphere,sphere),iterations = 1,dims = c(1),num_workers = 2)$test_statistics[[1]],0)
  expect_equal(permutation_test(list(circle,circle,circle),list(sphere,sphere,sphere),iterations = 1,dims = c(2),num_workers = 2)$test_statistics[[1]],0)
  expect_lte(abs(permutation_test(list(circle,circle,circle),list(sphere,circle,sphere),iterations = 1,dims = c(1),num_workers = 2)$test_statistics[[1]] - d^2/3),0.03)
  expect_lte(abs(permutation_test(list(circle,circle,circle),list(sphere,circle,circle),iterations = 1,dims = c(1),num_workers = 2)$test_statistics[[1]] - d^2/3),0.03)
  expect_lte(abs(permutation_test(list(circle,sphere,circle),list(sphere,circle,sphere),iterations = 1,dims = c(1),num_workers = 2)$test_statistics[[1]] - 2*d^2/3),0.03)
  expect_length(unique(permutation_test(list(circle,sphere,circle),list(circle,sphere,circle),paired = T,iterations = 3,dims = c(1),num_workers = 2)$permvals[[1]]),1)
  expect_length(unique(permutation_test(list(sphere,sphere,circle),list(sphere,sphere,circle),paired = T,iterations = 3,dims = c(1),num_workers = 2)$permvals[[1]]),1)

})

test_that("independence_test detects incorrect parameters correctly",{
  
  g1 <- lapply(X = 1:6,FUN = function(X){return(TDA::ripsDiag(X = TDA::circleUnif(n = 50,r = 1),maxdimension = 1,maxscale = 2))})
  g2 <- lapply(X = 1:6,FUN = function(X){return(TDA::ripsDiag(X = TDA::sphereUnif(n = 50,d = 2,r = 1),maxdimension = 1,maxscale = 2))})
  expect_error(independence_test(g1,g2,dims = c(0,1),sigma = 1,t = NA,num_workers = 2),"NA")
  expect_error(independence_test(g1,g2,dims = c(0,1),sigma = NaN,t = 1,num_workers = 2),"NaN")
  expect_error(independence_test(g1,g2,dims = c(0,1),sigma = c(1,2),t = 1,num_workers = 2),"single")
  expect_error(independence_test(g1,g2,dims = c(0,1.1),sigma = 1,t = 1,num_workers = 2),"whole")
  expect_error(independence_test(g1,g2[1:5],dims = c(0,1),sigma = 1,t = 1,num_workers = 2),"same length")
  expect_error(independence_test(g1[1:5],g2[1:5],dims = c(0,1),sigma = 1,t = 1,num_workers = 2),"6")
  expect_error(independence_test(list(g1[[1]],g1[[2]],g1[[3]],g1[[4]],g1[[5]],data.frame(dimension = numeric(),birth = numeric(),death = numeric())),g2,dims = c(0,1),sigma = 1,t = 1,num_workers = 2),"empty")
  expect_error(independence_test(g1,g2,dims = c(0,1),num_workers = 2.3),"num_workers")
  expect_error(independence_test(g1,g2,dims = c(0,1),num_workers = -2),"num_workers")

})

test_that("independence_test can accept inputs from TDA, TDAstats and diagram_to_df",{
  
  circle = TDA::ripsDiag(X = TDA::circleUnif(n = 20,r = 1),maxdimension = 1,maxscale = 2)
  sphere = TDAstats::calculate_homology(TDA::sphereUnif(n = 20,d = 2,r = 1),threshold = 2)
  expect_length(permutation_test(list(circle,circle,diagram_to_df(circle),circle,circle,circle),list(sphere,sphere,sphere,sphere,sphere,circle),iterations = 10,dims = c(1),num_workers = 2)$permvals[[1]],10)
  
})

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