
test_that("permutation_test detects incorrect parameters correctly",{
  
  circle = TDA::ripsDiag(X = TDA::circleUnif(n = 100,r = 1),maxdimension = 2,maxscale = 2)
  sphere = TDA::ripsDiag(X = TDA::sphereUnif(n = 100,d = 2,r = 1),maxdimension = 2,maxscale = 2)
  expect_error(permutation_test(list(circle,2,circle),list(sphere,sphere,sphere),iterations = 5,p = 2,q = 2,dims = c(0,1),paired = F,distance = "wasserstein",sigma = NULL,verbose = F),"Every diagram must")
  expect_error(permutation_test(list(circle,circle,circle),list(sphere,"2",sphere),iterations = 5,p = 2,q = 2,dims = c(0,1),paired = F,distance = "wasserstein",sigma = NULL,verbose = F),"Every diagram must")
  expect_error(permutation_test(list(circle,circle,circle),list(sphere,sphere,sphere),iterations = NA,p = 2,q = 2,dims = c(0,1),paired = F,distance = "wasserstein",sigma = NULL,verbose = F),"numeric")
  expect_error(permutation_test(list(circle,circle,circle),list(sphere,sphere,sphere),iterations = 5,p = NULL,q = 2,dims = c(0,1),paired = F,distance = "wasserstein",sigma = NULL,verbose = F),"NULL")
  expect_error(permutation_test(list(circle,circle,circle),list(sphere,sphere,sphere),iterations = 5,p = 2,q = NA,dims = c(0,1),paired = F,distance = "wasserstein",sigma = NULL,verbose = F),"numeric")
  expect_error(permutation_test(list(circle,circle,circle),list(sphere,sphere,sphere),iterations = 5,p = 2,q = 2,dims = c(-1,1),paired = F,distance = "wasserstein",sigma = NULL,verbose = F),"positive")
  expect_error(permutation_test(list(circle,circle),list(sphere,sphere,sphere),iterations = 5,p = 2,q = 2,dims = c(0,1),paired = T,distance = "wasserstein",sigma = NULL,verbose = F),"same number")
  expect_error(permutation_test(list(circle,circle,circle),list(sphere,sphere,sphere),iterations = 5,p = 2,q = 2,dims = c(0,1),paired = F,distance = "Wasserstein",sigma = NULL,verbose = F),"wasserstein")
  expect_error(permutation_test(list(circle,circle,circle),list(sphere,sphere,sphere),iterations = 5,p = 2,q = 2,dims = c(0,1),paired = F,distance = "fisher",sigma = NULL,verbose = F),"sigma")
  expect_error(permutation_test(list(circle,circle,circle),list(sphere,sphere),iterations = 5,p = 2,q = 2,dims = c(0,1),paired = T,distance = "fisher",sigma = 1,verbose = F),"paired")
  
})

test_that("permutation_test is working",{
  
  circle <- TDA::ripsDiag(X = TDA::circleUnif(n = 100,r = 1),maxdimension = 2,maxscale = 2)
  sphere <- TDA::ripsDiag(X = TDA::sphereUnif(n = 100,d = 2,r = 1),maxdimension = 2,maxscale = 2)
  d <- diagram_distance(circle,sphere,dim = 1) # wasserstein distance in dimension 1
  expect_length(permutation_test(list(circle,circle,circle),list(sphere,sphere,sphere),iterations = 10,dims = c(1))$permvals[[1]],10)
  expect_length(permutation_test(list(circle,circle,circle),list(sphere,sphere,sphere),iterations = 15,dims = c(1))$permvals[[1]],15)
  expect_equal(permutation_test(list(circle,circle,circle),list(sphere,sphere,sphere),iterations = 10,dims = c(1))$test_statistics[[1]],0)
  expect_equal(permutation_test(list(circle,circle,circle),list(sphere,sphere,sphere),iterations = 10,dims = c(2))$test_statistics[[1]],0)
  expect_equal(permutation_test(list(circle,circle,circle),list(sphere,circle,sphere),iterations = 10,dims = c(1))$test_statistics[[1]],d^2/3)
  expect_equal(permutation_test(list(circle,circle,circle),list(sphere,circle,circle),iterations = 10,dims = c(1))$test_statistics[[1]],d^2/3)
  expect_equal(permutation_test(list(circle,sphere,circle),list(sphere,circle,sphere),iterations = 10,dims = c(1))$test_statistics[[1]],2*d^2/3)
  expect_length(unique(permutation_test(list(circle,sphere,circle),list(circle,sphere,circle),paired = T,iterations = 10,dims = c(1))$permvals[[1]]),1)
  expect_length(unique(permutation_test(list(sphere,sphere,circle),list(sphere,sphere,circle),paired = T,iterations = 10,dims = c(1))$permvals[[1]]),1)
  
})

test_that("independence_test detects incorrect parameters correctly",{
  
  # (...,iterations = 100,p = 2,q = 2,dims = c(0,1),paired = F,distance = "wasserstein",sigma = NULL,verbose = FALSE)

})

test_that("independence_test is working",{
  
  # do, check for several types of datasets checking the test statistic, permutation values and p value
  
})