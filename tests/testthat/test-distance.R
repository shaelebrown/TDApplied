
test_that("diagram_distance detects incorrect parameters correctly",{
  
  D = data.frame(dimension = c(0),birth = c(0),death = c(1))
  expect_error(diagram_distance(D1 = NULL,D2 = D,dim = 1),"TDA computation or data frame")
  expect_error(diagram_distance(D1 = D,D2 = NULL,dim = 1),"TDA computation or data frame")
  expect_error(diagram_distance(D1 = D,D2 = D,dim = "2"),"whole number")
  expect_error(diagram_distance(D1 = D,D2 = D,dim = 1,p = "2"),"number")
  expect_error(diagram_distance(D1 = D,D2 = D,dim = 1,distance = "Wasserstein"),"distance must")
  expect_error(diagram_distance(D1 = D,D2 = D,dim = 1,distance = "fisher",sigma = NA),"sigma must")
  
})

test_that("diagram_distance can accept inputs from either TDA homology output or diagram_to_df function",{
  
  D1 = TDA::ripsDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxscale = 1,maxdimension = 1)
  D2 = TDA::alphaComplexDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxdimension = 1)
  expect_gte(diagram_distance(D1 = D1,D2 = D2,dim = 1),0)
  expect_gte(diagram_distance(D1 = diagram_to_df(D1),D2 = D2,dim = 1),0)
  expect_gte(diagram_distance(D1 = D1,D2 = diagram_to_df(D2),dim = 1),0)
  
})

test_that("diagram_distance is computing correctly",{
  
  D1 = data.frame(dimension = 0,birth = 2,death = 3)
  D2 = data.frame(dimension = 0,birth = c(2,5),death = c(3.1,6))
  expect_identical(diagram_distance(D1,D2,dim = 0,distance = "wasserstein",p = 2),sqrt(0.1^2+0.5^2))
  expect_identical(diagram_distance(D2,D1,dim = 0,distance = "wasserstein",p = 2),sqrt(0.1^2+0.5^2))
  expect_identical(diagram_distance(D1,D2,dim = 0,distance = "wasserstein",p = 3),(0.1^3+0.5^3)^(1/3))
  expect_equal(diagram_distance(D1 = D1,D2 = D2,distance = "fisher",dim = 0,sigma = 1),diagram_distance(D1 = D2,D2 = D1,distance = "fisher",dim = 0,sigma = 1))
  expect_identical(diagram_distance(D1 = D1,D2 = D2,p = Inf,distance = "wasserstein",dim = 0),0.5)
  expect_identical(diagram_distance(D1 = D2,D2 = D1,p = Inf,distance = "wasserstein",dim = 0),0.5)
  expect_identical(diagram_distance(D1 = D1,D2 = D1,p = Inf,distance = "wasserstein",dim = 0),0)
  expect_identical(diagram_distance(D1 = D1,D2 = D1,p = 2,distance = "wasserstein",dim = 0),0)
  expect_identical(diagram_distance(D1 = D1,D2 = D1,distance = "fisher",sigma = 1,dim = 0),0)
  
})