
test_that("diagram_distance detects incorrect parameters correctly",{
  
  D = data.frame(dimension = c(0),birth = c(0),death = c(1))
  expect_error(diagram_distance(D1 = NULL,D2 = D,dim = 1),"TDA/TDAstats")
  expect_error(diagram_distance(D1 = D,D2 = NULL,dim = 1),"TDA/TDAstats")
  expect_error(diagram_distance(D1 = D,D2 = D,dim = "2"),"numeric")
  expect_error(diagram_distance(D1 = D,D2 = D,dim = 1,p = "2"),"numeric")
  expect_error(diagram_distance(D1 = D,D2 = D,dim = 1,distance = "Wasserstein"),"distance must")
  expect_error(diagram_distance(D1 = D,D2 = D,dim = 1,distance = "fisher",sigma = NA),"sigma must")
  
})

test_that("diagram_distance can accept inputs from either TDA/TDAstats homology output or diagram_to_df function, with or without cycle location",{
  
  D1 = TDA::ripsDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxscale = 1,maxdimension = 1)
  D2 = TDA::alphaComplexDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxdimension = 1)
  D3 = TDA::ripsDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxscale = 1,maxdimension = 1,library = "dionysus",location = T)
  D4 = TDAstats::calculate_homology(data.frame(x = runif(50,0,1),y = runif(50,0,1)),threshold = 1)
  D5 = TDAstats::calculate_homology(data.frame(x = runif(50,0,1),y = runif(50,0,1)),threshold = 10,dim = 1)
  expect_gte(diagram_distance(D1 = D1,D2 = D2,dim = 1),0)
  expect_gte(diagram_distance(D1 = diagram_to_df(D1),D2 = D2,dim = 1),0)
  expect_gte(diagram_distance(D1 = D1,D2 = diagram_to_df(D2),dim = 1),0)
  expect_gte(diagram_distance(D1 = D3,D2 = diagram_to_df(D2),dim = 1),0)
  expect_gte(diagram_distance(D1 = D1,D2 = diagram_to_df(D3),dim = 1),0)
  expect_gte(diagram_distance(D1 = D1,D2 = D4,dim = 1),0)
  expect_error(diagram_distance(D1 = D1,D2 = D2,dim = 0),"Inf")
  
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

test_that("distance_matrix detects incorrect parameters correctly",{
  
  D1 <- data.frame(dimension = 0,birth = 2,death = 3)
  D2 <- data.frame(dimension = 0,birth = 2,death = 3.1)
  D3 <- data.frame(dimension = 0,birth = c(2,5),death = c(3.1,6))
  expect_error(distance_matrix(diagrams = list(D1,D2,D3),num_workers = 0),"num_workers")
  expect_error(distance_matrix(diagrams = list(D1,D2,D3),num_workers = "2"),"num_workers")
  expect_error(distance_matrix(diagrams = list(D1,D2,D3),num_workers = 1.1),"whole")
  
})

test_that("distance_matrix is computing correctly",{
  
  D1 <- data.frame(dimension = 0,birth = 2,death = 3)
  D2 <- data.frame(dimension = 0,birth = 2,death = 3.1)
  D3 <- data.frame(dimension = 0,birth = c(2,5),death = c(3.1,6))
  m1 <- matrix(data = c(0,diagram_distance(D1,D2,dim = 0,p = 2,distance = "wasserstein"),diagram_distance(D1,D2,dim = 0,p = 2,distance = "wasserstein"),0),byrow = T,nrow = 2,ncol = 2)
  m2 <- matrix(data = c(0,diagram_distance(D1,D2,dim = 0,p = 3,distance = "wasserstein"),diagram_distance(D1,D3,dim = 0,p = 3,distance = "wasserstein"),diagram_distance(D1,D2,dim = 0,p = 3,distance = "wasserstein"),0,diagram_distance(D2,D3,dim = 0,p = 3,distance = "wasserstein"),diagram_distance(D1,D3,dim = 0,p = 3,distance = "wasserstein"),diagram_distance(D3,D2,dim = 0,p = 3,distance = "wasserstein"),0),byrow = T,nrow = 3,ncol = 3)
  m3 <- matrix(data = c(0,diagram_distance(D1,D3,dim = 0,distance = "fisher",sigma = 1),diagram_distance(D1,D2,dim = 0,distance = "fisher",sigma = 1),diagram_distance(D3,D2,dim = 0,distance = "fisher",sigma = 1)),byrow = T,nrow = 2,ncol = 2)
  colnames(m3) <- c("result.1","result.2")
  expect_identical(distance_matrix(diagrams = list(D1,D2),dim = 0,distance = "wasserstein",p = 2,num_workers = 2),m1)
  expect_equal(distance_matrix(diagrams = list(D1,D2,D3),dim = 0,distance = "wasserstein",p = 3,num_workers = 2),m2)
  expect_equal(distance_matrix(diagrams = list(D1,D2),other_diagrams = list(D1,D3),dim = 0,distance = "fisher",sigma = 1,num_workers = 2),m3)
  
})
