
test_that("diagram_mds detects incorrect parameters correctly",{
  
  D <- data.frame(dimension = c(0),birth = c(0),death = c(1))
  expect_error(diagram_mds(diagrams = list(D,D,"D"),num_workers = 2),"Diagrams")
  expect_error(diagram_mds(diagrams = list(),num_workers = 2),"2")
  expect_error(diagram_mds(diagrams = list(D,D,D),distance = NaN,num_workers = 2),"distance")
  expect_error(diagram_mds(diagrams = list(D,D,D),distance = "fisher",sigma = NULL,num_workers = 2),"sigma")
  expect_error(diagram_mds(diagrams = list(D,D,D),p = NaN,num_workers = 2),"p")
  
})

test_that("diagram_mds is computing correctly",{
  
  D1 <- data.frame(dimension = 0,birth = 2,death = 3)
  D2 <- data.frame(dimension = 0,birth = 2,death = 3.1)
  D3 <- data.frame(dimension = 0,birth = c(2,5),death = c(3.1,6))
  d12 <- diagram_distance(D1,D2,dim = 0) # 2-wasserstein
  d13 <- diagram_distance(D1,D3,dim = 0)
  d23 <- diagram_distance(D2,D3,dim = 0)
  D <- matrix(data = c(0,d12,d13,d12,0,d23,d13,d23,0),byrow = T,nrow = 3,ncol = 3)^2
  D <- scale(D,center = T,scale = F)
  D <- t(scale(t(D),center = T,scale = F))
  S <- -D/2
  ev <- eigen(S)
  embedding <- -1*t(diag(sqrt(ev$values[1:2])) %*% t(ev$vectors[,1:2]))
  dimnames(embedding) <- list(NULL,NULL)
  dmds <- diagram_mds(diagrams = list(D1,D2,D3),num_workers = 2)
  if(embedding[1,1] < 0)
  {
    embedding[,1] <- embedding[,1]/-1
  }
  if(dmds[1,1] < 0)
  {
    dmds[,1] <- dmds[,1]/-1
  }
  if(embedding[1,2] < 0)
  {
    embedding[,2] <- embedding[,2]/-1
  }
  if(dmds[1,2] < 0)
  {
    dmds[,2] <- dmds[,2]/-1
  }
  expect_equal((abs(dmds[1,1])-abs(embedding[1,1]))+(abs(dmds[2,1])-abs(embedding[2,1]))+(abs(dmds[3,1])-abs(embedding[3,1])) + (abs(dmds[1,2])-abs(embedding[1,2]))+(abs(dmds[2,2])-abs(embedding[2,2]))+(abs(dmds[3,2])-abs(embedding[3,2])),0)
  
})

test_that("diagram_mds can accept inputs from TDA, TDAstats and diagram_to_df",{
  
  D1 = TDA::ripsDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxscale = 1,maxdimension = 1)
  D2 = TDA::alphaComplexDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxdimension = 1)
  D3 = TDA::ripsDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxscale = 1,maxdimension = 1,library = "dionysus",location = T)
  D4 = TDAstats::calculate_homology(data.frame(x = runif(50,0,1),y = runif(50,0,1)),threshold = 1)
  expect_type(diagram_mds(diagrams = list(D1,D2,D3,D4),dim = 1,num_workers = 2),"double")
  expect_error(diagram_mds(diagrams = list(D1,D2,D3,D4),dim = 0,num_workers = 2),"Inf")
  
})
