
test_that("diagram_MDS detects incorrect parameters correctly",{
  
  D <- data.frame(dimension = c(0),birth = c(0),death = c(1))
  expect_error(diagram_MDS(diagrams = list(D,D,"D"),num_workers = 2),"diagram")
  expect_error(diagram_MDS(diagrams = list(),num_workers = 2),"1")
  expect_error(diagram_MDS(diagrams = list(D,D,D),distance = NaN,num_workers = 2),"distance")
  expect_error(diagram_MDS(diagrams = list(D,D,D),distance = "fisher",sigma = NULL,num_workers = 2),"sigma")
  expect_error(diagram_MDS(diagrams = list(D,D,D),p = NaN,num_workers = 2),"p")
  
})

test_that("diagram_MDS is computing correctly",{
  
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
  expect_equal(diagram_MDS(diagrams = list(D1,D2,D3),num_workers = 2),embedding)
  
})



