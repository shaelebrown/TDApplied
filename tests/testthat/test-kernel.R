
test_that("diagram_kernel detects incorrect parameters correctly",{
  
  D <- data.frame(dimension = c(0),birth = c(0),death = c(1))
  expect_error(diagram_kernel(D1 = NULL,D2 = D,dim = 1),"TDA/TDAstats")
  expect_error(diagram_kernel(D1 = D,D2 = NULL,dim = 1),"TDA/TDAstats")
  expect_error(diagram_kernel(D1 = D,D2 = D,dim = "2"),"numeric")
  expect_error(diagram_kernel(D1 = D,D2 = D,sigma = "2"),"numeric")
  expect_error(diagram_kernel(D1 = D,D2 = D,t = NA),"NA")

})

test_that("diagram_kernel can accept inputs from either TDA/TDAstats homology output or diagram_to_df function, with or without cycle location",{
  
  D1 = TDA::ripsDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxscale = 1,maxdimension = 1)
  D2 = TDA::alphaComplexDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxdimension = 1)
  D3 = TDA::ripsDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxscale = 1,maxdimension = 1,library = "dionysus",location = T)
  D4 = TDAstats::calculate_homology(data.frame(x = runif(50,0,1),y = runif(50,0,1)),threshold = 1)
  expect_gte(diagram_kernel(D1 = D1,D2 = D2,dim = 1),0)
  expect_gte(diagram_kernel(D1 = diagram_to_df(D1),D2 = D2,dim = 1),0)
  expect_gte(diagram_kernel(D1 = D1,D2 = diagram_to_df(D2),dim = 1),0)
  expect_gte(diagram_kernel(D1 = D3,D2 = diagram_to_df(D2),dim = 1),0)
  expect_gte(diagram_kernel(D1 = D1,D2 = diagram_to_df(D3),dim = 1),0)
  expect_gte(diagram_kernel(D1 = D1,D2 = D4,dim = 1),0)
  expect_error(diagram_kernel(D1 = D1,D2 = D2,dim = 0),"Inf")
  
})

test_that("diagram_kernel is computing correctly",{
  
  D1 <- data.frame(dimension = 0,birth = 2,death = 3)
  D2 <- data.frame(dimension = 0,birth = 2,death = 3.1)
  D3 <- data.frame(dimension = 0,birth = c(2,5),death = c(3.1,6))
  sqrt_rho_1 <- function(sigma)
  {
    v <- (1/(2*pi*sigma^2))*c(exp(0)+exp(-(0.45^2+0.55^2)/(2*sigma^2)),exp(-(0.1^2)/(2*sigma^2))+exp(-(2*0.55^2)/(2*sigma^2)),exp(-(2*0.5^2)/(2*sigma^2)) + exp(-(2*0.05^2)/(2*sigma^2)),exp(-(0.45^2+0.55^2)/(2*sigma^2)) + exp(0))
    v <- v/sum(v)
    return(sqrt(v))
  }
  sqrt_rho_2 <- function(sigma)
  {
    v <- (1/(2*pi*sigma^2))*c(exp(-(0.1^2)/(2*sigma^2))+exp(-(2*0.5^2)/(2*sigma^2)),exp(0)+exp(-(0.5^2+0.6^2)/(2*sigma^2)),exp(-(0.5^2+0.6^2)/(2*sigma^2)) + exp(0),exp(-(2*0.55^2)/(2*sigma^2)) + exp(-(2*0.05^2)/(2*sigma^2)))
    v <- v/sum(v)
    return(sqrt(v))
  }
  v11 <- sqrt_rho_1(1)
  v21 <- sqrt_rho_2(1)
  v12 <- sqrt_rho_1(2)
  v22 <- sqrt_rho_2(2)
  norm_11 <- as.numeric(v11 %*% v21)
  norm_22 <- as.numeric(v12 %*% v22)
  if(norm_11 > 1)
  {
    norm_11 <- 1
  }
  if(norm_11 < -1)
  {
    norm_11 <- -1
  }
  if(norm_22 > 1)
  {
    norm_22 <- 1
  }
  if(norm_22 < -1)
  {
    norm_22 <- -1
  }
  val_1 <- acos(norm_11)
  val_2 <- acos(norm_22)
  expect_equal(diagram_kernel(D1,D2,dim = 0,sigma = 1,t = 1),exp(-1*val_1))
  expect_equal(diagram_kernel(D2,D1,dim = 0,sigma = 1,t = 1),exp(-1*val_1))
  expect_equal(diagram_kernel(D1,D2,dim = 0,sigma = 2,t = 1),exp(-1*val_2))
  expect_equal(diagram_kernel(D1 = D1,D2 = D2,dim = 0,sigma = 1,t = 2),exp(-2*val_1))
  expect_equal(diagram_kernel(D1 = D1,D2 = D2,sigma = 2,t = 2),exp(-2*val_2))
  expect_equal(diagram_kernel(D1 = D2,D2 = D1,sigma = 2,t = 2),exp(-2*val_2))
  expect_identical(diagram_kernel(D1,D1,sigma = 1,t = 1),1)
  
})

test_that("gram_matrix detect incorrect parameters correctly",{
  
  D1 <- data.frame(dimension = 0,birth = 2,death = 3)
  D2 <- data.frame(dimension = 0,birth = 2,death = 3.1)
  D3 <- data.frame(dimension = 0,birth = c(2,5),death = c(3.1,6))
  expect_error(gram_matrix(diagrams = list(D1,D2,D3),num_workers = NaN),"NaN")
  expect_error(gram_matrix(diagrams = list(D1,D2,D3),num_workers = "1"),"numeric")
  
})

test_that("gram_matrix is computing correctly",{
  
  D1 <- data.frame(dimension = 0,birth = 2,death = 3)
  D2 <- data.frame(dimension = 0,birth = 2,death = 3.1)
  D3 <- data.frame(dimension = 0,birth = c(2,5),death = c(3.1,6))
  m1 <- matrix(data = c(1,diagram_kernel(D1,D2,dim = 0,sigma = 1,t = 1),diagram_kernel(D1,D2,dim = 0,sigma = 1,t = 1),1),byrow = T,nrow = 2,ncol = 2)
  class(m1) <- "kernelMatrix"
  m2 <- matrix(data = c(1,diagram_kernel(D1,D2,dim = 0,sigma = 1,t = 1),diagram_kernel(D1,D3,dim = 0,sigma = 1,t = 1),diagram_kernel(D2,D1,dim = 0,sigma = 1,t = 1),1,diagram_kernel(D2,D3,dim = 0,sigma = 1,t = 1),diagram_kernel(D3,D1,dim = 0,sigma = 1,t = 1),diagram_kernel(D3,D2,dim = 0,sigma = 1,t = 1),1),byrow = T,nrow = 3,ncol = 3)
  class(m2) <- "kernelMatrix"
  m3 <- matrix(data = c(1,diagram_kernel(D1,D3,dim = 0,sigma = 1,t = 1),diagram_kernel(D1,D2,dim = 0,sigma = 1,t = 1),diagram_kernel(D2,D3,dim = 0,sigma = 1,t = 1)),byrow = T,nrow = 2,ncol = 2)
  class(m3) <- "kernelMatrix"
  colnames(m3) <- c("result.1","result.2")
  expect_identical(gram_matrix(diagrams = list(D1,D2),dim = 0,sigma = 1,t = 1,num_workers = 2),m1)
  expect_equal(gram_matrix(diagrams = list(D1,D2,D3),dim = 0,sigma = 1,t = 1,num_workers = 2),m2)
  expect_equal(gram_matrix(diagrams = list(D1,D2),other_diagrams = list(D1,D3),dim = 0,sigma = 1,t = 1,num_workers = 2),m3)
  
})


