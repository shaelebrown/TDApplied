
test_that("diagram_kpca detects incorrect parameters correctly",{
  
  D <- data.frame(dimension = c(0),birth = c(0),death = c(1))
  expect_error(diagram_kpca(diagrams = list(D,D,D[,1:2]),num_workers = 2),"three")
  expect_error(diagram_kpca(diagrams = list(D,D,D),t = -1,num_workers = 2),"t")
  expect_error(diagram_kpca(diagrams = list(D,D,D),sigma = 0,num_workers = 2),"sigma")
  expect_error(diagram_kpca(diagrams = list(D,D,D),dim = NULL,num_workers = 2),"dim")
  
})

test_that("diagram_kpca is computing correctly",{
  
  D1 <- data.frame(dimension = 0,birth = 2,death = 3)
  D2 <- data.frame(dimension = 0,birth = 2,death = 3.1)
  D3 <- data.frame(dimension = 0,birth = c(2,5),death = c(3.1,6))
  k12 <- diagram_kernel(D1,D2)
  k13 <- diagram_kernel(D1,D3)
  k23 <- diagram_kernel(D2,D3)
  K <- matrix(data = c(1,k12,k13,k12,1,k23,k13,k23,1),nrow = 3,ncol = 3,byrow = T)
  K <- scale(K,center = T,scale = F)
  K <- t(scale(t(K),center = T,scale = F))
  eig <- eigen(K)
  kpca <- diagram_kpca(diagrams = list(D1,D2,D3),features = 2,num_workers = 2)
  expect_equal(as.numeric(kpca$pca@pcv[,1]),(kpca$pca@pcv[1,1]/eig$vectors[1,1])*as.numeric(eig$vectors[,1]))
  expect_equal(as.numeric(kpca$pca@pcv[,2]),(kpca$pca@pcv[1,2]/eig$vectors[1,2])*as.numeric(eig$vectors[,2]))
  expect_equal(as.numeric(kpca$pca@eig)/sum(as.numeric(kpca$pca@eig)),eig$values[1:2]/sum(eig$values[1:2]))
  
})

test_that("diagram_kpca can accept inputs from TDA, TDAstats and diagram_to_df",{
  
  D1 = TDA::ripsDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxscale = 1,maxdimension = 1)
  D2 = TDA::alphaComplexDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxdimension = 1)
  D3 = TDA::ripsDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxscale = 1,maxdimension = 1,library = "dionysus",location = T)
  D4 = TDAstats::calculate_homology(data.frame(x = runif(50,0,1),y = runif(50,0,1)),threshold = 1,dim = 1)
  expect_error(diagram_kpca(diagrams = list(D1,D2,D3,D4),dim = 1,features = 2,num_workers = 2),"embedding")
  expect_error(diagram_kpca(diagrams = list(D1,D2,D3,D4),dim = 0,features = 2,num_workers = 2),"Inf")
  D1 <- data.frame(dimension = 0,birth = 2,death = 3)
  D2 <- data.frame(dimension = 0,birth = 2,death = 3.1)
  D3 <- data.frame(dimension = 0,birth = c(2,5),death = c(3.1,6))
  expect_s3_class(diagram_kpca(diagrams = list(D1,D2,D3),num_workers = 2),"diagram_kpca")
  
})

test_that("predict_diagram_kpca detects incorrect parameters correctly",{
  
  D1 <- data.frame(dimension = 0,birth = 2,death = 3)
  D2 <- data.frame(dimension = 0,birth = 2,death = 3.1)
  D3 <- data.frame(dimension = 0,birth = c(2,5),death = c(3.1,6))
  diagrams = list(D1,D2,D3)
  kpca <- diagram_kpca(diagrams = list(D1,D2,D3),features = 2,num_workers = 2)
  expect_error(predict_diagram_kpca(new_diagrams = list(),kpca,num_workers = 2),"1")
  expect_error(predict_diagram_kpca(new_diagrams = NA,kpca,num_workers = 2),"NA")
  expect_error(predict_diagram_kpca(new_diagrams = list(diagrams[[1]],1),kpca,num_workers = 2),"TDA/TDAstats")
  expect_error(predict_diagram_kpca(new_diagrams = list(D1,D2,D3),embedding = 2,num_workers = 2),"kpca")
  expect_error(predict_diagram_kpca(new_diagrams = list(D1,D2,D3),embedding = NULL,num_workers = 2),"supplied")
  
})

test_that("predict_diagram_kpca is computing correctly",{
  
  D1 <- data.frame(dimension = 0,birth = 2,death = 3)
  D2 <- data.frame(dimension = 0,birth = 2,death = 3.1)
  D3 <- data.frame(dimension = 0,birth = c(2,5),death = c(3.1,6))
  kpca <- diagram_kpca(diagrams = list(D1,D2,D3),features = 2,num_workers = 2)
  expect_equal(predict_diagram_kpca(new_diagrams = list(D1,D2,D3),embedding = kpca,num_workers = 2),kpca$pca@rotated)
  
})

test_that("predict_diagram_kpca can accept inputs from TDA, TDAstats and diagram_to_df",{
  
  skip_on_cran()
  
  D1 <- data.frame(dimension = 1,birth = 2,death = 3)
  D2 <- data.frame(dimension = 1,birth = 2,death = 3.1)
  D3 <- data.frame(dimension = 1,birth = c(2,5),death = c(3.1,6))
  kpca <- diagram_kpca(diagrams = list(D1,D2,D3),features = 2,num_workers = 2,dim = 1)
  D1 = TDA::ripsDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxscale = 1,maxdimension = 1)
  D2 = TDA::alphaComplexDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxdimension = 1)
  D3 = TDA::ripsDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxscale = 1,maxdimension = 1,library = "dionysus",location = T)
  D4 = TDAstats::calculate_homology(data.frame(x = runif(50,0,1),y = runif(50,0,1)),threshold = 1)
  expect_type(predict_diagram_kpca(new_diagrams = list(D1,D2,D3,D4),embedding = kpca,num_workers = 2),"double")
  D1 <- data.frame(dimension = 0,birth = 2,death = 3)
  D2 <- data.frame(dimension = 0,birth = 2,death = 3.1)
  D3 <- data.frame(dimension = 0,birth = c(2,5),death = c(3.1,6))
  kpca <- diagram_kpca(diagrams = list(D1,D2,D3),features = 2,num_workers = 2,dim = 0)
  D1 = TDA::ripsDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxscale = 1,maxdimension = 1)
  D2 = TDA::alphaComplexDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxdimension = 1)
  D3 = TDA::ripsDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxscale = 1,maxdimension = 1,library = "dionysus",location = T)
  D4 = TDAstats::calculate_homology(data.frame(x = runif(50,0,1),y = runif(50,0,1)),threshold = 1)
  expect_error(predict_diagram_kpca(new_diagrams = list(D1,D2,D3,D4),embedding = kpca,num_workers = 2),"Inf")
  
})

