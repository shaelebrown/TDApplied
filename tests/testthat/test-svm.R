
test_that("diagram_ksvm detects incorrect parameters correctly",{
  
  D1 <- data.frame(dimension = 0,birth = 2,death = 3)
  D2 <- data.frame(dimension = 0,birth = 2,death = 3.1)
  D3 <- data.frame(dimension = 0,birth = c(2,5),death = c(3.1,6))
  expect_error(diagram_ksvm(diagrams = list(D1,D2,NULL),y = c(0,1,2),num_workers = 2),"Diagrams")
  expect_error(diagram_ksvm(diagrams = list(D1,D2,D3),cv = NA,y = c(0,1,2),num_workers = 2),"cv")
  expect_error(diagram_ksvm(diagrams = list(D1,D2,D3),cv = 0,y = c(0,1,2),num_workers = 2),"cv")
  expect_error(diagram_ksvm(diagrams = list(D1,D2,D3),cv = 1.1,y = c(0,1,2),num_workers = 2),"cv")
  expect_error(diagram_ksvm(diagrams = list(D1,D2,D3),cv = c(1,2),y = c(0,1,2),num_workers = 2),"cv")
  expect_error(diagram_ksvm(diagrams = list(D1,D2,D3),t = -1,y = c(0,1,2),num_workers = 2),"t")
  expect_error(diagram_ksvm(diagrams = list(D1,D2,D3),dim = 0,sigma = 0,y = c(0,1,2),num_workers = 2),"sigma")
  expect_error(diagram_ksvm(diagrams = list(D1,D2,D3),dim = NaN,y = c(0,1,2),num_workers = 2),"dim")
  expect_error(diagram_ksvm(diagrams = list(D1,D2,D3),dim = 1,y = c(0,1),num_workers = 2),"number of elements")
  expect_error(diagram_ksvm(diagrams = list(D1,D2,D3),dim = 1,y = c("0","1","2"),num_workers = 2),"factor")
  
})

test_that("diagram_ksvm can accept inputs from TDA, TDAstats and diagram_to_df",{
  
  skip_on_cran()
  skip_if_not_installed("TDA")
  skip_if_not_installed("TDAstats")
  
  D1 = TDA::ripsDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxscale = 1,maxdimension = 1)
  D2 = TDA::alphaComplexDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxdimension = 1)
  D3 = TDA::ripsDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxscale = 1,maxdimension = 1,library = "dionysus",location = T)
  D4 = TDAstats::calculate_homology(data.frame(x = runif(50,0,1),y = runif(50,0,1)),threshold = 1)
  expect_s3_class(diagram_ksvm(diagrams = list(D1,D2,D3,D4),y = c(1,2,3,4),num_workers = 2,dim = c(1)),"diagram_ksvm")
  expect_error(diagram_ksvm(diagrams = list(D1,D2,D3,D4),y = c(1,2,3,4),num_workers = 2,dim = c(0)),"Inf")
  expect_error(diagram_ksvm(diagrams = list(D1,D2,D3,D4),y = c(1,2,3,4),num_workers = 2,cv = 2,dim = c(0)),"Inf")
  
})

test_that("predict_diagram_ksvm detects incorrect parameters correctly",{
  
  D1 <- data.frame(dimension = 0,birth = 2,death = 3)
  D2 <- data.frame(dimension = 0,birth = 2,death = 3.1)
  D3 <- data.frame(dimension = 0,birth = c(2,5),death = c(3.1,6))
  ksvm <- diagram_ksvm(diagrams = list(D1,D2,D3),dim = 0,y = c(1,2,3),t = c(1,2),num_workers = 2)
  expect_error(predict_diagram_ksvm(new_diagrams = list(),ksvm,num_workers = 2),"1")
  expect_error(predict_diagram_ksvm(new_diagrams = NULL,ksvm,num_workers = 2),"NULL")
  expect_error(predict_diagram_ksvm(new_diagrams = list(D1,"1"),ksvm,num_workers = 2),"Diagrams")
  expect_error(predict_diagram_ksvm(new_diagrams = list(D1,D2,D3),model = list(1,2,3),num_workers = 2),"ksvm")
  expect_error(predict_diagram_ksvm(new_diagrams = list(D1,D2,D3),model = NULL,num_workers = 2),"supplied")
  
})

test_that("predict_diagram_ksvm is computing correctly",{
  
  circle <- data.frame(dimension = c(0,1,2),birth = c(0,0,0),death = c(2,2,0))
  torus <- data.frame(dimension = c(0,1,1,2),birth = c(0,0,0,0),death = c(2,0.5,1.5,0.5))
  sphere <- data.frame(dimension = c(0,1,2),birth = c(0,0,0),death = c(2,0,2))
  
  circles <- lapply(X = 1:5,FUN = function(X){
    
    t <- circle
    t$death <- t$death + rnorm(nrow(t),mean = 0,sd = 0.01)
    t[which(t$death < 0),3L] <- 0.001
    return(t)
    
  })
  tori <- lapply(X = 1:5,FUN = function(X){
    
    t <- torus
    t$death <- t$death + rnorm(nrow(t),mean = 0,sd = 0.01)
    t[which(t$death < 0),3L] <- 0.001
    return(t)
    
  })
  spheres <- lapply(X = 1:5,FUN = function(X){
    
    t <- sphere
    t$death <- t$death + rnorm(nrow(t),mean = 0,sd = 0.01)
    t[which(t$death < 0),3L] <- 0.001
    return(t)
    
  })
  diagrams <- list(circles[[1]],circles[[2]],circles[[3]],circles[[4]],circles[[5]],
                   tori[[1]],tori[[2]],tori[[3]],tori[[4]],tori[[5]],
                   spheres[[1]],spheres[[2]],spheres[[3]],spheres[[4]],spheres[[5]])
  ksvm <- diagram_ksvm(diagrams = diagrams,dim = 1,y = as.factor(c(rep("circle",5),rep("torus",5),rep("sphere",5))),num_workers = 2)
  expect_equal(as.character(predict_diagram_ksvm(new_diagrams = diagrams,ksvm,num_workers = 2)),c(rep("circle",5),rep("torus",5),rep("sphere",5)))
  
})

test_that("diagram_ksvm can accept inputs from TDA, TDAstats and diagram_to_df",{
  
  skip_on_cran()
  skip_if_not_installed("TDA")
  skip_if_not_installed("TDAstats")
  
  D1 = TDA::ripsDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxscale = 1,maxdimension = 1)
  D2 = TDA::alphaComplexDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxdimension = 1)
  D3 = TDA::ripsDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxscale = 1,maxdimension = 1,library = "dionysus",location = T)
  D4 = TDAstats::calculate_homology(data.frame(x = runif(50,0,1),y = runif(50,0,1)),threshold = 1)
  ksvm = diagram_ksvm(diagrams = list(D1,D2,D3,D4),y = c(1,2,3,4),num_workers = 2,dim = c(1))
  expect_length(predict_diagram_ksvm(new_diagrams = list(D1,D2,D3,D4),model = ksvm,num_workers = 2),4)
  
})

