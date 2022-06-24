
test_that("diagram_ksvm detects incorrect parameters correctly",{
  
  D1 <- data.frame(dimension = 0,birth = 2,death = 3)
  D2 <- data.frame(dimension = 0,birth = 2,death = 3.1)
  D3 <- data.frame(dimension = 0,birth = c(2,5),death = c(3.1,6))
  expect_error(diagram_ksvm(diagrams = list(D1,D2,NULL),y = c(0,1,2)),"diagram")
  expect_error(diagram_ksvm(diagrams = list(D1,D2,D3),cv = NA,y = c(0,1,2)),"cv")
  expect_error(diagram_ksvm(diagrams = list(D1,D2,D3),cv = 0,y = c(0,1,2)),"cv")
  expect_error(diagram_ksvm(diagrams = list(D1,D2,D3),cv = 1.1,y = c(0,1,2)),"cv")
  expect_error(diagram_ksvm(diagrams = list(D1,D2,D3),cv = c(1,2),y = c(0,1,2)),"cv")
  expect_error(diagram_ksvm(diagrams = list(D1,D2,D3),t = -1,y = c(0,1,2)),"t")
  expect_error(diagram_ksvm(diagrams = list(D1,D2,D3),dim = 0,sigma = 0,y = c(0,1,2)),"sigma")
  expect_error(diagram_ksvm(diagrams = list(D1,D2,D3),dim = NaN,y = c(0,1,2)),"dim")
  expect_error(diagram_ksvm(diagrams = list(D1,D2,D3),dim = 1,y = c(0,1)),"number of elements")
  
})

test_that("predict_diagram_ksvm detects incorrect parameters correctly",{
  
  D1 <- data.frame(dimension = 0,birth = 2,death = 3)
  D2 <- data.frame(dimension = 0,birth = 2,death = 3.1)
  D3 <- data.frame(dimension = 0,birth = c(2,5),death = c(3.1,6))
  ksvm <- diagram_ksvm(diagrams = list(D1,D2,D3),dim = 0,y = c(1,2,3))
  expect_error(predict_diagram_ksvm(new_diagrams = list(),ksvm),"1")
  expect_error(predict_diagram_ksvm(new_diagrams = NULL,ksvm),"NULL")
  expect_error(predict_diagram_ksvm(new_diagrams = list(diagrams[[1]],"1"),ksvm),"diagram")
  expect_error(predict_diagram_ksvm(new_diagrams = list(D1,D2,D3),model = list(1,2,3)),"ksvm")
  
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
  ksvm <- diagram_ksvm(diagrams = diagrams,dim = 1,y = as.factor(c(rep("circle",5),rep("torus",5),rep("sphere",5))))
  expect_equal(as.character(predict_diagram_ksvm(new_diagrams = diagrams,ksvm)),c(rep("circle",5),rep("torus",5),rep("sphere",5)))
  
})

