
test_that("diagram_ksvm detects incorrect parameters correctly",{
  
  D1 <- data.frame(dimension = 0,birth = 2,death = 3)
  D2 <- data.frame(dimension = 0,birth = 2,death = 3.1)
  D3 <- data.frame(dimension = 0,birth = c(2,5),death = c(3.1,6))
  expect_error(diagram_ksvm(diagrams = list(D1,D2,NULL),y = c(0,1,2),num_workers = 2),"Diagrams")
  expect_error(diagram_ksvm(diagrams = list(D1,D2,D3),cv = NA,y = c(0,1,2),num_workers = 2),"cv")
  expect_error(diagram_ksvm(diagrams = list(D1,D2,D3),cv = 0,y = c(0,1,2),num_workers = 2),"cv")
  expect_error(diagram_ksvm(diagrams = list(D1,D2,D3),dim = 1,cv = 2,y = c(0,1,2),num_workers = 2),"cv")
  expect_error(diagram_ksvm(diagrams = list(D1,D2,D3),cv = 1.1,y = c(0,1,2),num_workers = 2,dim = 1),"cv")
  expect_error(diagram_ksvm(diagrams = list(D1,D2,D3),cv = c(1,2),y = c(0,1,2),num_workers = 2),"cv")
  expect_error(diagram_ksvm(diagrams = list(D1,D2,D3),t = -1,y = c(0,1,2),num_workers = 2),"t")
  expect_error(diagram_ksvm(diagrams = list(D1,D2,D3),dim = 0,sigma = 0,y = c(0,1,2),num_workers = 2),"sigma")
  expect_error(diagram_ksvm(diagrams = list(D1,D2,D3),dim = NaN,y = c(0,1,2),num_workers = 2),"dim")
  expect_error(diagram_ksvm(diagrams = list(D1,D2,D3),dim = 1,y = c(0,1),num_workers = 2),"number of elements")
  expect_error(diagram_ksvm(diagrams = list(D1,D2,D3),dim = 1,y = c("0","1","2"),num_workers = 2),"factor")
  expect_error(diagram_ksvm(diagrams = list(D1,D2,D3),dim = 1,y = as.factor(c("0","1","2")),cv = 2,num_workers = 2),"One class")
  
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

test_that("diagram_ksvm can accept precomputed distance matrices",{
  
  skip_on_cran()
  skip_if_not_installed("TDA")
  skip_if_not_installed("TDAstats")
  
  D1 = TDA::ripsDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxscale = 1,maxdimension = 1)
  D2 = TDA::ripsDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxscale = 1,maxdimension = 1,library = "dionysus",location = T)
  D3 = TDAstats::calculate_homology(data.frame(x = runif(50,0,1),y = runif(50,0,1)),threshold = 1)
  d0 = distance_matrix(diagrams = list(D1,D2,D3),dim = 0,num_workers = 2,distance = "fisher",sigma = 1)
  d1 = distance_matrix(diagrams = list(D1,D2,D3),dim = 1,num_workers = 2,distance = "fisher",sigma = 1)
  expect_s3_class(diagram_ksvm(diagrams = list(D1,D2,D3),y = c(1,2,3),distance_matrices = list(d0,d1),num_workers = 2,dim = c(0,1)),"diagram_ksvm")
  
  expect_error(diagram_ksvm(diagrams = list(D1,D2,D3),y = c(1,2,3),distance_matrices = list(d0,matrix(data = c(0,1,NA,0),nrow = 2,ncol = 2,byrow = T)),num_workers = 2,dim = c(0,1)),"missing")
  
  expect_error(diagram_ksvm(diagrams = list(D1,D2,D3),y = c(1,2,3),distance_matrices = list(d0,d1,d1),num_workers = 2,dim = c(0,1)),"expand.grid")
  
  expect_error(diagram_ksvm(diagrams = list(D1,D2,D3),y = c(1,2,3),distance_matrices = NA,num_workers = 2,dim = c(0,1)),"list")
  
})

test_that("diagram_ksvm can perform cross validation with any valid model type",{
  
  skip_on_cran()
  skip_if_not_installed("TDA")
  skip_if_not_installed("TDAstats")
  
  # create diags
  g <- lapply(X = 1:10,FUN = function(X){
    
    if(X <= 5)
    {
      return(TDAstats::calculate_homology(TDA::circleUnif(n = 50),threshold = 1,dim = 1))
    }
    return(TDAstats::calculate_homology(TDA::sphereUnif(n = 50,d = 2),threshold = 1,dim = 1))
    
  })
  
  # create models with CV
  expect_type(diagram_ksvm(diagrams = g,cv = 2,dim = 1,y = factor(c(rep("0",5),rep("1",5))),type = "C-svc"),"list")
  expect_type(diagram_ksvm(diagrams = g,cv = 2,dim = 1,y = factor(c(rep("0",5),rep("1",5))),type = "C-svc"),"list")
  expect_type(diagram_ksvm(diagrams = g,cv = 2,dim = 1,y = factor(c(rep("0",3),rep("1",3),rep("2",4))),type = "C-svc"),"list")
  expect_type(diagram_ksvm(diagrams = g,cv = 1,dim = 1,y = factor(c(rep("0",3),rep("1",3),rep("2",4))),type = "C-svc",prob.model = T),"list") # performs cv internally anyways
  expect_type(diagram_ksvm(diagrams = g,cv = 2,dim = 1,y = factor(c(rep("0",5),rep("1",5))),type = "nu-svc"),"list")
  expect_error(diagram_ksvm(diagrams = g,cv = 2,dim = 1,y = factor(c(rep("0",3),rep("1",3),rep("2",4))),type = "nu-svc"),"nu-svc")
  expect_error(diagram_ksvm(diagrams = g,cv = 2,dim = 1,y = factor(c(rep("0",5),rep("1",5))),type = "C-bsvc"),"type")
  expect_error(diagram_ksvm(diagrams = g,cv = 2,dim = 1,y = factor(c(rep("0",3),rep("1",3),rep("2",4))),type = "spoc-svc"),"type")
  expect_error(diagram_ksvm(diagrams = g,cv = 2,dim = 1,y = factor(c(rep("0",3),rep("1",3),rep("2",4))),type = "kbb-svc"),"type")
  expect_type(diagram_ksvm(diagrams = g,cv = 2,dim = 1,type = "one-svc"),"list")
  # expect_type(diagram_ksvm(diagrams = g,cv = 2,dim = 1,y = c(rep(1,5),rep(2,5)),type = "eps-svr"),"list")
  expect_type(diagram_ksvm(diagrams = g,cv = 2,dim = 1,y = c(rep(1,5),rep(2,5)),type = "nu-svr"),"list")
  expect_error(diagram_ksvm(diagrams = g,cv = 2,dim = 1,y = c(rep(1,5),rep(2,5)),type = "eps-bsvr"),"type")
  expect_type(diagram_ksvm(diagrams = g,cv = 2,dim = 1,y = factor(c(rep("0",3),rep("1",3),rep("2",4))),type = "C-svc"),"list")
})

test_that("diagram_ksvm can handle missing t",{
  
  skip_on_cran()
  skip_if_not_installed("TDA")
  skip_if_not_installed("TDAstats")
  
  # create diags
  g <- lapply(X = 1:10,FUN = function(X){
    
    if(X <= 5)
    {
      return(TDAstats::calculate_homology(TDA::circleUnif(n = 50),threshold = 1,dim = 1))
    }
    return(TDAstats::calculate_homology(TDA::sphereUnif(n = 50,d = 2),threshold = 1,dim = 1))
    
  })
  
  expect_error(diagram_ksvm(diagrams = g,cv = 1,dim = 1,y = factor(c(rep("0",5),rep("1",5))),t = NA),"t")
  expect_error(diagram_ksvm(diagrams = g,cv = 1,dim = 1,y = factor(c(rep("0",5),rep("1",5))),distance_matrices = list(matrix(data = 0,nrow = 10,ncol = 10)),t = NULL),"variance")
  expect_error(diagram_ksvm(diagrams = g,cv = 2,dim = 1,y = factor(c(rep("0",5),rep("1",5))),distance_matrices = list(matrix(data = 0,nrow = 10,ncol = 10)),t = NULL),"variance")
  # create models
  expect_type(diagram_ksvm(diagrams = g,cv = 1,dim = 1,y = factor(c(rep("0",5),rep("1",5))),t = NULL),"list")
  expect_type(diagram_ksvm(diagrams = g,cv = 2,dim = 1,y = factor(c(rep("0",5),rep("1",5))),t = NULL),"list")
  
})

test_that("diagram_ksvm can handle zero variance distances matrices",{
  
  skip_on_cran()
  skip_if_not_installed("TDA")
  skip_if_not_installed("TDAstats")
  
  # create diags and distance mats
  g <- lapply(X = 1:10,FUN = function(X){
    
    if(X <= 5)
    {
      return(TDAstats::calculate_homology(TDA::circleUnif(n = 20),threshold = 1,dim = 1))
    }
    return(TDAstats::calculate_homology(TDA::sphereUnif(n = 20,d = 2),threshold = 1,dim = 1))
    
  })
  D0 <- distance_matrix(diagrams = g,dim = 0,distance = "fisher",sigma = 1)
  D1 <- distance_matrix(diagrams = g,dim = 1,distance = "fisher",sigma = 1)
  D2 <- distance_matrix(diagrams = g,dim = 2,distance = "fisher",sigma = 1)
  D3 <- D2
  D3[1,2] <- 1
  D3[2,1] <- 1
  
  expect_error(diagram_ksvm(diagrams = g,cv = 1,dim = c(0,1,2),y = factor(c(rep("0",5),rep("1",5))),t = 1,distance_matrices = list(D2,D2,D2)),"0 variance")
  expect_error(diagram_ksvm(diagrams = g,cv = 2,dim = c(0,1,2),y = factor(c(rep("0",5),rep("1",5))),t = 1,distance_matrices = list(D2,D2,D2)),"one cv fold")
  expect_error(diagram_ksvm(diagrams = g,cv = 2,dim = c(0,1,2),y = factor(c(rep("0",5),rep("1",5))),t = 1,distance_matrices = list(D3,D2,D2)),"one cv fold")
  
  res <- diagram_ksvm(diagrams = g,cv = 2,dim = c(0,1,2),y = factor(c(rep("0",5),rep("1",5))),t = 1,distance_matrices = list(D0,D1,D2))
  expect_true(is.na(res$cv_results[3,7]))
  res <- diagram_ksvm(diagrams = g,cv = 2,dim = c(0,1,2),y = factor(c(rep("0",5),rep("1",5))),distance_matrices = list(D0,D1,D2))
  expect_true(is.na(res$cv_results[3,7]))
  res <- diagram_ksvm(diagrams = g,cv = 1,dim = c(0,1,2),y = factor(c(rep("0",5),rep("1",5))),distance_matrices = list(D0,D2,D1))
  expect_true(is.na(res$cv_results[3,7]))
  expect_true(res$cv_results[2,1] == 2)
  
  res <- diagram_ksvm(diagrams = g,cv = 2,dim = c(0,1,2),y = factor(c(rep("0",5),rep("1",5))),t = 1)
  expect_true(is.na(res$cv_results[3,7]))
  res <- diagram_ksvm(diagrams = g,cv = 2,dim = c(0,1,2),y = factor(c(rep("0",5),rep("1",5))))
  expect_true(is.na(res$cv_results[3,7]))
  res <- diagram_ksvm(diagrams = g,cv = 1,dim = c(0,1,2),y = factor(c(rep("0",5),rep("1",5))))
  expect_true(is.na(res$cv_results[3,7]))
  
})


test_that("predict_diagram_ksvm detects incorrect parameters correctly",{
  
  D1 <- data.frame(dimension = 0,birth = 2,death = 3)
  D2 <- data.frame(dimension = 0,birth = 2,death = 3.1)
  D3 <- data.frame(dimension = 0,birth = c(2,5),death = c(3.1,6))
  ksvm <- diagram_ksvm(diagrams = list(D1,D2,D3),dim = 0,y = c(1,2,3),t = c(1,2),num_workers = 2)
  expect_error(predict_diagram_ksvm(new_diagrams = list(),model = ksvm,num_workers = 2),"1")
  expect_error(predict_diagram_ksvm(new_diagrams = NULL,model = ksvm,num_workers = 2),"NULL")
  expect_error(predict_diagram_ksvm(new_diagrams = list(D1,"1"),model = ksvm,num_workers = 2),"Diagrams")
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
  expect_equal(as.character(predict_diagram_ksvm(new_diagrams = diagrams,model = ksvm,num_workers = 2)),c(rep("circle",5),rep("torus",5),rep("sphere",5)))
  
})

test_that("predict_diagram_ksvm can accept inputs from TDA, TDAstats and diagram_to_df",{
  
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

test_that("predict_diagram_ksvm can accept pre-computed Gram matrices",{
  
  skip_on_cran()
  skip_if_not_installed("TDA")
  skip_if_not_installed("TDAstats")
  
  D1 = TDA::ripsDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxscale = 1,maxdimension = 1)
  D2 = TDA::ripsDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxscale = 1,maxdimension = 1,library = "dionysus",location = T)
  D3 = TDAstats::calculate_homology(data.frame(x = runif(50,0,1),y = runif(50,0,1)),threshold = 1)
  d0 = distance_matrix(diagrams = list(D1,D2,D3),dim = 0,num_workers = 2,sigma = 1,distance = "fisher")
  d1 = distance_matrix(diagrams = list(D1,D2,D3),dim = 1,num_workers = 2,sigma = 1,distance = "fisher")
  model = diagram_ksvm(diagrams = list(D1,D2,D3),y = c(1,2,3),distance_matrices = list(d0,d1),num_workers = 2,dim = c(0,1))
  if(model$best_model$dim == 0)
  {
    K = exp(-1*d0)
  }else
  {
    K = exp(-1*d1)
  }
  K_small = K[1:2,1:2]
  class(K) = "kernelMatrix"
  class(K_small) = "kernelMatrix"
  expect_error(predict_diagram_ksvm(model = model,K = K_small,num_workers = 2),"number")
  expect_equal(predict_diagram_ksvm(model = model,K = K,num_workers = 2),predict_diagram_ksvm(new_diagrams = list(D1,D2,D3),model = model,num_workers = 2),tolerance = 1e-5)
  
})

