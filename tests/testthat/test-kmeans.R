
test_that("diagram_kkmeans detects incorrect parameters correctly",{
  
  D <- data.frame(dimension = c(0),birth = c(0),death = c(1))
  expect_error(diagram_kkmeans(diagrams = list(D,D,D[0,]),centers = 2,num_workers = 2),"empty")
  expect_error(diagram_kkmeans(diagrams = list(),centers = 2,num_workers = 2),"2")
  expect_error(diagram_kkmeans(diagrams = list(D,D,D),centers = 1,t = NaN,num_workers = 2),"t")
  expect_error(diagram_kkmeans(diagrams = list(D,D,D),centers = 1,sigma = NA,num_workers = 2),"sigma")
  expect_error(diagram_kkmeans(diagrams = list(D,D,D),dim = c(1,2),centers = 1,num_workers = 2),"single value")
  expect_error(diagram_kkmeans(diagrams = list(D,D,D),dim = c(1),centers = 4,num_workers = 2),"number of")
  
})

test_that("diagram_kkmeans is computing correctly",{
  
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
  
  cluster_labels_dim_1_circle_torus <- diagram_kkmeans(diagrams = diagrams[1:10],centers = 2,dim = 1,num_workers = 2)$clustering@.Data
  if(cluster_labels_dim_1_circle_torus[[1]] == 1)
  {
    expect_equal(cluster_labels_dim_1_circle_torus,c(1,1,1,1,1,2,2,2,2,2))
  }else
  {
    expect_equal(cluster_labels_dim_1_circle_torus,c(2,2,2,2,2,1,1,1,1,1))
  }
  cluster_labels_dim_2_circle_torus <- diagram_kkmeans(diagrams = diagrams[1:10],centers = 2,dim = 2,num_workers = 2)$clustering@.Data
  if(cluster_labels_dim_2_circle_torus[[1]] == 1)
  {
    expect_equal(cluster_labels_dim_2_circle_torus,c(1,1,1,1,1,2,2,2,2,2))
  }else
  {
    expect_equal(cluster_labels_dim_2_circle_torus,c(2,2,2,2,2,1,1,1,1,1))
  }
  
  cluster_labels_dim_1_circle_sphere <- diagram_kkmeans(diagrams = diagrams[c(1:5,11:15)],centers = 2,dim = 1,num_workers = 2)$clustering@.Data
  if(cluster_labels_dim_1_circle_sphere[[1]] == 1)
  {
    expect_equal(cluster_labels_dim_1_circle_sphere,c(1,1,1,1,1,2,2,2,2,2))
  }else
  {
    expect_equal(cluster_labels_dim_1_circle_sphere,c(2,2,2,2,2,1,1,1,1,1))
  }
  cluster_labels_dim_2_circle_sphere <- diagram_kkmeans(diagrams = diagrams[c(1:5,11:15)],centers = 2,dim = 2,num_workers = 2)$clustering@.Data
  if(cluster_labels_dim_2_circle_sphere[[1]] == 1)
  {
    expect_equal(cluster_labels_dim_2_circle_sphere,c(1,1,1,1,1,2,2,2,2,2))
  }else
  {
    expect_equal(cluster_labels_dim_2_circle_sphere,c(2,2,2,2,2,1,1,1,1,1))
  }
  
  cluster_labels_dim_1_sphere_torus <- diagram_kkmeans(diagrams = diagrams[6:15],centers = 2,dim = 1,num_workers = 2)$clustering@.Data
  if(cluster_labels_dim_1_sphere_torus[[1]] == 1)
  {
    expect_equal(cluster_labels_dim_1_sphere_torus,c(1,1,1,1,1,2,2,2,2,2))
  }else
  {
    expect_equal(cluster_labels_dim_1_sphere_torus,c(2,2,2,2,2,1,1,1,1,1))
  }
  cluster_labels_dim_2_sphere_torus <- diagram_kkmeans(diagrams = diagrams[6:15],centers = 2,dim = 2,num_workers = 2)$clustering@.Data
  if(cluster_labels_dim_2_sphere_torus[[1]] == 1)
  {
    expect_equal(cluster_labels_dim_2_sphere_torus,c(1,1,1,1,1,2,2,2,2,2))
  }else
  {
    expect_equal(cluster_labels_dim_2_sphere_torus,c(2,2,2,2,2,1,1,1,1,1))
  }
  
})

test_that("diagram_kkmeans can accept inputs from TDA, TDAstats and diagram_to_df",{
  
  skip_on_cran()
  skip_if_not_installed("TDA")
  skip_if_not_installed("TDAstats")
  
  D1 = TDA::ripsDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxscale = 1,maxdimension = 1)
  D2 = TDA::alphaComplexDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxdimension = 1)
  D3 = TDA::ripsDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxscale = 1,maxdimension = 1,library = "dionysus",location = T)
  D4 = TDAstats::calculate_homology(data.frame(x = runif(50,0,1),y = runif(50,0,1)),threshold = 1)
  expect_length(diagram_kkmeans(diagrams = list(D1,D2,D3,D4,D1,D1,D1,D2,D2,D2),centers = 2,dim = 1,num_workers = 2)$clustering@.Data,10)
  expect_error(diagram_kkmeans(diagrams = list(D1,D2,D3,D4,D1,D1,D1,D2,D2,D2),centers = 2,dim = 0,num_workers = 2),"Inf")
  
})

test_that("diagram_kkmeans can accept a precomputed Gram matrix",{
  
  circle <- data.frame(dimension = c(0,1),birth = c(0,0),death = c(2,2))
  torus <- data.frame(dimension = c(0,1,1,2),birth = c(0,0,0,0),death = c(2,0.5,1.5,0.5))
  sphere <- data.frame(dimension = c(0,2),birth = c(0,0),death = c(2,2))
  
  K = gram_matrix(list(circle,torus,sphere),dim = 1,num_workers = 2)
  
  expect_error(diagram_kkmeans(diagrams = list(circle,torus,sphere),K = K,centers = 2,dim = 1,num_workers = 2),"empty")
  expect_error(diagram_kkmeans(diagrams = list(circle,torus,sphere),K = K,centers = 4,dim = 1,num_workers = 2),"column")
  expect_error(diagram_kkmeans(diagrams = list(circle,torus),K = K,centers = 2,dim = 1,num_workers = 2),"length")
  
  K = gram_matrix(diagrams = list(circle,torus,sphere,circle,torus,sphere,circle,torus,sphere),dim = 1,num_workers = 2)
  
  expect_type(diagram_kkmeans(diagrams = list(circle,torus,sphere,circle,torus,sphere,circle,torus,sphere),K = K,centers = 2,dim = 1,num_workers = 2),"list")
  
})

test_that("predict_diagram_kkmeans detects incorrect parameters correctly",{
  
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
  dkk <- diagram_kkmeans(diagrams = diagrams,centers = 2,dim = 1,num_workers = 2)
  expect_error(predict_diagram_kkmeans(new_diagrams = list(),clustering = dkk,num_workers = 2),"1")
  expect_error(predict_diagram_kkmeans(new_diagrams = "D",clustering = dkk,num_workers = 2),"list")
  # expect_error(predict_diagram_kkmeans(new_diagrams = list(diagrams[[1]],diagrams[[2]][0,]),dkk,num_workers = 2),"empty")
  expect_error(predict_diagram_kkmeans(new_diagrams = diagrams,clustering = diagrams,num_workers = 2),"kkmeans")
  expect_error(predict_diagram_kkmeans(new_diagrams = diagrams,clustering = NULL,num_workers = 2),"NULL")
  
})

test_that("predict_diagram_kkmeans is computing correctly",{
  
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
  dkk <- diagram_kkmeans(diagrams = diagrams,centers = 2,dim = 1,num_workers = 2)
  expect_equal(predict_diagram_kkmeans(new_diagrams = diagrams,clustering = dkk,num_workers = 2),dkk$clustering@.Data)
  
})

test_that("predict_diagram_kkmeans can accept inputs from TDA, TDAstats and diagram_to_df",{
  
  skip_on_cran()
  skip_if_not_installed("TDA")
  skip_if_not_installed("TDAstats")
  
  D1 = TDA::ripsDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxscale = 1,maxdimension = 1)
  D2 = TDA::alphaComplexDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxdimension = 1)
  D3 = TDA::ripsDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxscale = 1,maxdimension = 1,library = "dionysus",location = T)
  D4 = TDAstats::calculate_homology(data.frame(x = runif(50,0,1),y = runif(50,0,1)),threshold = 1)
  dkk <- diagram_kkmeans(diagrams = list(D1,D1,D1,D2,D2,D2,D3,D3,D3,D4,D4,D4),centers = 2,dim = 1,num_workers = 2)
  expect_length(predict_diagram_kkmeans(new_diagrams = list(D1,D2,D3,D4),clustering = dkk,num_workers = 2),4)

})

test_that("predict_diagram_kkmeans can accept a precomputed Gram matrix",{
  
  skip_on_cran()
  skip_if_not_installed("TDA")
  skip_if_not_installed("TDAstats")
  
  D1 = TDA::ripsDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxscale = 1,maxdimension = 1)
  D2 = TDA::alphaComplexDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxdimension = 1)
  D3 = TDA::ripsDiag(data.frame(x = runif(50,0,1),y = runif(50,0,1)),maxscale = 1,maxdimension = 1,library = "dionysus",location = T)
  D4 = TDAstats::calculate_homology(data.frame(x = runif(50,0,1),y = runif(50,0,1)),threshold = 1)
  K = gram_matrix(diagrams = list(D1,D2,D3,D4),other_diagrams = list(D1,D1,D1,D2,D2,D2,D3,D3,D3,D4,D4,D4),dim = 1,num_workers = 2)
  dkk <- diagram_kkmeans(diagrams = list(D1,D1,D1,D2,D2,D2,D3,D3,D3,D4,D4,D4),centers = 2,dim = 1,num_workers = 2)
  expect_length(predict_diagram_kkmeans(K = K,clustering = dkk,num_workers = 2),4)
  expect_length(predict_diagram_kkmeans(new_diagrams = list(D1),clustering = dkk,num_workers = 2),1)
  expect_identical(predict_diagram_kkmeans(K = K,clustering = dkk,num_workers = 2),predict_diagram_kkmeans(new_diagrams = list(D1,D2,D3,D4),clustering = dkk,num_workers = 2))
  
  K = matrix(data = c(1,2,3),nrow = 3)
  class(K) = "kernelMatrix"
  expect_error(predict_diagram_kkmeans(K = K,clustering = dkk,num_workers = 2),"rows")
  
})
