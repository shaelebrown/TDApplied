
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
  
  skip_if_not_installed("TDA")
  skip_if_not_installed("TDAstats")
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
  expect_identical(diagram_distance(D1 = D1,D2 = D2,dim = 1),0)
  
  D1$dimension = 1
  expect_identical(diagram_distance(D1 = D1,D2 = D2,dim = 1,p = 2,distance = "wasserstein"),sqrt(0.5^2))
  expect_identical(diagram_distance(D1 = D1,D2 = D2,dim = 1,p = Inf,distance = "wasserstein"),0.5)
  expect_identical(diagram_distance(D1 = D2,D2 = D1,dim = 1,p = 2,distance = "wasserstein"),sqrt(0.5^2))
  expect_identical(diagram_distance(D1 = D2,D2 = D1,dim = 1,p = Inf,distance = "wasserstein"),0.5)

  # this example was picked the TDA function wasserstein disagrees with the actual minimum values
  # for p = 2,3, but diagram_distance gets the correct answer
  D1 = data.frame(dimension = c(0,0),birth = c(0,0),death = c(0.9640122,1.3467424))
  D2 = data.frame(dimension = c(0,0),birth = c(0,0),death = c(1.233867,1.398447))
  phom1 = D1
  phom2 = D2
  
  D1_subset <- D1[,2:3]
  D2_subset <- D2[,2:3]
  diag1 <- D1_subset[0,]
  diag2 <- D2_subset[0,]
  D1_subset <- D1_subset[which(D1_subset[,1] != D1_subset[,2]),]
  D2_subset <- D2_subset[which(D2_subset[,1] != D2_subset[,2]),]
  if(nrow(D1_subset) > 0)
  {
    for(i in 1:nrow(D1_subset))
    {
      proj_diag <- mean(as.numeric(D1_subset[i,]))
      diag1 <- rbind(diag1,data.frame(birth = proj_diag,death = proj_diag))
    }
  }
  if(nrow(D2_subset) > 0)
  {
    for(i in 1:nrow(D2_subset))
    {
      proj_diag <- mean(as.numeric(D2_subset[i,]))
      diag2 <- rbind(diag2,data.frame(birth = proj_diag,death = proj_diag))
    }
  }
  D1_subset <- rbind(D1_subset,diag2)
  D2_subset <- rbind(D2_subset,diag1)

  dist_mat_bottleneck <- as.matrix(rdist::cdist(D1_subset,D2_subset,metric = "maximum"))
  dist_mat_2 <- dist_mat_bottleneck^2
  dist_mat_3 <- dist_mat_bottleneck^3
  
  min_bottleneck = Inf
  min_wass_2 = Inf
  min_wass_3 = Inf
  
  perms = matrix(data = c(1,2,3,4,1,2,4,3,1,3,2,4,1,3,4,2,1,4,2,3,1,4,3,2,2,1,3,4,2,1,4,3,2,3,1,4,2,3,4,1,2,4,1,3,2,4,3,1,3,1,2,4,3,1,4,2,3,2,1,4,3,2,4,1,3,4,1,2,3,4,2,1,4,1,2,3,4,1,3,2,4,2,1,3,4,2,3,1,4,3,1,2,4,3,2,1),nrow = 24,ncol = 4,byrow = T)
  class(perms) <- c("matrix","array")
  
  for(i in 1:nrow(perms))
  {
    for(j in 1:nrow(perms))
    {
      if(i != j)
      {
        temp = cbind(data.frame(x = perms[i,]),data.frame(y = perms[j,]))
        temp = as.matrix(temp[which(temp[,1] <= (nrow(D1_subset) - nrow(diag2)) | temp[,2] <= (nrow(D2_subset) - nrow(diag1))),])
        if(max(dist_mat_bottleneck[temp]) < min_bottleneck)
        {
          min_bottleneck <- max(dist_mat_bottleneck[temp])
        }
        if(sqrt(sum(dist_mat_2[temp])) < min_wass_2)
        {
          min_wass_2 <- sqrt(sum(dist_mat_2[temp]))
        }
        if((sum(dist_mat_3[temp]))^(1/3) < min_wass_3)
        {
          min_wass_3 <- (sum(dist_mat_3[temp]))^(1/3)
        }
      }
    }
  }
  
  expect_equal(diagram_distance(phom1,phom2,p = 2),min_wass_2)
  expect_equal(diagram_distance(phom1,phom2,p = 3),min_wass_3)
  expect_equal(diagram_distance(phom1,phom2,p = Inf),min_bottleneck)
  
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
