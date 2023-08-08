
test_that("analyze_representatives can detect incorrect parameters correctly",{
  
  D = data.frame(dimension = c(0),birth = c(0),death = c(1))
  expect_error(analyze_representatives(diagrams = list(),dim = 1,num_points = 10),"2")
  expect_error(analyze_representatives(diagrams = 5,dim = 1,num_points = 10),"list")
  expect_error(analyze_representatives(diagrams = list(D,D,D),dim = -1,num_points = 10),"at least one")
  expect_error(analyze_representatives(diagrams = list(D,D,D),dim = 1,num_points = 3),"4")
  
  expect_error(analyze_representatives(diagrams = list(D,D,D),dim = 1,num_points = 4,return_contributions = "F"),"logical")
  expect_error(analyze_representatives(diagrams = list(D,D,D),dim = 1,num_points = 4,return_contributions = c(F,T)),"single")
  expect_error(analyze_representatives(diagrams = list(D,D,D),dim = 1,num_points = 4,return_contributions = NULL),"NULL")
  expect_error(analyze_representatives(diagrams = list(D,D,D),dim = 1,num_points = 4,return_contributions = NA),"NA")
  
  expect_error(analyze_representatives(diagrams = list(D,D,D),dim = 1,num_points = 4,plot_heatmap = "F"),"logical")
  expect_error(analyze_representatives(diagrams = list(D,D,D),dim = 1,num_points = 4,plot_heatmap = c(F,T)),"single")
  expect_error(analyze_representatives(diagrams = list(D,D,D),dim = 1,num_points = 4,plot_heatmap = NULL),"NULL")
  expect_error(analyze_representatives(diagrams = list(D,D,D),dim = 1,num_points = 4,plot_heatmap = NA),"NA")
  
  expect_error(analyze_representatives(diagrams = list(D,D,2),dim = 1,num_points = 4),"diagram")
  expect_error(analyze_representatives(diagrams = list(D,D,D),dim = 1,num_points = 4),"Representative")
  expect_error(analyze_representatives(diagrams = list(list(D,representatives = 2),D,D),dim = 1,num_points = 4),"list")
  expect_error(analyze_representatives(diagrams = list(list(D,representatives = list(matrix(data = c(1.1,2),nrow = 1))),D,D),dim = 1,num_points = 4),"integer")
  expect_error(analyze_representatives(diagrams = list(list(D,representatives = list(matrix(data = c(1,5),nrow = 1))),D,D),dim = 1,num_points = 4),"num_points")
  expect_error(analyze_representatives(diagrams = list(list(D,representatives = list(matrix(data = c(1,4),nrow = 1))),D,D),dim = 1,num_points = 4),"diagram")

})

test_that("analyze_representatives is computing properly",{
  
  skip_on_cran()
  skip_if_not_installed("TDA")
  skip_if(T)
  
  ripser <- import_ripser()
  
  circs_ripsDiag <- lapply(X = 1:10,FUN = function(X){
    
    return(TDA::ripsDiag(X = as.matrix(dist(TDA::circleUnif(n = 25))),maxdimension = 1,maxscale = 1,library = "dionysus",location = T,dist = "arbitrary"))
    
  })
  circs_PyH <- lapply(X = 1:10,FUN = function(X){
    
    return(PyH(X = as.matrix(dist(TDA::circleUnif(n = 25))),maxdim = 1,thresh = 1,ripser = ripser,calculate_representatives = T,distance_mat = T))
    
  })
  circs_ripsDiag_boot <- lapply(X = 1:10,FUN = function(X){
    
    return(bootstrap_persistence_thresholds(X = as.matrix(dist(TDA::circleUnif(n = 25))),maxdim = 1,thresh = 2,distance_mat = T,calculate_representatives = T,FUN = "ripsDiag",return_subsetted = T))
    
  })
  circs_PyH_boot <- lapply(X = 1:10,FUN = function(X){
    
    return(bootstrap_persistence_thresholds(X = as.matrix(dist(TDA::circleUnif(n = 25))),maxdim = 1,thresh = 2,distance_mat = T,calculate_representatives = T,FUN = "PyH",ripser = ripser,return_subsetted = T))
    
  })
  
  check1 <- analyze_representatives(circs_ripsDiag,dim = 1,num_points = 25,plot_heatmap = F)
  expect_identical(dim(check1)[[2]],25L)
  
  check2 <- analyze_representatives(circs_PyH,dim = 1,num_points = 25,plot_heatmap = F)
  expect_identical(dim(check2)[[2]],25L)
  
  check3 <- analyze_representatives(circs_ripsDiag_boot,dim = 1,num_points = 25,plot_heatmap = F)
  expect_identical(dim(check3)[[2]],25L)
  
  check4 <- analyze_representatives(circs_PyH_boot,dim = 1,num_points = 25,plot_heatmap = F)
  expect_identical(dim(check4)[[2]],25L)
  
  check5 <- analyze_representatives(circs_PyH_boot,dim = 1,num_points = 25,plot_heatmap = F,return_contributions = T)
  expect_identical(length(check5$contributions),25L)
  
  expect_error(analyze_representatives(circs_PyH_boot,dim = 2,num_points = 25,plot_heatmap = F,return_contributions = T),"dimension")

})


