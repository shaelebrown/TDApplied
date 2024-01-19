
# test_that("analyze_representatives can detect incorrect parameters correctly",{
#   
#   skip_if_not_installed("TDA")
#   
#   D = data.frame(dimension = c(0),birth = c(0),death = c(1))
#   expect_error(analyze_representatives(diagrams = list(),dim = 1,num_points = 10),"2")
#   expect_error(analyze_representatives(diagrams = 5,dim = 1,num_points = 10),"list")
#   expect_error(analyze_representatives(diagrams = list(D,D,D),dim = -1,num_points = 10),"at least one")
#   expect_error(analyze_representatives(diagrams = list(D,D,D),dim = 1,num_points = 3),"4")
#   
#   expect_error(analyze_representatives(diagrams = list(D,D,D),dim = 1,num_points = 4,return_contributions = "F"),"logical")
#   expect_error(analyze_representatives(diagrams = list(D,D,D),dim = 1,num_points = 4,return_contributions = c(F,T)),"single")
#   expect_error(analyze_representatives(diagrams = list(D,D,D),dim = 1,num_points = 4,return_contributions = NULL),"NULL")
#   expect_error(analyze_representatives(diagrams = list(D,D,D),dim = 1,num_points = 4,return_contributions = NA),"NA")
#   
#   expect_error(analyze_representatives(diagrams = list(D,D,D),dim = 1,num_points = 4,plot_heatmap = "F"),"logical")
#   expect_error(analyze_representatives(diagrams = list(D,D,D),dim = 1,num_points = 4,plot_heatmap = c(F,T)),"single")
#   expect_error(analyze_representatives(diagrams = list(D,D,D),dim = 1,num_points = 4,plot_heatmap = NULL),"NULL")
#   expect_error(analyze_representatives(diagrams = list(D,D,D),dim = 1,num_points = 4,plot_heatmap = NA),"NA")
#   
#   expect_error(analyze_representatives(diagrams = list(D,D,2),dim = 1,num_points = 4),"diagram")
#   expect_error(analyze_representatives(diagrams = list(D,D,D),dim = 1,num_points = 4),"Representative")
#   expect_error(analyze_representatives(diagrams = list(list(D,representatives = 2),D,D),dim = 1,num_points = 4),"list")
#   expect_error(analyze_representatives(diagrams = list(list(D,representatives = list(matrix(data = c(1.1,2),nrow = 1))),D,D),dim = 1,num_points = 4),"integer")
#   expect_error(analyze_representatives(diagrams = list(list(D,representatives = list(matrix(data = c(1,5),nrow = 1))),D,D),dim = 1,num_points = 4),"num_points")
#   expect_error(analyze_representatives(diagrams = list(list(D,representatives = list(matrix(data = c(1,4),nrow = 1))),D,D),dim = 1,num_points = 4),"diagram")
#   
#   circs_ripsDiag <- lapply(X = 1:10,FUN = function(X){
#     
#     return(TDA::ripsDiag(X = as.matrix(dist(TDA::circleUnif(n = 25))),maxdimension = 1,maxscale = 2,library = "dionysus",location = T,dist = "arbitrary"))
#     
#   })
#   expect_error(analyze_representatives(diagrams = circs_ripsDiag,dim = 1,num_points = 25,boxed_reps = NA),"boxed_reps")
#   expect_error(analyze_representatives(diagrams = circs_ripsDiag,dim = 1,num_points = 25,boxed_reps = data.frame()),"two")
#   expect_error(analyze_representatives(diagrams = circs_ripsDiag,dim = 1,num_points = 25,boxed_reps = data.frame(diagrams = c(1),reps = c(1))),"rep")
#   
#   expect_error(analyze_representatives(diagrams = circs_ripsDiag,dim = 1,num_points = 25,lwd = "2",boxed_reps = data.frame(diagram = c(1),rep = c(1))),"numeric")
#   expect_error(analyze_representatives(diagrams = circs_ripsDiag,dim = 1,num_points = 25,lwd = c(1,2),boxed_reps = data.frame(diagram = c(1),rep = c(1))),"single")
#   expect_error(analyze_representatives(diagrams = circs_ripsDiag,dim = 1,num_points = 25,lwd = 0,boxed_reps = data.frame(diagram = c(1),rep = c(1))),"positive")
#   expect_error(analyze_representatives(diagrams = circs_ripsDiag,dim = 1,num_points = 25,lwd = NA,boxed_reps = data.frame(diagram = c(1),rep = c(1))),"NA")
#   expect_error(analyze_representatives(diagrams = circs_ripsDiag,dim = 1,num_points = 25,lwd = Inf,boxed_reps = data.frame(diagram = c(1),rep = c(1))),"finite")
#   
#   expect_error(analyze_representatives(diagrams = circs_ripsDiag,dim = 1,num_points = 25,title = NaN),"NaN")
#   expect_error(analyze_representatives(diagrams = circs_ripsDiag,dim = 1,num_points = 25,title = 1),"character")
#   expect_error(analyze_representatives(diagrams = circs_ripsDiag,dim = 1,num_points = 25,title = c("1","2")),"single")
#   
#   expect_error(analyze_representatives(diagrams = list(D,D,D),dim = 1,num_points = 4,return_clust = "F"),"logical")
#   expect_error(analyze_representatives(diagrams = list(D,D,D),dim = 1,num_points = 4,return_clust = c(F,T)),"single")
#   expect_error(analyze_representatives(diagrams = list(D,D,D),dim = 1,num_points = 4,return_clust = NULL),"NULL")
#   expect_error(analyze_representatives(diagrams = list(D,D,D),dim = 1,num_points = 4,return_clust = NA),"NA")
#   
#   expect_error(analyze_representatives(diagrams = circs_ripsDiag,dim = 1,num_points = 25,d = matrix(data = 0,nrow = 2,ncol = 2)),"dist")
#   expect_error(analyze_representatives(diagrams = circs_ripsDiag,dim = 1,num_points = 25,d = stats::dist(matrix(data = 0,nrow = 2,ncol = 2))),"rows")
# 
# })

# test_that("analyze_representatives is computing properly",{
#   
#   skip_on_cran()
#   skip_if_not_installed("TDA")
#   skip_if(T)
#   
#   ripser <- import_ripser()
#   
#   circs_ripsDiag <- lapply(X = 1:10,FUN = function(X){
#     
#     return(TDA::ripsDiag(X = as.matrix(dist(TDA::circleUnif(n = 25))),maxdimension = 1,maxscale = 2,library = "dionysus",location = T,dist = "arbitrary"))
#     
#   })
#   circs_PyH <- lapply(X = 1:10,FUN = function(X){
#     
#     return(PyH(X = as.matrix(dist(TDA::circleUnif(n = 25))),maxdim = 1,thresh = 2,ripser = ripser,calculate_representatives = T,distance_mat = T))
#     
#   })
#   circs_ripsDiag_boot <- lapply(X = 1:10,FUN = function(X){
#     
#     return(bootstrap_persistence_thresholds(X = as.matrix(dist(TDA::circleUnif(n = 25))),maxdim = 1,thresh = 2,distance_mat = T,calculate_representatives = T,FUN_diag = "ripsDiag",return_subsetted = T))
#     
#   })
#   circs_PyH_boot <- lapply(X = 1:10,FUN = function(X){
#     
#     return(bootstrap_persistence_thresholds(X = as.matrix(dist(TDA::circleUnif(n = 25))),maxdim = 1,thresh = 2,distance_mat = T,calculate_representatives = T,FUN_diag = "PyH",ripser = ripser,return_subsetted = T))
#     
#   })
#   
#   tryCatch(expr = {check1 <- analyze_representatives(circs_ripsDiag,dim = 1,num_points = 25,plot_heatmap = F)},
#            error = function(err){
#              a <- 1
#            },
#            finally = {
#              if(exists("check1")){expect_identical(dim(check1)[[2]],25L)}
#                })
#   
#   tryCatch(expr = {check2 <- analyze_representatives(circs_PyH,dim = 1,num_points = 25,plot_heatmap = F)},
#            error = function(err){
#              a <- 1
#            },
#            finally = {
#              if(exists("check2")){expect_identical(dim(check2)[[2]],25L)}
#            })
#   
#   tryCatch(expr = {check3 <- analyze_representatives(circs_ripsDiag_boot,dim = 1,num_points = 25,plot_heatmap = F)},
#            error = function(err){
#              a <- 1
#            },
#            finally = {
#              if(exists("check3")){expect_identical(dim(check3)[[2]],25L)}
#            })
#   
#   tryCatch(expr = {check4 <- analyze_representatives(circs_PyH_boot,dim = 1,num_points = 25,plot_heatmap = F)},
#            error = function(err){
#              a <- 1
#            },
#            finally = {
#              if(exists("check4")){expect_identical(dim(check4)[[2]],25L)}
#            })
#   
#   tryCatch(expr = {check5 <- analyze_representatives(circs_PyH_boot,dim = 1,num_points = 25,plot_heatmap = F,return_contributions = T)},
#            error = function(err){
#              a <- 1
#            },
#            finally = {
#              if(exists("check5")){expect_identical(length(check5$contributions),25L)}
#            })
#   
#   expect_error(analyze_representatives(circs_PyH_boot,dim = 2,num_points = 25,plot_heatmap = F,return_contributions = T),"dimension")
#   
#   circs <- lapply(X = 1:10,FUN = function(X){
#     
#     return(TDA::circleUnif(n = 25,r = 1))
#     
#   })
#   circs_PH <- lapply(X = circs,FUN = function(X){
#     
#     TDA::ripsDiag(X = as.matrix(dist(X)),maxdimension = 1,maxscale = 2,library = "dionysus",location = T,dist = "arbitrary")
#     
#   })
#   d <- matrix(data = 0,nrow = length(circs),ncol = length(circs))
#   r <- abs(stats::rnorm(45))
#   d[which(upper.tri(d),arr.ind = T)] <- r
#   d[which(upper.tri(d),arr.ind = T)[,c("col","row")]] <- r
#   d <- stats::as.dist(d)
#   check6 <- analyze_representatives(circs_PH,dim = 1,num_points = 25,plot_heatmap = T,d = d)
#   expect_identical(ncol(check6),25L)
# 
# })


