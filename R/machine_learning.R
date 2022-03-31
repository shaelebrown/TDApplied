
#### Multidimensional scaling ####
#' Calculate the metric multidimensional scaling of a group of persistence diagrams
#'
#' Returns the output of cmdscale on the desired distance matrix of a group of persistence diagrams
#' in a particular dimension.
#'
#' The `diagrams` parameter should be a list of persistence diagrams computed from TDA.
#' The `distance` parameter is a string representing which determines which distance metric to use.
#' The `dim` parameter should be a positive finite integer.
#' The `sigma` parameter is the positive bandwith for the Fisher information metric, and
#' `p` is the wasserstein power parameter. `k` is the desired dimension of the embedding, and
#' `eig`, `add`, `x.ret` and `list.` are cmdscale parameters providing optional additional
#' information to be returned.
#'
#' @param diagrams a list of persistence diagrams, as the output of a TDA calculation.
#' @param distance a string representing the desired distance metric to be used, either 'wasserstein' (default) or 'fisher'.
#' @param dim the homological dimension in which the distance is to be computed.
#' @param p the wasserstein power, a number at least 1 (infinity for the bottleneck distance), default 2.
#' @param sigma a positive number representing the bandwith for the Fisher information metric, default NULL.
#' @param k the maximum dimension of the space which the data are to be represented in; must be in {1,2,...,n-1}.
#' @param eig indicates whether eigenvalues should be returned.
#' @param add logical indicating if an additive constant c* should be computed, and added to the non-diagonal dissimilarities such that the modified dissimilarities are Euclidean.
#' @param x.ret indicates whether the doubly centered symmetric distance matrix should be returned.
#' @param list. local indicating if a list should be returned or just the n*k matrix.
#'
#' @return the output of cmdscale on the diagram distance matrix, either just the embedding matrix or a list.
#' @importFrom stats cmdscale
#' @importFrom parallel makeCluster stopCluster clusterExport clusterEvalQ
#' @importFrom parallelly availableCores
#' @importFrom doParallel registerDoParallel
#' @export
#' @examples
#'
#' # create ten diagrams with package TDA based on 2D Gaussians
#' g <- lapply(X = 1:10,FUN = function(X){
#'
#' diag <- TDA::ripsDiag(data.frame(x = rnorm(100,mean = 0,sd = 1),
#' y = rnorm(100,mean = 0,sd = 1)),
#' maxscale = 1,
#' maxdimension = 1)
#' df <- diagram_to_df(d = diag)
#' return(df)
#'
#' })
#'
#' # calculate their 2D mds embedding in dimension 1 with the bottleneck distance
#' embedding <- diagram_MDS(diagrams = g,dim = 1,p = Inf,k = 2)

diagram_MDS <- function(diagrams,distance = "wasserstein",dim = 0,p = 2,sigma = NULL,k = 2,eig = FALSE,add = FALSE,x.ret = FALSE,list. = eig || add || x.ret){
  
  # error check diagrams argument
  if(is.null(diagrams))
  {
    stop("diagrams must be a list of persistence diagrams.")
  }
  if(!is.list(diagrams) | length(diagrams) < 2)
  {
    stop("diagrams must be a list of persistence diagrams of length at least 2.")
  }
  diagrams <- all_diagrams(diagram_groups = list(diagrams),inference = "independence")[[1]]
  
  # error check other distance parameters
  check_params(iterations = 10,p = p,q = 2,dims = c(dim),paired = F,distance = distance,sigma = sigma)
  
  # compute distance matrix
  d <- distance_matrix(diagrams = diagrams,dim = dim,distance = distance,p = p,sigma = sigma)

  # return metric multidimensional scaling with d as input
  return(stats::cmdscale(d = d,k = k,eig = eig,add = add,x.ret = x.ret,list. = list.))
  
}

#### KERNEL PCA ####
#' Calculate the kernel PCA embedding of a group of persistence diagrams
#'
#' Returns the output of cmdscale on the desired distance matrix of a group of persistence diagrams
#' in a particular dimension.
#'
#' The `diagrams` parameter should be a list of persistence diagrams computed from TDA.
#' The `dim` parameter should be a positive finite integer.
#' The `sigma` and `t` parameters are the positive bandwith for the Fisher information metric and
#' the positive scale for the persistence Fisher kernel respectively.
#' `features` is the number of desired features (principal components) in the embedding, and
#' `...` are additional parameters to eh kpca kernlab function.
#'
#' @param diagrams a list of persistence diagrams, as the output of a TDA calculation.
#' @param dim the homological dimension in which the distance is to be computed.
#' @param t the positive scale for the persistence Fisher kernel, default 1.
#' @param sigma a positive number representing the bandwith for the Fisher information metric, default 1
#' @param features number of features (principal components) to return, default 1.
#' @param ... additional parameters.
#'
#' @return a list containing the output of cmdscale on the diagram distance matrix, either just the embedding matrix or a list, the diagram groups, dimension, t, sigma and features. The class of this object is 'diagram_kpca'.
#' @export
#' @examples
#'
#' # create ten diagrams with package TDA based on 2D Gaussians
#' g <- lapply(X = 1:10,FUN = function(X){
#'
#' diag <- TDA::ripsDiag(data.frame(x = rnorm(100,mean = 0,sd = 1),
#' y = rnorm(100,mean = 0,sd = 1)),
#' maxscale = 1,
#' maxdimension = 1)
#' df <- diagram_to_df(d = diag)
#' return(df)
#'
#' })
#'
#' # calculate their 2D PCA embedding in dimension 1 with sigma = t = 2
#' embedding <- diagram_kpca(diagrams = g,dim = 1,t = 2,sigma = 2,features = 2)

diagram_kpca <- function(diagrams,dim = 0,t = 1,sigma = 1,features = 1,...){
  
  # error check diagrams argument
  if(is.null(diagrams))
  {
    stop("diagrams must be a list of persistence diagrams.")
  }
  if(!is.list(diagrams) | length(diagrams) < 2)
  {
    stop("diagrams must be a list of persistence diagrams of length at least 2.")
  }
  diagrams <- all_diagrams(diagram_groups = list(diagrams),inference = "independence")[[1]]
  
  # error check sigma and dim parameters
  check_params(iterations = 10,p = 2,q = 2,dims = c(dim),paired = T,distance = "fisher",sigma)
  
  # error check t parameter
  if(is.null(t))
  {
    stop("t must not be NULL.")
  }
  if(!is.numeric(t) | length(t) > 1 | is.na(t) | is.nan(t) | t <= 0)
  {
    stop("t must be a positive number.")
  }
  
  # compute Gram matrix
  K <- gram_matrix(diagrams = diagrams,t = t,sigma = sigma,dim = dim)
  
  # return kernlab computation
  ret_list <- list(pca = kernlab::kpca(x = K,features = features,...),diagrams = diagrams,t = t,sigma = sigma,dim = dim)
  class(ret_list) <- "diagram_kpca"
  return(ret_list)
  
}

#### PREDICT WITH KERNEL PCA OBJECT ####
#' Predict the kernel PCA embedding of new persistence diagrams using a precomputed diagram_kpca object
#'
#' Returns the embedding matrix for the new points.
#'
#' The `new_diagrams` parameter should be a list of persistence diagrams computed from TDA.
#' The `embedding` parameter is the diagram_kpca embedding object to be used for embedding
#' the new diagrams.
#'
#' @param new_diagrams a list of persistence diagrams, as the output of a TDA calculation.
#' @param embedding the output to a diagram_kpca function call.
#'
#' @return the embedding matrix.
#' @export
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom parallel makeCluster stopCluster clusterExport clusterEvalQ
#' @importFrom parallelly availableCores
#' @importFrom doParallel registerDoParallel
#' @examples
#'
#' # create ten diagrams with package TDA based on 2D Gaussians
#' g <- lapply(X = 1:10,FUN = function(X){
#'
#' diag <- TDA::ripsDiag(data.frame(x = rnorm(100,mean = 0,sd = 1),
#' y = rnorm(100,mean = 0,sd = 1)),
#' maxscale = 1,
#' maxdimension = 1)
#' df <- diagram_to_df(d = diag)
#' return(df)
#'
#' })
#'
#' # calculate their 2D PCA embedding in dimension 1 with sigma = t = 2
#' embedding <- diagram_kpca(diagrams = g,dim = 1,t = 2,sigma = 2,features = 2)
#' 
#' # create ten new diagrams with package TDA based on 2D Gaussians
#' g_new <- lapply(X = 1:10,FUN = function(X){
#'
#' diag <- TDA::ripsDiag(data.frame(x = rnorm(100,mean = 0,sd = 1),
#' y = rnorm(100,mean = 0,sd = 1)),
#' maxscale = 1,
#' maxdimension = 1)
#' df <- diagram_to_df(d = diag)
#' return(df)
#'
#' })
#'
#' # calculate their 2D PCA embedding using the embedding object
#' embedding_new <- predict_diagram_kpca(new_diagrams = g_new,embedding = embedding)

predict_diagram_kpca <- function(new_diagrams,embedding){
  
  # set internal variables to NULL to avoid build issues
  r <- NULL
  X <- NULL
  
  # error check new_diagrams argument
  if(is.null(new_diagrams))
  {
    stop("new_diagrams must be a list of persistence diagrams.")
  }
  if(!is.list(new_diagrams) | length(new_diagrams) < 1)
  {
    stop("new_diagrams must be a list of persistence diagrams of length at least 1.")
  }
  new_diagrams <- all_diagrams(diagram_groups = list(new_diagrams),inference = "independence")[[1]]
  
  # error check embedding argument
  if(is.null(embedding))
  {
    stop("embedding must be supplied.")
  }
  if(class(embedding) != "diagram_kpca")
  {
    stop("embedding must be the output of a diagram_kpca function call.")
  }
  
  # compute kernel matrix, storing the value of each kernel computation between the new diagrams and the old ones
  num_workers <- parallelly::availableCores(omit = 1)
  cl <- parallel::makeCluster(num_workers)
  doParallel::registerDoParallel(cl)
  parallel::clusterExport(cl,c("diagram_distance","diagram_kernel"))
  force(new_diagrams) # required for parallel computation in this environment
  force(embedding)
  force(check_diagram)
  
  # compute in parallel along the diagram list (new or old) which has more elements
  if(length(new_diagrams) > length(embedding$diagrams))
  {
    K <- foreach::`%dopar%`(obj = foreach::foreach(r = 1:length(new_diagrams),.combine = rbind),ex = {
      
      return(unlist(lapply(X = 1:length(embedding$diagrams),FUN = function(X){return(diagram_kernel(D1 = new_diagrams[[r]],D2 = embedding$diagrams[[X]],dim = embedding$dim,sigma = embedding$sigma,t = embedding$t))})))
      
    })
  }else
  {
    K <- foreach::`%do%`(obj = foreach::foreach(r = 1:length(new_diagrams),.combine = rbind),ex = {
      
      return(foreach::`%dopar%`(obj = foreach::foreach(X = 1:length(new_diagrams),.combine = c),ex = {
        
        return(diagram_kernel(D1 = new_diagrams[[r]],D2 = embedding$diagrams[[X]],dim = embedding$dim,sigma = embedding$sigma,t = embedding$t))
        
      }))

    })
  }
  
  parallel::stopCluster(cl)
  
  # project the new diagrams into the embedding space
  return(K %*% embedding$pca@pcv)
  
}





