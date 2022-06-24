
#### Multidimensional scaling ####
#' Calculate the metric multidimensional scaling of a group of persistence diagrams
#'
#' Returns the output of cmdscale on the desired distance matrix of a group of persistence diagrams
#' in a particular dimension.
#'
#' The `diagrams` parameter should be a list of persistence diagrams computed from a TDA calculation like \code{\link[TDA]{ripsDiag}} or from \code{\link{diagram_to_df}}.
#' The `distance` parameter is a string representing which determines which distance metric to use.
#' The `dim` parameter should be a positive finite integer.
#' The `sigma` parameter is the positive bandwith for the Fisher information metric, and
#' `p` is the wasserstein power parameter. `k` is the desired dimension of the embedding, and
#' `eig`, `add`, `x.ret` and `list.` are cmdscale parameters providing optional additional
#' information to be returned.
#'
#' @param diagrams a list of persistence diagrams, as the output of a TDA calculation like \code{\link[TDA]{ripsDiag}} or a \code{\link{diagram_to_df}} function call.
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
#' @return the output of \code{\link[stats]{cmdscale}} on the diagram distance matrix. If `list.` is false (as per default),
#' a matrix with `k` columns whose rows give the coordinates of the points chosen to represent the dissimilarities.
#' 
#' Otherwise, a list containing the following components.
#' 
#' \describe{
#' 
#'  \item{points}{a matrix with up to `k` columns whose rows give the coordinates of the points chosen to represent the dissimilarities.}
#' 
#'  \item{eig}{the \eqn{n} eigenvalues computed during the scaling process if `eig` is true.}
#'  
#'  \item{x}{the doubly centered distance matrix if `x.ret` is true.}
#'  
#'  \item{ac}{the additive constant \eqn{c^*}, 0 if `add` = FALSE.}
#'  
#'  \item{GOF}{the numeric vector of length 2, representing the sum of all the eigenvalues divided by the sum of their absolute values (first vector element) or by the sum of the max of each eigenvalue and 0 (second vector element).}
#' 
#' }
#' 
#' @importFrom stats cmdscale
#' @importFrom parallel makeCluster stopCluster clusterExport clusterEvalQ
#' @importFrom parallelly availableCores
#' @importFrom doParallel registerDoParallel
#' @export
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
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
  check_param("diagrams",diagrams,numeric = F,multiple = T)
  diagrams <- all_diagrams(diagram_groups = list(diagrams),inference = "independence")[[1]]

  # compute distance matrix
  d <- distance_matrix(diagrams = diagrams,dim = dim,distance = distance,p = p,sigma = sigma)

  # return metric multidimensional scaling with d as input
  return(stats::cmdscale(d = d,k = k,eig = eig,add = add,x.ret = x.ret,list. = list.))
  
}

#### KERNEL KMEANS ####
#' Cluster a group of persistence diagrams using kernel k-means
#'
#' Returns the output of kkmeans on the desired Gram matrix of a group of persistence diagrams
#' in a particular dimension.
#'
#' The `diagrams` parameter should be a list of persistence diagrams computed from TDA.
#' The `dim` parameter should be a positive finite integer.
#' The `sigma` and `t` parameters are the positive bandwith for the Fisher information metric and
#' the positive scale for the persistence Fisher kernel respectively.
#' `centers` is the number of desired clusters, and
#' `...` are additional parameters to the kkmeans kernlab function.
#'
#' @param diagrams a list of persistence diagrams, as the output of a TDA calculation.
#' @param dim the homological dimension in which the distance is to be computed.
#' @param t the positive scale for the persistence Fisher kernel, default 1.
#' @param sigma a positive number representing the bandwith for the Fisher information metric, default 1
#' @param centers number of clusters to initialize.
#' @param ... additional parameters.
#'
#' @return a 'diagram_kkmeans' object containing the output of \code{\link[kernlab]{kkmeans}} on the diagram distance matrix, i.e. a list containing
#' 
#' \describe{
#' 
#' \item{clustering}{an S4 object of class specc, the output of a \code{\link[kernlab]{kkmeans}} function call. The .Data slot of this object contains cluster memberships, withinss contains the within-cluster sum of squares for each cluster, etc.}
#' 
#' \item{diagrams}{the input `diagrams` argument.}
#' 
#' \item{dim}{the input `dim` argument.}
#' 
#' \item{t}{the input `t` argument.}
#' 
#' \item{sigma}{the input `sigma` argument.}
#' 
#' }
#' 
#' @export
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @importFrom kernlab kkmeans
#' @seealso \code{\link{diagram_nearest_clusters}} for predicting clusters of new diagrams.
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
#' # calculate kmeans clusters with centers = 2 in dimension 1 with sigma = t = 2
#' clusters <- diagram_kkmeans(diagrams = g,centers = 2,dim = 1,t = 2,sigma = 2)

diagram_kkmeans <- function(diagrams,centers,dim = 0,t = 1,sigma = 1,...){
  
  # error check arguments
  check_param("diagrams",diagrams,numeric = F,multiple = T)
  diagrams <- all_diagrams(diagram_groups = list(diagrams),inference = "independence")[[1]]
  check_param("centers",centers,whole_numbers = T,at_least_one = T)
  
  # compute Gram matrix
  K <- gram_matrix(diagrams = diagrams,dim = dim,sigma = sigma,t = t)
  
  # return kernlab calculation, saving in a list of class diagram_kkmeans for later interface with prediction calculations
  ret_list <- list(clustering = kernlab::kkmeans(x = K,centers = centers,...),diagrams = diagrams,dim = dim,sigma = sigma,t = t)
  
  class(ret_list) <- "diagram_kkmeans"
  
  return(ret_list)
  
}

#### PREDICT KERNEL KMEANS ####
#' Find the nearest kkmeans cluster center to a list of new diagrams to return approximate cluster labels
#'
#' Returns the nearest kkmeans cluster center labels on the desired Gram matrix of a group of persistence diagrams
#' in a particular dimension.
#'
#' The `new_diagrams` parameter should be a list of persistence diagrams computed from a TDA calculation like \code{\link[TDA]{ripsDiag}} or from a 
#' \code{\link{diagram_to_df}} function call, and the
#' `clustering` parameter should be the output of a \code{\link{diagram_kkmeans}} function call.
#'
#' @param new_diagrams a list of persistence diagrams, as the output of a TDA calculation like \code{\link[TDA]{ripsDiag}} or \code{\link{diagram_to_df}}.
#' @param clustering the output of a \code{link{diagram_kkmeans}} function call.
#'
#' @return a vector of the predicted cluster labels for the new diagrams.
#' @export
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @seealso \code{\link{diagram_kkmeans}} for clustering persistence diagrams.
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
#' # calculate kmeans clusters with centers = 2 in dimension 1 with sigma = t = 2
#' clusters <- diagram_kkmeans(diagrams = g,centers = 2,dim = 1,t = 2,sigma = 2)
#' 
#' # predict the nearest cluster for all diagrams in g
#' nearest_clusters <- diagram_nearest_clusters(new_diagrams = g,clustering = clusters)

diagram_nearest_clusters <- function(new_diagrams,clustering){
  
  # set internal variables to NULL to avoid build issues
  X <- NULL
  
  # error check diagrams argument
  check_param("new_diagrams",new_diagrams,numeric = F,multiple = T)
  new_diagrams <- all_diagrams(diagram_groups = list(new_diagrams),inference = "independence")[[1]]
  
  # error check clustering argument
  if(is.null(clustering))
  {
    stop("clustering must not be NULL.")
  }
  if(class(clustering) != "diagram_kkmeans")
  {
    stop("clustering object must be the output of a diagram_kkmeans function call.")
  }
  
  # compute cross Gram matrix
  K = gram_matrix(diagrams = new_diagrams,other_diagrams = clustering$diagrams,dim = clustering$dim,sigma = clustering$sigma,t = clustering$t)
  
  # return predicted class for each new diagram
  predicted_clusters <- unlist(lapply(X = 1:length(new_diagrams),FUN = function(X){
    
    return(clustering$clustering@.Data[[as.numeric(which.max(K[X,]))]])
    
  }))
  
  return(predicted_clusters)
  
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
#' `...` are additional parameters to the kpca kernlab function.
#'
#' @param diagrams a list of persistence diagrams, as the output of a TDA calculation.
#' @param dim the homological dimension in which the distance is to be computed.
#' @param t the positive scale for the persistence Fisher kernel, default 1.
#' @param sigma a positive number representing the bandwith for the Fisher information metric, default 1
#' @param features number of features (principal components) to return, default 1.
#' @param ... additional parameters.
#'
#' @return a list containing 
#' 
#' \describe{
#' 
#' \item{pca}{the output of \code{\link[kernlab]{kpca}} on the Gram matrix, an S4 object containing the slots pcv (a matrix containing the principal component vectors (column wise)), eig (the corresponding eigenvalues), rotated (the original data projected (rotated) on the principal components) and xmatrix (the original data matrix).}
#' 
#' \item{diagrams}{the input `diagrams` argument.}
#' 
#' \item{t}{the input `t` argument.}
#' 
#' \item{sigma}{the input `sigma` argument.}
#' 
#' \item{dim}{the input `dim` argument.}
#' 
#' }
#' 
#' @export
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @importFrom kernlab kpca
#' @seealso \code{\link{predict_diagram_kpca}} for predicting embedding coordinates of new diagrams.
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
  check_param("diagrams",diagrams,numeric = F,multiple = T)
  diagrams <- all_diagrams(diagram_groups = list(diagrams),inference = "independence")[[1]]
  
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
#' The `new_diagrams` parameter should be a list of persistence diagrams computed from a TDA calculation like ripsDiag
#' or \code{\link{diagram_to_df}}.
#' The `embedding` parameter is the diagram_kpca embedding object to be used for embedding
#' the new diagrams.
#'
#' @param new_diagrams a list of persistence diagrams, as the output of a TDA calculation like \code{\link[TDA]{ripsDiag}} or \code{\link{diagram_to_df}}.
#' @param embedding the output to a diagram_kpca function call.
#'
#' @return the data projection (rotation), stored as a numeric matrix. Each row corresponds to the same-index diagram in `new_diagrams`.
#' @export
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom parallel makeCluster stopCluster clusterExport clusterEvalQ
#' @importFrom parallelly availableCores
#' @importFrom doParallel registerDoParallel
#' @seealso \code{\link{diagram_kpca}} for embedding persistence diagrams into a low-dimensional space.
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
  check_param("new_diagrams",new_diagrams,numeric = F,multiple = T)
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
  
  # compute cross-Gram matrix and scale twice
  K <- gram_matrix(diagrams = new_diagrams,other_diagrams = embedding$diagrams,dim = embedding$dim,sigma = embedding$sigma,t = embedding$t)
  K <- scale(K,center = T,scale = F)
  K <- t(scale(t(K),center = T,scale = F))
  
  # project the new diagrams into the embedding space
  return(K %*% embedding$pca@pcv)
  
}

#### KERNEL SVM ####
#' Fit a support vector machine model where the training set is a list of persistence diagrams
#'
#' Returns the output of ksvm on the Gram matrix of a group of persistence diagrams
#' in a particular dimension. Cross validation is carried out in parallel, using a trick
#' noted in \url{https://doi.org/10.1007/s41468-017-0008-7}.
#'
#' The `diagrams` parameter should be a list of persistence diagrams computed from a TDA calculation like \code{\link[TDA]{ripsDiag}} or \code{\link{diagram_to_df}}.
#' The `dim` parameter should be a positive finite integer.
#' The `sigma` and `t` parameters are the positive bandwith for the Fisher information metric and
#' the positive scale for the persistence Fisher kernel respectively.
#' `type`, `C`, `nu`, `epsilon`, `prob.model`, `class.weights`, `cross`, `fit`, `cache`, `tol`, and `shrinking` 
#' are additional parameters to the ksvm kernlab function.
#'
#' @param diagrams a list of persistence diagrams, as the output of a TDA calculation.
#' @param cv a positive number at most the length of `diagrams` which determines the number of cross validation splits to be performed (default 1, aka no cross-validation).
#' @param dim the homological dimension in which the distance is to be computed.
#' @param t the positive scale for the persistence Fisher kernel, default 1.
#' @param sigma a positive number representing the bandwith for the Fisher information metric, default 1
#' @param y a response vector with one label for each persistence diagram.
#' @param type type of task to be performed.
#' @param C cost of contraints violation (default 1) this is the 'C'-constant of the regularization term in the Lagrange formulation.
#' @param nu parameter needed for nu-svc, one-svc and nu-svr. The `nu` parameter sets the upper bound on the training error and the lower bound on the fraction of data points to become Support Vector (default 0.2).
#' @param epsilon epsilon in the insensitive-loss function used for eps-svr, nu-svr and eps-bsvm (default 0.1).
#' @param fit indicates whether the fitted values should be computed and included in the model or not (default TRUE).
#' @param prob.model if set to TRUE builds a model for calculating class probabilities or in case of regression, calculates the scaling parameter of the Laplacian distribution fitted on the residuals. Fitting is done on output data created by performing a 3-fold cross-validation on the training data. For details see references (default FALSE).
#' @param class.weights a named vector of weights for the different classes, used for asymmetric class sizes. Not all factor levels have to be supplied (default weight: 1). All components have to be named.
#' @param cache cache memory in MB (default 40).
#' @param tol tolerance of termination criteria (default 0.001).
#' @param shrinking option whether to use the shrinking-heuristics (default TRUE).
#' @return a list containing 
#' 
#' \describe{
#' 
#' \item{models}{the cross-validation results - a matrix storing the parameters for each model in the tuning grid and its mean cross-validation error over all split.}
#' 
#' \item{best_model}{the output of \code{\link[kernlab]{ksvm}} run on the whole dataset with the optimal model parameters found during cross-validation. See the help page for \code{\link[kernlab]{ksvm}} for more details about this object.}
#' 
#' \item{diagrams}{the diagrams which were support vectors in the best_model. These are used for downstream prediction.}
#' 
#' \item{dim}{the input `dim` argument.}
#' 
#' \item{t}{the input `t` argument.}
#' 
#' \item{sigma}{the input `sigma` argument.}
#' 
#' }
#' 
#' @export
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @seealso \code{\link{predict_diagram_ksvm}} for predicting labels of new diagrams.
#' @importFrom kernlab ksvm
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel makeCluster stopCluster clusterExport clusterEvalQ
#' @importFrom parallelly availableCores
#' @importFrom doParallel registerDoParallel
#' @importFrom iterators iter
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
#' # create random response vector
#' y <- stats::runif(10,min = 0,max = 1)
#'
#' # calculate models over a grid with 5-fold CV
#' model_svm <- diagram_ksvm(diagrams = g,cv = 5,dim = c(1,2),y = y,sigma = c(1,0.1))

diagram_ksvm <- function(diagrams,cv = 1,dim,t = 1,sigma = 1,y,type = NULL,C = 1,nu = 0.2,epsilon = 0.1,prob.model = F,class.weights = NULL,fit = T,cache = 40,tol = 0.001,shrinking = T){
  
  # set internal variables to NULL to avoid build issues
  r <- NULL
  s <- NULL
  X <- NULL
  
  # error check diagrams argument
  check_param("diagrams",diagrams,numeric = F,multiple = T)
  diagrams <- all_diagrams(diagram_groups = list(diagrams),inference = "independence")[[1]]
  
  # error check cv parameters
  check_param("cv",cv,whole_numbers = T)
  if(cv < 1 | cv > length(diagrams))
  {
    stop("cv must be at least one and at most the number of diagrams.")
  }
  
  # basic error check for y
  if(length(y) != length(diagrams))
  {
    stop("y must be a vector with the same number of elements as diagrams.")
  }
  
  # set defaults for type argument
  if(is.null(type))
  {
    if(is.factor(y))
    {
      type <- "C-svc"
    }else
    {
      type <- "eps-svr"
    }
  }
  
  # expand parameter grid
  params <- expand.grid(dim = dim,t = t,sigma = sigma,C = C,nu = nu,epsilon = epsilon)
  
  # make vector of row memberships
  v <- rep(1:cv,floor(length(diagrams)/cv))
  if(length(v) != length(diagrams))
  {
    v <- c(v,1:(length(diagrams) %% floor(length(diagrams)/cv)))
  }
  v <- sample(v,size = length(v),replace = F)
  
  # split diagrams by membership vector
  diagrams_split <- lapply(X = 1:cv,FUN = function(X){
    
    return(which(v == X))
    
  })
  
  # for each pair of sigma and dim values compute the Fisher distance matrix to avoid recomputing Gram matrices
  dim_and_sigma <- expand.grid(dim = dim,sigma = sigma)
  distance_matrices <- lapply(X = 1:nrow(dim_and_sigma),FUN = function(X){
    
    return(distance_matrix(diagrams = diagrams,dim = dim_and_sigma[X,1],distance = "fisher",sigma = dim_and_sigma[X,2]))
    
  })
  names(distance_matrices) <- paste(dim_and_sigma$dim,dim_and_sigma$sigma,sep = "_")
  
  # set up for parallel computation
  num_workers <- parallelly::availableCores(omit = 1)
  cl <- parallel::makeCluster(num_workers)
  doParallel::registerDoParallel(cl)
  parallel::clusterEvalQ(cl,c(library(clue),library(rdist)))
  parallel::clusterExport(cl,c("y","predict_diagram_ksvm","all_diagrams","check_diagram","gram_matrix","diagram_distance","diagram_kernel","check_param"),envir = environment())
  force(diagrams)
  force(distance_matrices)
  force(diagrams_split)
  force(type)
  force(prob.model)
  force(class.weights)
  force(fit)
  force(cache)
  force(tol)
  force(shrinking)
  
  # for each model (combination of parameters), train the model on each subset and get the
  # prediction error on the hold out set
  if(cv < nrow(params))
  {
    model_errors <- foreach(r = iterators::iter(params,by = "row"),.combine = c) %dopar%
      {
        K <- exp(-1*r[[2]]*distance_matrices[[paste0(r[[1]],"_",r[[3]])]])
        return(mean(foreach(s = 1:cv,.combine = c) %do%
                      {
                        if(cv > 1)
                        {
                          training_indices <- unlist(diagrams_split[setdiff(1:cv,s)])
                          training_indices <- training_indices[order(training_indices)]
                          test_indices <- setdiff(1:length(diagrams),training_indices)
                          m <- length(test_indices)
                          K_subset <- K[training_indices,training_indices]
                          class(K_subset) <- "kernelMatrix"
                          model <- list(model = kernlab::ksvm(x = K_subset,y = y[training_indices],type = type,C = r[[4]],nu = r[[5]],epsilon = r[[6]],prob.model = prob.model,class.weights = class.weights,fit = fit,cache = cache,tol = tol,shrinking = shrinking),
                                        dim = r[[1]],
                                        sigma = r[[3]],
                                        t = r[[2]],
                                        CV = NULL)
                          model$diagrams <- diagrams[model$model@SVindex]
                          model_list <- list(cv_results = NULL,best_model = model)
                          class(model_list) <- "diagram_ksvm"
                          predictions <- predict_diagram_ksvm(new_diagrams = diagrams[test_indices],model = model)
                          
                          if(type == "C-svc" || type == "nu-svc" || type == "spoc-svc" || type == "kbb-svc" || type == "C-bsvc")
                          {
                            tab <- as.matrix(table(y[test_indices],as.integer(predictions)))
                            error  <- 1 - sum(diag(tab))/sum(as.numeric(tab))
                          }
                          if(type == "one-svc")
                          {
                            error <- sum(!predictions)/m
                          }
                          if(type == "eps-svr" || type == "nu-svr" || type == "eps-bsvr")
                          {
                            error <- drop(crossprod(predictions - y[setdiff(1:length(diagrams),training_indices)])/m)
                          }
                          
                          return(error)
                        }else
                        {
                          class(K) <- "kernelMatrix"
                          model <- kernlab::ksvm(x = K,y = y,type = type,C = r[[4]],nu = r[[5]],epsilon = r[[6]],prob.model = prob.model,class.weights = class.weights,fit = fit,cache = cache,tol = tol,shrinking = shrinking)
                          return(model@error)
                          
                        }
                        
                      }))
      }
  }else
  {
    model_errors <- foreach(r = iterators::iter(params,by = "row"),.combine = c) %do%
      {
        K <- exp(-1*r[[2]]*distance_matrices[[paste0(r[[1]],"_",r[[3]])]])
        return(mean(foreach(s = 1:cv,.combine = c) %dopar%
                      {
                        if(cv > 1)
                        {
                          training_indices <- unlist(diagrams_split[setdiff(1:cv,s)])
                          training_indices <- training_indices[order(training_indices)]
                          test_indices <- setdiff(1:length(diagrams),training_indices)
                          m <- length(test_indices)
                          K_subset <- K[training_indices,training_indices]
                          class(K_subset) <- "kernelMatrix"
                          model <- list(model = kernlab::ksvm(x = K_subset,y = y[training_indices],type = type,C = r[[4]],nu = r[[5]],epsilon = r[[6]],prob.model = prob.model,class.weights = class.weights,fit = fit,cache = cache,tol = tol,shrinking = shrinking),
                                        dim = r[[1]],
                                        sigma = r[[3]],
                                        t = r[[2]],
                                        CV = NULL)
                          model$diagrams <- diagrams[model$model@SVindex]
                          model_list <- list(cv_results = NULL,best_model = model)
                          class(model_list) <- "diagram_ksvm"
                          predictions <- predict_diagram_ksvm(new_diagrams = diagrams[test_indices],model = model_list)
                          
                          if(type == "C-svc" || type == "nu-svc" || type == "spoc-svc" || type == "kbb-svc" || type == "C-bsvc")
                          {
                            tab <- as.matrix(table(y[test_indices],as.integer(predictions)))
                            error  <- 1 - sum(diag(tab))/sum(as.numeric(tab))
                          }
                          if(type == "one-svc")
                          {
                            error <- sum(!predictions)/m
                          }
                          if(type == "eps-svr" || type == "nu-svr" || type == "eps-bsvr")
                          {
                            error <- drop(crossprod(predictions - y[setdiff(1:length(diagrams),training_indices)])/m)
                          }
                          
                          return(error)
                        }else
                        {
                          class(K) <- "kernelMatrix"
                          model <- kernlab::ksvm(x = K,y = y,type = type,C = r[[4]],nu = r[[5]],epsilon = r[[6]],prob.model = prob.model,class.weights = class.weights,fit = fit,cache = cache,tol = tol,shrinking = shrinking)
                          return(model@error)
                          
                        }
                        
                      }))
      }
  }

  parallel::stopCluster(cl)
  
  # get best parameters
  params$error <- model_errors
  best_params <- params[which.min(model_errors),1:(ncol(params) - 1)]
  
  # fit full model using those parameters
  K <- exp(-1*best_params[[2]]*distance_matrices[[paste0(best_params[[1]],"_",best_params[[3]])]])
  class(K) <- "kernelMatrix"
  
  # return kernlab calculation
  best_model <- list(model = kernlab::ksvm(x = K,y = y,type = type,C = best_params[[4]],nu = best_params[[5]],epsilon = best_params[[6]],prob.model = prob.model,class.weights = class.weights,fit = fit,cache = cache,tol = tol,shrinking = shrinking),
                   dim = best_params[[1]],
                   sigma = best_params[[3]],
                   t = best_params[[2]],
                   CV = params)

  best_model$diagrams <- diagrams[best_model$model@SVindex]

  ret_list <- list(cv_results = params,best_model = best_model)
  class(ret_list) <- "diagram_ksvm"
  
  return(ret_list)
  
}

#### PREDICT KERNEL SVM ####
#' Predict the response to a new group of persistence diagrams from a computed diagram_ksvm model
#'
#' Returns the predicted response vector of the model on the new diagrams.
#'
#' The `new_diagrams` parameter should be a list of persistence diagrams computed from a TDA calculation like \code{\link[TDA]{ripsDiag}} or \code{\link{diagram_to_df}}.
#' The `model` parameter should be the output from a diagram_ksvm function call.
#'
#' @param new_diagrams a list of new persistence diagrams, as the output of a TDA calculation.
#' @param model the diagram_ksvm model to be used for prediction.
#' @return a vector containing the output of \code{\link[kernlab]{predict.ksvm}} on the cross Gram matrix of the new diagrams and the support vector diagrams stored in the model.
#' @export
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @seealso \code{\link{diagram_ksvm}} for training a SVM model on a training set of persistence diagrams.
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom parallel makeCluster stopCluster clusterExport clusterEvalQ
#' @importFrom parallelly availableCores
#' @importFrom doParallel registerDoParallel
#' @importFrom kernlab predict as.kernelMatrix
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
#' # create random response vector
#' y <- stats::runif(10,min = 0,max = 1)
#'
#' # calculate model in dimension 1
#' model_svm <- diagram_ksvm(diagrams = g,dim = 1,y = y)
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
#' # predict responses
#' y_pred <- predict_diagram_ksvm(new_diagrams = g_new,model = model_svm)

predict_diagram_ksvm <- function(new_diagrams,model){
  
  # set internal variables to NULL to avoid build issues
  r <- NULL
  X <- NULL
  
  # error check new_diagrams argument
  check_param("new_diagrams",new_diagrams,numeric = F,multiple = T)
  new_diagrams <- all_diagrams(diagram_groups = list(new_diagrams),inference = "independence")[[1]]
  
  # error check model argument
  if(is.null(model))
  {
    stop("model must be supplied.")
  }
  if(class(model) != "diagram_ksvm")
  {
    stop("model must be the output of a diagram_ksvm function call.")
  }
  
  # compute kernel matrix, storing the value of each kernel computation between the new diagrams and the old ones
  K <- gram_matrix(diagrams = new_diagrams,other_diagrams = model$best_model$diagrams,dim = model$best_model$dim,sigma = model$best_model$sigma,t = model$best_model$t)

  return(kernlab::predict(object = model$best_model$model,kernlab::as.kernelMatrix(K)))
  
}

