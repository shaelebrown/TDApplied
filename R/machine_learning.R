
#### Multidimensional scaling ####
#' Dimension reduction of a group of persistence diagrams via metric multidimensional scaling.
#'
#' Projects a group of persistence diagrams into a low-dimensional embedding space via metric multidimensional
#' scaling. Such a projection can be used for visualization of data, or a static analysis of the embedding
#' dimensions.
#' 
#' Returns the output of \code{\link[stats]{cmdscale}} on the desired distance matrix of a group of persistence diagrams
#' in a particular dimension. If `distance` is "fisher" then `sigma` must not be NULL.
#'
#' @param diagrams a list of n>=2 persistence diagrams which are either the output of a TDA/TDAstats calculation like \code{\link[TDA]{ripsDiag}}/\code{\link[TDAstats]{calculate_homology}}, or \code{\link{diagram_to_df}}.
#' @param k the dimension of the space which the data are to be represented in; must be in {1,2,...,n-1}.
#' @param distance a string representing the desired distance metric to be used, either 'wasserstein' (default) or 'fisher'.
#' @param dim the non-negative integer homological dimension in which the distance is to be computed, default 0.
#' @param p a positive number representing the wasserstein power, a number at least 1 (infinity for the bottleneck distance), default 2.
#' @param sigma a positive number representing the bandwidth for the Fisher information metric, default NULL.
#' @param eig a boolean indicating whether the eigenvalues should be returned.
#' @param add a boolean indicating if an additive constant c* should be computed, and added to the non-diagonal dissimilarities such that the modified dissimilarities are Euclidean.
#' @param x.ret a boolean indicating whether the doubly centered symmetric distance matrix should be returned.
#' @param list. a boolean indicating if a list should be returned or just the n*k matrix.
#' @param num_workers the number of cores used for parallel computation, default is one less than the number of cores on the machine.
#'
#' @return the output of \code{\link[stats]{cmdscale}} on the diagram distance matrix. If `list.` is false (as per default),
#' a matrix with `k` columns whose rows give the coordinates of the points chosen to represent the dissimilarities.
#' 
#' Otherwise, a list containing the following components.
#' 
#' \describe{
#' 
#'  \item{points}{a matrix with `k` columns whose rows give the coordinates of the points chosen to represent the dissimilarities.}
#' 
#'  \item{eig}{the \eqn{n} eigenvalues computed during the scaling process if `eig` is true.}
#'  
#'  \item{x}{the doubly centered distance matrix if `x.ret` is true.}
#'  
#'  \item{ac}{the additive constant \eqn{c*}, 0 if `add` = FALSE.}
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
#' @references 
#' Cox M and Cox F (2008). "Multidimensional Scaling." \doi{10.1007/978-3-540-33037-0_14}.
#' 
#' @examples
#'
#' # create two diagrams
#' D1 <- TDAstats::calculate_homology(TDA::circleUnif(n = 10,r = 1),
#'                                    dim = 0,threshold = 2)
#' D2 <- TDAstats::calculate_homology(TDA::circleUnif(n = 10,r = 1),
#'                                    dim = 0,threshold = 2)
#' g <- list(D1,D2)
#' 
#' # calculate their 1D MDS embedding in dimension 0 with the bottleneck distance
#' mds <- diagram_mds(diagrams = g,k = 1,dim = 0,p = Inf,num_workers = 2)

diagram_mds <- function(diagrams,k = 2,distance = "wasserstein",dim = 0,p = 2,sigma = NULL,eig = FALSE,add = FALSE,x.ret = FALSE,list. = eig || add || x.ret,num_workers = parallelly::availableCores(omit = 1)){
  
  # error check diagrams argument
  check_param("diagrams",diagrams,numeric = F,multiple = T)
  diagrams <- all_diagrams(diagram_groups = list(diagrams),inference = "independence")[[1]]

  # compute distance matrix
  d <- distance_matrix(diagrams = diagrams,dim = dim,distance = distance,p = p,sigma = sigma,num_workers = num_workers)

  # return metric multidimensional scaling with d as input
  return(stats::cmdscale(d = d,k = k,eig = eig,add = add,x.ret = x.ret,list. = list.))
  
}

#### KERNEL KMEANS ####
#' Cluster a group of persistence diagrams using kernel k-means.
#' 
#' Finds latent cluster labels for a group of persistence diagrams, using a kernelized version
#' of the popular k-means algorithm. An optimal number of clusters may be determined by analyzing
#' the withinss field of the clustering object over several values of k.
#'
#' Returns the output of \code{\link[kernlab]{kkmeans}} on the desired Gram matrix of a group of persistence diagrams
#' in a particular dimension. The additional list elements stored in the output are needed
#' to estimate cluster labels for new persistence diagrams in the `predict_diagram_kkmeans`
#' function.
#'
#' @param diagrams a list of persistence diagrams which are either the output of a TDA/TDAstats calculation like \code{\link[TDA]{ripsDiag}}/\code{\link[TDAstats]{calculate_homology}}, or the \code{\link{diagram_to_df}} function.
#' @param dim the non-negative integer homological dimension in which the distance is to be computed, default 0.
#' @param t a positive number representing the scale for the persistence Fisher kernel, default 1.
#' @param sigma a positive number representing the bandwidth for the Fisher information metric, default 1
#' @param centers number of clusters to initialize, no more than the number of diagrams although smaller values are recommended.
#' @param num_workers the number of cores used for parallel computation, default is one less than the number of cores on the machine.
#' @param ... additional parameters for the \code{\link[kernlab]{kkmeans}} kernlab function.
#'
#' @return a 'diagram_kkmeans' S3 object containing the output of \code{\link[kernlab]{kkmeans}} on the Gram matrix, i.e. a list containing the elements
#' 
#' \describe{
#' 
#' \item{clustering}{an S4 object of class specc, the output of a \code{\link[kernlab]{kkmeans}} function call. The `.Data` slot of this object contains cluster memberships, `withinss` contains the within-cluster sum of squares for each cluster, etc.}
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
#' @seealso \code{\link{predict_diagram_kkmeans}} for predicting cluster labels of new diagrams.
#' @references 
#' Dhillon, I and Guan, Y and Kulis, B (2004). "A Unified View of Kernel k-means , Spectral Clustering and Graph Cuts." \url{https://people.bu.edu/bkulis/pubs/spectral_techreport.pdf}.
#' @examples
#'
#' # create two diagrams
#' D1 <- TDAstats::calculate_homology(TDA::circleUnif(n = 10,r = 1),
#'                                    dim = 0,threshold = 2)
#' D2 <- TDAstats::calculate_homology(TDA::circleUnif(n = 10,r = 1),
#'                                    dim = 0,threshold = 2)
#' g <- list(D1,D1,D2,D2)
#' 
#' # calculate kmeans clusters with centers = 2, and sigma = t = 2 in dimension 0
#' clust <- diagram_kkmeans(diagrams = g,centers = 2,dim = 0,t = 2,sigma = 2,num_workers = 2)

diagram_kkmeans <- function(diagrams,centers,dim = 0,t = 1,sigma = 1,num_workers = parallelly::availableCores(omit = 1),...){
  
  # error check arguments
  check_param("diagrams",diagrams,numeric = F,multiple = T)
  diagrams <- all_diagrams(diagram_groups = list(diagrams),inference = "independence")[[1]]
  check_param("centers",centers,whole_numbers = T,at_least_one = T)
  
  # make sure there aren't more centers than data point - this tends to lead to kernlab errors
  if(centers > length(diagrams))
  {
    stop("centers must be at most the number of diagrams, although we recommend choosing smaller values.")
  }
  
  # compute Gram matrix
  K <- gram_matrix(diagrams = diagrams,dim = dim,sigma = sigma,t = t,num_workers = num_workers)
  
  # return kernlab calculation, saving in a list of class diagram_kkmeans for later interface with prediction calculations
  max_iters <- 20
  iter <- 1
  while(T)
  {
    # there are some kernlab errors that can't always be avoided and caused by randomness
    # when certain clusters are empty or unclosed connections
    # so we catch those errors and rerun if necessary, up to max_iters many times
    # if a different error occurs or we rerun max_iters many times then stop with the error.
    tryCatch(expr = {clustering <- kernlab::kkmeans(x = K,centers = centers,...)},
             error = function(e){
               
               if(grepl(pattern = "sum\\(abs\\(dc\\)\\)",x = e) == F & grepl(pattern = "\'x\' must be an array of at least two dimensions",x = e) == F)
               {
                 stop(e)
               }
               
             },
             warning = function(w){
               
               if(grepl(pattern = "closing unused connection",w) == F)
               {
                 message(w)
               }
               
             })
    iter <- iter + 1
    if(exists("clustering"))
    {
      break
    }
    if(iter > max_iters)
    {
      stop("One or more clusters was empty - try decreasing the number of clusters.")
    }
    
  }
  ret_list <- list(clustering = clustering,diagrams = diagrams,dim = dim,sigma = sigma,t = t)
  
  class(ret_list) <- "diagram_kkmeans"
  
  return(ret_list)
  
}

#### PREDICT KERNEL KMEANS ####
#' Predict the cluster labels for new persistence diagrams using a pre-computed clustering.
#'
#' Returns the nearest (highest kernel value) \code{\link[kernlab]{kkmeans}} cluster center label for new persistence diagrams.
#' This allows for reusing old cluster models for new tasks, or to perform cross validation.
#'
#' @param new_diagrams a list of persistence diagrams which are either the output of a TDA/TDAstats calculation like \code{\link[TDA]{ripsDiag}}/\code{\link[TDAstats]{calculate_homology}}, or \code{\link{diagram_to_df}}.
#' @param clustering the output of a \code{\link{diagram_kkmeans}} function call.
#' @param num_workers the number of cores used for parallel computation, default is one less than the number of cores on the machine.
#'
#' @return a vector of the predicted cluster labels for the new diagrams.
#' @importFrom methods is
#' @export
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @seealso \code{\link{diagram_kkmeans}} for clustering persistence diagrams.
#' @examples
#'
#' # create two diagrams
#' D1 <- TDAstats::calculate_homology(TDA::circleUnif(n = 10,r = 1),
#'                                    dim = 0,threshold = 2)
#' D2 <- TDAstats::calculate_homology(TDA::circleUnif(n = 10,r = 1),
#'                                    dim = 0,threshold = 2)
#' g <- list(D1,D1,D2,D2)
#' 
#' # calculate kmeans clusters with centers = 2, and sigma = t = 2 in dimension 0
#' clust <- diagram_kkmeans(diagrams = g,centers = 2,dim = 0,t = 2,sigma = 2,num_workers = 2)
#' 
#' # create two new diagrams
#' D4 <- TDAstats::calculate_homology(TDA::circleUnif(n = 10,r = 1),
#'                                    dim = 0,threshold = 2)
#' D5 <- TDAstats::calculate_homology(TDA::circleUnif(n = 10,r = 1),
#'                                    dim = 0,threshold = 2)
#' g_new <- list(D4,D5)
#' 
#' # predict cluster labels
#' predict_diagram_kkmeans(new_diagrams = g_new,clustering = clust,num_workers = 2)

predict_diagram_kkmeans <- function(new_diagrams,clustering,num_workers = parallelly::availableCores(omit = 1)){
  
  # set internal variables to NULL to avoid build issues
  X <- NULL
  
  # error check diagrams argument
  check_param("new_diagrams",new_diagrams,numeric = F,multiple = T,min_length = 1)
  new_diagrams <- all_diagrams(diagram_groups = list(new_diagrams),inference = "independence")[[1]]
  
  # error check clustering argument
  if(is.null(clustering))
  {
    stop("clustering must not be NULL.")
  }
  if(!methods::is(clustering,"diagram_kkmeans"))
  {
    stop("clustering object must be the output of a diagram_kkmeans function call.")
  }
  
  # compute cross Gram matrix
  K = gram_matrix(diagrams = new_diagrams,other_diagrams = clustering$diagrams,dim = clustering$dim,sigma = clustering$sigma,t = clustering$t,num_workers = num_workers)
  
  # return predicted class for each new diagram
  predicted_clusters <- unlist(lapply(X = 1:length(new_diagrams),FUN = function(X){
    
    return(clustering$clustering@.Data[[as.numeric(which.max(K[X,]))]])
    
  }))
  
  return(predicted_clusters)
  
}

#### KERNEL PCA ####
#' Calculate the kernel PCA embedding of a group of persistence diagrams.
#'
#' Project a group of persistence diagrams into a low-dimensional embedding space using
#' a kernelized version of the popular PCA algorithm. 
#'
#' Returns the output of kernlab's \code{\link[kernlab]{kpca}} function on the desired Gram matrix of a group of persistence diagrams
#' in a particular dimension. The prediction function \code{\link{predict_diagram_kpca}} can be used to 
#' project new persistence diagrams using an old embedding, and this could be one practical
#' advantage of using \code{\link{diagram_kpca}} over \code{\link{diagram_mds}}. The embedding coordinates can also
#' be used for further analysis, or simply as a data visualization tool for persistence diagrams.
#'
#' @param diagrams a list of persistence diagrams which are either the output of a TDA/TDAstats calculation like \code{\link[TDA]{ripsDiag}}/\code{\link[TDAstats]{calculate_homology}}, or \code{\link{diagram_to_df}}.
#' @param dim the non-negative integer homological dimension in which the distance is to be computed, default 0.
#' @param t a positive number representing the scale for the persistence Fisher kernel, default 1.
#' @param sigma a positive number representing the bandwidth for the Fisher information metric, default 1
#' @param features number of features (principal components) to return, default 1.
#' @param num_workers the number of cores used for parallel computation, default is one less than the number of cores on the machine.
#' @param th the threshold value under which principal components are ignored (default 0.0001).
#'
#' @return a list containing the elements
#' 
#' \describe{
#' 
#' \item{pca}{the output of kernlab's \code{\link[kernlab]{kpca}} function on the Gram matrix: an S4 object containing the slots `pcv` (a matrix containing the principal component vectors (column wise)), `eig` (the corresponding eigenvalues), `rotated` (the original data projected (rotated) on the principal components) and `xmatrix` (the original data matrix).}
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
#' @references 
#' Scholkopf, B and Smola, A and Muller, K (1998). "Nonlinear Component Analysis as a Kernel Eigenvalue Problem." \url{https://www.mlpack.org/papers/kpca.pdf}.
#' @examples
#'
#' # create six diagrams
#' D1 <- TDAstats::calculate_homology(TDA::circleUnif(n = 50,r = 1),
#'                                    dim = 1,threshold = 2)
#' D2 <- TDAstats::calculate_homology(TDA::sphereUnif(n = 50,d = 2,r = 1),
#'                                    dim = 1,threshold = 2)
#' D3 <- TDAstats::calculate_homology(TDA::torusUnif(n = 50,a = 0.25,c = 0.75),
#'                                    dim = 1,threshold = 2)
#' D4 <- TDAstats::calculate_homology(TDA::circleUnif(n = 50,r = 1),
#'                                    dim = 1,threshold = 2)
#' D5 <- TDAstats::calculate_homology(TDA::sphereUnif(n = 50,d = 2,r = 1),
#'                                    dim = 1,threshold = 2)
#' D6 <- TDAstats::calculate_homology(TDA::torusUnif(n = 50,a = 0.25,c = 0.75),
#'                                    dim = 1,threshold = 2)
#' g <- list(D1,D2,D3,D4,D5,D6)
#' 
#' # calculate their 2D PCA embedding with sigma = t = 2 in dimension 1
#' pca <- diagram_kpca(diagrams = g,dim = 1,t = 2,sigma = 2,features = 2,num_workers = 2)

diagram_kpca <- function(diagrams,dim = 0,t = 1,sigma = 1,features = 1,num_workers = parallelly::availableCores(omit = 1),th = 1e-4){
  
  # error check diagrams argument
  check_param("diagrams",diagrams,numeric = F,multiple = T)
  diagrams <- all_diagrams(diagram_groups = list(diagrams),inference = "independence")[[1]]
  
  # compute Gram matrix
  K <- gram_matrix(diagrams = diagrams,t = t,sigma = sigma,dim = dim,num_workers = num_workers)
  
  # return kernlab computation, ignore unhelpful kernlab warnings
  tryCatch(expr = {ret_list <- list(pca = kernlab::kpca(x = K,features = features),diagrams = diagrams,t = t,sigma = sigma,dim = dim)},
           error = function(e){stop(e)},
           warning = function(w){
             
             if(grepl(pattern = "closing unused connection",w) == F & grepl(pattern = "below threshold!",w) == F)
             {
               message(w)
             }
             
           },
           finally = {
             
             if(!exists(x = "ret_list"))
             {
               stop(paste0("The embedding could not be computed - try decreasing the th parameter (input value was ",th,")."))
             }
             
           })
  
  class(ret_list) <- "diagram_kpca"
  return(ret_list)
  
}

#### PREDICT WITH KERNEL PCA OBJECT ####
#' Project persistence diagrams into a low-dimensional space via a pre-computed kernel PCA embedding.
#'
#' Compute the location in low-dimensional space of each element of a list of new persistence diagrams using a
#' previously-computed kernel PCA embedding (from the \code{\link{diagram_kpca}} function).
#'
#' @param new_diagrams a list of persistence diagrams which are either the output of a TDA/TDAstats calculation like \code{\link[TDA]{ripsDiag}}/\code{\link[TDAstats]{calculate_homology}}, or \code{\link{diagram_to_df}}.
#' @param embedding the output of a \code{\link{diagram_kpca}} function call.
#' @param num_workers the number of cores used for parallel computation, default is one less than the number of cores on the machine.
#'
#' @return the data projection (rotation), stored as a numeric matrix. Each row corresponds to the same-index diagram in `new_diagrams`.
#' @export
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom parallel makeCluster stopCluster clusterExport clusterEvalQ
#' @importFrom parallelly availableCores
#' @importFrom doParallel registerDoParallel
#' @importFrom methods is
#' @seealso \code{\link{diagram_kpca}} for embedding persistence diagrams into a low-dimensional space.
#' @examples
#'
#' # create six diagrams
#' D1 <- TDAstats::calculate_homology(TDA::circleUnif(n = 50,r = 1),
#'                                    dim = 1,threshold = 2)
#' D2 <- TDAstats::calculate_homology(TDA::sphereUnif(n = 50,d = 2,r = 1),
#'                                    dim = 1,threshold = 2)
#' D3 <- TDAstats::calculate_homology(TDA::torusUnif(n = 50,a = 0.25,c = 0.75),
#'                                    dim = 1,threshold = 2)
#' D4 <- TDAstats::calculate_homology(TDA::circleUnif(n = 50,r = 1),
#'                                    dim = 1,threshold = 2)
#' D5 <- TDAstats::calculate_homology(TDA::sphereUnif(n = 50,d = 2,r = 1),
#'                                    dim = 1,threshold = 2)
#' D6 <- TDAstats::calculate_homology(TDA::torusUnif(n = 50,a = 0.25,c = 0.75),
#'                                    dim = 1,threshold = 2)
#' g <- list(D1,D2,D3,D4,D5,D6)
#' 
#' # calculate their 2D PCA embedding with sigma = t = 2 in dimension 0
#' pca <- diagram_kpca(diagrams = g,dim = 1,t = 2,sigma = 2,features = 2,num_workers = 2)
#' 
#' # project two new diagrams onto old model
#' D7 <- TDAstats::calculate_homology(TDA::circleUnif(n = 50,r = 1),
#'                                    dim = 0,threshold = 2)
#' D8 <- TDAstats::calculate_homology(TDA::circleUnif(n = 50,r = 1),
#'                                    dim = 0,threshold = 2)
#' g_new <- list(D4,D5)
#' 
#' # calculate new embedding coordinates
#' new_pca <- predict_diagram_kpca(new_diagrams = g_new,embedding = pca,num_workers = 2)

predict_diagram_kpca <- function(new_diagrams,embedding,num_workers = parallelly::availableCores(omit = 1)){
  
  # set internal variables to NULL to avoid build issues
  r <- NULL
  X <- NULL
  
  # error check new_diagrams argument
  check_param("new_diagrams",new_diagrams,numeric = F,multiple = T,min_length = 1)
  new_diagrams <- all_diagrams(diagram_groups = list(new_diagrams),inference = "independence")[[1]]
  
  # error check embedding argument
  if(is.null(embedding))
  {
    stop("embedding must be supplied.")
  }
  if(!methods::is(embedding,"diagram_kpca"))
  {
    stop("embedding must be the output of a diagram_kpca function call.")
  }
  
  # compute cross-Gram matrix and scale twice
  K <- gram_matrix(diagrams = new_diagrams,other_diagrams = embedding$diagrams,dim = embedding$dim,sigma = embedding$sigma,t = embedding$t,num_workers = num_workers)
  K <- scale(K,center = T,scale = F)
  K <- t(scale(t(K),center = T,scale = F))
  
  # project the new diagrams into the embedding space
  return(K %*% embedding$pca@pcv)
  
}

#### KERNEL SVM ####
#' Fit a support vector machine model where each training set instance is a persistence diagram.
#'
#' Returns the output of kernlab's \code{\link{ksvm}} function on the Gram matrix of the list of persistence diagrams
#' in a particular dimension.
#' 
#' Cross validation is carried out in parallel, using a trick
#' noted in \doi{10.1007/s41468-017-0008-7} - since the persistence Fisher kernel can be
#' written as \eqn{d_{PF}(D_1,D_2)=exp(t*d_{FIM}(D_1,D_2))=exp(d_{FIM}(D_1,D_2))^t}, we can
#' store the Fisher information metric distance matrix for each sigma value in the parameter grid to avoid
#' recomputing distances. Parallelization occurs either over CV folds or over the parameter combinations,
#' whichever has more elements. Note that the response parameter `y` must be a factor for classification - 
#' a character vector for instance will throw an error.
#'
#' @param diagrams a list of persistence diagrams which are either the output of a TDA/TDAstats calculation like \code{\link[TDA]{ripsDiag}}/\code{\link[TDAstats]{calculate_homology}}, or \code{\link{diagram_to_df}}.
#' @param cv a positive number at most the length of `diagrams` which determines the number of cross validation splits to be performed (default 1, aka no cross-validation).
#' @param dim a non-negative integer vector of homological dimensions in which the model is to be fit.
#' @param t a vector of positive numbers representing the grid of values for the scale of the persistence Fisher kernel, default 1.
#' @param sigma a vector of positive numbers representing the grid of values for the bandwidth of the Fisher information metric, default 1
#' @param y a response vector with one label for each persistence diagram. Must be either numeric or factor.
#' @param type a string representing the type of task to be performed.
#' @param C a number representing the cost of constraints violation (default 1) this is the 'C'-constant of the regularization term in the Lagrange formulation.
#' @param nu numeric parameter needed for nu-svc, one-svc and nu-svr. The `nu` parameter sets the upper bound on the training error and the lower bound on the fraction of data points to become Support Vector (default 0.2).
#' @param epsilon epsilon in the insensitive-loss function used for eps-svr, nu-svr and eps-bsvm (default 0.1).
#' @param fit indicates whether the fitted values should be computed and included in the model or not (default TRUE).
#' @param prob.model if set to TRUE builds a model for calculating class probabilities or in case of regression, calculates the scaling parameter of the Laplacian distribution fitted on the residuals. Fitting is done on output data created by performing a 3-fold cross-validation on the training data. For details see references (default FALSE).
#' @param class.weights a named vector of weights for the different classes, used for asymmetric class sizes. Not all factor levels have to be supplied (default weight: 1). All components have to be named.
#' @param cache cache memory in MB (default 40).
#' @param tol tolerance of termination criteria (default 0.001).
#' @param shrinking option whether to use the shrinking-heuristics (default TRUE).
#' @param num_workers the number of cores used for parallel computation, default is one less the number of cores on the machine.
#' @return a list containing the elements
#' 
#' \describe{
#' 
#' \item{models}{the cross-validation results - a matrix storing the parameters for each model in the tuning grid and its mean cross-validation error over all splits.}
#' 
#' \item{best_model}{the output of \code{\link[kernlab]{ksvm}} run on the whole dataset with the optimal model parameters found during cross-validation. See the help page for \code{\link[kernlab]{ksvm}} for more details about this object.}
#' 
#' \item{diagrams}{the diagrams which were support vectors in the `best_model`. These are used for downstream prediction.}
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
#' @references 
#' Murphy, K. "Machine learning: a probabilistic perspective." MIT press (2012).
#' @examples
#'
#' # create four diagrams
#' D1 <- TDAstats::calculate_homology(TDA::circleUnif(n = 10,r = 1),
#'                                    dim = 0,threshold = 2)
#' D2 <- TDAstats::calculate_homology(TDA::circleUnif(n = 10,r = 1),
#'                                    dim = 0,threshold = 2)
#' D3 <- TDAstats::calculate_homology(TDA::circleUnif(n = 10,r = 1),
#'                                    dim = 0,threshold = 2)
#' D4 <- TDAstats::calculate_homology(TDA::circleUnif(n = 10,r = 1),
#'                                    dim = 0,threshold = 2)
#' g <- list(D1,D2,D3,D4)
#' 
#' # create response vector
#' y <- as.factor(c("circle","sphere","circle","sphere"))
#' 
#' # fit model without cross validation
#' model_svm <- diagram_ksvm(diagrams = g,cv = 1,dim = c(0),
#'                           y = y,sigma = c(1),t = c(1),
#'                           num_workers = 2)
                          
diagram_ksvm <- function(diagrams,cv = 1,dim,t = 1,sigma = 1,y,type = NULL,C = 1,nu = 0.2,epsilon = 0.1,prob.model = FALSE,class.weights = NULL,fit = TRUE,cache = 40,tol = 0.001,shrinking = TRUE,num_workers = parallelly::availableCores(omit = 1)){
  
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
  
  # basic error checks for y
  if(length(y) != length(diagrams))
  {
    stop("y must be a vector with the same number of elements as diagrams.")
  }
  if(class(y) %in% c("numeric","factor") == F)
  {
    stop("y should be either a numeric or factor vector.")
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
    
    return(distance_matrix(diagrams = diagrams,dim = dim_and_sigma[X,1],distance = "fisher",sigma = dim_and_sigma[X,2],num_workers = num_workers))
    
  })
  names(distance_matrices) <- paste(dim_and_sigma$dim,dim_and_sigma$sigma,sep = "_")
  
  # set up for parallel computation
  cl <- parallel::makeCluster(num_workers)
  doParallel::registerDoParallel(cl)
  parallel::clusterEvalQ(cl,c(library(clue),library(rdist)))
  parallel::clusterExport(cl,c("y","predict_diagram_ksvm","all_diagrams","check_diagram","gram_matrix","diagram_distance","diagram_kernel","check_param","diagram_to_df","%do%","%dopar%"),envir = environment())
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
    # parallelize over parameter combinations
    model_errors <- foreach::`%dopar%`(obj = foreach::foreach(r = iterators::iter(params,by = "row"),.combine = c),ex = 
      {
        # use precomputed distance matrices to calculate Gram matrix from parameters
        K <- exp(-1*r[[2]]*distance_matrices[[paste0(r[[1]],"_",r[[3]])]])
        return(mean(foreach::`%do%`(obj = foreach::foreach(s = 1:cv,.combine = c),ex =
                      {
                        # foreach but not in parallel for cv folds
                        if(cv > 1)
                        {
                          # split into training and test set based on cv parameter
                          training_indices <- unlist(diagrams_split[setdiff(1:cv,s)])
                          training_indices <- training_indices[order(training_indices)]
                          test_indices <- setdiff(1:length(diagrams),training_indices)
                          m <- length(test_indices)
                          K_subset <- K[training_indices,training_indices]
                          class(K_subset) <- "kernelMatrix"
                          # suppress unhelpful kernlab warnings when fitting model
                          tryCatch(expr = {model = kernlab::ksvm(x = K_subset,y = y[training_indices],type = type,C = r[[4]],nu = r[[5]],epsilon = r[[6]],prob.model = prob.model,class.weights = class.weights,fit = fit,cache = cache,tol = tol,shrinking = shrinking)},
                                   error = function(e){stop(e)},
                                   warning = function(w){
                                     
                                     if(grepl(pattern = "closing unused connection",w) == F)
                                     {
                                       message(w)
                                     }
                                     
                                   })
                          
                          return(model@error)
                        }else
                        {
                          class(K) <- "kernelMatrix"
                          # suppress unhelpful kernlab warnings when fitting model
                          tryCatch(expr = {model <- kernlab::ksvm(x = K,y = y,type = type,C = r[[4]],nu = r[[5]],epsilon = r[[6]],prob.model = prob.model,class.weights = class.weights,fit = fit,cache = cache,tol = tol,shrinking = shrinking)},
                                   error = function(e){stop(e)},
                                   warning = function(w){
                                     
                                     if(grepl(pattern = "closing unused connection",w) == F)
                                     {
                                       message(w)
                                     }
                                     
                                   })
                          return(model@error)
                          
                        }
                        
                      })))
      })
  }else
  {
    # parallelize over cv folds
    model_errors <- foreach::`%do%`(obj = foreach::foreach(r = iterators::iter(params,by = "row"),.combine = c),ex =
      {
        # use precomputed distance matrices to calculate Gram matrix from parameters
        K <- exp(-1*r[[2]]*distance_matrices[[paste0(r[[1]],"_",r[[3]])]])
        return(mean(foreach::`%dopar%`(obj = foreach::foreach(s = 1:cv,.combine = c),ex = 
                      {
                        if(cv > 1)
                        {
                          # split into training and test set based on cv parameter
                          training_indices <- unlist(diagrams_split[setdiff(1:cv,s)])
                          training_indices <- training_indices[order(training_indices)]
                          test_indices <- setdiff(1:length(diagrams),training_indices)
                          m <- length(test_indices)
                          K_subset <- K[training_indices,training_indices]
                          class(K_subset) <- "kernelMatrix"
                          # suppress unhelpful kernlab warnings when fitting model
                          tryCatch(expr = {model = kernlab::ksvm(x = K_subset,y = y[training_indices],type = type,C = r[[4]],nu = r[[5]],epsilon = r[[6]],prob.model = prob.model,class.weights = class.weights,fit = fit,cache = cache,tol = tol,shrinking = shrinking)},
                                   error = function(e){stop(e)},
                                   warning = function(w){
                                     
                                     if(grepl(pattern = "closing unused connection",w) == F)
                                     {
                                       message(w)
                                     }
                                     
                                   })
                          
                          return(model@error)
                        }else
                        {
                          class(K) <- "kernelMatrix"
                          # suppress unhelpful kernlab warnings when fitting model
                          tryCatch(expr = {model <- kernlab::ksvm(x = K,y = y,type = type,C = r[[4]],nu = r[[5]],epsilon = r[[6]],prob.model = prob.model,class.weights = class.weights,fit = fit,cache = cache,tol = tol,shrinking = shrinking)},
                                   error = function(e){stop(e)},
                                   warning = function(w){
                                     
                                     if(grepl(pattern = "closing unused connection",w) == F)
                                     {
                                       message(w)
                                     }
                                     
                                   })
                          return(model@error)
                          
                        }
                        
                      })))
      })
  }

  parallel::stopCluster(cl)
  
  # get best parameters
  params$error <- model_errors
  best_params <- params[which.min(model_errors),1:(ncol(params) - 1)]
  
  # fit full model using those parameters
  K <- exp(-1*best_params[[2]]*distance_matrices[[paste0(best_params[[1]],"_",best_params[[3]])]])
  class(K) <- "kernelMatrix"
  
  # return kernlab calculation, suppressing unhelpful warnings
  tryCatch(expr = {best_model <- list(model = kernlab::ksvm(x = K,y = y,type = type,C = best_params[[4]],nu = best_params[[5]],epsilon = best_params[[6]],prob.model = prob.model,class.weights = class.weights,fit = fit,cache = cache,tol = tol,shrinking = shrinking),
                                      dim = best_params[[1]],
                                      sigma = best_params[[3]],
                                      t = best_params[[2]],
                                      CV = params)},
           error = function(e){stop(e)},
           warning = function(w){
             
             if(grepl(pattern = "closing unused connection",w) == F)
             {
               message(w)
             }
             
           })
  
  best_model$diagrams <- diagrams[best_model$model@SVindex]

  ret_list <- list(cv_results = params,best_model = best_model)
  class(ret_list) <- "diagram_ksvm"
  
  return(ret_list)
  
}

#### PREDICT KERNEL SVM ####
#' Predict the outcome labels for a list of persistence diagrams using a pre-trained diagram ksvm model.
#'
#' Returns the predicted response vector of the model on the new diagrams.
#' 
#' This function is a wrapper of the kernlab \code{\link{predict}} function.
#'
#'
#' @param new_diagrams a list of persistence diagrams which are either the output of a TDA/TDAstats calculation like \code{\link[TDA]{ripsDiag}}/\code{\link[TDAstats]{calculate_homology}}, or \code{\link{diagram_to_df}}.
#' @param model the output of a \code{\link{diagram_ksvm}} function call.
#' @param num_workers the number of cores used for parallel computation, default is one less than the number of cores on the machine.
#' @return a vector containing the output of \code{\link[kernlab]{predict.ksvm}} on the cross Gram matrix of the new diagrams and the support vector diagrams stored in the model.
#' @export
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @seealso \code{\link{diagram_ksvm}} for training a SVM model on a training set of persistence diagrams.
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom parallel makeCluster stopCluster clusterExport clusterEvalQ
#' @importFrom parallelly availableCores
#' @importFrom doParallel registerDoParallel
#' @importFrom kernlab predict as.kernelMatrix
#' @importFrom methods is
#' @examples
#'
#' # create four diagrams
#' D1 <- TDAstats::calculate_homology(TDA::circleUnif(n = 10,r = 1),
#'                                    dim = 0,threshold = 2)
#' D2 <- TDAstats::calculate_homology(TDA::circleUnif(n = 10,r = 1),
#'                                    dim = 0,threshold = 2)
#' D3 <- TDAstats::calculate_homology(TDA::circleUnif(n = 10,r = 1),
#'                                    dim = 0,threshold = 2)
#' D4 <- TDAstats::calculate_homology(TDA::circleUnif(n = 10,r = 1),
#'                                    dim = 0,threshold = 2)
#' g <- list(D1,D2,D3,D4)
#' 
#' # create response vector
#' y <- as.factor(c("circle","sphere","circle","sphere"))
#' 
#' # fit model without cross validation
#' model_svm <- diagram_ksvm(diagrams = g,cv = 1,dim = c(0),
#'                           y = y,sigma = c(1),t = c(1),
#'                           num_workers = 2)
#'
#' # create two new diagrams
#' D5 <- TDAstats::calculate_homology(TDA::circleUnif(n = 10,r = 1),
#'                                    dim = 0,threshold = 2)
#' D6 <- TDAstats::calculate_homology(TDA::circleUnif(n = 10,r = 1),
#'                                    dim = 0,threshold = 2)
#' g_new <- list(D5,D6)
#' 
#' # predict
#' predict_diagram_ksvm(new_diagrams = g_new,model = model_svm,num_workers = 2)

predict_diagram_ksvm <- function(new_diagrams,model,num_workers = parallelly::availableCores(omit = 1)){
  
  # set internal variables to NULL to avoid build issues
  r <- NULL
  X <- NULL
  
  # error check new_diagrams argument
  check_param("new_diagrams",new_diagrams,numeric = F,multiple = T,min_length = 1)
  new_diagrams <- all_diagrams(diagram_groups = list(new_diagrams),inference = "independence")[[1]]
  
  # error check model argument
  if(is.null(model))
  {
    stop("model must be supplied.")
  }
  if(!methods::is(model,"diagram_ksvm"))
  {
    stop("model must be the output of a diagram_ksvm function call.")
  }
  
  # compute kernel matrix, storing the value of each kernel computation between the new diagrams and the old ones
  K <- gram_matrix(diagrams = new_diagrams,other_diagrams = model$best_model$diagrams,dim = model$best_model$dim,sigma = model$best_model$sigma,t = model$best_model$t,num_workers = num_workers)

  # suppress unhelpful kernlab warnings when predicting
  tryCatch(expr = {predictions <- kernlab::predict(object = model$best_model$model,kernlab::as.kernelMatrix(K))},
           error = function(e){stop(e)},
           warning = function(w){
             
             if(grepl(pattern = "closing unused connection",w) == F)
             {
               message(w)
             }
             
           })
  
  return(predictions)
  
}
