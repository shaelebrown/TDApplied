
#### Multidimensional scaling ####
#' Dimension reduction of a group of persistence diagrams via metric multidimensional scaling.
#'
#' Projects a group of persistence diagrams (or a precomputed distance matrix of diagrams) into a low-dimensional 
#' embedding space via metric multidimensional scaling. Such a projection can be used for visualization of data, 
#' or a static analysis of the embedding dimensions.
#' 
#' Returns the output of \code{\link[stats]{cmdscale}} on the desired distance matrix of a group of persistence diagrams
#' in a particular dimension. If `distance` is "fisher" then `sigma` must not be NULL.
#'
#' @param diagrams a list of n>=2 persistence diagrams which are either the output of a persistent homology calculation like ripsDiag/\code{\link[TDAstats]{calculate_homology}}/\code{\link{PyH}}, or \code{\link{diagram_to_df}}. Only one of `diagrams` and `D` need to be supplied.
#' @param D an optional precomputed distance matrix of persistence diagrams, default NULL. If not NULL then `diagrams` parameter does not need to be supplied.
#' @param k the dimension of the space which the data are to be represented in; must be in \{1,2,...,n-1\}.
#' @param distance a string representing the desired distance metric to be used, either 'wasserstein' (default) or 'fisher'.
#' @param dim the non-negative integer homological dimension in which the distance is to be computed, default 0.
#' @param p a positive number representing the wasserstein power, a number at least 1 (infinity for the bottleneck distance), default 2.
#' @param sigma a positive number representing the bandwidth for the Fisher information metric, default NULL.
#' @param rho an optional positive number representing the heuristic for Fisher information metric approximation, see \code{\link{diagram_distance}}. Default NULL. If supplied, distance matrix calculation is sequential.
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
#' @export
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @references 
#' Cox M and Cox F (2008). "Multidimensional Scaling." \doi{10.1007/978-3-540-33037-0_14}.
#' 
#' @examples
#'
#' if(require("TDAstats"))
#' {
#'   # create two diagrams
#'   D1 <- TDAstats::calculate_homology(TDAstats::circle2d[sample(1:100,10),],
#'                       dim = 1,threshold = 2)
#'   D2 <- TDAstats::calculate_homology(TDAstats::circle2d[sample(1:100,10),],
#'                       dim = 1,threshold = 2)
#'   g <- list(D1,D2)
#' 
#'   # calculate their 1D MDS embedding in dimension 0 with the bottleneck distance
#'   mds <- diagram_mds(diagrams = g,k = 1,dim = 0,p = Inf,num_workers = 2)
#'   
#'   # repeat but with a precomputed distance matrix, gives same result just much faster
#'   Dmat <- distance_matrix(diagrams = list(D1,D2),dim = 0,p = Inf,num_workers = 2)
#'   mds <- diagram_mds(D = Dmat,k = 1)
#'   
#' }

diagram_mds <- function(diagrams,D = NULL,k = 2,distance = "wasserstein",dim = 0,p = 2,sigma = NULL,rho = NULL,eig = FALSE,add = FALSE,x.ret = FALSE,list. = eig || add || x.ret,num_workers = parallelly::availableCores(omit = 1)){
  
  # function to embed a group of persistence diagrams in low dimensions
  # using MDS
  
  if(is.null(D))
  {
    # distance matrix not supplied
    
    # error check diagrams argument
    check_param("diagrams",diagrams,min_length = 2)
    diagrams <- all_diagrams(diagram_groups = list(diagrams),inference = "independence")[[1]]
    
    # compute distance matrix
    D <- distance_matrix(diagrams = diagrams,dim = dim,distance = distance,p = p,sigma = sigma,rho = rho,num_workers = num_workers)
  }else
  {
    # check supplied distance matrix
    check_matrix(M = D,name = "D",type = "matrix")
  }
  
  # check k parameter (embedding dimension)
  check_param("k",k,whole_numbers = T,at_least_one = T,numeric = T,finite = T,multiple = F)
  
  # return metric multidimensional scaling with d as input
  return(stats::cmdscale(d = D,k = k,eig = eig,add = add,x.ret = x.ret,list. = list.))
  
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
#' @param diagrams a list of n>=2 persistence diagrams which are either the output of a persistent homology calculation like ripsDiag/\code{\link[TDAstats]{calculate_homology}}/\code{\link{PyH}}, or the \code{\link{diagram_to_df}} function.
#' @param K an optional precomputed Gram matrix of persistence diagrams, default NULL.
#' @param dim the non-negative integer homological dimension in which the distance is to be computed, default 0.
#' @param t a positive number representing the scale for the persistence Fisher kernel, default 1.
#' @param sigma a positive number representing the bandwidth for the Fisher information metric, default 1.
#' @param rho an optional positive number representing the heuristic for Fisher information metric approximation, see \code{\link{diagram_distance}}. Default NULL. If supplied, Gram matrix calculation is sequential.
#' @param centers number of clusters to initialize, no more than the number of diagrams although smaller values are recommended.
#' @param num_workers the number of cores used for parallel computation, default is one less than the number of cores on the machine.
#' @param ... additional parameters for the \code{\link[kernlab]{kkmeans}} kernlab function.
#'
#' @return a list of class 'diagram_kkmeans' containing the output of \code{\link[kernlab]{kkmeans}} on the Gram matrix, i.e. a list containing the elements
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
#' if(require("TDAstats"))
#' {
#'   # create two diagrams
#'   D1 <- TDAstats::calculate_homology(TDAstats::circle2d[sample(1:100,20),],
#'                       dim = 1,threshold = 2)
#'   D2 <- TDAstats::calculate_homology(TDAstats::circle2d[sample(1:100,20),],
#'                       dim = 1,threshold = 2)
#'   g <- list(D1,D1,D2,D2)
#' 
#'   # calculate kmeans clusters with centers = 2, and sigma = t = 2 in dimension 0
#'   clust <- diagram_kkmeans(diagrams = g,centers = 2,dim = 0,t = 2,sigma = 2,num_workers = 2)
#'   
#'   # repeat with precomputed Gram matrix, gives the same result just much faster
#'   K <- gram_matrix(diagrams = g,num_workers = 2,t = 2,sigma = 2)
#'   cluster <- diagram_kkmeans(diagrams = g,K = K,centers = 2,dim = 0,sigma = 2,t = 2)
#'   
#' }

diagram_kkmeans <- function(diagrams,K = NULL,centers,dim = 0,t = 1,sigma = 1,rho = NULL,num_workers = parallelly::availableCores(omit = 1),...){
  
  # function to cluster a group of persistence diagrams
  
  # error check centers (number of clusters) parameter
  check_param("centers",centers,whole_numbers = T,at_least_one = T,numeric = T,finite = T,multiple = F)
  
  if(is.null(K))
  {
    # Gram matrix not supplied
    
    # error check diagrams
    check_param("diagrams",diagrams,min_length = 2)
    diagrams <- all_diagrams(diagram_groups = list(diagrams),inference = "independence")[[1]]
    
    # make sure there aren't more centers than data points - this tends to lead to kernlab errors
    if(centers > length(diagrams))
    {
      stop("centers must be at most the number of diagrams, although we recommend choosing smaller values.")
    }
    
    # compute Gram matrix
    K <- gram_matrix(diagrams = diagrams,dim = dim,sigma = sigma,t = t,rho = rho,num_workers = num_workers) 
  }else
  {
    check_matrix(M = K,name = "K")
    # make sure there aren't more centers than data points - this tends to lead to kernlab errors
    if(centers > ncol(K))
    {
      stop("centers must be at most the number of columns of K, although we recommend choosing smaller values.")
    }
    
    if(ncol(K) != length(diagrams))
    {
      stop("K must have the same number of rows and columns as the length of diagrams.")
    }
    
  }
  
  # return kernlab calculation, saving in a list of class diagram_kkmeans for later interface with prediction calculations
  
  # sometimes the algorithm doesn't converge by fluke (randomness), so we will
  # run up to 20 iterations to see if we have successful convergence
  # otherwise an error will be thrown
  max_iters <- 20
  iter <- 1
  while(T)
  {
    # there are some kernlab errors that can't always be avoided and caused by randomness
    # when certain clusters are empty or unclosed connections
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
  
  # set up return list
  ret_list <- list(clustering = clustering,diagrams = diagrams,dim = dim,sigma = sigma,t = t,rho = rho)

  # set class for easy recognition in prediction method
  class(ret_list) <- "diagram_kkmeans"

  return(ret_list)
  
}

#### PREDICT KERNEL KMEANS ####
#' Predict the cluster labels for new persistence diagrams using a pre-computed clustering.
#'
#' Returns the nearest (highest kernel value) \code{\link[kernlab]{kkmeans}} cluster center label for new persistence diagrams.
#' This allows for reusing old cluster models for new tasks, or to perform cross validation.
#'
#' @param new_diagrams a list of persistence diagrams which are either the output of a persistent homology calculation like ripsDiag/\code{\link[TDAstats]{calculate_homology}}/\code{\link{PyH}}, or \code{\link{diagram_to_df}}. Only one of `new_diagrams` and `K` need to be supplied.
#' @param K an optional precomputed cross Gram matrix of the new diagrams and the diagrams used in `clustering`, default NULL. If not NULL then `new_diagrams` does not need to be supplied.
#' @param clustering the output of a \code{\link{diagram_kkmeans}} function call, of class 'diagram_kkmeans'.
#' @param num_workers the number of cores used for parallel computation, default is one less than the number of cores on the machine.
#'
#' @return a vector of the predicted cluster labels for the new diagrams.
#' @export
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @seealso \code{\link{diagram_kkmeans}} for clustering persistence diagrams.
#' @examples
#'
#' if(require("TDAstats"))
#' {
#'   # create two diagrams
#'   D1 <- TDAstats::calculate_homology(TDAstats::circle2d[sample(1:100,20),],
#'                       dim = 1,threshold = 2)
#'   D2 <- TDAstats::calculate_homology(TDAstats::circle2d[sample(1:100,20),],
#'                       dim = 1,threshold = 2)
#'   g <- list(D1,D1,D2,D2)
#' 
#'   # calculate kmeans clusters with centers = 2, and sigma = t = 2 in dimension 0
#'   clust <- diagram_kkmeans(diagrams = g,centers = 2,dim = 0,t = 2,sigma = 2,num_workers = 2)
#' 
#'   # create two new diagrams
#'   D3 <- TDAstats::calculate_homology(TDAstats::circle2d[sample(1:100,20),],
#'                       dim = 1,threshold = 2)
#'   D4 <- TDAstats::calculate_homology(TDAstats::circle2d[sample(1:100,20),],
#'                       dim = 1,threshold = 2)
#'   g_new <- list(D3,D4)
#' 
#'   # predict cluster labels
#'   predict_diagram_kkmeans(new_diagrams = g_new,clustering = clust,num_workers = 2)
#'   
#'   # predict cluster labels with precomputed Gram matrix, gives same result but
#'   # much faster
#'   K <- gram_matrix(diagrams = g_new,other_diagrams = clust$diagrams,
#'                    dim = clust$dim,t = clust$t,sigma = clust$sigma,
#'                    num_workers = 2)
#'   predict_diagram_kkmeans(K = K,clustering = clust)
#'   
#' }

predict_diagram_kkmeans <- function(new_diagrams,K = NULL,clustering,num_workers = parallelly::availableCores(omit = 1)){
  
  # function to predict the cluster labels of new persistence diagrams
  # with a pre-fit clustering model
  
  # set internal variables to NULL to avoid build issues
  X <- NULL
  
  # error check clustering argument
  if(is.null(clustering))
  {
    stop("clustering must not be NULL.")
  }
  if(!inherits(clustering,"diagram_kkmeans"))
  {
    stop("clustering object must be the output of a diagram_kkmeans function call.")
  }
  
  if(is.null(K))
  {
    # cross Gram matrix not supplied
    
    # error check diagrams argument
    check_param("new_diagrams",new_diagrams,min_length = 1)
    new_diagrams <- all_diagrams(diagram_groups = list(new_diagrams),inference = "independence")[[1]]
    
    # compute cross Gram matrix
    K = gram_matrix(diagrams = new_diagrams,other_diagrams = clustering$diagrams,dim = clustering$dim,sigma = clustering$sigma,t = clustering$t,rho = clustering$rho,num_workers = num_workers)
  }else
  {
    # cross Gram matrix supplied
    if(missing(new_diagrams))
    {
      new_diagrams <- rep(0,nrow(K))
    }
    
    # error check supplied matrix
    check_matrix(M = K,name = "K",symmetric = F)
    if(nrow(K) != length(new_diagrams) | ncol(K) != length(clustering$diagrams))
    {
      stop("K must have the same number of rows as the length of new_diagrams and the same number of columns as the length of diagrams in clustering.")
    }
  }
  
  # return predicted class for each new diagram
  predicted_clusters <- unlist(lapply(X = 1:nrow(K),FUN = function(X){
    
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
#' @param diagrams a list of persistence diagrams which are either the output of a persistent homology calculation like ripsDiag/\code{\link[TDAstats]{calculate_homology}}/\code{\link{PyH}}, or \code{\link{diagram_to_df}}.
#' @param K an optional precomputed Gram matrix of the persistence diagrams in `diagrams`, default NULL.
#' @param dim the non-negative integer homological dimension in which the distance is to be computed, default 0.
#' @param t a positive number representing the scale for the persistence Fisher kernel, default 1.
#' @param sigma a positive number representing the bandwidth for the Fisher information metric, default 1.
#' @param rho an optional positive number representing the heuristic for Fisher information metric approximation, see \code{\link{diagram_distance}}. Default NULL. If supplied, Gram matrix calculation is sequential.
#' @param features number of features (principal components) to return, default 1.
#' @param num_workers the number of cores used for parallel computation, default is one less than the number of cores on the machine.
#' @param th the threshold value under which principal components are ignored (default 0.0001).
#'
#' @return a list of class 'diagram_kpca' containing the elements
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
#' if(require("TDAstats"))
#' {
#'   # create six diagrams
#'   D1 <- TDAstats::calculate_homology(TDAstats::circle2d[sample(1:100,20),],
#'                       dim = 1,threshold = 2)
#'   D2 <- TDAstats::calculate_homology(TDAstats::circle2d[sample(1:100,20),],
#'                       dim = 1,threshold = 2)
#'   D3 <- TDAstats::calculate_homology(TDAstats::sphere3d[sample(1:100,20),],
#'                       dim = 1,threshold = 2)
#'   D4 <- TDAstats::calculate_homology(TDAstats::sphere3d[sample(1:100,20),],
#'                       dim = 1,threshold = 2)
#'   D5 <- TDAstats::calculate_homology(TDAstats::sphere3d[sample(1:100,20),],
#'                       dim = 1,threshold = 2)
#'   D6 <- TDAstats::calculate_homology(TDAstats::sphere3d[sample(1:100,20),],
#'                       dim = 1,threshold = 2)
#'   g <- list(D1,D2,D3,D4,D5,D6)
#' 
#'   # calculate their 2D PCA embedding with sigma = t = 2 in dimension 1
#'   pca <- diagram_kpca(diagrams = g,dim = 1,t = 2,sigma = 2,features = 2,num_workers = 2,th = 1e-6)
#'   
#'   # repeat with precomputed Gram matrix, gives same result but much faster
#'   K <- gram_matrix(diagrams = g,dim = 1,t = 2,sigma = 2,num_workers = 2)
#'   pca <- diagram_kpca(diagrams = g,K = K,dim = 1,t = 2,sigma = 2,features = 2,th = 1e-6)
#'   
#' }

diagram_kpca <- function(diagrams,K = NULL,dim = 0,t = 1,sigma = 1,rho = NULL,features = 1,num_workers = parallelly::availableCores(omit = 1),th = 1e-4){
  
  # function to embed a group of diagrams into low dimensions
  # using kernel PCA
  
  # error check diagrams argument
  check_param("diagrams",diagrams,min_length = 2)
  diagrams <- all_diagrams(diagram_groups = list(diagrams),inference = "independence")[[1]]
  
  if(is.null(K))
  {
    # compute Gram matrix
    K <- gram_matrix(diagrams = diagrams,t = t,sigma = sigma,rho = rho,dim = dim,num_workers = num_workers)
  }else
  {
    # error check Gram matrix
    check_matrix(M = K,name = "K")
    if(length(diagrams) != nrow(K))
    {
      stop("K must have the same number of rows and columns as the length of diagrams.")
    }
  }
  
  # kernlab computation, ignore unhelpful kernlab warnings
  tryCatch(expr = {ret_list <- list(pca = kernlab::kpca(x = K,features = features,th = th),diagrams = diagrams,t = t,sigma = sigma,dim = dim,rho = rho)},
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
  
  # set class for easy recognition in prediction method
  class(ret_list) <- "diagram_kpca"
  return(ret_list)
  
}

#### PREDICT WITH KERNEL PCA OBJECT ####
#' Project persistence diagrams into a low-dimensional space via a pre-computed kernel PCA embedding.
#'
#' Compute the location in low-dimensional space of each element of a list of new persistence diagrams using a
#' previously-computed kernel PCA embedding (from the \code{\link{diagram_kpca}} function).
#'
#' @param new_diagrams a list of persistence diagrams which are either the output of a persistent homology calculation like ripsDiag/\code{\link[TDAstats]{calculate_homology}}/\code{\link{PyH}}, or \code{\link{diagram_to_df}}. Only one of `new_diagrams` and `K` need to be supplied.
#' @param K an optional precomputed cross-Gram matrix of the new diagrams and the ones used in `embedding`, default NULL. If not NULL then `new_diagrams` does not need to be supplied.
#' @param embedding the output of a \code{\link{diagram_kpca}} function call, of class 'diagram_kpca'.
#' @param num_workers the number of cores used for parallel computation, default is one less than the number of cores on the machine.
#'
#' @return the data projection (rotation), stored as a numeric matrix. Each row corresponds to the same-index diagram in `new_diagrams`.
#' @export
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @seealso \code{\link{diagram_kpca}} for embedding persistence diagrams into a low-dimensional space.
#' @examples
#'
#' if(require("TDAstats"))
#' {
#'   # create six diagrams
#'   D1 <- TDAstats::calculate_homology(TDAstats::circle2d[sample(1:100,20),],
#'                       dim = 1,threshold = 2)
#'   D2 <- TDAstats::calculate_homology(TDAstats::circle2d[sample(1:100,20),],
#'                       dim = 1,threshold = 2)
#'   D3 <- TDAstats::calculate_homology(TDAstats::sphere3d[sample(1:100,20),],
#'                       dim = 1,threshold = 2)
#'   D4 <- TDAstats::calculate_homology(TDAstats::sphere3d[sample(1:100,20),],
#'                       dim = 1,threshold = 2)
#'   D5 <- TDAstats::calculate_homology(TDAstats::sphere3d[sample(1:100,20),],
#'                       dim = 1,threshold = 2)
#'   D6 <- TDAstats::calculate_homology(TDAstats::sphere3d[sample(1:100,20),],
#'                       dim = 1,threshold = 2)
#'   g <- list(D1,D2,D3,D4,D5,D6)
#' 
#'   # calculate their 2D PCA embedding with sigma = t = 2 in dimension 0
#'   pca <- diagram_kpca(diagrams = g,dim = 1,t = 2,sigma = 2,
#'                       features = 2,num_workers = 2,th = 1e-6)
#' 
#'   # project two new diagrams onto old model
#'   D7 <- TDAstats::calculate_homology(TDAstats::circle2d[sample(1:100,50),],
#'                                      dim = 0,threshold = 2)
#'   D8 <- TDAstats::calculate_homology(TDAstats::circle2d[sample(1:100,50),],
#'                                      dim = 0,threshold = 2)
#'   g_new <- list(D7,D8)
#' 
#'   # calculate new embedding coordinates
#'   new_pca <- predict_diagram_kpca(new_diagrams = g_new,embedding = pca,num_workers = 2)
#'   
#'   # repeat with precomputed Gram matrix, gives same result but much faster
#'   K <- gram_matrix(diagrams = g_new,other_diagrams = pca$diagrams,dim = pca$dim,
#'                    t = pca$t,sigma = pca$sigma,num_workers = 2)
#'   new_pca <- predict_diagram_kpca(K = K,embedding = pca,num_workers = 2)
#' }

predict_diagram_kpca <- function(new_diagrams,K = NULL,embedding,num_workers = parallelly::availableCores(omit = 1)){
  
  # function to predict embedding coordinates of new diagrams
  # given a pre-fit kPCA model
  
  # set internal variables to NULL to avoid build issues
  r <- NULL
  X <- NULL
  
  # error check embedding argument
  if(is.null(embedding))
  {
    stop("embedding must be supplied.")
  }
  if(!inherits(embedding,"diagram_kpca"))
  {
    stop("embedding must be the output of a diagram_kpca function call.")
  }
  
  # compute cross-Gram matrix and scale twice
  if(is.null(K))
  {
    # cross Gram matrix not supplied
    
    # error check new_diagrams argument
    check_param("new_diagrams",new_diagrams,min_length = 1)
    new_diagrams <- all_diagrams(diagram_groups = list(new_diagrams),inference = "independence")[[1]]
    
    # compute cross Gram matrix
    K <- gram_matrix(diagrams = new_diagrams,other_diagrams = embedding$diagrams,dim = embedding$dim,sigma = embedding$sigma,t = embedding$t,rho = embedding$rho,num_workers = num_workers)
  }else
  {
   # cross Gram matrix supplied
   if(missing(new_diagrams))
   {
     new_diagrams <- rep(0,nrow(K))
   }
    
   # check matrix
   check_matrix(M = K,name = "K",symmetric = F) 
   if(nrow(K) != length(new_diagrams) | ncol(K) != length(embedding$diagrams))
   {
     stop("K must have the same number of rows as the length of new_diagrams and the same number of columns as the length of diagrams in embedding.")
   }
  }
  
  # scale and double center Gram matrix
  K <- scale(K,center = T,scale = F)
  K <- t(scale(t(K),center = T,scale = F))
  
  # project the new diagrams into the embedding space
  return(K %*% embedding$pca@pcv)
  
}

#### KERNEL SVM ####
#' Fit a support vector machine model where each training set instance is a persistence diagram.
#'
#' Returns the output of kernlab's \code{\link[kernlab]{ksvm}} function on the Gram matrix of the list of persistence diagrams
#' in a particular dimension.
#' 
#' Cross validation is carried out in parallel, using a trick
#' noted in \doi{10.1007/s41468-017-0008-7} - since the persistence Fisher kernel can be
#' written as \eqn{d_{PF}(D_1,D_2)=exp(t*d_{FIM}(D_1,D_2))=exp(d_{FIM}(D_1,D_2))^t}, we can
#' store the Fisher information metric distance matrix for each sigma value in the parameter grid to avoid
#' recomputing distances, and cross validation is therefore performed in parallel. 
#' Note that the response parameter `y` must be a factor for classification - 
#' a character vector for instance will throw an error. If `t` is NULL then 1/`t` is selected as
#' the 1,2,5,10,20,50 percentiles of the upper triangle of the distance matrix of its training sample (per fold in the case of cross-validation). 
#' This is the process suggested in the persistence Fisher kernel paper. If
#' any of these values would divide by 0 (i.e. if the training set is small) then the minimum non-zero element
#' is taken as the denominator (and hence the returned parameters may have duplicate rows except for differing error values). If
#' cross-validation is performed then the mean error across folds is still recorded, but the best `t` parameter
#' across all folds is recorded in the cv results table.
#'
#' @param diagrams a list of persistence diagrams which are either the output of a persistent homology calculation like ripsDiag/\code{\link[TDAstats]{calculate_homology}}/\code{\link{PyH}}, or \code{\link{diagram_to_df}}.
#' @param cv a positive number at most the length of `diagrams` which determines the number of cross validation splits to be performed (default 1, aka no cross-validation). If `prob.model` is TRUE then cv is set to 1 since kernlab performs 3-fold CV internally in this case. When performing classification, classes are balanced within each cv fold.
#' @param dim a non-negative integer vector of homological dimensions in which the model is to be fit.
#' @param t either a vector of positive numbers representing the grid of values for the scale of the persistence Fisher kernel or NULL, default 1. If NULL then t is selected automatically, see details.
#' @param sigma a vector of positive numbers representing the grid of values for the bandwidth of the Fisher information metric, default 1.
#' @param rho an optional positive number representing the heuristic for Fisher information metric approximation, see \code{\link{diagram_distance}}. Default NULL. If supplied, distance matrix calculations are sequential.
#' @param y a response vector with one label for each persistence diagram. Must be either numeric or factor, but doesn't need to be supplied when `type` is "one-svc".
#' @param type a string representing the type of task to be performed. Can be any one of "C-svc","nu-svc","one-svc","eps-svr","nu-svr" - default for regression is "eps-svr" and for classification is "C-svc". See \code{\link[kernlab]{ksvm}} for details.
#' @param distance_matrices an optional list of precomputed Fisher distance matrices, corresponding to the rows in `expand.grid(dim = dim,sigma = sigma)`, default NULL.
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
#' @return a list of class 'diagram_ksvm' containing the elements
#' 
#' \describe{
#' 
#' \item{cv_results}{the cross-validation results - a matrix storing the parameters for each model in the tuning grid and its mean cross-validation error over all splits.}
#' 
#' \item{best_model}{a list containing the output of \code{\link[kernlab]{ksvm}} run on the whole dataset with the optimal model parameters found during cross-validation, as well as the optimal kernel parameters for the model.}
#' 
#' \item{diagrams}{the diagrams which were supplied in the function call.}
#' 
#' }
#' 
#' @export
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @seealso \code{\link{predict_diagram_ksvm}} for predicting labels of new diagrams.
#' @importFrom kernlab ksvm
#' @importFrom foreach foreach %dopar% %:%
#' @importFrom parallel makeCluster stopCluster clusterExport
#' @importFrom parallelly availableCores
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom iterators iter
#' @references 
#' Murphy, K. "Machine learning: a probabilistic perspective." MIT press (2012).
#' @examples
#'
#' if(require("TDAstats"))
#' {
#'   # create four diagrams
#'   D1 <- TDAstats::calculate_homology(TDAstats::circle2d[sample(1:100,20),],
#'                       dim = 1,threshold = 2)
#'   D2 <- TDAstats::calculate_homology(TDAstats::circle2d[sample(1:100,20),],
#'                       dim = 1,threshold = 2)
#'   D3 <- TDAstats::calculate_homology(TDAstats::sphere3d[sample(1:100,20),],
#'                       dim = 1,threshold = 2)
#'   D4 <- TDAstats::calculate_homology(TDAstats::sphere3d[sample(1:100,20),],
#'                       dim = 1,threshold = 2)
#'   g <- list(D1,D2,D3,D4)
#' 
#'   # create response vector
#'   y <- as.factor(c("circle","circle","sphere","sphere"))
#' 
#'   # fit model without cross validation
#'   model_svm <- diagram_ksvm(diagrams = g,cv = 1,dim = c(0),
#'                             y = y,sigma = c(1),t = c(1),
#'                             num_workers = 2)
#' }
                          
diagram_ksvm <- function(diagrams,cv = 1,dim,t = 1,sigma = 1,rho = NULL,y,type = NULL,distance_matrices = NULL,C = 1,nu = 0.2,epsilon = 0.1,prob.model = FALSE,class.weights = NULL,fit = TRUE,cache = 40,tol = 0.001,shrinking = TRUE,num_workers = parallelly::availableCores(omit = 1)){
  
  # fit a SVM model to a group of persistence diagrams and a set of labels
  
  # set internal variables to NULL to avoid build issues
  r <- NULL
  s <- NULL
  X <- NULL
  
  # error check diagrams argument
  check_param("diagrams",diagrams,min_length = 2)
  diagrams <- all_diagrams(diagram_groups = list(diagrams),inference = "independence")[[1]]
  
  # error check cv parameters
  check_param("cv",cv,whole_numbers = T,finite = T,numeric = T,multiple = F)
  if(cv < 1 | cv > length(diagrams))
  {
    stop("cv must be at least one and at most the number of diagrams.")
  }
  
  # basic error checks for y
  if(missing(y))
  {
    if(!is.null(type))
    {
      if(type == "one-svc")
      {
        y <- rep(0,length(diagrams))
      }else
      {
        stop("y must be supplied.")
      }
    }else
    {
      stop("y must be supplied.")
    }
  }
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
  }else
  {
    if(type %in% c("C-svc","nu-svc","one-svc","eps-svr","nu-svr") == F)
    {
      stop("type must be one of \'C-svc\', \'nu-svc\', \'one-svc\', \'eps-svr\' or \'nu-svr\'.")
    }
  }
  
  if(type == "nu-svc" & is.factor(y) & length(levels(y)) > 2)
  {
    stop("Currently there are errors for multiclass classification with nu-svc. Try using C-svc instead.")
  }
  
  # if prob.model is TRUE set cv = 1
  if(prob.model == T)
  {
    cv <- 1
  }
  
  # error check t
  if(!is.null(t))
  {
    check_param(param_name = "t",param = t,finite = T,numeric = T,multiple = T,positive = T) 
  }
  
  # expand parameter grid
  if(!is.null(t))
  {
    params <- expand.grid(dim = dim,t = t,sigma = sigma,C = C,nu = nu,epsilon = epsilon)
    estimate_t <- F
  }else
  {
    params <- expand.grid(dim = dim,t = 0.01*c(1,2,5,10,20,50),sigma = sigma,C = C,nu = nu,epsilon = epsilon)
    estimate_t <- T
  }
  
  # make vector of row memberships for cv
  # balanced by class for C-svc and nu-svc
  if(type %in% c("C-svc","nu-svc") == F)
  {
    v <- rep(1:cv,floor(length(diagrams)/cv))
    if(length(v) != length(diagrams))
    {
      if(length(diagrams) - length(v) == 1)
      {
        v <- c(v,1)
      }else
      {
        v <- c(v,sample(1:cv,size = length(diagrams) - length(v))) 
      }
    }
    v <- sample(v,size = length(v),replace = F)
  }else
  {
    v <- rep(1,length(diagrams))
    classes <- levels(y)
    for(cl in classes)
    {
      class_inds <- which(y == cl)
      if(length(class_inds) < cv)
      {
        stop("One class of y has fewer training examples than the number of cv folds - try decreasing cv parameter.")
      }
      v_class <- rep(1:cv,floor(length(class_inds)/cv))
      if(length(v_class) != length(class_inds))
      {
        v_class <- c(v_class,sample(1:cv,size = length(class_inds) - length(v_class))) 
      }
      v_class <- sample(v_class,size = length(v_class),replace = F)
      v[class_inds] <- v_class
    }
  }
  
  if(cv > 1 & sum(as.numeric(table(v)) == 1) > 0)
  {
    stop(paste0("Too few training examples to perform cv with ",cv," folds. Try decreasing cv parameter or fitting the model without cv."))
  } 
  
  # split diagrams by membership vector
  diagrams_split <- lapply(X = 1:cv,FUN = function(X){
    
    return(which(v == X))
    
  })
  
  # for each pair of sigma and dim values compute the Fisher distance matrix to avoid recomputing Gram matrices
  # if these matrices have not been supplied
  dim_and_sigma <- expand.grid(dim = dim,sigma = sigma)
  if(is.null(distance_matrices))
  {
    # distance matrices not supplied, compute them
    distance_matrices <- lapply(X = 1:nrow(dim_and_sigma),FUN = function(X){
      
      return(distance_matrix(diagrams = diagrams,dim = dim_and_sigma[X,1],distance = "fisher",sigma = dim_and_sigma[X,2],rho = rho,num_workers = num_workers))
      
    })
  }else
  {
    # error check supplied distance matrices
    if(!is.list(distance_matrices))
    {
      stop("distance_matrices must be a list.")
    }
    lapply(X = distance_matrices,FUN = function(X){
      
      check_matrix(M = X,name = "distance matrices",type = "matrix")
      
    })
    if(length(distance_matrices)!=nrow(dim_and_sigma))
    {
      stop("distance_matrices must have one entry per row of expand.grid(dim = dim,sigma = sigma).")
    }
  }
  names(distance_matrices) <- paste(dim_and_sigma$dim,dim_and_sigma$sigma,sep = "_")
  
  # check if distance matrices (within cv folds) have 0 variance
  zero_var_inds <- which(unlist(lapply(X = distance_matrices,FUN = function(X){
    
    X_temp <- X
    prod_vars <- unlist(lapply(X = 1:cv,FUN = function(X){
      
      cv_inds <- which(v == X)
      return(stats::var(as.vector(X_temp[cv_inds,cv_inds])))
      
    }))
    prod_vars <- prod(prod_vars)
    
    return(prod_vars == 0)
    
  })))
  
  # throw error if all distance matrices have 0 variances
  if(length(zero_var_inds) == length(distance_matrices))
  {
    if(cv == 1)
    {
      stop("All distance matrices have 0 variance.")
    }else
    {
      stop("All distance matrices have 0 variance in at least one cv fold.")
    }
  }
  
  # remove "zero variance parameters" from model fitting
  if(length(zero_var_inds) > 0)
  {
    zero_var_names <- names(zero_var_inds)
    zero_var_inds <- unlist(lapply(X = zero_var_names,FUN = function(X){
      
      # split name to get dim and sigma
      s <- strsplit(X,split = "_")[[1]]
      d <- s[[1]]
      sig <- s[[2]]
      return(which(params$dim == d & params$sigma == sig))
      
    }))
    zero_var_params <- params[zero_var_inds,]
    zero_var_params$error <- NA
  }else
  {
    zero_var_params <- params[0,]
    zero_var_params$error <- numeric()
  }
  params <- params[setdiff(1:nrow(params),zero_var_inds),]
  
  # set up for parallel computation
  cl <- parallel::makeCluster(num_workers)
  doParallel::registerDoParallel(cl)
  parallel::clusterExport(cl,c("all_diagrams","check_diagram","gram_matrix","diagram_distance","diagram_kernel","check_param","diagram_to_df",".classAgreement"),envir = environment())
  
  # for each model (combination of parameters), train the model on each subset and get the
  # prediction error on the hold out set, wrapped in tryCatch to ensure cluster is stopped
  if(!estimate_t)
  {
    tryCatch(expr = {model_errors <- foreach::`%dopar%`(foreach::`%:%`(foreach::foreach(r = iter(params,by = "row"),.combine = c,.noexport = c("diagrams")),foreach::foreach(s = 1:cv,.combine = "+",.packages = c("clue","rdist","kernlab","stats"),.noexport = c("diagrams"))),ex = {
      
      # use precomputed distance matrices to calculate Gram matrix from parameters
      K <- exp(-1*r[[2]]*distance_matrices[[paste0(r[[1]],"_",r[[3]])]])
      
      if(cv > 1)
      {
        # split into training and test set based on cv parameter
        training_indices <- unlist(diagrams_split[setdiff(1:cv,s)])
        training_indices <- training_indices[order(training_indices)]
        test_indices <- setdiff(1:nrow(K),training_indices)
        
        # fit model on training folds
        K_subset <- K[training_indices,training_indices]
        class(K_subset) <- "kernelMatrix"
        model = kernlab::ksvm(x = K_subset,y = y[training_indices],type = type,C = r[[4]],nu = r[[5]],epsilon = r[[6]],prob.model = prob.model,class.weights = class.weights,fit = fit,cache = cache,tol = tol,shrinking = shrinking)
        
        # predict model on single test fold
        K_cross <- K[test_indices,training_indices[model@SVindex]]
        class(K_cross) <- "kernelMatrix"
        predictions <- kernlab::predict(model,K_cross)
        
        # calculate error depending on model type
        if(type %in% c("C-svc","nu-svc"))
        {
          return((1 - .classAgreement(table(y[test_indices],as.character(predictions)))))
        }
        if(type == "one-svc")
        {
          return((1 - sum(predictions)/length(predictions)))
        }
        if(type %in% c("eps-svr","nu-svr"))
        {
          return(drop(crossprod(predictions - y[test_indices])/nrow(K)))
        }
        
      }else
      {
        class(K) <- "kernelMatrix"
        model <- kernlab::ksvm(x = K,y = y,type = type,C = r[[4]],nu = r[[5]],epsilon = r[[6]],prob.model = prob.model,class.weights = class.weights,fit = fit,cache = cache,tol = tol,shrinking = shrinking)
        return(model@error)
      }
    })/cv},
    error = function(e){stop(e)},
    warning = function(w){warning(w)},
    finally = {
      # clean up
      doParallel::stopImplicitCluster()
      parallel::stopCluster(cl)
      
    }) 
  }else
  {
    tryCatch(expr = {model_errors <- foreach::`%dopar%`(foreach::`%:%`(foreach::foreach(r = iter(params,by = "row"),.noexport = c("diagrams")),foreach::foreach(s = 1:cv,.noexport = c("diagrams"))),ex = {
      
      # use precomputed distance matrices to calculate Gram matrix from parameters
      D <- distance_matrices[[paste0(r[[1]],"_",r[[3]])]]
      
      if(cv > 1)
      {
        # split into training and test set based on cv parameter
        training_indices <- unlist(diagrams_split[setdiff(1:cv,s)])
        training_indices <- training_indices[order(training_indices)]
        test_indices <- setdiff(1:nrow(D),training_indices)
        
        # fit model on training folds
        D_subset <- D[training_indices,training_indices]
        t <- as.numeric(stats::quantile(as.vector(D_subset[upper.tri(D_subset)]),probs = c(r[[2]])))
        if(t == 0)
        {
          non_zero_inds <- which(D_subset != 0,arr.ind = T)
          if(nrow(non_zero_inds) > 0)
          {
            t <- min(D_subset[non_zero_inds])
          }else
          {
            stop("Zero variance found in a cv-fold. Try decreasing the cv parameter or specifying values for t.") 
          }
        }
        t <- 1/t
        K_subset <- exp(-1*t*D_subset)
        class(K_subset) <- "kernelMatrix"
        model = kernlab::ksvm(x = K_subset,y = y[training_indices],type = type,C = r[[4]],nu = r[[5]],epsilon = r[[6]],prob.model = prob.model,class.weights = class.weights,fit = fit,cache = cache,tol = tol,shrinking = shrinking)
        
        # predict model on single test fold
        K_cross <- exp(-1*t*D[test_indices,training_indices[model@SVindex]])
        class(K_cross) <- "kernelMatrix"
        predictions <- kernlab::predict(model,K_cross)
        
        # calculate error depending on model type
        if(type %in% c("C-svc","nu-svc"))
        {
          return((1 - .classAgreement(table(y[test_indices],as.character(predictions)))))
        }
        if(type == "one-svc")
        {
          return((1 - sum(predictions)/length(predictions)))
        }
        if(type %in% c("eps-svr","nu-svr"))
        {
          return(drop(crossprod(predictions - y[test_indices])/nrow(K)))
        }
        
      }else
      {
        t <- as.numeric(stats::quantile(as.vector(D[upper.tri(D)]),probs = c(r[[2]])))
        if(t == 0)
        {
          non_zero_inds <- which(D != 0,arr.ind = T)
          if(nrow(non_zero_inds) > 0)
          {
            t <- min(D[non_zero_inds])
          }else
          {
            stop("Zero variance found in a cv-fold. Try decreasing the cv parameter or specifying values for t.") 
          }
        }
        t <- 1/t
        K <- exp(-1*t*D)
        class(K) <- "kernelMatrix"
        model <- kernlab::ksvm(x = K,y = y,type = type,C = r[[4]],nu = r[[5]],epsilon = r[[6]],prob.model = prob.model,class.weights = class.weights,fit = fit,cache = cache,tol = tol,shrinking = shrinking)
        return(model@error)
      }
    })},
    error = function(e){stop(e)},
    warning = function(w){warning(w)},
    finally = {
      # clean up
      doParallel::stopImplicitCluster()
      parallel::stopCluster(cl)
      
    })
    
    if(cv > 1)
    {
      model_errors <- lapply(X = model_errors,FUN = function(X){
        
        errors <- unlist(lapply(X,"[[",1))
        return(list(mean(unlist(X))))
        
      })
    }
    
    # get mean model errors
    model_errors <- unlist(lapply(model_errors,"[[",1))
  }
  
  # get best parameters
  params$error <- model_errors
  best_params <- params[which.min(model_errors),c("dim","t","sigma","C","nu","epsilon")]
  params <- rbind(params,zero_var_params)
  # if t estimated then get estimate from whole dataset (for predictions)
  if(estimate_t)
  {
    D <- distance_matrices[[paste0(best_params$dim,"_",best_params$sigma)]]^2
    t <- as.numeric(stats::quantile(as.vector(D[upper.tri(D)]),probs = c(best_params$t)))
    if(t == 0)
    {
      t <- min(D[which(D != 0,arr.ind = T)])
    }
    t <- 1/t
    best_params$t <- t
  }
  
  # fit full model using those parameters
  K <- exp(-1*best_params[[2]]*distance_matrices[[paste0(best_params[[1]],"_",best_params[[3]])]])
  class(K) <- "kernelMatrix"
  
  # kernlab calculation, suppressing unhelpful warnings
  tryCatch(expr = {best_model <- list(model = kernlab::ksvm(x = K,y = y,type = type,C = best_params[[4]],nu = best_params[[5]],epsilon = best_params[[6]],prob.model = prob.model,class.weights = class.weights,fit = fit,cache = cache,tol = tol,shrinking = shrinking),
                                      dim = best_params[[1]],
                                      sigma = best_params[[3]],
                                      t = best_params[[2]],
                                      rho = rho)},
           error = function(e){stop(e)},
           warning = function(w){
             
             if(grepl(pattern = "closing unused connection",w) == F)
             {
               message(w)
             }
             
           })

  ret_list <- list(cv_results = params,best_model = best_model,diagrams = diagrams)
  
  # set class for easy recognition in prediction function
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
#' @param new_diagrams a list of persistence diagrams which are either the output of a persistent homology calculation like ripsDiag/\code{\link[TDAstats]{calculate_homology}}/\code{\link{PyH}}, or \code{\link{diagram_to_df}}. Only one of `new_diagrams` and `K` need to be supplied.
#' @param model the output of a \code{\link{diagram_ksvm}} function call, of class 'diagram_ksvm'.
#' @param K an optional cross-Gram matrix of the new diagrams and the diagrams in `model`, default NULL. If not NULL then `new_diagrams` does not need to be supplied.
#' @param num_workers the number of cores used for parallel computation, default is one less than the number of cores on the machine.
#' @return a vector containing the output of \code{\link[kernlab]{predict.ksvm}} on the cross Gram matrix of the new diagrams and the support vector diagrams stored in the model.
#' @export
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @seealso \code{\link{diagram_ksvm}} for training a SVM model on a training set of persistence diagrams and labels.
#' @importFrom kernlab predict as.kernelMatrix
#' @importFrom methods is
#' @examples
#'
#' if(require("TDAstats"))
#' {
#'   # create four diagrams
#'   D1 <- TDAstats::calculate_homology(TDAstats::circle2d[sample(1:100,20),],
#'                       dim = 1,threshold = 2)
#'   D2 <- TDAstats::calculate_homology(TDAstats::circle2d[sample(1:100,20),],
#'                       dim = 1,threshold = 2)
#'   D3 <- TDAstats::calculate_homology(TDAstats::sphere3d[sample(1:100,20),],
#'                       dim = 1,threshold = 2)
#'   D4 <- TDAstats::calculate_homology(TDAstats::sphere3d[sample(1:100,20),],
#'                       dim = 1,threshold = 2)
#'   g <- list(D1,D2,D3,D4)
#' 
#'   # create response vector
#'   y <- as.factor(c("circle","circle","sphere","sphere"))
#' 
#'   # fit model without cross validation
#'   model_svm <- diagram_ksvm(diagrams = g,cv = 1,dim = c(0),
#'                             y = y,sigma = c(1),t = c(1),
#'                             num_workers = 2)
#'
#'   # create two new diagrams
#'   D5 <- TDAstats::calculate_homology(TDAstats::circle2d[sample(1:100,20),],
#'                       dim = 1,threshold = 2)
#'   D6 <- TDAstats::calculate_homology(TDAstats::sphere3d[sample(1:100,20),],
#'                       dim = 1,threshold = 2)
#'   g_new <- list(D5,D6)
#' 
#'   # predict with precomputed Gram matrix
#'   K <- gram_matrix(diagrams = g_new,other_diagrams = model_svm$diagrams,
#'                    dim = model_svm$best_model$dim,sigma = model_svm$best_model$sigma,
#'                    t = model_svm$best_model$t,num_workers = 2)
#'   predict_diagram_ksvm(K = K,model = model_svm,num_workers = 2)
#' }

predict_diagram_ksvm <- function(new_diagrams,model,K = NULL,num_workers = parallelly::availableCores(omit = 1)){
  
  # function to predict labels of new diagrams
  # from pre-fit SVM model
  
  # set internal variables to NULL to avoid build issues
  r <- NULL
  X <- NULL
  
  # error check model argument
  if(is.null(model))
  {
    stop("model must be supplied.")
  }
  if(!inherits(model,"diagram_ksvm"))
  {
    stop("model must be the output of a diagram_ksvm function call.")
  }
  
  if(is.null(K))
  {
    # cross Gram matrix not supplied
    
    # error check new_diagrams argument
    check_param("new_diagrams",new_diagrams,min_length = 1)
    new_diagrams <- all_diagrams(diagram_groups = list(new_diagrams),inference = "independence")[[1]]
    
    # compute kernel matrix, storing the value of each kernel computation between the new diagrams and the old ones
    K <- gram_matrix(diagrams = new_diagrams,other_diagrams = model$diagrams[model$best_model$model@SVindex],dim = model$best_model$dim,sigma = model$best_model$sigma,t = model$best_model$t,rho = model$best_model$rho,num_workers = num_workers)
  }else
  {
    # work with precomputed matrix
    new_diagrams <- NULL
    check_matrix(M = K,name = "K",symmetric = F)
    if(ncol(K) != length(model$diagrams))
    {
      stop("K must have the same number of columns as the number of the diagrams used in the model.")
    }
    
    # subset by support vector indices
    class(K) <- "matrix"
    K <- K[,model$best_model$model@SVindex]
    class(K) <- "kernelMatrix"
  }
  
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
