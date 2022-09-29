#### PYTHON SETUP ####
#' Configure R to interface with python for fast persistent homology calculations.
#'
#' Installs the "reticulate" package, python and the python module "ripser". 
#'
#' The "ripser" python module has faster persistent homology calculations than the R wrapper for ripser, the TDAstats
#' package, even when using the reticulate package as an interface. This function allows
#' access to the added python speedup, and only needs to be called one time after TDApplied has been
#' installed (not at the beginning of every R session).
#' 
#' @importFrom reticulate py_available install_python py_list_packages py_install
#' @importFrom utils install.packages
#' @export
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @examples
#'
#' # set up python for fast PH
#' PyH_setup()
PyH_setup <- function(){
  
  # first make sure that reticulate is installed
  if(requireNamespace("reticulate",quietly = T) == F)
  {
    message("The reticulate package could not be found, trying to download now...")
    utils::install.packages("reticulate")
  }
  
  # then make sure that python is installed
  if(reticulate::py_available() == F)
  {
    message("python has not yet been installed, trying now...")
    reticulate::install_python()
  }
  
  modules <- reticulate::py_list_packages()[,1]
  
  # finally make sure that ripser module has been downloaded
  if("ripser" %in% modules == F)
  {
    message("ripser has not yet been installed, trying now...")
    reticulate::py_install("ripser")
  }
  # 
  # # finally make sure that persim module has been downloaded
  # if("persim" %in% modules == F)
  # {
  #   message("persim has not yet been installed, trying now...")
  #   reticulate::py_install("persim")
  # }
  
}

#### CHECK PYTHON SETUP ####
#' Make sure that python has been configured correctly for persistent homology calculations.
#' 
#' Ensures that the reticulate package, python and the python module "ripser" have all been installed.
#' 
#' @importFrom reticulate py_available py_list_packages
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
check_PyH_setup <- function(){
  
  # first make sure that reticulate is installed
  if(requireNamespace("reticulate",quietly = T) == F)
  {
    stop("reticulate package must be installed. Try using PyH_setup() function.")
  }
  
  # then make sure that python is installed
  if(reticulate::py_available() == F)
  {
    stop("python must be installed. Try using PyH_setup() function.")
  }
  
  # finally make sure that ripser module has been downloaded
  modules <- reticulate::py_list_packages()[,1]
  if("ripser" %in% modules == F)
  {
    stop("ripser must be installed. Try using PyH_setup() function.")
  }
  # 
  # # finally make sure that persim module has been downloaded
  # if("persim" %in% modules == F)
  # {
  #   stop("persim must be installed. Try using PyH_setup() function.")
  # }
  
}

#### IMPORT RIPSER MODULE ####
#' Import the python module ripser for fast persistent homology calculations.
#' 
#' Same as "reticulate::import("ripser")".
#' 
#' @importFrom reticulate import
#' @export
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @examples
#' # set up python for fast PH
#' PyH_setup()
#' 
#' # import ripser
#' ripser <- import_ripser()
import_ripser <- function(){
  
  # check python configuration
  check_PyH_setup()
  
  # if all good, import
  return(reticulate::import("ripser"))
  
}

#### CHECK RIPSER MODULE ####
#' Verify an imported ripser module.
#' 
#' @importFrom methods is
#' @param ripser the ripser module object.
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
check_ripser <- function(ripser){
  
  e <- "ripser object should be created using the import_ripser() function."
  
  # check class
  if((methods::is(ripser,"python.builtin.module")) == F | (methods::is(ripser,"python.builtin.object")) == F)
  {
    stop(e)
  }
  
  # check function names
  if(length(which(names(ripser)%in% c("lower_star_img","Rips","ripser"))) < 3)
  {
    stop(e)
  }
  
}

#### PYTHON PERSISTENT HOMOLOGY ####
#' Fast persistent homology calculations with python.
#'
#' This function is a wrapper of the python wrapper of the ripser engine, but is still faster than using the
#' R package TDAstats (see the TDApplied package vignette for details). 
#' 
#' If `distance_mat` is `TRUE` then `X` must be a square matrix. The `ripser` parameter should be the
#' result of an `import_ripser` function call, but since that function is slow the ripser object should
#' be explicitly created before a PyH function call (see examples).
#' 
#' @param X either a matrix or dataframe, representing either point cloud data or a distance matrix. In either case there
#' must be at least two rows and 1 column.
#' @param maxdim the non-negative integer maximum dimension for persistent homology, default 1.
#' @param thresh the non-negative numeric radius threshold for the Vietoris-Rips filtration.
#' @param distance_mat a boolean representing whether the input X is a distance matrix or not, default FALSE.
#' @param ripser the ripser python module.
#' @param ignore_infinite_clusters a boolean representing whether to remove clusters (0 dimensional cycles) which
#' die at the threshold value. Default is TRUE as this is the default for TDAstats homology calculations, but can be set to
#' FALSE which is the default for python ripser. 
#' @importFrom methods is
#' @export
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @examples
#'
#' # set up python for fast PH
#' PyH_setup()
#' 
#' # create sample data
#' df <- data.frame(x = 1:10,y = 1:10)
#' 
#' # import the ripser module
#' ripser <- import_ripser()
#' 
#' # calculate persistent homology up to dimension 1 with a maximum
#' # radius of 5
#' phom <- PyH(X = df,thresh = 5,ripser = ripser)
PyH <- function(X,maxdim = 1,thresh,distance_mat = FALSE,ripser,ignore_infinite_clusters = TRUE){
  
  # error check parameters
  if(is.null(distance_mat))
  {
    stop("distance_mat must not be NULL.")
  }
  if(length(distance_mat) > 1 | !methods::is(distance_mat,"logical"))
  {
    stop("distance_mat must be a single logical (i.e. T or F).")
  }
  if(is.na(distance_mat) | is.nan(distance_mat) )
  {
    stop("distance_mat must not be NA/NAN.")
  }
  
  if(is.null(ignore_infinite_clusters))
  {
    stop("ignore_infinite_clusters must not be NULL.")
  }
  if(length(ignore_infinite_clusters) > 1 | !methods::is(ignore_infinite_clusters,"logical"))
  {
    stop("ignore_infinite_clusters must be a single logical (i.e. T or F).")
  }
  if(is.na(ignore_infinite_clusters) | is.nan(ignore_infinite_clusters) )
  {
    stop("ignore_infinite_clusters must not be NA/NAN.")
  }
  
  check_param(param = maxdim,param_name = "maxdim",numeric = T,whole_numbers = T,multiple = F,finite = T,non_negative = T)
  check_param(param = thresh,param_name = "thresh",numeric = T,whole_numbers = T,multiple = F,finite = T,non_negative = T)
  check_ripser(ripser)
  if(!methods::is(X,"data.frame") & !methods::is(X,"matrix"))
  {
    stop("X must either be a dataframe or a matrix.")
  }
  if(nrow(X) < 2 | ncol(X) < 1)
  {
    stop("X must have at least two rows and one column.")
  }
  if(length(which(stats::complete.cases(X) == F)) > 0)
  {
    stop("X must not contain any missing values.")
  }
  if(distance_mat == T & (ncol(X) != nrow(X) | !methods::is(X,"matrix")))
  {
    stop("if distance_mat is TRUE then X must be a square matrix.")
  }
  if((methods::is(X,"matrix") & !methods::is(X[1,1],"numeric")) | (methods::is(X,"data.frame") & length(which(unlist(lapply(X,is.numeric)))) < ncol(X)))
  {
    stop("X must have only numeric entries.")
  }
  
  # calculate persistent homology
  PH <- ripser$ripser(X = X,maxdim = maxdim,thresh = thresh,distance_matrix = distance_mat)
  
  # format output
  PH <- do.call(rbind,lapply(X = 1:(maxdim + 1),FUN = function(X){
    
    d <- as.data.frame(PH$dgms[[X]])
    if(nrow(d) == 0)
    {
      return(data.frame(dimension = numeric(),birth = numeric(),death = numeric()))
    }
    d$D <- X - 1
    d <- d[,c(3,1,2)]
    colnames(d) = c("dimension","birth","death")
    d[which(d$death == Inf),3] <- thresh
    if(ignore_infinite_clusters == T & X == 1)
    {
      d <- d[which(d$birth > 0 | d$death < thresh),]
    }
    return(d)
    
  }))
  
  PH <- as.data.frame(PH)
  
  return(PH)
  
}

# old functions which can use python to calculate distances...

#### IMPORT PERSIM MODULE ####
#' Import the python persim module for fast persistent homology distance calculations.
#' 
#' @importFrom reticulate import
# @export
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @examples
#' # set up python for fast PH
#' PyH_setup()
#' 
#' # import persim
#' persim <- import_persim()
# import_persim <- function(){
#   
#   # check python configuration
#   check_PyH_setup()
#   
#   # if all good, import
#   return(reticulate::import("persim"))
#   
# }

#### CHECK PERSIM MODULE ####
#' Verify an imported persim module.
#' 
#' @importFrom methods is
#' @param ripser the persim module object.
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
# check_persim <- function(persim){
#   
#   e <- "persim object should be created using the import_persim() function."
#   
#   # check class
#   if((methods::is(persim,"python.builtin.module")) == F | (methods::is(persim,"python.builtin.object")) == F)
#   {
#     stop(e)
#   }
#   
#   # check function names
#   if(length(which(names(persim)%in% c("bottleneck","bottleneck_matching","gromov_hausdorff","heat","images","images_kernels","images_weights","PersImage","PersistenceImager","plot_diagrams","sliced_wasserstein","visuals","wasserstein","wasserstein_matching"))) < 14)
#   {
#     stop(e)
#   }
#   
# }

#### CONVERT DF TO PERSIM INPUT ####
#' 
#' @param ripser the data frame persistence diagram.
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
# convert_df_to_persim_input <- function(df){
#   
#   # check diagram first
#   check_diagram(df,ret = F)
#   
#   maxdim <- max(df$dimension)
#   
#   return(lapply(X = 0:maxdim,FUN = function(X){
#     
#     d = df[which(df$dimension == X),2:3]
#     if(nrow(d) == 0)
#     {
#       return(data.frame(dimension = numeric(),birth = numeric(),death = numeric()))
#     }
#     colnames(d) <- NULL
#     return(as.matrix(d))
#     
#   }))
#   
# }

#### PYTHON PERSISTENT HOMOLOGY DISTANCES ####
#' Fast persistent homology distance calculations with python.
#'
#' This function is a wrapper of the python module persim, which provides a very fast wasserstein and bottleneck
#' distance calculation.
#' 
#' The persim bottleneck function behaves the same as the diagram_distance bottleneck function, but the wasserstein
#' functions are slightly different. The formulas being used are different.
#' 
#' @param D1 the first persistence diagram, either computed from a TDA/TDAstats function like \code{\link[TDA]{ripsDiag}}/\code{\link[TDAstats]{calculate_homology}}, or such an object converted to a data frame with \code{\link{diagram_to_df}}.
#' @param D2 the second persistence diagram, either computed from TDA/TDAstats function like \code{\link[TDA]{ripsDiag}}/\code{\link[TDAstats]{calculate_homology}}, or such an object converted to a data frame with \code{\link{diagram_to_df}}.
#' @param dim the non-negative integer homological dimension in which to calculate distances, default 0.
#' @param p a number representing the wasserstein power parameter, at least 1 and default 2.
#' @param ripser the ripser python module.
#' @importFrom methods is
# @export
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @examples
#'
#' # set up python for fast PH
#' PyH_setup()
#' 
#' # create sample data
#' df <- data.frame(x = 1:10,y = 1:10)
#' 
#' # import the ripser module
#' ripser <- import_ripser()
#' 
#' # calculate persistent homology up to dimension 1 with a maximum
#' # radius of 5
#' phom <- PyH(X = df,thresh = 5,ripser = ripser)
# py_diagram_distance <- function(D1,D2,dim = 0,p = 2,persim){
#   
#   # error check persim parameter
#   check_persim(persim = persim)
#   
#   # convert diagrams to persim input
#   D1 <- convert_df_to_persim_input(D1)
#   D2 <- convert_df_to_persim_input(D2)
#   
#   # calculate optimal matching
#   dist_and_match <- persim$bottleneck(D1[[dim + 1]],D2[[dim + 1]],matching = TRUE)
#   
#   # return correct distance based on p
#   if(is.infinite(p))
#   {
#     return(dist_and_match[[1]])
#   }
#   return((sum(dist_and_match[[2]][,3]^p))^(1/p))
#   
# }


