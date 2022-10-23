#### BOOTSTRAPPING FOR PERSISTENCE DIAGRAMS ####
#' Estimate persistence threshold(s) for topological features in a data set using bootstrapping.
#'
#' Bootstrapping is used to find a conservative estimate of a "confidence interval" around
#' each point in the persistence diagram of the data set, and points whose intervals do not
#' overlap with the diagonal (birth = death) would be considered "significant" or "real".
#' Either one threshold can be computed across all desired dimensions or one threshold for
#' each dimension.
#' 
#' The thresholds are determined by calculating the 1-alpha percentile of the bottleneck
#' distance values between the real persistence diagram and other diagrams obtained
#' by bootstrap resampling the data. Note that since \code{\link[TDAstats]{calculate_homology}} and
#' \code{\link{PyH}} (when `ignore_infinite_cluster` is TRUE as is the default) can ignore the longest-lived
#' cluster, fewer "real" clusters may be found when using these functions. To avoid this possibility
#' either set `FUN` equal to 'ripsDiag' or to 'PyH' and set `ignore_infinite_cluster` to FALSE.
#'
#' @param X the input dataset, must either be a matrix or data frame.
#' @param FUN a string representing the persistent homology function to use, either
#' 'calculate_homology', 'ripsDiag' or 'PyH' (default).
#' @param maxdim the integer maximum homological dimension for persistent homology, default 0.
#' @param thresh the positive numeric maximum radius of the Vietoris-Rips filtration.
#' @param distance_mat a boolean representing if `X` is a distance matrix (TRUE) or not (FALSE, default).
#' @param global_threshold a boolean representing if one threshold should be calculated for all 
#' dimensions together (TRUE, the default) or if one threshold should be calculated for each dimension separately (FALSE).
#' @param ripser either NULL (default) or the ripser module for `PyH` calculations, as imported by the `import_ripser` function.
#' @param ignore_infinite_cluster a boolean representing whether to ignore the infinitely persisting cluster when using the `PyH` function, default TRUE.
#' @param calculate_representatives a boolean representing whether to calculate representative (co)cycles, default FALSE. Note that representatives cant be
#' calculated when using the 'calculate_homology' function.
#' @param num_samples the positive integer number of bootstrap samples, default 30.
#' @param alpha the type-1 error threshold, default 0.05.
#' @param return_diag a boolean representing whether or not to return the calculated persistence diagram, default TRUE.
#' @param return_subsetted a boolean representing whether or not to return the subsetted persistence diagram (with or without representatives), default TRUE.
#' @param num_workers the integer number of cores used for parallelizing (over bootstrap samples), default one less the maximum amount of cores on the machine.
#' @return either one numeric threshold value (if `global_threshold` is TRUE) or a vector of such values with one for each
#' dimension 0..`maxdim` (in that order).
#' @export
#' @importFrom methods is
#' @importFrom stats quantile
#' @importFrom foreach foreach
#' @importFrom future plan
#' @importFrom doFuture registerDoFuture
#' @importFrom parallel makeCluster stopCluster
#' @importFrom parallelly availableCores
#' @importFrom doRNG %dorng%
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @examples
#'
#' # create a persistence diagram from a sample of the unit circle
#' df = TDA::circleUnif(n = 50)
#'
#' # calculate persistence thresholds for alpha = 0.05 
#' # and return the calculated diagram as well as the subsetted diagram
#' bootstrapped_diagram <- bootstrap_persistence_thresholds(X = df,
#' FUN = "calculate_homology",maxdim = 1,thresh = 2)

bootstrap_persistence_thresholds <- function(X,FUN = "PyH",maxdim = 0,thresh,distance_mat = FALSE,global_threshold = TRUE,ripser = NULL,ignore_infinite_cluster = TRUE,calculate_representatives = FALSE,num_samples = 30,alpha = 0.05,return_subsetted = TRUE,return_diag = TRUE,num_workers = parallelly::availableCores(omit = 1)){

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
  
  if(is.null(return_subsetted))
  {
    stop("return_subsetted must not be NULL.")
  }
  if(length(return_subsetted) > 1 | !methods::is(return_subsetted,"logical"))
  {
    stop("return_subsetted must be a single logical (i.e. T or F).")
  }
  if(is.na(return_subsetted) | is.nan(return_subsetted) )
  {
    stop("return_subsetted must not be NA/NAN.")
  }
  
  if(is.null(return_diag))
  {
    stop("return_diag must not be NULL.")
  }
  if(length(return_diag) > 1 | !methods::is(return_diag,"logical"))
  {
    stop("return_diag must be a single logical (i.e. T or F).")
  }
  if(is.na(return_diag) | is.nan(return_diag) )
  {
    stop("return_diag must not be NA/NAN.")
  }
  
  check_param(param = maxdim,param_name = "maxdim",numeric = T,whole_numbers = T,multiple = F,finite = T,non_negative = T)
  check_param(param = thresh,param_name = "thresh",numeric = T,whole_numbers = F,multiple = F,finite = T,non_negative = T,positive = T)
  if(!is.null(ripser))
  {
    check_ripser(ripser)
  }
  
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
  
  if(is.null(FUN))
  {
    stop("FUN must not be NULL.")
  }
  if(length(FUN) > 1 | !methods::is(FUN,"character"))
  {
    stop("FUN must be a single string.")
  }
  if(is.na(FUN) | is.nan(FUN) )
  {
    stop("FUN must not be NA/NAN.")
  }
  if(FUN %in% c("PyH","calculate_homology","ripsDiag") == F)
  {
    stop("FUN must be either \'PyH\', \'calculate_homology\' or \'ripsDiag\'.")
  }
  if(FUN == 'calculate_homology')
  {
    if(requireNamespace("TDAstats",quietly = T) == F)
    {
      stop("To use the \'calculate_homology\' function the package TDAstats must be installed.")
    }
  }
  if(FUN == 'ripsDiag')
  {
    if(requireNamespace("TDA",quietly = T) == F)
    {
      stop("To use the \'ripsDiag\' function the package TDA must be installed.")
    }
  }
  
  if(is.null(global_threshold))
  {
    stop("global_threshold must not be NULL.")
  }
  if(length(global_threshold) > 1 | !methods::is(global_threshold,"logical"))
  {
    stop("global_threshold must be a single boolean value.")
  }
  if(is.na(global_threshold) | is.nan(global_threshold) )
  {
    stop("global_threshold must not be NA/NAN.")
  }
  
  if(is.null(ignore_infinite_cluster))
  {
    stop("ignore_infinite_cluster must not be NULL.")
  }
  if(length(ignore_infinite_cluster) > 1 | !methods::is(ignore_infinite_cluster,"logical"))
  {
    stop("ignore_infinite_cluster must be a single boolean value.")
  }
  if(is.na(ignore_infinite_cluster) | is.nan(ignore_infinite_cluster) )
  {
    stop("ignore_infinite_cluster must not be NA/NAN.")
  }
  
  if(is.null(calculate_representatives))
  {
    stop("calculate_representatives must not be NULL.")
  }
  if(length(calculate_representatives) > 1 | !methods::is(calculate_representatives,"logical"))
  {
    stop("calculate_representatives must be a single boolean value.")
  }
  if(is.na(calculate_representatives) | is.nan(calculate_representatives) )
  {
    stop("calculate_representatives must not be NA/NAN.")
  }
  check_param(param_name = "num_samples",param = num_samples,numeric = T,whole_numbers = T,multiple = F,finite = T,at_least_one = T)
  check_param(param_name = "alpha",param = alpha,numeric = T,whole_numbers = F,multiple = F,finite = T,non_negative = T)
  if(alpha <= 0 | alpha >= 1)
  {
    stop("alpha must be between 0 and 1 (non-inclusive).")
  }
  check_param("num_workers",num_workers,whole_numbers = T,at_least_one = T,finite = T,numeric = T,multiple = F)
  if(num_workers > parallelly::availableCores())
  {
    warning("num_workers is greater than the number of available cores - setting to maximum value less one.")
    num_workers <- parallelly::availableCores(omit = 1)
  }
  
  # first calculate the "real" persistence diagram, storing representative cycles if need be
  if(FUN == "PyH")
  {
    diag <- PyH(X = X,maxdim = maxdim,thresh = thresh,distance_mat = distance_mat,ripser = ripser,ignore_infinite_cluster = ignore_infinite_cluster,calculate_representatives = calculate_representatives)
    if(calculate_representatives == T)
    {
      representatives <- diag$representatives
      diag <- diagram_to_df(diag$diagram)
    }else
    {
      representatives <- NULL
    }
  }
  if(FUN == "calculate_homology")
  {
    diag <- diagram_to_df(TDAstats::calculate_homology(mat = X,dim = maxdim,threshold = thresh,format = ifelse(test = distance_mat == T,yes = "distmat",no = "cloud")))
    representatives <- NULL
  }
  if(FUN == "ripsDiag")
  {
    diag <- TDA::ripsDiag(X = X,maxdimension = maxdim,maxscale = thresh,dist = ifelse(test = distance_mat == F,yes = "euclidean",no = "arbitrary"),library = "dionysus",location = calculate_representatives,printProgress = F)
    if(calculate_representatives == T)
    {
      representatives <- diag$cycleLocation
    }
    diag <- diagram_to_df(diag)
  }
  
  # compute distance matrix in parallel
  doFuture::registerDoFuture()
  cl = parallel::makeCluster(num_workers)
  future::plan(strategy = "cluster",workers = cl)

  # perform bootstrapping in parallel
  tryCatch(expr = {
    
    bootstrap_values <- doRNG::`%dorng%`(foreach::foreach(N = 1:num_samples,.combine = c),ex = {
      
      # sample data points with replacement
      s <- sample(1:nrow(X),size = nrow(X),replace = T)
      if(distance_mat == F)
      {
        X_sample <- X[s,]
      }else
      {
        # also subset columns for a distance matrix
        X_sample <- X[s,s]
      }
      
      # calculate diagram (without representatives)
      if(FUN == "PyH")
      {
        # calculate persistent homology
        bootstrap_diag <- ripser$ripser(X = X_sample,maxdim = maxdim,thresh = thresh,distance_matrix = distance_mat,do_cocycles = F)
        
        # format output
        bootstrap_diag <- do.call(rbind,lapply(X = 1:(maxdim + 1),FUN = function(X){
          
          d <- as.data.frame(bootstrap_diag$dgms[[X]])
          if(nrow(d) == 0)
          {
            return(data.frame(dimension = numeric(),birth = numeric(),death = numeric()))
          }
          d$D <- X - 1
          d <- d[,c(3,1,2)]
          colnames(d) = c("dimension","birth","death")
          d[which(d$death == Inf),3] <- thresh
          if(ignore_infinite_cluster == T & X == 1)
          {
            d <- d[which(d$birth > 0 | d$death < thresh),]
          }
          return(d)
          
        }))
        bootstrap_diag <- as.data.frame(bootstrap_diag)
        
      }
      if(FUN == "calculate_homology")
      {
        bootstrap_diag <- diagram_to_df(TDAstats::calculate_homology(mat = X_sample,dim = maxdim,threshold = thresh,format = ifelse(test = distance_mat == T,yes = "distmat",no = "cloud")))
      }
      if(FUN == "ripsDiag")
      {
        bootstrap_diag <- diagram_to_df(TDA::ripsDiag(X = X_sample,maxdimension = maxdim,maxscale = thresh,dist = ifelse(test = distance_mat == F,yes = "euclidean",no = "arbitrary"),library = "dionysus",location = F,printProgress = F))
      }
      
      if(global_threshold == T)
      {
        # return bottleneck distance with original diagram considering each point to have the same dimension
        bootstrap_diag$dimension = rep(0,nrow(bootstrap_diag))
        return(diagram_distance(D1 = data.frame(dimension = rep(0,nrow(diag)),birth = diag$birth,death = diag$death),D2 = bootstrap_diag,p = Inf))
      }else
      {
        # return bottleneck distance with original diagram in each dimension
        ret_vec <- list()
        for(d in 0:maxdim)
        {
          ret_vec[[length(ret_vec) + 1]] <- diagram_distance(D1 = diag,D2 = bootstrap_diag,p = Inf,dim = d)
        }
        return(ret_vec)
      }
    })
    
  },
  warning = function(w){
    
    if(grepl(pattern = "Every diagram must be non-empty.",w))
    {
      warning("A diagram in some dimension was empty. Try decreasing the maxdim parameter.")
    }else
    {
      warning(w)
    }},
  error = function(e){stop(e)},
  finally = {
             
    parallel::stopCluster(cl)
             
  })
  
  # convert distance threshold(s) into persistence threshold(s)
  if(length(bootstrap_values) == num_samples)
  {
    thresholds <- 2*stats::quantile(unlist(bootstrap_values),probs = c(1-alpha))[[1]]
  }else
  {
    thresholds <- unlist(lapply(X = 0:maxdim,FUN = function(X){
      
      return(2*stats::quantile(unlist(bootstrap_values[seq(X + 1,length(bootstrap_values),maxdim + 1)]),probs = c(1-alpha))[[1]])
      
    }))
  }
  
  # make return list
  ret_list <- list(thresholds = thresholds)
  if(return_diag == T)
  {
    ret_list$diag <- diag
  }
  if(calculate_representatives == T)
  {
    ret_list$representatives <- representatives
  }
  if(return_subsetted == T)
  {
    if(length(thresholds) == 1)
    {
      # subset globally over all dimensions
      inds <- which(diag$death - diag$birth > thresholds)
      ret_list$subsetted_diag <- diag[inds,]
    }else
    {
      # subset in each dimension
      inds <- c()
      for(d in 0:maxdim)
      {
        inds <- c(inds,which(diag$dimension == d & diag$death - diag$birth > thresholds[[d + 1]]))
      }
      ret_list$subsetted_diag <- diag[inds,]
    }
    if(calculate_representatives == T & FUN == "ripsDiag")
    {
      ret_list$subsetted_representatives = representatives[inds]
    }
    if(calculate_representatives == T & FUN == "PyH")
    {
      ret_list$subsetted_representatives = representatives
      if(maxdim > 0)
      {
        for(d in 1:maxdim)
        {
          min_ind <- min(which(diag$dimension == d))
          max_ind <- max(which(diag$dimension == d))
          ret_list$subsetted_representatives[[d]] <- ret_list$subsetted_representatives[[d]][inds[which(inds %in% c(min_ind:max_ind))] - min_ind + 1]
        }
      }
    }
  }
  
  return(ret_list)
  
}

