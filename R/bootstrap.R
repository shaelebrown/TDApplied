#### BOOTSTRAPPING FOR PERSISTENCE DIAGRAMS ####
#' Estimate persistence threshold(s) for topological features in a data set using bootstrapping.
#'
#' Bootstrapping is used to find a conservative estimate of a 1-`alpha` percent "confidence interval" around
#' each point in the persistence diagram of the data set, and points whose intervals do not
#' touch the diagonal (birth == death) would be considered "significant" or "real".
#' One threshold is computed for each dimension in the diagram.
#' 
#' The thresholds are then determined by calculating the 1-`alpha'` percentile of the bottleneck
#' distance values between the real persistence diagram and other diagrams obtained
#' by bootstrap resampling the data. Since `ripsDiag` is the slowest homology engine but is the
#' only engine which calculates representative cycles (as opposed to co-cycles with `PyH`), two
#' homology engines are input to this function - one to calculate the actual persistence diagram, `FUN_diag`
#' (possibly with representative (co)cycles) and one to calculate the bootstrap diagrams, `FUN_boot` (this should be
#' a faster engine, like `calculate_homology` or `PyH`).
#' p-values can be calculated for any feature which survives the thresholding if both `return_subsetted` and `return_pvals` are `TRUE`, 
#' however these values may be larger than the original `alpha` value in some cases. Note that this is not part of the original bootstrap procedure.
#' If stricter thresholding is desired,
#' or the p-values must be less than `alpha`, set `p_less_than_alpha` to `TRUE`. The minimum
#' possible p-value is always 1/(`num_samples` + 1).
#' Note that since \code{\link[TDAstats]{calculate_homology}} 
#' can ignore the longest-lived cluster, fewer "real" clusters may be found. To avoid this possibility
#' try setting `FUN_diag` equal to 'ripsDiag'. Please note that due to the TDA package no longer being available on CRAN,
#' if `FUN_diag` or `FUN_boot` are 'ripsDiag' then `bootstrap_persistence_thresholds` will look for the ripsDiag function in the global environment, 
#' so the TDA package should be attached with `library("TDA")` prior to use.
#'
#' @param X the input dataset, must either be a matrix or data frame.
#' @param FUN_diag a string representing the persistent homology function to use for calculating the full persistence diagram, either
#' 'calculate_homology' (the default), 'PyH' or 'ripsDiag'.
#' @param FUN_boot a string representing the persistent homology function to use for calculating the bootstrapped persistence diagrams, either
#' 'calculate_homology' (the default), 'PyH' or 'ripsDiag'.
#' @param maxdim the integer maximum homological dimension for persistent homology, default 0.
#' @param thresh the positive numeric maximum radius of the Vietoris-Rips filtration.
#' @param distance_mat a boolean representing if `X` is a distance matrix (TRUE) or not (FALSE, default).
#' dimensions together (TRUE, the default) or if one threshold should be calculated for each dimension separately (FALSE).
#' @param ripser the imported ripser module when `FUN_diag` or `FUN_boot` is `PyH`.
#' @param ignore_infinite_cluster a boolean indicating whether or not to ignore the infinitely lived cluster when `FUN_diag` or `FUN_boot` is `PyH`.
#' @param calculate_representatives a boolean representing whether to calculate representative (co)cycles, default FALSE. Note that representatives cant be
#' calculated when using the 'calculate_homology' function.
#' @param num_samples the positive integer number of bootstrap samples, default 30.
#' @param alpha the type-1 error threshold, default 0.05.
#' @param return_diag a boolean representing whether or not to return the calculated persistence diagram, default TRUE.
#' @param return_subsetted a boolean representing whether or not to return the subsetted persistence diagram (with or without representatives), default FALSE.
#' @param return_pvals a boolean representing whether or not to return p-values for features in the subsetted diagram, default FALSE.
#' @param p_less_than_alpha a boolean representing whether or not subset further and return only feature whose p-values are strictly less than `alpha`, default `FALSE`. Note that this is not part of the original bootstrap procedure.
#' @param num_workers the integer number of cores used for parallelizing (over bootstrap samples), default one less the maximum amount of cores on the machine.
#' @return either a numeric vector of threshold values, with one for each dimension 0..`maxdim` (in that order), or a list containing those thresholds and elements (if desired) 
#' @export
#' @importFrom stats quantile
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom parallel makeCluster stopCluster clusterExport clusterEvalQ
#' @importFrom parallelly availableCores
#' @importFrom doParallel registerDoParallel
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @references
#' Chazal F et al (2017). "Robust Topological Inference: Distance to a Measure and Kernel Distance." \url{https://www.jmlr.org/papers/volume18/15-484/15-484.pdf}.
#' @examples
#'
#' if(require("TDAstats"))
#' {
#'   # create a persistence diagram from a sample of the unit circle
#'   df <- TDAstats::circle2d[sample(1:100,size = 50),]
#'
#'   # calculate persistence thresholds for alpha = 0.05 
#'   # and return the calculated diagram as well as the subsetted diagram
#'   bootstrapped_diagram <- bootstrap_persistence_thresholds(X = df,
#'   maxdim = 1,thresh = 2,num_workers = 2)
#' }

bootstrap_persistence_thresholds <- function(X,FUN_diag = "calculate_homology",FUN_boot = "calculate_homology",maxdim = 0,thresh,distance_mat = FALSE,ripser = NULL,ignore_infinite_cluster = TRUE,calculate_representatives = FALSE,num_samples = 30,alpha = 0.05,return_subsetted = FALSE,return_pvals = FALSE,return_diag = TRUE,num_workers = parallelly::availableCores(omit = 1),p_less_than_alpha = FALSE){
  
  # function for thresholding the points computed in a diagram
  # based on persistence
  
  # error check parameters
  if(is.null(distance_mat))
  {
    stop("distance_mat must not be NULL.")
  }
  if(length(distance_mat) > 1 | !inherits(distance_mat,"logical"))
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
  if(length(return_subsetted) > 1 | !inherits(return_subsetted,"logical"))
  {
    stop("return_subsetted must be a single logical (i.e. T or F).")
  }
  if(is.na(return_subsetted) | is.nan(return_subsetted) )
  {
    stop("return_subsetted must not be NA/NAN.")
  }
  
  if(is.null(return_pvals))
  {
    stop("return_pvals must not be NULL.")
  }
  if(length(return_pvals) > 1 | !inherits(return_pvals,"logical"))
  {
    stop("return_pvals must be a single logical (i.e. T or F).")
  }
  if(is.na(return_pvals) | is.nan(return_pvals) )
  {
    stop("return_pvals must not be NA/NAN.")
  }
  
  if(is.null(p_less_than_alpha))
  {
    stop("p_less_than_alpha must not be NULL.")
  }
  if(length(p_less_than_alpha) > 1 | !inherits(p_less_than_alpha,"logical"))
  {
    stop("p_less_than_alpha must be a single logical (i.e. T or F).")
  }
  if(is.na(p_less_than_alpha) | is.nan(p_less_than_alpha) )
  {
    stop("p_less_than_alpha must not be NA/NAN.")
  }
  
  if(is.null(return_diag))
  {
    stop("return_diag must not be NULL.")
  }
  if(length(return_diag) > 1 | !inherits(return_diag,"logical"))
  {
    stop("return_diag must be a single logical (i.e. T or F).")
  }
  if(is.na(return_diag) | is.nan(return_diag) )
  {
    stop("return_diag must not be NA/NAN.")
  }
  
  check_param(param = maxdim,param_name = "maxdim",numeric = T,whole_numbers = T,multiple = F,finite = T,non_negative = T)
  check_param(param = thresh,param_name = "thresh",numeric = T,whole_numbers = F,multiple = F,finite = T,non_negative = T,positive = T)
  
  if(!inherits(X,"data.frame") & !inherits(X,"matrix"))
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
  if(distance_mat == T & (ncol(X) != nrow(X) | !inherits(X,"matrix")))
  {
    stop("if distance_mat is TRUE then X must be a square matrix.")
  }
  if((inherits(X,"matrix") & !inherits(X[1,1],"numeric")) | (inherits(X,"data.frame") & length(which(unlist(lapply(X,is.numeric)))) < ncol(X)))
  {
    stop("X must have only numeric entries.")
  }
  
  if(is.null(FUN_diag))
  {
    stop("FUN_diag must not be NULL.")
  }
  if(length(FUN_diag) > 1 | !inherits(FUN_diag,"character"))
  {
    stop("FUN_diag must be a single string.")
  }
  if(is.na(FUN_diag) | is.nan(FUN_diag) )
  {
    stop("FUN_diag must not be NA/NAN.")
  }
  if(FUN_diag %in% c("calculate_homology","ripsDiag","PyH") == F)
  {
    stop("FUN_diag must be either \'calculate_homology\', \'PyH\' or \'ripsDiag\'.")
  }
  if(FUN_diag == 'calculate_homology')
  {
    if(requireNamespace("TDAstats",quietly = T) == F)
    {
      stop("To use the \'calculate_homology\' function the package TDAstats must be installed.")
    }
  }
  if(FUN_diag == 'ripsDiag')
  {
    tryCatch(expr = {ripsDiag <- get("ripsDiag",envir = globalenv())},
             error = function(e){
               
               stop("To use the \'ripsDiag\' function the package TDA must be attached in the global environment.")
               
             })
  }
  
  if(FUN_diag == "PyH")
  {
    
    if(is.null(ignore_infinite_cluster))
    {
      stop("ignore_infinite_cluster must not be NULL.")
    }
    if(length(ignore_infinite_cluster) > 1 | !inherits(ignore_infinite_cluster,"logical"))
    {
      stop("ignore_infinite_cluster must be a single logical (i.e. T or F).")
    }
    if(is.na(ignore_infinite_cluster) | is.nan(ignore_infinite_cluster) )
    {
      stop("ignore_infinite_cluster must not be NA/NAN.")
    }
    
    check_ripser(ripser)
    
  }
  
  if(is.null(FUN_boot))
  {
    stop("FUN_boot must not be NULL.")
  }
  if(length(FUN_boot) > 1 | !inherits(FUN_boot,"character"))
  {
    stop("FUN_boot must be a single string.")
  }
  if(is.na(FUN_boot) | is.nan(FUN_boot) )
  {
    stop("FUN_boot must not be NA/NAN.")
  }
  if(FUN_boot %in% c("calculate_homology","ripsDiag","PyH") == F)
  {
    stop("FUN_boot must be either \'calculate_homology\', \'PyH\' or \'ripsDiag\'.")
  }
  if(FUN_boot == 'calculate_homology')
  {
    if(requireNamespace("TDAstats",quietly = T) == F)
    {
      stop("To use the \'calculate_homology\' function the package TDAstats must be installed.")
    }
  }
  if(FUN_boot == 'ripsDiag')
  {
    tryCatch(expr = {ripsDiag <- get("ripsDiag",envir = globalenv())},
             error = function(e){
               
               stop("To use the \'ripsDiag\' function the package TDA must be attached in the global environment.")
               
             })
  }
  
  if(FUN_boot == "PyH")
  {
    
    if(is.null(ignore_infinite_cluster))
    {
      stop("ignore_infinite_cluster must not be NULL.")
    }
    if(length(ignore_infinite_cluster) > 1 | !inherits(ignore_infinite_cluster,"logical"))
    {
      stop("ignore_infinite_cluster must be a single logical (i.e. T or F).")
    }
    if(is.na(ignore_infinite_cluster) | is.nan(ignore_infinite_cluster) )
    {
      stop("ignore_infinite_cluster must not be NA/NAN.")
    }
    
    check_ripser(ripser)
    
  }
  
  if(FUN_diag == "ripsDiag" & FUN_boot == "PyH")
  {
    if(ignore_infinite_cluster == T)
    {
      message("Setting ignore_infinite_cluster to F because FUN_diag is 'ripsDiag' and FUN_boot is 'PyH'.")
      ignore_infinite_cluster <- F
    }
  }
  
  if(is.null(calculate_representatives))
  {
    stop("calculate_representatives must not be NULL.")
  }
  if(length(calculate_representatives) > 1 | !inherits(calculate_representatives,"logical"))
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
  if(FUN_diag == "PyH")
  {
    diag <- PyH(X = X,maxdim = maxdim,thresh = thresh,distance_mat = distance_mat,ripser = ripser,ignore_infinite_cluster = ignore_infinite_cluster,calculate_representatives = calculate_representatives)
    if(calculate_representatives == T)
    {
      representatives = diag$representatives
      diag <- diag$diagram
    }
  }
  if(FUN_diag == "calculate_homology")
  {
    diag <- diagram_to_df(TDAstats::calculate_homology(mat = X,dim = maxdim,threshold = thresh,format = ifelse(test = distance_mat == T,yes = "distmat",no = "cloud")))
    representatives <- NULL
  }
  if(FUN_diag == "ripsDiag")
  {
    diag <- ripsDiag(X = X,maxdimension = maxdim,maxscale = thresh,dist = ifelse(test = distance_mat == F,yes = "euclidean",no = "arbitrary"),library = "dionysus",location = calculate_representatives,printProgress = F)
    if(calculate_representatives == T)
    {
      representatives <- diag$cycleLocation
    }
    diag <- diagram_to_df(diag)
  }
  
  # compute in parallel if FUN != "PyH"
  if(FUN_boot != "PyH")
  {
    cl <- parallel::makeCluster(num_workers)
    doParallel::registerDoParallel(cl)
    parallel::clusterEvalQ(cl,c(library("clue"),library("rdist")))
    parallel::clusterExport(cl,varlist = c("diagram_distance","check_diagram","check_param","diagram_to_df"),envir = environment())
    foreach_func <- foreach::`%dopar%`
  }else
  {
    foreach_func <- foreach::`%do%`
  }
  
  # perform bootstrapping in parallel
  tryCatch(expr = {
    
    bootstrap_values <- foreach_func(foreach::foreach(N = 1:num_samples,.combine = c),ex = {
      
      # sample data points with replacement
      s <- unique(sample(1:nrow(X),size = nrow(X),replace = T))
      if(distance_mat == F)
      {
        X_sample <- X[s,]
      }else
      {
        # also subset columns for a distance matrix
        X_sample <- X[s,s]
      }
      
      # calculate diagram (without representatives)
      if(FUN_boot == "PyH")
      {
        bootstrap_diag <- PyH(X = X_sample,maxdim = maxdim,thresh = thresh,distance_mat = distance_mat,ripser = ripser,ignore_infinite_cluster = T,calculate_representatives = F)
      }
      if(FUN_boot == "calculate_homology")
      {
        bootstrap_diag <- diagram_to_df(TDAstats::calculate_homology(mat = X_sample,dim = maxdim,threshold = thresh,format = ifelse(test = distance_mat == T,yes = "distmat",no = "cloud")))
      }
      if(FUN_boot == "ripsDiag")
      {
        bootstrap_diag <- diagram_to_df(ripsDiag(X = X_sample,maxdimension = maxdim,maxscale = thresh,dist = ifelse(test = distance_mat == F,yes = "euclidean",no = "arbitrary"),library = "dionysus",location = F,printProgress = F))
      }
      
      # return bottleneck distance with original diagram in each dimension
      ret_vec <- list()
      for(d in 0:maxdim)
      {
        ret_vec[[length(ret_vec) + 1]] <- diagram_distance(D1 = diag,D2 = bootstrap_diag,p = Inf,dim = d)
      }
      return(ret_vec)
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
           
    # make sure to close cluster  
    if(FUN_boot != "PyH")
    {
      parallel::stopCluster(cl)
    }
             
  })
  
  # convert distance threshold(s) into persistence threshold(s)
  thresholds <- unlist(lapply(X = 0:maxdim,FUN = function(X){
    
    return(2*stats::quantile(unlist(bootstrap_values[seq(X + 1,length(bootstrap_values),maxdim + 1)]),probs = c(1-alpha))[[1]])
    
  }))
  
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
    # subset in each dimension
    inds <- c()
    for(d in 0:maxdim)
    {
      inds <- c(inds,which(diag$dimension == d & diag$death - diag$birth > thresholds[[d + 1]])) 
    }
    ret_list$subsetted_diag <- diag[inds,]
    
    # deal with representatives
    if(calculate_representatives == T & FUN_diag == "ripsDiag")
    {
      ret_list$subsetted_representatives = representatives[inds]
    }
    if(calculate_representatives == T & FUN_diag == "PyH")
    {
      ret_list$subsetted_representatives = representatives
      if(maxdim > 0)
      {
        for(d in 1:maxdim)
        {
          if(length(which(diag$dimension == d)) > 0)
          {
            min_ind <- min(which(diag$dimension == d))
            max_ind <- max(which(diag$dimension == d))
            ret_list$subsetted_representatives[[d + 1]] <- ret_list$subsetted_representatives[[d + 1]][inds[which(inds %in% c(min_ind:max_ind))] - min_ind + 1]
          }else
          {
            ret_list$subsetted_representatives[[d + 1]] <- list()
          }
        }
      }
    }
    
    if(return_pvals)
    {
      # compute p-values for each remaining feature
      # p-value is (Z + 1)/(N + 1) where N is the number of samples and
      # Z is the number of bootstrap samples which do not have smaller distances to the diagonal
      pvals <- c()
      ret_list$pvals <- pvals
      if(nrow(ret_list$subsetted_diag) > 0)
      {
        pvals <- unlist(lapply(X = 1:nrow(ret_list$subsetted_diag),FUN = function(X){
          
          pers <- ret_list$subsetted_diag[X,3L] - ret_list$subsetted_diag[X,2L]
          d <- ret_list$subsetted_diag[X,1L]
          Z <- length(which(2*unlist(bootstrap_values[seq(d + 1,length(bootstrap_values),maxdim + 1)]) >= pers))
          return((Z + 1)/(num_samples + 1))
          
        }))
        
        ret_list$pvals <- pvals
        
        if(p_less_than_alpha)
        {
          # perform stricter thresholding
          inds <- which(pvals < alpha)
          ret_list$pvals <- ret_list$pvals[inds]
          if(calculate_representatives == T & FUN_diag == "ripsDiag")
          {
            ret_list$subsetted_representatives <- ret_list$subsetted_representatives[inds]
          }
          if(calculate_representatives == T & FUN_diag == "PyH")
          {
            if(maxdim > 0)
            {
              for(d in 1:maxdim)
              {
                min_ind <- min(which(ret_list$subsetted_diag$dimension == d))
                max_ind <- max(which(ret_list$subsetted_diag$dimension == d))
                ret_list$subsetted_representatives[[d + 1]] <- ret_list$subsetted_representatives[[d + 1]][inds[which(inds %in% c(min_ind:max_ind))] - min_ind + 1]
              }
            }
          }
          ret_list$subsetted_diag <- ret_list$subsetted_diag[rownames(ret_list$subsetted_diag)[inds],]
        }
        
      }
      
    }
  }
  
  return(ret_list)
  
}

