
#' @importFrom stats complete.cases

# error checks for function parameters

check_diagram <- function(d,ret){

  # error checks for a diagram d stored as a data frame, and conversion
  if((is.list(d) && ((length(d) == 1 && all(names(d) %in% "diagram") && (inherits(d$diagram,"diagram")) || inherits(d$diagram,"data.frame")) || ((length(d) == 4 && all(names(d) %in% c("diagram","birthLocation","deathLocation","cycleLocation")) && inherits(d$diagram,"diagram"))))) || (inherits(d,"matrix") && inherits(d,"array") & all(colnames(d) %in% c("dimension","birth","death"))))
  {
    # d is the output from a TDA/TDAstats calculation
    d <- diagram_to_df(d)
  }else
  {
    if(!inherits(d,"data.frame"))
    {
      stop("Diagrams must either be the output of a TDA/TDAstats computation or a data frame.")
    }
  }

  # if(nrow(d) == 0)
  # {
  #   stop("Every diagram must be non-empty.")
  # }

  if(ncol(d) != 3)
  {
    stop("Every diagram must have three columns.")
  }

  if(!is.numeric(d[,1L]) | !is.numeric(d[,2L]) | !is.numeric(d[,3L]))
  {
    stop("Diagrams must have numeric columns.")
  }

  if(length(which(d[,1L] != round(d[,1L]))) > 0)
  {
    stop("Homology dimensions must be whole numbers.")
  }

  if(length(which(d[,1L] < 0)) > 0)
  {
    stop("Homology dimensions must be >= 0.")
  }

  if(length(which(d[,2L] < 0)) > 0 | length(which(d[,3L] < 0)) > 0)
  {
    stop("Birth and death radii must be >= 0.")
  }

  if(length(which(stats::complete.cases(d) == T)) != nrow(d))
  {
    stop("Diagrams can't have missing values.")
  }
  
  if(length(which(d[,3L] <= d[,2L])) > 0)
  {
    stop("Death values must always be larger than birth values.")
  }
  
  # if(length(which(is.infinite(d[,3L]))) > 0)
  # {
  #   stop("Death values must be finite.")
  # }
  
  if(ret == T)
  {
    return(d) 
  }

}

all_diagrams <- function(diagram_groups,inference){

  # function to make sure all diagram groups are lists or vectors of diagrams,
  # to convert the diagrams to data frames and to error check each diagram.
  # diagram_groups is a vector or list of vectors or lists of diagrams
  # inference is a string, either 'difference' for the permutation test or
  # 'independence' for the independence test
  
  if(inference == "difference")
  {
    # compute cumulative sums of groups lengths in order to correctly compute diagram indices
    csum_group_sizes <- cumsum(unlist(lapply(diagram_groups,FUN = length)))
    csum_group_sizes <- c(0,csum_group_sizes)
  }
  
  # loop through all diagram groups
  for(g in 1:length(diagram_groups))
  {
    # loop through each diagram in each group
    for(diag in 1:length(diagram_groups[[g]]))
    {
      # check to make sure each diagram is actually the output of some TDA/TDAstats computation or a data frame
      check_diagram(diagram_groups[[g]][[diag]],ret = F) # will stop here if there is an error
      diagram_groups[[g]][[diag]] <- check_diagram(diagram_groups[[g]][[diag]],ret = T)
      # if of the right form, format into a data frame and store diagram index for difference inference
      if(inference == "difference")
      {
        diagram_groups[[g]][[diag]] <- list(diag = diagram_groups[[g]][[diag]],
                                            ind = csum_group_sizes[g] + diag)
      }
    }
  }
  
  # return diagram groups with reformatted diagrams
  return(diagram_groups)
  
}

check_param <- function(param_name,param,...){
  
  extra_params <- list(...)
  
  # check if a single parameter satisfies certain constraints
  if(!is.list(param) & (!is.vector(param) | length(param) == 1))
  {
    if(is.null(param))
    {
      stop(paste0(param_name," must not be NULL."))
    }
    if(is.na(param))
    {
      stop(paste0(param_name," must not be NA/NaN."))
    }
  }
  
  if(param_name == "diagrams" | param_name == "other_diagrams" | param_name == "diagram_groups" | param_name == "new_diagrams")
  {
    if(!is.list(param) | length(param) < extra_params$min_length)
    {
      stop(paste0(param_name," must be a list of persistence diagrams of length at least ",extra_params$min_length,"."))
    }
  }
  
  if("multiple" %in% names(extra_params))
  {
    if(extra_params$multiple == F & length(param) > 1)
    {
      stop(paste0(param_name," must be a single value."))
    }
  }
  
  if(param_name == "distance")
  {
    if(is.null(param) | length(param) > 1 | (param %in% c("wasserstein","fisher")) == F)
    {
      stop("distance must either be \'wasserstein\' or \'fisher\'.")
    }
  }
  
  if("numeric" %in% names(extra_params))
  {
    if(extra_params$numeric == T)
    {
      if(is.numeric(param) == F)
      {
        stop(paste0(param_name," must be numeric."))
      }
    }else
    {
      if(is.logical(param) == F)
      {
        stop(paste0(param_name," must be T or F."))
      }
      
    }
  }
  
  if("numeric" %in% names(extra_params) & "whole_numbers" %in% names(extra_params))
  {
    if(extra_params$numeric == T & extra_params$whole_numbers == T & length(which(floor(param) != param)) > 0)
    {
      stop(paste0(param_name," must be whole numbers."))
    }
  }
  
  if("finite" %in% names(extra_params))
  {
    if(extra_params$finite == T & length(which(!is.finite(param))) > 0)
    {
      stop(paste0(param_name," must be finite."))
    }
  }
  
  if("non_negative" %in% names(extra_params))
  {
    if(extra_params$non_negative == T & length(which(param < 0)) > 0)
    {
      stop(paste0(param_name," must be non-negative"))
    }
  }
  
  if("positive" %in% names(extra_params))
  {
    if(extra_params$positive == T & length(which(param <= 0)) > 0)
    {
      stop(paste0(param_name," must be positive."))
    }
  }
  
  if("at_least_one" %in% names(extra_params))
  {
    if(extra_params$at_least_one == T & length(which(param < 1)) > 0)
    {
      stop(paste0(param_name," must be at least one."))
    }
  }
  
}

#' @importFrom stats complete.cases

check_matrix <- function(M,name,type = "kernel",symmetric = T){
  
  if(type == "kernel")
  {
    if(inherits(M,"kernelMatrix") == F)
    {
      stop(paste0(name," must be of type kernelMatrix."))
    }
  }else
  {
    if(inherits(M,"matrix") == F)
    {
      stop(paste0(name," must be of type matrix."))
    }
  }
  
  if(nrow(M)*ncol(M) == 0)
  {
    stop(paste0(name," must have at least one row and column."))
  }

  if(length(which(stats::complete.cases(M) == F)) > 0)
  {
    stop(paste0(name," must not have missing values."))
  }
  
  if(symmetric == T)
  {
    if(type == "kernel")
    {
      if(length(which(diag(M) != rep(1,ncol(M)))) > 0)
      {
        stop(paste0(name," must have 1's on its diagonal."))
      }
    }else
    {
      if(length(which(diag(M) != rep(0,ncol(M)))) > 0)
      {
        stop(paste0(name," must have 0's on its diagonal."))
      }
    }
    
    if(nrow(M) != ncol(M))
    {
      stop(paste0(name," must have the same number of rows and columns."))
    }
    
    if(type == "kernel")
    {
      class(M) = "matrix"
    }
    if(isSymmetric(M) == F)
    {
      stop(paste0(name," must be symmetric."))
    }
    if(type == "kernel")
    {
      class(M) = "kernelMatrix"
    } 
  }
  
}

# function for cv in diagram_ksvm
# copied from kernlab package v0.9-29
.classAgreement <- function (tab) {
  n <- sum(tab)
  if (!is.null(dimnames(tab))) {
    lev <- intersect(colnames(tab), rownames(tab))
    p0 <- sum(diag(tab[lev, lev])) / n
  } else {
    m <- min(dim(tab))
    p0 <- sum(diag(tab[1:m, 1:m])) / n
  }
  return(p0)
}


