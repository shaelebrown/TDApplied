
#' @importFrom stats complete.cases

# error checks for permutation_test function parameters

check_diagram <- function(d,ret){

  # error checks for a diagram d stored as a data frame, and conversion
  if(is.list(d) && length(d) == 1 && names(d) == "diagram" && class(d$diagram) == "diagram")
  {
    # d is the output from a TDA calculation
    d <- diagram_to_df(d)
  }else
  {
    if(class(d) != "data.frame")
    {
      stop("Diagrams must either be the output of a TDA computation or data frame.")
    }
  }

  if(nrow(d) == 0)
  {
    stop("Every diagram must be non-empty.")
  }

  if(ncol(d) != 3)
  {
    stop("Every diagram must have three columns.")
  }

  if(class(d[,1]) != "numeric" | class(d[,2]) != "numeric" | class(d[,3]) != "numeric")
  {
    stop("Diagrams must have numeric columns.")
  }

  if(!all.equal(d[,1],as.integer(d[,1])))
  {
    stop("Homology dimensions must be whole numbers.")
  }

  if(length(which(d[,1] < 0)) > 0)
  {
    stop("Homology dimensions must be >= 0.")
  }

  if(length(which(d[,2] < 0)) > 0 | length(which(d[,3] < 0)) > 0)
  {
    stop("Birth and death radii must be >= 0.")
  }

  if(length(which(stats::complete.cases(d))) != nrow(d))
  {
    stop("Diagrams can't have missing values.")
  }
  
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
      # check to make sure each diagram is actually the output of some TDA computation or a data frame
      if(class(diagram_groups[[g]][[diag]]) != "data.frame" && !(is.list(diagram_groups[[g]][[diag]]) && length(diagram_groups[[g]][[diag]]) == 1 && names(diagram_groups[[g]][[diag]]) == "diagram" && class(diagram_groups[[g]][[diag]]$diagram) == "diagram"))
      {
        stop(paste0("Every diagram must be the output from a homology calculation from TDA or diagram_to_df."))
      }else
      {
        # if of the right form, format into a data frame and store diagram index
        if(class(diagram_groups[[g]][[diag]]) == "data.frame")
        {
          if(inference == "difference")
          {
            diagram_groups[[g]][[diag]] <- list(diag = diagram_groups[[g]][[diag]],ind = csum_group_sizes[g] + diag)
          }
        }else
        {
          if(inference == "difference")
          {
            diagram_groups[[g]][[diag]] <- list(diag = diagram_to_df(diagram_groups[[g]][[diag]]),ind = csum_group_sizes[g] + diag)
          }
        }
      }
      
      # make sure the converted diagram has appropriate attributes for further use
      if(inference == "difference")
      {
        check_diagram(diagram_groups[[g]][[diag]]$diag,ret = F)
      }else
      {
        check_diagram(diagram_groups[[g]][[diag]],ret = F) 
      }
      
    }
  }
  
  # return diagram groups with reformatted diagrams
  return(diagram_groups)
  

}

check_param <- function(param_name,param,numeric = T,multiple = F,whole_numbers = F,finite = T,at_least_one = F,positive = T,min_length = 1){
  
  # check if a single parameter satisfies certain constraints
  if(!is.list(param) & !is.vector(param))
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
    if(!is.list(param) | length(param) < min_length)
    {
      stop(paste0(param_name," must be a list of persistence diagrams of length at least ",min_length,"."))
    }
    return()
  }
  
  if(multiple == F & length(param) > 1)
  {
    stop(paste0(param_name," must be a single value."))
  }
  
  if(param_name == "distance")
  {
    if(is.null(param) | length(param) > 1 | (param %in% c("wasserstein","fisher")) == F)
    {
      stop("distance must either be \'wasserstein\' or \'fisher\'.")
    }
    return()
  }
  
  if(numeric == T)
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
    
    return()
  }
  
  if(numeric == T & whole_numbers == T & length(which(floor(param) != param)) > 0)
  {
    stop("dims must be whole numbers.")
  }
  
  if(finite == T & length(which(!is.finite(param))) > 0)
  {
    stop(paste0(param_name," must be finite."))
  }
  
  if(positive == T & length(which(param < 0)) > 0)
  {
    stop(paste0(param_name," must be positive."))
  }
  
  if(at_least_one == T & length(which(param < 1)) > 0)
  {
    stop(paste0(param_name," must be at least one."))
  }
  
}
