
#' @importFrom stats complete.cases

# error checks for permutation_test function parameters

check_diagram <- function(d){

  # error checks for a diagram d stored as a data frame

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

}

all_diagrams <- function(diagram_groups){

  # function to make sure all diagram groups are lists or vectors of diagrams,
  # to convert the diagrams to data frames and to error check each diagram.
  # diagram_groups is a vector or list of vectors or lists of diagrams

  # compute cumulative sums of groups lengths in order to correctly compute diagram indices
  csum_group_sizes <- cumsum(unlist(lapply(diagram_groups,FUN = length)))
  csum_group_sizes <- c(0,csum_group_sizes)

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
          diagram_groups[[g]][[diag]] <- list(diag = diagram_groups[[g]][[diag]],ind = csum_group_sizes[g] + diag)
        }else
        {
          diagram_groups[[g]][[diag]] <- list(diag = diagram_to_df(diagram_groups[[g]][[diag]]),ind = csum_group_sizes[g] + diag)
        }
      }

      # make sure the converted diagram has appropriate attributes for further use
      check_diagram(diagram_groups[[g]][[diag]]$diag)

    }
  }

  # return diagram groups with reformatted diagrams
  return(diagram_groups)

}

check_params <- function(iterations,p,q,dims,paired,distance){

  # error checks on the parameters for the permutation_test function

  if(!is.numeric(iterations))
  {
    stop("Iterations must be numeric.")
  }

  if(iterations != as.integer(iterations))
  {
    stop("Iterations must be a whole number.")
  }

  if(iterations <= 1 | is.infinite(iterations))
  {
    stop("Iterations must be a finite number at least 1.")
  }

  if(!as.numeric(p))
  {
    stop("p must be numeric.")
  }

  if(p < 1)
  {
    stop("p must be at least 1.")
  }

  if(!is.numeric(q))
  {
    stop("q must be numeric.")
  }

  if(q < 1 | is.infinite(q))
  {
    stop("q must be a finite number at least 1.")
  }

  if(!is.numeric(dims))
  {
    stop("dims must be a numeric vector or value.")
  }

  if(!all.equal(dims,as.integer(dims)))
  {
    stop("dims must be whole numbers.")
  }

  if(length(which(dims < 0 | is.infinite(dims))) > 0)
  {
    stop("dims must be positive and finite.")
  }

  if(!is.logical(paired))
  {
    stop("paired must be T or F.")
  }

  if(!is.character(distance) | distance %in% c("wasserstein","Turner") == F)
  {
    stop("distance must be a single character, either wasserstein or Turner.")
  }

}
