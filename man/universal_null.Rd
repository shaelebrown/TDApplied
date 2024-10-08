% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inference.R
\name{universal_null}
\alias{universal_null}
\title{Filtering topological features with the universal null distribution.}
\usage{
universal_null(
  X,
  FUN_diag = "calculate_homology",
  maxdim = 1,
  thresh,
  distance_mat = FALSE,
  ripser = NULL,
  ignore_infinite_cluster = TRUE,
  calculate_representatives = FALSE,
  alpha = 0.05,
  return_pvals = FALSE,
  infinite_cycle_inference = FALSE
)
}
\arguments{
\item{X}{the input dataset, must either be a matrix or data frame.}

\item{FUN_diag}{a string representing the persistent homology function to use for calculating the full persistence diagram, either
'calculate_homology' (the default), 'PyH' or 'ripsDiag'.}

\item{maxdim}{the integer maximum homological dimension for persistent homology, default 0.}

\item{thresh}{the positive numeric maximum radius of the Vietoris-Rips filtration.}

\item{distance_mat}{a boolean representing if `X` is a distance matrix (TRUE) or not (FALSE, default).
dimensions together (TRUE, the default) or if one threshold should be calculated for each dimension separately (FALSE).}

\item{ripser}{the imported ripser module when `FUN_diag` is `PyH`.}

\item{ignore_infinite_cluster}{a boolean indicating whether or not to ignore the infinitely lived cluster when `FUN_diag` is `PyH`.}

\item{calculate_representatives}{a boolean representing whether to calculate representative (co)cycles, default FALSE. Note that representatives cant be
calculated when using the 'calculate_homology' function.}

\item{alpha}{the type-1 error threshold, default 0.05.}

\item{return_pvals}{a boolean representing whether or not to return p-values for features in the subsetted diagram as well as a list of p-value thresholds, default FALSE.}

\item{infinite_cycles}{a boolean representing whether or not to perform inference for features with infinite (i.e. `thresh`) death values, default FALSE.}
}
\value{
a list containing the full persistence diagram, the subsetted diagram, representatives and/or subsetted representatives if desired, the p-values of subsetted features and the Bonferroni p-value thresholds in each dimension if desired.
}
\description{
An inference procedure to determine which topological features (if any) of a datasets are likely signal (i.e. significant)
vs noise (not).
}
\details{
For each feature in a diagram we compute its persistence ratio $\pi = \frac{death radius}{birth radius}$, and a
test statistic $A \log \log \pi + B$ (where $A$ and $B$ are constants). This statistic is compared to a left-skewed Gumbel distribution
to get a p-value. A Bonferroni correction is applied to all the p-values across all features, so when `return_pvals` is TRUE a list of 
p-value thresholds is also returned, one for each dimension, which is `alpha` divided by the number of features in that dimension. 
This function is significantly faster than the \code{\link{bootstrap_persistence_thresholds}} function. Note that the `calculate_homology`
function does not seem to store infinite cycles (i.e. cycles that have death value equal to `thresh`).
}
\examples{

if(require("TDAstats"))
{
  # create two datasets
  D1 <- TDAstats::calculate_homology(TDAstats::circle2d[sample(1:100,10),],
                                     dim = 0,threshold = 2)
  D2 <- TDAstats::calculate_homology(TDAstats::circle2d[sample(1:100,10),],
                                     dim = 0,threshold = 2)

  # do model inference test with 1 iteration (for speed, more
  # iterations should be used in practice)
  model_test <- permutation_model_inference(D1, D2, iterations = 1,
                                            thresh = 1.75,num_samples = 3,
                                            num_workers = 2L)
  
  # with more iterations, p-values show a difference in the 
  # clustering of points but not in the arrangement of loops
  model_test$p_values
  
  # to supply samp, when we believe there is a correspondence between
  # the rows in D1 and the rows in D2
  # note that the number of entries of samp (3 in this case) must
  # match the num_samples parameter to the function call
  samp <- lapply(X = 1:3,FUN = function(X){

           return(unique(sample(1:nrow(D1),size = nrow(D1),replace = TRUE)))

          })
  
  # model inference will theoretically have higher power now for a
  # paired test 
  model_test2 <- permutation_model_inference(D1, D2, iterations = 1,
                                             thresh = 1.75,num_samples = 3,
                                             paired = TRUE,samp = samp,
                                             num_workers = 2L)
  model_test2$p_values
}
}
\references{
Bobrowski O, Skraba P (2023). "A universal null-distribution for topological data analysis." \url{https://www.nature.com/articles/s41598-023-37842-2}.
}
\author{
Shael Brown - \email{shaelebrown@gmail.com}
}
