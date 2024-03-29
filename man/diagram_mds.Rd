% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/machine_learning.R
\name{diagram_mds}
\alias{diagram_mds}
\title{Dimension reduction of a group of persistence diagrams via metric multidimensional scaling.}
\usage{
diagram_mds(
  diagrams,
  D = NULL,
  k = 2,
  distance = "wasserstein",
  dim = 0,
  p = 2,
  sigma = NULL,
  rho = NULL,
  eig = FALSE,
  add = FALSE,
  x.ret = FALSE,
  list. = eig || add || x.ret,
  num_workers = parallelly::availableCores(omit = 1)
)
}
\arguments{
\item{diagrams}{a list of n>=2 persistence diagrams which are either the output of a persistent homology calculation like ripsDiag/\code{\link[TDAstats]{calculate_homology}}/\code{\link{PyH}}, or \code{\link{diagram_to_df}}. Only one of `diagrams` and `D` need to be supplied.}

\item{D}{an optional precomputed distance matrix of persistence diagrams, default NULL. If not NULL then `diagrams` parameter does not need to be supplied.}

\item{k}{the dimension of the space which the data are to be represented in; must be in \{1,2,...,n-1\}.}

\item{distance}{a string representing the desired distance metric to be used, either 'wasserstein' (default) or 'fisher'.}

\item{dim}{the non-negative integer homological dimension in which the distance is to be computed, default 0.}

\item{p}{a positive number representing the wasserstein power, a number at least 1 (infinity for the bottleneck distance), default 2.}

\item{sigma}{a positive number representing the bandwidth for the Fisher information metric, default NULL.}

\item{rho}{an optional positive number representing the heuristic for Fisher information metric approximation, see \code{\link{diagram_distance}}. Default NULL. If supplied, distance matrix calculation is sequential.}

\item{eig}{a boolean indicating whether the eigenvalues should be returned.}

\item{add}{a boolean indicating if an additive constant c* should be computed, and added to the non-diagonal dissimilarities such that the modified dissimilarities are Euclidean.}

\item{x.ret}{a boolean indicating whether the doubly centered symmetric distance matrix should be returned.}

\item{list.}{a boolean indicating if a list should be returned or just the n*k matrix.}

\item{num_workers}{the number of cores used for parallel computation, default is one less than the number of cores on the machine.}
}
\value{
the output of \code{\link[stats]{cmdscale}} on the diagram distance matrix. If `list.` is false (as per default),
a matrix with `k` columns whose rows give the coordinates of the points chosen to represent the dissimilarities.

Otherwise, a list containing the following components.

\describe{

 \item{points}{a matrix with `k` columns whose rows give the coordinates of the points chosen to represent the dissimilarities.}

 \item{eig}{the \eqn{n} eigenvalues computed during the scaling process if `eig` is true.}
 
 \item{x}{the doubly centered distance matrix if `x.ret` is true.}
 
 \item{ac}{the additive constant \eqn{c*}, 0 if `add` = FALSE.}
 
 \item{GOF}{the numeric vector of length 2, representing the sum of all the eigenvalues divided by the sum of their absolute values (first vector element) or by the sum of the max of each eigenvalue and 0 (second vector element).}

}
}
\description{
Projects a group of persistence diagrams (or a precomputed distance matrix of diagrams) into a low-dimensional 
embedding space via metric multidimensional scaling. Such a projection can be used for visualization of data, 
or a static analysis of the embedding dimensions.
}
\details{
Returns the output of \code{\link[stats]{cmdscale}} on the desired distance matrix of a group of persistence diagrams
in a particular dimension. If `distance` is "fisher" then `sigma` must not be NULL.
}
\examples{

if(require("TDAstats"))
{
  # create two diagrams
  D1 <- TDAstats::calculate_homology(TDAstats::circle2d[sample(1:100,10),],
                      dim = 1,threshold = 2)
  D2 <- TDAstats::calculate_homology(TDAstats::circle2d[sample(1:100,10),],
                      dim = 1,threshold = 2)
  g <- list(D1,D2)

  # calculate their 1D MDS embedding in dimension 0 with the bottleneck distance
  mds <- diagram_mds(diagrams = g,k = 1,dim = 0,p = Inf,num_workers = 2)
  
  # repeat but with a precomputed distance matrix, gives same result just much faster
  Dmat <- distance_matrix(diagrams = list(D1,D2),dim = 0,p = Inf,num_workers = 2)
  mds <- diagram_mds(D = Dmat,k = 1)
  
}
}
\references{
Cox M and Cox F (2008). "Multidimensional Scaling." \doi{10.1007/978-3-540-33037-0_14}.
}
\author{
Shael Brown - \email{shaelebrown@gmail.com}
}
