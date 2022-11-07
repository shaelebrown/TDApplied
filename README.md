
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TDApplied

<!-- badges: start -->
<!-- badges: end -->

TDApplied is an R package for applied topological data analysis using
machine learning and statistical inference, and uses the output of
persistent homology calculations from the R packages TDA/TDAstats as
input to its methods.

R package TDA:

> Fasy, Brittany T., Jisu Kim, Fabrizio Lecci, Clement Maria, David L.
> Millman, and Vincent Rouvreau. 2021. TDA: Statistical Tools for
> Topological Data Analysis. <https://CRAN.R-project.org/package=TDA>.

R package TDAstats:

> Wadhwa, Raoul R., Drew R. K. Williamson, Andrew Dhawan, and Jacob G.
> Scott. 2018. TDAstats: R pipeline for computing persistent homology in
> topological data analysis.
> <https://CRAN.R-project.org/package=TDAstats>.

## Installation

To install the latest version of this R package directly from github:

    install.packages("devtools")
    library(devtools)
    devtools::install_github("shaelebrown/TDApplied")
    library(TDApplied)

To install from Github you might need:

-   **Windows:** Rtools
    (<https://cran.r-project.org/bin/windows/Rtools/>)
-   **OS X:** xcode (from the app store)
-   **Linux:** apt-get install r-base-dev (or similar).

To install the stable version of this R package from CRAN:

    install.packages("TDApplied")

## Examples

These are basic examples which show you how to use the package:

``` r
library(TDApplied)
```

For these examples we will use three base persistence diagrams:

``` r
D1 = data.frame(dimension = c(0),birth = c(2),death = c(3))
D2 = data.frame(dimension = c(0),birth = c(2,0),death = c(3.3,0.5))
D3 = data.frame(dimension = c(0),birth = c(0),death = c(0.5))
```

Plotting a diagram:

``` r
plot_diagram(D1,title = "D1")
```

Computing distances between persistence diagrams:

``` r
# calculate 2-wasserstein distance between D1 and D2
diagram_distance(D1,D2,dim = 0,p = 2,distance = "wasserstein")
#> [1] 0.3905125

# calculate bottleneck distance between D1 and D3
diagram_distance(D1,D3,dim = 0,p = Inf,distance = "wasserstein")
#> [1] 0.5

# Fisher information metric calculation between D1 and D2 for sigma = 1
diagram_distance(D1,D2,dim = 0,distance = "fisher",sigma = 1)
#> [1] 0.02354779

# Fisher information metric calculation between D1 and D3 for sigma = 2
diagram_distance(D1,D3,dim = 0,distance = "fisher",sigma = 2)
#> [1] 0.01485812
```

Computing kernel values between persistence diagrams:

``` r
# calculate the kernel value between D1 and D2 with sigma = 2, t = 2
diagram_kernel(D1,D2,dim = 0,sigma = 2,t = 2)
#> [1] 0.9872455
# calculate the kernel value between D1 and D3 with sigma = 2, t = 2
diagram_kernel(D1,D3,dim = 0,sigma = 2,t = 2)
#> [1] 0.9707209
```

Computing a MDS projection of persistence diagrams:

``` r
# create three diagrams
D1 <- TDA::ripsDiag(TDA::circleUnif(n = 20,r = 1),
                   maxdimension = 1,maxscale = 2)
D2 <- TDA::ripsDiag(TDA::sphereUnif(n = 20,d = 2,r = 1),
                   maxdimension = 1,maxscale = 2)
D3 <- TDA::ripsDiag(TDA::torusUnif(n = 20,a = 0.25,c = 0.75),
                   maxdimension = 1,maxscale = 2)
g <- list(D1,D2,D3)

# calculate their 2D MDS embedding in dimension 1 with the bottleneck distance
mds <- diagram_mds(diagrams = g,dim = 0,p = Inf,k = 2,num_workers = 2)
```

Looking for group differences in groups of persistence diagrams:

``` r
# create two groups of diagrams
D1 <- TDA::ripsDiag(TDA::circleUnif(n = 20,r = 1),
                    maxdimension = 1,maxscale = 2)
D2 <- TDA::ripsDiag(TDA::circleUnif(n = 20,r = 1),
                    maxdimension = 1,maxscale = 2)
g1 <- list(D1,D2)
g2 <- list(D1,D2)
perm_test <- permutation_test(g1,g2,
                              num_workers = 2,
                              dims = c(0))
```

Clustering persistence diagrams with kernel k-means:

``` r
# create three diagrams
D1 <- TDA::ripsDiag(TDA::circleUnif(n = 20,r = 1),
                    maxdimension = 1,maxscale = 2)
D2 <- TDA::ripsDiag(TDA::sphereUnif(n = 20,d = 2,r = 1),
                    maxdimension = 1,maxscale = 2)
D3 <- TDA::ripsDiag(TDA::torusUnif(n = 20,a = 0.25,c = 0.75),
                    maxdimension = 1,maxscale = 2)
g <- list(D1,D1,D1,D2,D2,D2,D3,D3,D3)

# calculate kmeans clusters with centers = 3, and sigma = t = 2
clust <- diagram_kkmeans(diagrams = g,centers = 3,dim = 0,t = 2,sigma = 2,num_workers = 2)
```

Predicting new cluster labels:

``` r
# create three new diagrams
D4 <- TDA::ripsDiag(TDA::circleUnif(n = 20,r = 1),
                    maxdimension = 1,maxscale = 2)
D5 <- TDA::ripsDiag(TDA::sphereUnif(n = 20,d = 2,r = 1),
                    maxdimension = 1,maxscale = 2)
D6 <- TDA::ripsDiag(TDA::torusUnif(n = 20,a = 0.25,c = 0.75),
                    maxdimension = 1,maxscale = 2)
g_new <- list(D4,D5,D6)

# predict cluster labels
predict_diagram_kkmeans(new_diagrams = g_new,clustering = clust,num_workers = 2)
#> [1] 1 2 3
```

Computing a kernel PCA embedding of persistence diagrams:

``` r
# create three diagrams
D1 <- TDA::ripsDiag(TDA::circleUnif(n = 20,r = 1),
                    maxdimension = 1,maxscale = 2)
D2 <- TDA::ripsDiag(TDA::sphereUnif(n = 20,d = 2,r = 1),
                    maxdimension = 1,maxscale = 2)
D3 <- TDA::ripsDiag(TDA::torusUnif(n = 20,a = 0.25,c = 0.75),
                    maxdimension = 1,maxscale = 2)
g <- list(D1,D2,D3)

# calculate their 2D PCA embedding with sigma = t = 2
pca <- diagram_kpca(diagrams = g,dim = 0,t = 2,sigma = 2,features = 2,num_workers = 2)
```

Project new persistence diagrams into a kernel PCA embedding:

``` r
# project new diagrams onto old model
D4 <- TDA::ripsDiag(TDA::circleUnif(n = 20,r = 1),
                    maxdimension = 1,maxscale = 2)
D5 <- TDA::ripsDiag(TDA::sphereUnif(n = 20,d = 2,r = 1),
                    maxdimension = 1,maxscale = 2)
g_new <- list(D4,D5)

# predict cluster labels
new_pca <- predict_diagram_kpca(new_diagrams = g_new,embedding = pca,num_workers = 2)
```

Fit a kernel SVM model on persistence diagrams:

``` r
# create four diagrams
D1 <- TDA::ripsDiag(TDA::circleUnif(n = 20,r = 1),
                    maxdimension = 1,maxscale = 2)
D2 <- TDA::ripsDiag(TDA::sphereUnif(n = 20,d = 2,r = 1),
                    maxdimension = 1,maxscale = 2)
D3 <- TDA::ripsDiag(TDA::circleUnif(n = 20,r = 1),
                    maxdimension = 1,maxscale = 2)
D4 <- TDA::ripsDiag(TDA::sphereUnif(n = 20,d = 2,r = 1),
                    maxdimension = 1,maxscale = 2)
g <- list(D1,D2,D3,D4)

# create response vector
y <- as.factor(c("circle","sphere","circle","sphere"))

# fit model without cross validation
model_svm <- diagram_ksvm(diagrams = g,cv = 1,dim = c(0),
                          y = y,sigma = c(1),t = c(1),
                          num_workers = 2)
```

Predict labels for new persistence diagrams:

``` r
# create new diagrams
D5 <- TDA::ripsDiag(TDA::circleUnif(n = 20,r = 1),
                    maxdimension = 1,maxscale = 2)
D6 <- TDA::ripsDiag(TDA::sphereUnif(n = 20,d = 2,r = 1),
                    maxdimension = 1,maxscale = 2)
g_new <- list(D5,D6)

# predict
predict_diagram_ksvm(new_diagrams = g_new,model = model_svm,num_workers = 2)
#> [1] circle sphere
#> Levels: circle sphere
```

Check if two groups of persistence diagrams are independent or not:

``` r
# create two independent groups of diagrams of length 6, which
# is the minimum length
D1 <- TDA::ripsDiag(TDA::circleUnif(n = 20,r = 1),
                    maxdimension = 1,maxscale = 2)
D2 <- TDA::ripsDiag(TDA::circleUnif(n = 20,r = 1),
                    maxdimension = 1,maxscale = 2)
g1 <- list(D1,D2,D2,D2,D2,D2)
g2 <- list(D2,D1,D1,D1,D1,D1)

# do independence test with sigma = t = 1 in dimension 1
indep_test <- independence_test(g1,g2,dims = c(1),num_workers = 2)
```

Performing fast persistent homology with python:

``` r
# uniformly sample from a unit circle
circ <- TDA::circleUnif(n = 50,r = 1)

# import the ripser python module
ripser <- import_ripser()

# run persistent homology
diagram <- PyH(circ,maxdim = 1,thresh = 1,ripser = ripser)
```

Finding real topological features in data:

```r
# uniformly sample from a unit circle
circ <- TDA::circleUnif(n = 50,r = 1)

# find real topological features
boot <- bootstrap_persistence_thresholds(X = circ,FUN = "calculate_homology",maxdim = 1,thresh = 2)
```

Plot a diagram with persistence thresholds:

```r
diag <- TDAstats::calculate_homology(circ,dim = 1)
plot_diagram(diag,thresholds = boot)
```