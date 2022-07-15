
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TDAML

<!-- badges: start -->
<!-- badges: end -->

TDAML is an R package for applied topological data analysis using
machine learning and statistical inference, and uses the output of
persistent homology calculations from the R package TDA as input to its
methods.

R package TDA:

> Fasy, Brittany T., Jisu Kim, Fabrizio Lecci, Clement Maria, David L.
> Millman, and Vincent Rouvreau. 2021.TDA: Statistical Tools for
> Topological Data Analysis. <https://CRAN.R-project.org/package=TDA>.

## Installation

To install the latest version of this R package directly from github:

    install.packages("devtools")
    library(devtools)
    devtools::install_github("shaelebrown/TDAML")
    library(TDAML)

To install from Github you might need:

-   **Windows:** Rtools
    (<https://cran.r-project.org/bin/windows/Rtools/>)
-   **OS X:** xcode (from the app store)
-   **Linux:** apt-get install r-base-dev (or similar).

To install the stable version of this R package from CRAN:

    install.packages("TDAML")

## Examples

These are basic examples which show you how to use the package:

``` r
library(TDAML)
```

For these examples we will use three base persistence diagrams:

``` r
D1 = generate_TDAML_test_data(1,0,0)
D2 = generate_TDAML_test_data(0,1,0)
D3 = generate_TDAML_test_data(0,0,1)
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
g <- generate_TDAML_test_data(3,3,3)

# calculate their 2D MDS embedding in dimension 1 with the bottleneck distance
mds <- diagram_mds(diagrams = g,dim = 0,p = Inf,k = 2,num_workers = 2)
```

Looking for group differences in groups of persistence diagrams:

``` r
g1 <- generate_TDAML_test_data(3,0,0)
g2 <- generate_TDAML_test_data(0,3,0)
g3 <- generate_TDAML_test_data(0,0,3)
perm_test <- permutation_test(g1,g2,g3,
                              num_workers = 2,
                              dims = c(0))
```

Clustering persistence diagrams with kernel k-means:

``` r
g <- generate_TDAML_test_data(3,3,3)

# calculate kmeans clusters with centers = 3, and sigma = t = 2
clust <- diagram_kkmeans(diagrams = g,centers = 3,dim = 0,t = 2,sigma = 2,num_workers = 2)
```

Predicting new cluster labels:

``` r
# create nine new diagrams
g_new <- generate_TDAML_test_data(3,3,3)

# predict cluster labels
predict_diagram_kkmeans(new_diagrams = g_new,clustering = clust,num_workers = 2)
#> [1] 3 3 3 1 1 1 2 2 2
```

Computing a kernel PCA embedding of persistence diagrams:

``` r
g <- generate_TDAML_test_data(3,3,3)

# calculate their 2D PCA embedding with sigma = t = 2
pca <- diagram_kpca(diagrams = g,dim = 0,t = 2,sigma = 2,features = 2,num_workers = 2)
```

Project new persistence diagrams into a kernel PCA embedding:

``` r
# project new diagrams onto old model
g_new <- generate_TDAML_test_data(3,3,3)

# predict cluster labels
new_pca <- predict_diagram_kpca(new_diagrams = g_new,embedding = pca,num_workers = 2)
```

Fit a kernel SVM model on persistence diagrams:

``` r
# create thirty diagrams
g <- generate_TDAML_test_data(10,10,10)

# create response vector
y <- as.factor(rep(c("D1","D2","D3"),each = 10))

# fit model with cross validation
model_svm <- diagram_ksvm(diagrams = g,cv = 2,dim = c(0),
                          y = y,sigma = c(1,0.1),t = c(1,2),
                          num_workers = 2)
```

Predict labels for new persistence diagrams:

``` r
# create nine new diagrams
g_new <- generate_TDAML_test_data(3,3,3)

# predict
predict_diagram_ksvm(new_diagrams = g_new,model = model_svm,num_workers = 2)
#> [1] D1 D1 D1 D2 D2 D2 D3 D3 D3
#> Levels: D1 D2 D3
```

Check if two groups of persistence diagrams are independent or not:

``` r
# create copies of D1 and D2
g1 <- generate_TDAML_test_data(10,0,0)
g2 <- generate_TDAML_test_data(0,10,0)

# do independence test with sigma = t = 1
indep_test <- independence_test(g1,g2,dims = c(0),num_workers = 2)
```
