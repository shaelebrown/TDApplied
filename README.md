# Machine learning and statistical inference of persistence diagrams with TDAML

## Description

TDAML is an R package for applied topological data analysis using machine learning and statistical inference, and uses the output of persistent homology calculations from the R package TDA as input to its methods.

R package TDA:

> Fasy, Brittany T., Jisu Kim, Fabrizio Lecci, Clement Maria, David L. Millman, and Vincent Rouvreau. 2021.TDA: Statistical Tools for Topological Data Analysis. https://CRAN.R-project.org/package=TDA.

## Installation

To install the latest version of this R package directly from github:

    install.packages("devtools")
    library(devtools)
    devtools::install_github("shaelebrown/TDAML")
    library(TDAML)

To install from Github you might need: 

- **Windows:** Rtools (http://cran.r-project.org/bin/windows/Rtools/)
- **OS X:** xcode (from the app store)
- **Linux:** apt-get install r-base-dev (or similar).

To install the stable version of this R package from CRAN:

    install.packages("TDAML")

## Examples

Computing distances between persistence diagrams:

```R
# create two diagrams with package TDA based on torus and sphere
torus <- TDA::ripsDiag(TDA::torusUnif(n = 100,a = 1,c = 2),
maxscale = 2,
maxdimension = 2)
sphere <- TDA::ripsDiag(TDA::sphereUnif(n = 100,d = 2,r = 1),
maxscale = 2,
maxdimension = 2)

# calculate their wasserstein distance in dimension 1
diagram_distance(D1 = torus,D2 = sphere,dim = 1,p = 2,distance = "wasserstein")

# calculate their bottleneck distance in dimension 2
diagram_distance(D1 = torus,D2 = sphere,dim = 2,p = Inf,distance = "wasserstein")

# Fisher information metric calculation in dimension 1
diagram_distance(D1 = torus,D2 = sphere,dim = 1,distance = "fisher",sigma = 1)

# Fisher information metric calculation in dimension 2
diagram_distance(D1 = torus,D2 = sphere,dim = 2,distance = "fisher",sigma = 1)
```

Computing kernel values between persistence diagrams:
```R
# create two diagrams with package TDA based on torus and sphere
torus <- TDA::ripsDiag(TDA::torusUnif(n = 100,a = 1,c = 2),
maxscale = 2,
maxdimension = 2)
sphere <- TDA::ripsDiag(TDA::sphereUnif(n = 100,d = 2,r = 1),
maxscale = 2,
maxdimension = 2)

# calculate their kernel value in dimension 2 with sigma = 2, t = 2
diagram_kernel(D1 = torus,D2 = sphere,dim = 2,sigma = 2,t = 2)
```

Computing a MDS projection of persistence diagrams:

```R
# create 9 diagrams with package TDA based on spheres, circles and tori
g <- lapply(X = 1:9,FUN = function(X){

  if(X %% 3 == 0)
  {
    df <- TDA::ripsDiag(TDA::circleUnif(n = 100,r = 1),
                        maxdimension = 1,
                        maxscale = 2)
  }
  if(X %% 3 == 1)
  {
    df <- TDA::ripsDiag(TDA::sphereUnif(n = 100,d = 2,r = 1),
                        maxdimension = 1,
                        maxscale = 2)
  }
  if(X %% 3 == 2)
  {
    df <- TDA::ripsDiag(TDA::torusUnif(n = 100,a = 1,c = 2),
                        maxdimension = 1,
                        maxscale = 2)
  }
  df <- diagram_to_df(d = df)
  return(df)

})

# calculate their 2D MDS embedding in dimension 1 with the bottleneck distance
mds <- diagram_MDS(diagrams = g,dim = 1,p = Inf,k = 2)
```

Looking for group differences in groups of persistence diagrams:

```R
# create two groups of persistence diagrams from Diagram 1 and Diagram 2 respectively
g1 <- lapply(X = 1:10,FUN = function(X){
  
  t <- D1
  t$dimension <- as.numeric(as.character(t$dimension))
  t$birth <- t$birth + stats::rnorm(n = 2,mean = 0,sd = 0.05)
  t[which(t$birth < 0),2] <- 0
  t$death <- t$death + stats::rnorm(n = 2,mean = 0,sd = 0.05)
  return(t)
  
})

g2 <- lapply(X = 1:10,FUN = function(X){
  
  t <- D2
  t$dimension <- as.numeric(as.character(t$dimension))
  t$birth <- t$birth + stats::rnorm(n = 2,mean = 0,sd = 0.05)
  t[which(t$birth < 0),2] <- 0
  t$death <- t$death + stats::rnorm(n = 2,mean = 0,sd = 0.05)
  return(t)
  
})

# do permutation test with 30 iterations, p,q = 2, in dimensions 0 and 1, with
# no pairing using wasserstein distance 
permutation_test(g1,g2,
                 iterations = 30,
                 dims = c(0,1))$p_values
```

Clustering persistence diagrams with kernel k-means:

```R
# concatenate the two lists based on diagram 1 and diagram 2, 
# they should be clearly separable in dimension 1
g <- append(g1,g2)

# calculate kmeans clusters with centers = 2 in dimension 1 with sigma = t = 2
clust <- diagram_kkmeans(diagrams = g,centers = 2,dim = 1,t = 2,sigma = 2)

# display cluster labels
clust$clustering@.Data
```

Predicting new cluster labels:

```R
# create ten new diagrams with package TDA based on spheres and circles
g_new <- lapply(X = 1:10,FUN = function(X){
  
  if(X <= 5)
  {
    t <- D1
    t$dimension <- as.numeric(as.character(t$dimension))
    t$birth <- t$birth + stats::rnorm(n = 2,mean = 0,sd = 0.05)
    t[which(t$birth < 0),2] <- 0
    t$death <- t$death + stats::rnorm(n = 2,mean = 0,sd = 0.05)
  }else
  {
    t <- D2
    t$dimension <- as.numeric(as.character(t$dimension))
    t$birth <- t$birth + stats::rnorm(n = 2,mean = 0,sd = 0.05)
    t[which(t$birth < 0),2] <- 0
    t$death <- t$death + stats::rnorm(n = 2,mean = 0,sd = 0.05)
  }

  return(t)

})

# predict cluster labels
diagram_nearest_clusters(new_diagrams = g_new,clustering = clust)
```

Computing a kernel PCA embedding of persistence diagrams:

```R
# create ten diagrams with package TDA based on spheres and circles
g <- lapply(X = 1:10,FUN = function(X){
  
  if(X <= 5)
  {
    diag <- TDA::ripsDiag(TDA::circleUnif(n = 100,r = 1),
                          maxscale = 2,
                          maxdimension = 1)
  }else
  {
    diag <- TDA::ripsDiag(TDA::sphereUnif(n = 100,d = 2,r = 1),
                          maxscale = 2,
                          maxdimension = 1)
  }
  
  df <- diagram_to_df(d = diag)
  return(df)

})

# calculate their 2D PCA embedding in dimension 1 with sigma = t = 2
pca <- diagram_kpca(diagrams = g,dim = 1,t = 2,sigma = 2,features = 2)
```

Project new persistence diagrams into a kernel PCA embedding:

```R
# create ten new diagrams with package TDA based on spheres and circles
g_new <- lapply(X = 1:10,FUN = function(X){
  
  if(X <= 5)
  {
    diag <- TDA::ripsDiag(TDA::circleUnif(n = 100,r = 1),
                          maxscale = 2,
                          maxdimension = 1)
  }else
  {
    diag <- TDA::ripsDiag(TDA::sphereUnif(n = 100,d = 2,r = 1),
                          maxscale = 2,
                          maxdimension = 1)
  }
  
  df <- diagram_to_df(d = diag)
  return(df)

})

# project new diagrams onto old model
new_pca <- predict_diagram_kpca(new_diagrams = g_new,embedding = pca)
```

Fit a kernel SVM model on persistence diagrams:

```R
# create thirty diagrams with package TDA based on the circle and torus
g <- lapply(X = 1:30,FUN = function(X){
  
  if(X <= 15)
  {
    diag <- TDA::ripsDiag(TDA::torusUnif(n = 100,a = 1,c = 2),
                          maxscale = 2,
                          maxdimension = 2)
  }else
  {
    diag <- TDA::ripsDiag(TDA::circleUnif(n = 100,r = 1),
                          maxscale = 2,
                          maxdimension = 2)
  }
  
  df <- diagram_to_df(d = diag)
  return(df)

})

# create response vector
y <- as.factor(c(rep("torus",15),rep("circle",15)))

# fit model with cross validation
model_svm <- diagram_ksvm(diagrams = g,cv = 2,dim = c(1,2),y = y,sigma = c(1,0.1))
```

Predict labels for new persistence diagrams:

```R
# create ten new diagrams with package TDA based on spheres and circles
g_new <- lapply(X = 1:10,FUN = function(X){
  
  if(X <= 5)
  {
    diag <- TDA::ripsDiag(TDA::circleUnif(n = 100,r = 1),
                          maxscale = 2,
                          maxdimension = 2)
  }else
  {
    diag <- TDA::ripsDiag(TDA::torusUnif(n = 100,a = 1,c = 2),
                          maxscale = 2,
                          maxdimension = 2)
  }
  
  df <- diagram_to_df(d = diag)
  return(df)

})

# answer should be c(rep("circle",5),rep("torus",5))
predict_diagram_ksvm(new_diagrams = g_new,model = model_svm)
```

Check if two groups of persistence diagrams are independent or not:

```R
# create two groups of persistence diagrams from the circle and sphere
circles <- lapply(X = 1:10,FUN = function(X){

  diag <- TDA::ripsDiag(TDA::circleUnif(n = 100,r = 1),
  maxscale = 2,
  maxdimension = 1)
  df <- diagram_to_df(d = diag)
  return(df)

})

spheres <- lapply(X = 1:10,FUN = function(X){

  diag <- TDA::ripsDiag(TDA::sphereUnif(n = 100,d = 2,r = 1),
  maxscale = 2,
  maxdimension = 1)
  df <- diagram_to_df(d = diag)
  return(df)

})

# do independence test with sigma = 1, t = 1, in dimensions 0 and 1
independence_test(circles,spheres,verbose = TRUE)
```



