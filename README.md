
<!-- README.md is generated from README.Rmd. Please edit that file -->

# **TDApplied**

<!-- badges: start -->
<!-- badges: end -->

**TDApplied** is an R package for analyzing persistence diagrams using
machine learning and statistical inference, and is designed to interface
with persistent (co)homology calculations from the R packages **TDA**
and **TDAstats**.

R package **TDA**:

> Fasy, Brittany T., Jisu Kim, Fabrizio Lecci, Clement Maria, David L.
> Millman, and Vincent Rouvreau. 2021. TDA: Statistical Tools for
> Topological Data Analysis. <https://CRAN.R-project.org/package=TDA>.

R package **TDAstats**:

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

## Citation

To cite this package in publication please use the BibTex entry:

@Manual{TDApplied, title = {TDApplied: Machine Learning and Inference
for Topological Data Analysis}, author = {Shael Brown and Dr. Reza
Farivar}, note = {R package version 3.0.0}, url =
{<https://github.com/shaelebrown/TDApplied>}, }

If you wish to cite a particular method used in **TDApplied** see the
REFERENCES.bib file in the vignette directory.

## Functionality

**TDApplied** has three major modules:

1.  Computing and interpreting persistence diagrams. The `PyH` function
    connects with python creating a fast persistent (co)homology engine
    compared to alternatives. The `plot_diagram` function can be used to
    plot diagrams computed from `PyH` or the **TDA** and **TDAstats**
    packages. The `rips_graphs` and `plot_rips_graphs` functions can be
    used to visualize dataset structure at the scale of particular
    topological features. The `bootstrap_persistence_thresholds`
    function can be used to identify statistically significant
    topological features in a dataset.
2.  Machine learning. The functions `diagram_mds`, `diagram_kpca` and
    `predict_diagram_kpca` can be used to project a group of diagrams
    into a low dimensional space (i.e. dimension reduction). The
    functions `diagram_kkmeans` and `predict_diagram_kkmeans` can be
    used to cluster a group of diagrams. The functions `diagram_ksvm`
    and `predict_diagram_ksvm` can be used to link, through a prediction
    function, persistence diagrams and an outcome (i.e. dependent)
    variable.
3.  Statistics. The `permutation_test` function acts like an ANOVA test
    for identifying group differences of persistence diagrams. The
    `independence_test` function can determine if two groups of paired
    persistence diagrams are likely independent or not.

Not only does **TDApplied** provide methods for the applied analysis of
persistence diagrams which were previously unavailable, but an emphasis
on speed and scalability through parallelization, C code, avoiding
redundant slow computations, etc., makes **TDApplied** a powerful tool
for carrying out applied analyses of persistence diagrams.

## Example Code

This example creates nine persistence diagrams, plots one and projects
all nine into 2D space using multidimensional scaling (MDS) to
demonstrate **TDApplied** functionalities.

``` r
library(TDApplied)

# create 9 persistence diagrams
# 3 from circles, 3 from tori and 3 from spheres
circ1 <- TDAstats::calculate_homology(TDA::circleUnif(n = 50),dim = 1,threshold = 1)
circ2 <- TDAstats::calculate_homology(TDA::circleUnif(n = 50),dim = 1,threshold = 1)
circ3 <- TDAstats::calculate_homology(TDA::circleUnif(n = 50),dim = 1,threshold = 1)
torus1 <- TDAstats::calculate_homology(TDA::torusUnif(n = 50,a = 0.25,c = 1),dim = 1,threshold = 1)
torus2 <- TDAstats::calculate_homology(TDA::torusUnif(n = 50,a = 0.25,c = 1),dim = 1,threshold = 1)
torus3 <- TDAstats::calculate_homology(TDA::torusUnif(n = 50,a = 0.25,c = 1),dim = 1,threshold = 1)
sphere1 <- TDAstats::calculate_homology(TDA::sphereUnif(n = 50,d = 2),dim = 1,threshold = 1)
sphere2 <- TDAstats::calculate_homology(TDA::sphereUnif(n = 50,d = 2),dim = 1,threshold = 1)
sphere3 <- TDAstats::calculate_homology(TDA::sphereUnif(n = 50,d = 2),dim = 1,threshold = 1)

# plot a diagram
plot_diagram(circ1,title = "Circle 1")

# project into 2D and plot
proj_2D <- diagram_mds(list(circ1,circ2,circ3,torus1,torus2,torus3,sphere1,sphere2,sphere3),dim = 1,k = 2)
plot(x = proj_2D[,1],y = proj_2D[,2])
```

## Documentation

**TDApplied** has five major vignettes:

1.  “TDApplied Theory and Practice”, which documents the background
    theory and practical usage of all functions (on simple simulated
    data).
2.  “Human Connectome Project Analysis”, which provides a sample
    analysis of real neurological data using **TDApplied**.
3.  “Benchmarking and Speedups”, which describes all implemented
    optimizations of **TDApplied** functions and compares the runtime of
    **TDApplied** functions with functions from other packages.
4.  “Personalized Analyses with TDApplied”, which demonstrates how
    machine learning (or statistical) models and pipelines, other than
    those implemented in **TDApplied**, can be fit to persistence
    diagrams.
5.  “Comparing Distance Calculations”, which accounts for differences in
    distance functions of persistence diagrams across R packages.
