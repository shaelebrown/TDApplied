
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TDApplied

<!-- badges: start -->
<!-- badges: end -->

TDApplied is an R package for analyzing persistence diagrams using
machine learning and statistical inference, and is designed to interface
with persistent (co)homology calculations from the R packages TDA and
TDAstats.

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

## Citation

To cite this package in publication please use the BibTex entry:

@Manual{TDApplied, title = {TDApplied: Machine Learning and Inference
for Topological Data Analysis}, author = {Shael Brown and Dr. Reza
Farivar}, note = {R package version 3.0.0}, url =
{<https://github.com/shaelebrown/TDApplied>}, }

If you wish to cite a particular method used in TDApplied see the
REFERENCES.bib file in the vignette directory.

## Functionality

TDApplied has five major goals:

1.  Deliver a fast engine for calculating persistence diagrams: the
    `PyH` function connects with python creating a fast persistent
    (co)homology engine compared to alternatives.
2.  Convert persistence diagrams computed using the R packages TDA and
    TDAstats into a commonly-used data type for data analyses (a data
    frame): `diagram_to_df` performs this conversion.
3.  Quantify differences and similarities between pairs of persistence
    diagrams with high speed: the `diagram_distance` and
    `diagram_kernel` functions allow for fast distance and kernel
    calculations respectively.
4.  Contribute tools for interpreting persistence diagrams: the
    `plot_diagram` and `bootstrap_persistence_thresholds` functions can
    be used to plot and threshold persistence diagrams (for
    “significant” features).
5.  Provide parallelized methods for machine learning and inference for
    persistence diagrams: `diagram_mds`, `diagram_kpca`,
    `diagram_kkmeans` and `diagram_ksvm` perform machine learning with
    persistence diagrams (multidimensional scaling, kernel principal
    components analysis, kernel k-means clustering and support vector
    machines respectively) and `permutation_test` and
    `independence_test` perform inference on groups of diagrams (a
    permutation test for group differences and a two-sample independence
    test respectively).

Through accomplishing these goals TDApplied is a fast and scalable
one-stop-shop for all things persistent homology.

## Example Code

This example creates nine persistence diagrams, plots one and projects
all nine into 2D space using MDS to demonstrate TDApplied
functionalities.

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

TDApplied has three major vignettes:

1.  “Introduction to TDApplied”, which documents the background theory
    and practical usage of all functions (on simple simulated data),
2.  “Benchmarking and Speedups”, which describes all implemented
    optimizations of TDApplied functions and compares the runtime of
    TDApplied functions with functions from other packages, and
3.  “Practical Example - Human Connectome Project Analysis”, which
    provides a sample analysis of real neurological data using
    TDApplied.
