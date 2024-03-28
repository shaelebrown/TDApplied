
<!-- README.md is generated from README.Rmd. Please edit that file -->

# **TDApplied**

<!-- badges: start -->
<!-- badges: end -->

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CRAN version](http://www.r-pkg.org/badges/version/TDApplied)](https://CRAN.R-project.org/package=TDApplied)
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/grand-total/TDApplied)](https://CRAN.R-project.org/package=TDApplied)

[![JOSS DOI](https://joss.theoj.org/papers/10.21105/joss.06321/status.svg)](https://doi.org/10.21105/joss.06321)
[![Zenodo DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10814141.svg)](https://doi.org/10.5281/zenodo.10814141)

## Overview

**TDApplied** is an R package for analyzing persistence diagrams using
machine learning and statistical inference, and is designed to interface
with persistent (co)homology calculations from the R packages **TDA**
and **TDAstats**. Please note that during the development of
**TDApplied**, **TDA** was available on CRAN and therefore included in
package examples and tests, however since that is presently not the case
the dependency on **TDA** has been removed (and therefore some examples
and tests have been modified) but **TDApplied** will still work with
**TDA** computed persistence diagrams and **TDA** functions if a user
already has a working version installed.

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

To install the latest version of this R package directly from GitHub:

    install.packages("devtools")
    library(devtools)
    devtools::install_github("shaelebrown/TDApplied")
    library(TDApplied)

To install from GitHub you might need:

-   **Windows:** Rtools
    (<https://cran.r-project.org/bin/windows/Rtools/>)
-   **OS X:** xcode (from the app store)
-   **Linux:** apt-get install r-base-dev (or similar).

To install the stable version of this R package from CRAN:

    install.packages("TDApplied")

## Citation

If you use TDApplied, please consider citing as:

- Brown et al., (2024). TDApplied: An R package for machine learning and inference with persistence diagrams. Journal of Open Source Software, 9(95), 6321, https://doi.org/10.21105/joss.06321

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

This example creates six persistence diagrams, plots one and projects
all six into 2D space using multidimensional scaling (MDS) to
demonstrate **TDApplied** functionalities.

``` r
library(TDApplied)

# create 6 persistence diagrams
# 3 from circles and 3 from spheres
circ1 <- TDAstats::calculate_homology(TDAstats::circle2d[sample(1:100,size = 50),],dim = 1,threshold = 2)
circ2 <- TDAstats::calculate_homology(TDAstats::circle2d[sample(1:100,size = 50),],dim = 1,threshold = 2)
circ3 <- TDAstats::calculate_homology(TDAstats::circle2d[sample(1:100,size = 50),],dim = 1,threshold = 2)
sphere1 <- TDAstats::calculate_homology(TDAstats::sphere3d[sample(1:100,size = 50),],dim = 1,threshold = 2)
sphere2 <- TDAstats::calculate_homology(TDAstats::sphere3d[sample(1:100,size = 50),],dim = 1,threshold = 2)
sphere3 <- TDAstats::calculate_homology(TDAstats::sphere3d[sample(1:100,size = 50),],dim = 1,threshold = 2)

# plot a diagram
plot_diagram(circ1,title = "Circle 1")

# project into 2D and plot
proj_2D <- diagram_mds(list(circ1,circ2,circ3,sphere1,sphere2,sphere3),dim = 1,k = 2)
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

## Contribute

To contribute to **TDApplied** you can create issues for any
bugs/suggestions on the [issues page](https://github.com/shaelebrown/TDApplied/issues). You can also fork the **TDApplied**
repository and create pull requests to add features you think will be
useful for users.

## Published applications

- Shael Brown and Reza Farivar. The topology of representational geometry. bioRxiv, 2024.
- Yashbir Singh, Colleen M. Farrelly, Quincy A. Hathaway, Tim Leiner, Jaidip Jagtap, Gunnar E. Carlsson, and Bradley J. Erickson. Topological data analysis in medical imaging:
current state of the art. Insights into Imaging, 14(1):58, 2023.
- Rui Dong. Linguistics from a topological viewpoint. arXiv, 2024.
