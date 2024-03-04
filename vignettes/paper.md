
---
title: 'TDApplied: An R package for machine learning and inference with persistence diagrams'
tags:
  - R
  - topological data analysis
  - persistent homology
authors:
  - name: Shael Brown
    orcid: 0000-0001-8868-2867
    affiliation: 1
  - name: Reza Farivar-Mohseni
    orcid: 0000-0002-3123-2627
    affiliation: 2
affiliations:
 - name: Department of Quantitative Life Sciences, McGill University, Montreal Canada.
   index: 1
 - name: McGill Vision Research, Department of Opthamology, McGill University, Montreal Canada.
   index: 2
date: 24 January 2024
bibliography: paper.bib

---

# Summary

Topological data analysis is a collection of tools, based on the mathematical fields of topology and geometry, for finding structure in whole datasets. Its main tool, persistent homology [@PHoriginal;@ComputingPH], computes a shape descriptor of a dataset called a persistence diagram which encodes information about holes that exist in the dataset (example applications span a variety of areas, see for example [@TDA_ADHD;@word_embeddings;@TDA_chemistry]). These types of features cannot be identified by other methods, making persistence diagrams a unique and valuable data science object for studying and comparing datasets. The two most popular data science tools for analyzing multiple objects are machine learning and inference, but to date there has been no open source implementation of published methods for machine learning and inference of persistence diagrams.

# Statement of need

`TDApplied` is the first R package for machine learning and inference of persistence diagrams, building on the main R packages for the calculation of persistence diagrams `TDA` [@R-TDA] and `TDAstats` [@R-TDAstats;@TDAstats2018] and publications of applied analysis methods for persistence diagrams [@Robinson_Turner;@persistence_fisher]. `TDApplied` is intended to be used by academic researchers and industry professionals wanting to integrate persistence diagrams into their analysis workflows. An example `TDApplied` workflow, in which the topological differences between three datasets are visualized in 2D using multidimensional scaling (MDS) [@Cox2008], is visualized in figure \autoref{fig:software}: 

![An example `TDApplied` workflow. A dataset (D1, left) contains one loop (yellow) and two clusters (the loop forms one cluster and the three points on the bottom are another cluster, and clusters are denoted by the color red). These topological features are captured with persistent homology in a persistence diagram PD1 (middle top), and two other data sets, D2 and D3 (not shown), have their persistence diagrams, PD2 and PD3, computed (middle center and middle bottom). PD1 and PD2 are not very topologically different in terms of their loops, with both containing a loop with similar birth and death values, and this is represented by a dashed-line relationship. On the other hand, PD2 and PD3 are topologically different in terms of their loops because PD3 does not contain a loop, and this is represented by a dotted-line relationship. `TDApplied` can quantify these topological differences and use MDS to project the persistence diagrams into three points in a 2D embedding space (right) where interpoint distances reflect the topological differences between the persistence diagrams. \label{fig:software}](software.pdf){width=100%}

The `TDApplied` package is built on three main pillars:

1. User-friendly -- internal preprocessing of persistence diagrams that would normally be left to R users to figure out ad hoc, and functions designed to easily flow from input diagrams to output metrics.
2. Efficient -- parallelization, C code, computational tricks and storage of reusable and cumbersome calculations significantly increases the feasibility of topological analyses (compared to existing R packages).
3. Flexible -- ability to interface with other data science packages to create  personalized analyses.

`TDApplied` has already been featured in a [conference workshop](https://github.com/WoComtoQC/wocomtoqc.github.io/blob/main/abstract.md) and a [conference tutorial](https://www.ihcisociety.org/program/tutorial-lecture), utilized in a journal publication [@Yash] and downloaded over 4400 times. Therefore, we propose `TDApplied` as a user-friendly, efficient and flexible R package for the analysis of multiple datasets using machine learning and inference via topological data analysis.

# Project Management

*Installation and availability:* `TDApplied` can be installed directly from CRAN using the command `install.packages("TDApplied")`, or from GitHub using the `devtools` package [@R-devtools]. `TDApplied` is distributed under the GPL-3 license.

*Code quality:* Code has been tested using the `testthat` package [@R-testthat], with 91.45\% coverage of R code when not skipping tests involving Python code (or 88.44\% coverage when skipping the Python tests).

*Documentation:* `TDApplied` contains five main vignettes: 

1. "TDApplied Theory and Practice" provides example function usage on simulated data as well as mathematical background and intuition, 
2. "Human Connectome Project Analysis" demonstrates an applied example analysis of neurological data, 
3. "Benchmarking and Speedups" outlines the package's optimization strategies and highlights performance gains compared to other packages, 
4. "Personalized Analyses with TDApplied" demonstrates how to interface `TDApplied` with other data science packages, and 
5. "Comparing Distance Calculations" accounts for differences in computed distance values between persistence diagrams across comparable packages.

# Acknowledgements

We acknowledge funding from the CIHR 2016 grant for cortical mechanisms of 3-D scene and object recognition in the primate brain.

# References