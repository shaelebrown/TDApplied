> All changes to TDApplied are documented here.

> Additions referenced with relevant [GitHub Issue](https://github.com/shaelebrown/TDApplied/issues) or
[Pull Request](https://github.com/shaelebrown/TDApplied/pulls) number.
Please see those for more details.

# 3.0.2
- all CRAN issues for this update were caused by the rho parameter, which invokes external C++ code. We therefore fixed these issues by removing the rho parameter from the predict_diagram_kpca and diagram_distance examples, from all tests and in the ML_and_inference.Rmd file. This parameter has been kept and is still tested, just not tested on CRAN
- sped up the independence_test example by only showing the Gram-matrix approach
- removed warnings from benchmarking plots in Speed.Rmd
- removed dependency on package TDA which is currently unavailable on CRAN

# 3.0.1
- same updates as 3.0.0 but with more efficient vignette building

# 3.0.0
- added ability to precompute distance/Gram matrices for ML and inference functions
- added fast approximation to Fisher information metric
- added vignettes for speedups, HCP analysis, personalized analyses and distance calculation comparisons (and removed those parts from the main vignette)
- fixed issues with cv model fitting in diagram_ksvm
- added automatic calculation of t parameters in diagram_ksvm
- decreased memory load on parallel functions (except for permutation test loss function)
- added checks for 0 variance distance matrices in diagram_ksvm
- added comparisons against package rgudhi
- updated DESCRIPTION
- added interpretations tools for vr graphs and multiple representative (co) cycles
- improved HCP analysis
- resolved some distance 0 cases in diagram_distance

# 2.0.4
- fixed build issues related to use of suggested packages in tests, examples and vignettes

# 2.0.3
- fixed bootstrap reference in vignette

# 2.0.2
- set seed in vignette for reproducibility (which is reset at the end)
- added more examples of TDA applications in publications

# 2.0.1
- increased testing coverage
- fixed issue with th parameter in diagram_kpca
- fixed issue with gamma distribution in independence_test
- added applied analysis of TDApplied on HCP data to package vignette

# 2.0.0

- added PyH function for fast persistence diagram calculations with python
- added bootstrap_persistence_thresholds for finding "real" topological features in a data set
- added plot_diagram function for plotting persistence diagrams, with or without persistence thresholds
- fixed problem with diagram_distance in which one of the two diagrams was empty in the
desired dimension

# 0.1.3

- fixed small bug with computing mean cv model error for svm
- added tryCatch's around parallelized code to ensure that clusters are closed even when errors occur

# 0.1.2

- fixed bug with mds test and properly cleaned up parallelization clusters

# 0.1.1

- Fixed bug with one diagram_mds test, although code was working properly

# 0.1.0

- Initial version