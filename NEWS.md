> All changes to TDApplied are documented here.

> Additions referenced with relevant [GitHub Issue](https://github.com/shaelebrown/TDApplied/issues) or
[Pull Request](https://github.com/shaelebrown/TDApplied/pulls) number.
Please see those for more details.

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