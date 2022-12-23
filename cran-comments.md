
## Test environments
* local Mac OS X install, R 4.1.2
* win-builder (devel and release)
* rhub windows server 2022, R-devel, 64 bit
* rhub Ubuntu Linux 20.04.1 LTS, R-release GCC
* rhub Fedora Linux, R-devel, clang, gfortran

## R CMD check results

0 errors | 0 warnings | 0 notes

## NOTES

* on Windows I get an error "No suitable spell-checker program found" on rhub, although the build gets marked as a success
* some of the examples run for over 5s, however these examples have been made as small and fast as possible without throwing errors
* on Fedora Linux when checking the rebuilding of package vignettes I got the following error: pandoc-citeproc: Error in $: Incompatible API versions: encoded with [1,20] but attempted to decode with [1,17,0,4]. But I'm not sure how/if I need to fix this
* there are domain-specific words and author names in ML_and_Inference.Rmd which were flagged by devtools::check_spelling() but to the author's knowledge they are all spelled correctly


