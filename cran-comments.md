
## Test environments
* local Mac OS X install, R 4.1.2
* win-builder (devel and release)
* rhub windows server 2022, R-devel, 64 bit
* rhub Ubuntu Linux 20.04.1 LTS, R-release GCC
* rhub Fedora Linux, R-devel, clang, gfortran

## R CMD check results

0 errors | 0 warnings | 0 notes

## NOTES

* only on the Windows check for rhub I got strange notes that a 'NULL' directory was found, and a lastMiKTeXException exception, but after checking online these may be due to issues with rhub.
* some of the examples run for over 5s, however these examples have been made as small and fast as possible without throwing errors
* there are domain-specific words and author names in ML_and_Inference.Rmd which were flagged by devtools::check_spelling() but to the author's knowledge they are all spelled correctly


