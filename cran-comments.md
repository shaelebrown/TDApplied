
## Test environments
* local Mac OS X install, R 4.1.2
* win-builder (devel and release)
* rhub windows server 2022, R-devel, 64 bit
* rhub Ubuntu Linux 20.04.1 LTS, R-release GCC
* rhub Fedora Linux, R-devel, clang, gfortran

## R CMD check results

0 errors | 0 warnings | 0 notes

## NOTES

* a note about how this package used to be on CRAN but was removed due to not fixing problems in time is always there
* when using rhub to check for cran, a note is given that says a detritus was found in the temp directory, but then nothing is listed 
* a 'lastMiKTeXException' was found when looking for a detritus on windows when checking on rhub, but I've read that this shouldn't be an issue (and no other detritus was found in other checks)
* some of the examples run for over 5s, however these examples have been made as small and fast as possible without throwing errors
* on Ubuntu there was a PREP error, however in the log there were only notes (like the ones above) and the build was successful
* there are domain-specific words and author names in ML_and_Inference.Rmd which were flagged by devtools::check_spelling() but to the author's knowledge they are all spelled correctly


