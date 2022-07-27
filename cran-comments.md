
## Test environments
* local Mac OS X install, R 4.1.2
* win-builder (devel and release)
* rhub windows server 2022, R-devel, 64 bit
* rhub Ubuntu Linux 20.04.1 LTS, R-release GCC
* rhub Fedora Linux, R-devel, clang, gfortran

## R CMD check results

0 errors | 0 warnings | 0 notes

## Notes

* This is a first submission.
* The words "TDA", "TDAML", and "scalable", in the DESCRIPTION file are all spelled correctly. 
* There are other domain-specific words and author names in ML_and_Inference.Rmd which were flagged by devtools::check_spelling() but to the author's knowledge they are all spelled correctly (this did not result in a Note, however).
* On the rhub windows server there are five functions with examples taking longer than 5s in elapsed time, but it is not possible for the examples for these functions to be any smaller as it will throw an error in each case.
* Occasionally a Note is found like "checking for detritus in the temp directory ... NOTE
  Found the following files/directories:..." however there doesn't seem to be a clear fix for this.
