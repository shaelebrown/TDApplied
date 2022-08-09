
## Test environments
* local Mac OS X install, R 4.1.2
* win-builder (devel and release)
* rhub windows server 2022, R-devel, 64 bit
* rhub Ubuntu Linux 20.04.1 LTS, R-release GCC
* rhub Fedora Linux, R-devel, clang, gfortran

## R CMD check results

0 errors | 0 warnings | 0 notes

## Notes

* This is a third submission for this package (the first two were rejected). The first submission contained two failed tests, and this has been fixed. The second submission had one bad url, which has been fixed.
* A previous version of this package was submitted under the name "TDAML", however that submission should be ignored.
* The words "TDApplied" and "scalable", in the DESCRIPTION file are all spelled correctly. 
* There are other domain-specific words and author names in ML_and_Inference.Rmd which were flagged by devtools::check_spelling() but to the author's knowledge they are all spelled correctly (this did not result in a Note, however).
* On the various platforms there are up to eleven functions with examples taking longer than 5s in elapsed time. These examples were made to be as small as possible while not throwing an error and not being unrealistic. The functions in this long examples are often intended to be run with more cores than two, however two cores were used to avoid causing issues on the server.
* Occasionally a Note is found like "checking for detritus in the temp directory ... NOTE
  Found the following files/directories:..." however there doesn't seem to be a clear fix for this.
* There is a PREPERROR when checking on rhub Ubuntu Linux, however the log files say the build was successful.
