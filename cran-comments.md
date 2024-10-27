
## Test environments
* local Mac OS X install, R 4.1.2
* win-builder (devel and release)
* rhub windows virtual machine
* rhub macos virtual machine
* rhub linux virtual machine
* rhub ubuntu-release, valgrind, ubuntu-clang, clang19 and atlas containers

## R CMD check results

0 errors | 0 warnings | 1 note

## NOTES

* the note on R CMD check is for large sub directory size (necessary for the extensive documentation needed for journal publication).
* on rhub there are build errors for gcc14 (Fedora Linux R devel) and macos-arm64, seemingly because some of the package dependencies are not available on those platforms.
* some of the examples run for over 5s, however these examples have been made as small and fast as possible without throwing errors.
* there are domain-specific words and author names in ML_and_Inference.Rmd which were flagged by devtools::check_spelling() but to the author's knowledge they are all spelled correctly.


