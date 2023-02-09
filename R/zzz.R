
# unload C++ DLL for proper cleanup
.onUnload <- function (libpath) {
  library.dynam.unload("TDApplied", libpath)
}