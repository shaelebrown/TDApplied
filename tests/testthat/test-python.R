# all python tests are skipped to avoid build errors, even though they succeed locally
# to run the tests the reticulate package must be installed, correctly hooked up to
# python, and the ripser module must be downloaded.

test_that("ripser can be imported and verified.",{
  
  skip_if(T)
  ripser <- import_ripser()
  expect_invisible(check_ripser(ripser))
  expect_error(check_ripser(2),"ripser object")
  expect_error(check_ripser(NULL),"ripser object")
  np <- reticulate::import("numpy")
  expect_error(check_ripser(np),"ripser object")
  
})

test_that("PyH can detect bad input parameters.",{
  
  skip_if(T)
  ripser <- import_ripser()
  expect_error(PyH(X = data.frame(),maxdim = 1,thresh = 1,distance_mat = F,ripser = ripser),"two rows")
  expect_error(PyH(X = NULL,maxdim = 1,thresh = 1,distance_mat = F,ripser = ripser),"dataframe")
  expect_error(PyH(X = data.frame(x = 1:2,y = c("1","2")),maxdim = 1,thresh = 1,distance_mat = F,ripser = ripser),"numeric")
  expect_error(PyH(X = data.frame(x = c(1,NA,2)),maxdim = 1,thresh = 1,distance_mat = F,ripser = ripser),"missing")
  expect_error(PyH(X = data.frame(x = 1:3),maxdim = 1,thresh = 1,distance_mat = T,ripser = ripser),"square")
  expect_error(PyH(X = data.frame(x = 1:3,y = 1:3,z = 1:3),maxdim = 1,thresh = 1,distance_mat = T,ripser = ripser),"matrix")
  expect_error(PyH(X = data.frame(x = 1:3,y = 1:3,z = 1:3),maxdim = NA,thresh = 1,distance_mat = F,ripser = ripser),"maxdim")
  expect_error(PyH(X = data.frame(x = 1:3,y = 1:3,z = 1:3),maxdim = -1,thresh = 1,distance_mat = F,ripser = ripser),"maxdim")
  expect_error(PyH(X = data.frame(x = 1:3,y = 1:3,z = 1:3),maxdim = 1,thresh = 1.2,distance_mat = F,ripser = ripser),"whole")
  expect_error(PyH(X = data.frame(x = 1:3,y = 1:3,z = 1:3),maxdim = 1,thresh = 1,distance_mat = NULL,ripser = ripser),"NULL")
  expect_error(PyH(X = data.frame(x = 1:3,y = 1:3,z = 1:3),maxdim = 1,thresh = 1,distance_mat = NA,ripser = ripser),"NA")
  expect_error(PyH(X = data.frame(x = 1:3,y = 1:3,z = 1:3),maxdim = 1,thresh = 1,distance_mat = c(T,F),ripser = ripser),"logical")
  
})

test_that("PyH is computing correctly.",{
  
  skip_if(T)
  D1 <- data.frame(x = stats::rnorm(20),y = stats::rnorm(20))
  D2 <- data.frame(x = stats::rnorm(20),y = stats::rnorm(20))
  D3 <- data.frame(x = stats::rnorm(20),y = stats::rnorm(20))
  
  phom_TDA_1 <- diagram_to_df(TDAstats::calculate_homology(D1,threshold = 5))
  phom_TDA_2 <- diagram_to_df(TDAstats::calculate_homology(D2,threshold = 5))
  phom_TDA_3 <- diagram_to_df(TDAstats::calculate_homology(D3,threshold = 5))
  
  ripser <- import_ripser()
  
  phom_py_1 <- PyH(D1,thresh = 5,ripser = ripser)
  phom_py_2 <- PyH(D2,thresh = 5,ripser = ripser)
  phom_py_3 <- PyH(D3,thresh = 5,ripser = ripser)
  
  expect_equal(phom_TDA_1,phom_py_1,tolerance = 0.00001)
  expect_equal(phom_TDA_2,phom_py_2,tolerance = 0.00001)
  expect_equal(phom_TDA_3,phom_py_3,tolerance = 0.000001)
  
  phom_with_extra_cluster <- PyH(D1,thresh = 5,ripser = ripser,ignore_infinite_clusters = F)
  
  expect_length(which(phom_with_extra_cluster$dimension == 0),20)
  
})