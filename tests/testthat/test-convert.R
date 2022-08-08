
test_that("diagram_to_df can accept the right kinds of input",{
  
  D <- TDA::circleUnif(n = 20,r = 1)
  phom_TDA <- TDA::ripsDiag(X = D,maxdimension = 1,maxscale = 2)
  phom_TDAstats <- TDAstats::calculate_homology(mat = D,threshold = 2)
  expect_s3_class(diagram_to_df(phom_TDA),"data.frame")
  expect_s3_class(diagram_to_df(phom_TDAstats),"data.frame")
  expect_s3_class(diagram_to_df(diagram_to_df(phom_TDA)),"data.frame")

})