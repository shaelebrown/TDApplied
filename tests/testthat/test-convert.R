
# test_that("diagram_to_df can accept the right kinds of input",{
# 
#   skip_if_not_installed("TDA")
#   skip_if_not_installed("TDAstats")
#   D <- TDA::circleUnif(n = 20,r = 1)
#   phom_TDA <- TDA::ripsDiag(X = D,maxdimension = 1,maxscale = 2)
#   phom_TDAstats <- TDAstats::calculate_homology(mat = D,threshold = 2)
#   simulated_PyH_phom <- list(diagram = diagram_to_df(phom_TDA),representatives = list())
#   expect_s3_class(diagram_to_df(phom_TDA),"data.frame")
#   expect_s3_class(diagram_to_df(phom_TDAstats),"data.frame")
#   expect_s3_class(diagram_to_df(diagram_to_df(phom_TDA)),"data.frame")
#   expect_s3_class(diagram_to_df(simulated_PyH_phom),"data.frame")
# 
# })

test_that("diagram_to_df can detect incorrect parameters properly",{
  
  expect_error(diagram_to_df(2),"computation")
  
})