
test_that("plot_diagram can detect incorrect parameters",{
  
  expect_error(plot_diagram(D = data.frame(dimension = c(0:13),birth = rep(0,14),death = rep(1,14))),"12")
  expect_error(plot_diagram(D = data.frame(dimension = c(0,1,2,3),birth = c(0,0,0,0),death = c(0,0,0,Inf))),"finite")
  expect_error(plot_diagram(D = data.frame(dimension = c(0,1,2,3),birth = c(0,0,0,0),death = c(0,0,0,2)),title = NA),"NA")
  expect_error(plot_diagram(D = data.frame(dimension = c(0,1,2,3),birth = c(0,0,0,0),death = c(0,0,0,2)),max_radius = NA),"numeric")
  expect_error(plot_diagram(D = data.frame(dimension = c(0,1,2,3),birth = c(0,0,0,0),death = c(0,0,0,2)),max_radius = -1),"positive")
  expect_error(plot_diagram(D = data.frame(dimension = c(0,1,2,3),birth = c(0,0,0,0),death = c(0,0,0,2)),thresholds = c(0,1,2,NA)),"NA")
  
})