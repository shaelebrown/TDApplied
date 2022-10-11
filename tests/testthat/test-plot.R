
test_that("plot_diagram can detect incorrect parameters",{
  
  expect_error(plot_diagram(D = data.frame(dimension = c(0,1,2,3,4),birth = c(0,0,0,0,0),death = c(0,0,0,0,0))),"3")
  expect_error(plot_diagram(D = data.frame(dimension = c(0,1,2,3),birth = c(0,0,0,0),death = c(0,0,0,Inf))),"finite")
  expect_error(plot_diagram(D = data.frame(dimension = c(0,1,2,3),birth = c(0,0,0,0),death = c(0,0,0,2)),title = NA),"NA")
  
})