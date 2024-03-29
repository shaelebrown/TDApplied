
test_that("plot_diagram can detect incorrect parameters",{
  
  expect_error(plot_diagram(D = data.frame(dimension = c(0:13),birth = rep(0,14),death = rep(1,14))),"12")
  expect_error(plot_diagram(D = data.frame(dimension = c(0,1,2,3),birth = c(0,0,0,0),death = c(1,1,1,Inf))),"finite")
  expect_error(plot_diagram(D = data.frame(dimension = c(0,1,2,3),birth = c(0,0,0,0),death = c(1,1,1,2)),title = NA),"NA")
  expect_error(plot_diagram(D = data.frame(dimension = c(0,1,2,3),birth = c(0,0,0,0),death = c(1,1,1,2)),title = 2),"character")
  expect_error(plot_diagram(D = data.frame(dimension = c(0,1,2,3),birth = c(0,0,0,0),death = c(1,1,1,2)),max_radius = NA),"numeric")
  expect_error(plot_diagram(D = data.frame(dimension = c(0,1,2,3),birth = c(0,0,0,0),death = c(1,1,1,2)),max_radius = c(1,2)),"single")
  expect_error(plot_diagram(D = data.frame(dimension = c(0,1,2,3),birth = c(0,0,0,0),death = c(1,1,1,2)),max_radius = Inf),"finite")
  expect_error(plot_diagram(D = data.frame(dimension = c(0,1,2,3),birth = c(0,0,0,0),death = c(1,1,1,2)),max_radius = -1),"positive")
  expect_error(plot_diagram(D = data.frame(dimension = c(0,1,2,3),birth = c(0,0,0,0),death = c(1,1,1,2)),legend = NULL),"NULL")
  expect_error(plot_diagram(D = data.frame(dimension = c(0,1,2,3),birth = c(0,0,0,0),death = c(1,1,1,2)),legend = c(T,F)),"single")
  expect_error(plot_diagram(D = data.frame(dimension = c(0,1,2,3),birth = c(0,0,0,0),death = c(1,1,1,2)),thresholds = c(0,1,2,NA)),"NA")
  expect_error(plot_diagram(D = data.frame(dimension = c(0,1,2,3),birth = c(0,0,0,0),death = c(1,1,1,2)),thresholds = list(thresholds = c(0,1,2,NA))),"NA")
  expect_error(plot_diagram(D = data.frame(dimension = c(0,1,2,3),birth = c(0,0,0,0),death = c(1,1,1,2)),thresholds = list(foo = c(1,2,3))),"list element")
  expect_error(plot_diagram(D = data.frame(dimension = c(0,1,2,3),birth = c(0,0,0,0),death = c(1,1,1,2)),thresholds = list(thresholds = c(1,2,3))),"element")
  expect_error(plot_diagram(D = data.frame(dimension = c(0,1,2,3),birth = c(0,0,0,0),death = c(1,1,1,2)),thresholds = c(1,2,3)),"element")
  expect_error(plot_diagram(D = data.frame(dimension = c(0,1,2,3),birth = c(0,0,0,0),death = c(1,1,1,2)),thresholds = c(1,2,3,"5")),"numeric")
  expect_error(plot_diagram(D = data.frame(dimension = c(0,1,2,3),birth = c(0,0,0,0),death = c(1,1,1,2)),thresholds = list(thresholds = c(1,2,3,"5"))),"numeric")
  expect_error(plot_diagram(D = data.frame(dimension = c(0,1,2,3),birth = c(0,0,0,0),death = c(1,1,1,2)),thresholds = c(1,2,3,NA)),"NA")
  
})

test_that("plot_diagram is working correctly",{
  
  expect_identical(plot_diagram(D = data.frame(dimension = numeric(),birth = numeric(),death = numeric())),NULL)
  expect_identical(plot_diagram(D = data.frame(dimension = c(0),birth = c(0),death = c(1))),NULL)
  
})

