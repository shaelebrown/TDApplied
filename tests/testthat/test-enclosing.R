
test_that("enclosing_radius can detect incorrect inputs",{
  
  expect_error(enclosing_radius(NULL, NULL), "distance_mat")
  expect_error(enclosing_radius(NULL, c(T,F)), "single")
  expect_error(enclosing_radius(NULL, NA), "NA")
  expect_error(enclosing_radius(NULL, T), "X")
  expect_error(enclosing_radius(data.frame(),T),"X")
  expect_error(enclosing_radius(data.frame(x = 1),T),"X")
  expect_error(enclosing_radius(data.frame(x = c(1,2)),T),"X")
  expect_error(enclosing_radius(X = NULL,T),"X")
  expect_error(enclosing_radius(X = data.frame(x = c(1,NA)),T),"missing")
  expect_error(enclosing_radius(data.frame(x = c(1),y = c(2)),T),"two")
  expect_error(enclosing_radius(data.frame(x = c(1,2,3),y = c(2,1,2)),T),"square")
  
})

test_that("enclosing_radius is computing properly",{
  
  X <- data.frame(x = c(1:10),y = c(1:10))
  dist_X <- as.matrix(dist(X))
  expect_equal(enclosing_radius(X, F), dist_X[1,6])
  expect_equal(enclosing_radius(dist_X, T), dist_X[1,6])
  
  theta <- runif(n = 100,min = 0,max = 2*pi)
  x <- cos(theta)
  y <- sin(theta)
  df <- data.frame(x = x,y = y)
  dist_df <- as.matrix(dist(df))
  expect_equal(enclosing_radius(df, F),enclosing_radius(dist_df, T))
    
})
