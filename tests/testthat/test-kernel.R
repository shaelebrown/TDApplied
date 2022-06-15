
test_that("diagram_kernel detects incorrect parameters correctly",{
  
  D = data.frame(dimension = c(0),birth = c(0),death = c(1))
  expect_error(diagram_kernel(D1 = NULL,D2 = D,dim = 1),"TDA computation or data frame")
  expect_error(diagram_kernel(D1 = D,D2 = NULL,dim = 1),"TDA computation or data frame")
  expect_error(diagram_kernel(D1 = D,D2 = D,dim = "2"),"whole number")
  expect_error(diagram_kernel(D1 = D,D2 = D,sigma = "2"),"number")
  expect_error(diagram_kernel(D1 = D,D2 = D,t = NA),"number")

})

test_that("diagram_kernel is computing correctly",{
  
  D1 = data.frame(dimension = 0,birth = 2,death = 3)
  D2 = data.frame(dimension = 0,birth = c(2,5),death = c(3.1,6))
  #expect_identical(diagram_kernel(D1,D2,dim = 0,sigma = 1,t = 1),sqrt(0.1^2+0.5^2))
  #expect_identical(diagram_kernel(D1,D2,dim = 0,sigma = 2,t = 1),(0.1^3+0.5^3)^(1/3))
  #expect_equal(diagram_kernel(D1 = D1,D2 = D2,dim = 0,sigma = 1,t = 2),0.9984164)
  #expect_identical(diagram_kernel(D1 = D1,D2 = D2,sigma = 2,t = 2),0.5)
  
})


