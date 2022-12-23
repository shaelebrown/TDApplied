
test_that("utilities are working properly",{
  
  expect_error(check_diagram(data.frame(dimension = c(1,2,3),birth = c("1","2","3"),death = c(1,2,3))),"numeric")
  expect_error(check_diagram(data.frame(dimension = c(1.1,2,3),birth = c(1,2,3),death = c(1,2,3))),"whole")
  expect_error(check_diagram(data.frame(dimension = c(-1,2,3),birth = c(1,2,3),death = c(1,2,3))),">= 0")
  expect_error(check_diagram(data.frame(dimension = c(1,2,3),birth = c(1,-2,3),death = c(1,2,3))),">= 0")
  expect_error(check_diagram(data.frame(dimension = c(1,2,3),birth = c(1,2,3),death = c(1,2,NA))),"missing")
  expect_error(check_diagram(data.frame(dimension = c(1,2,3),birth = c(1,2,3),death = c(1,2,2.9))),"at least")
  expect_error(check_param(param_name = "test",param = "T",numeric = F),"T or F")
  
})