
test_that("utilities are working properly",{
  
  expect_error(check_diagram(data.frame(dimension = c(1,2,3),birth = c("1","2","3"),death = c(1,2,3))),"numeric")
  expect_error(check_diagram(data.frame(dimension = c(1.1,2,3),birth = c(1,2,3),death = c(1,2,3))),"whole")
  expect_error(check_diagram(data.frame(dimension = c(-1,2,3),birth = c(1,2,3),death = c(1,2,3))),">= 0")
  expect_error(check_diagram(data.frame(dimension = c(1,2,3),birth = c(1,-2,3),death = c(1,2,3))),">= 0")
  expect_error(check_diagram(data.frame(dimension = c(1,2,3),birth = c(1,2,3),death = c(1,2,NA))),"missing")
  expect_error(check_diagram(data.frame(dimension = c(1,2,3),birth = c(1,2,3),death = c(1,2,2.9))),"larger")
  expect_error(check_param(param_name = "test",param = "T",numeric = F),"T or F")
  
})

test_that("check_matrix works",{
  
  d1 = data.frame(dimension = rep(0,5),birth = 1:5,death = 1:5 + 0.1)
  d2 = data.frame(dimension = rep(0,5),birth = 1:5,death = 1:5 + 0.2)
  D = distance_matrix(list(d1,d2),dim = 0,num_workers = 2)
  K = gram_matrix(list(d1,d2),dim = 0,num_workers = 2)
  expect_error(check_matrix(D,"D"),"kernel")
  expect_error(check_matrix(K,"K","matrix"),"matrix")
  expect_error(check_matrix(rbind(D,c(1,2)),"D","matrix"),"rows")
  D[1,2] = NA
  D[2,1] = NaN
  expect_error(check_matrix(D,"D","matrix"),"missing")
  D = distance_matrix(list(d1,d2),dim = 0,num_workers = 2)
  D[1,1] = 1
  expect_error(check_matrix(D,"D","matrix"),"0's")
  D[1,1] = 0
  K[1,1] = 0
  expect_error(check_matrix(K,"K"),"1's")
  K[1,1] = 1
  K[1,2] = 1
  expect_error(check_matrix(K,"K"),"symmetric")
  D[1,2] = 0
  expect_error(check_matrix(D,"D","matrix"),"symmetric")
  expect_silent(check_matrix(D,"D",type = "matrix",symmetric = F))
  expect_error(check_matrix(D[0,],"D",type = "matrix"),"at least")
  
})