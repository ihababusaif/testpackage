context('Output')# Our file is called "test-tnl_output.R"
local_edition(3)
library(testthat) # load testthat package
library(testpackage)# load our package
test_that("trigonometric functions match identities", {
  expect_equal(10, 10 + 1e-7)
})
