test_that("simLong generates coherent binary longitudinal data", {
  sim <- boostmtree::simLong(n = 5, ntest = 2, N = 2, rho = 0.3,
                             model = 2, family = "Binary", q = 3)
  expect_equal(length(unique(sim$dta$id)), 7)
  expect_equal(ncol(sim$dtaL$features), 7) # 4 signal + q noise
  expect_true(all(sim$dta$y %in% c(0, 1)))
  expect_equal(length(sim$trn), sum(sim$dta$id <= 5))
  expect_match(sim$f.true, "time\\^2\\s*\\*\\s*x2\\^2")
})

test_that("simLong continuous responses retain numeric scale", {
  sim <- boostmtree::simLong(n = 4, ntest = 0, N = 3, rho = 0.1,
                             model = 1, family = "Continuous", q = 0)
  expect_type(sim$dta$y, "double")
  expect_gt(sd(sim$dta$y), 0)
})
