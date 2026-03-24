## ----------------------------------------------------------------------------
## simLong tests
## ----------------------------------------------------------------------------

test_that("simLong generates coherent binary longitudinal data", {
  sim <- boostmtree::simLong(n = 5, ntest = 2, N = 2, rho = 0.3,
                             model = 2, family = "Binary", q = 3)
  expect_equal(length(unique(sim$dta$id)), 7)
  expect_equal(ncol(sim$dtaL$features), 7)   # 4 signal + 3 noise
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

test_that("simLong trn indices are within the correct range", {
  sim <- boostmtree::simLong(n = 10, ntest = 5, N = 2,
                             model = 0, family = "Continuous", q = 0)
  trn_ids <- unique(sim$dta$id[sim$trn])
  expect_true(all(trn_ids <= 10))
  expect_true(all(trn_ids >= 1))
})

test_that("simLong model 0 formula contains only main effects", {
  sim <- boostmtree::simLong(n = 5, ntest = 0, N = 2,
                             model = 0, family = "Continuous")
  expect_false(grepl("time", sim$f.true))
})

test_that("simLong model 3 formula has three-way interaction", {
  sim <- boostmtree::simLong(n = 5, ntest = 0, N = 2,
                             model = 3, family = "Continuous")
  expect_match(sim$f.true, "time\\^2.*x2\\^2.*x3")
})

test_that("simLong stops on invalid family", {
  expect_error(
    boostmtree::simLong(n = 5, family = "Poisson"),
    regexp = "family"
  )
})

test_that("simLong with q = 0 returns exactly 4 features", {
  sim <- boostmtree::simLong(n = 5, ntest = 0, N = 2, q = 0,
                             model = 0, family = "Continuous")
  expect_equal(ncol(sim$dtaL$features), 4L)
})

test_that("simLong dtaL and dta are consistent", {
  sim <- boostmtree::simLong(n = 6, ntest = 0, N = 3,
                             model = 1, family = "Binary", q = 1)
  expect_equal(nrow(sim$dta), nrow(sim$dtaL$features))
  expect_equal(length(sim$dtaL$y), nrow(sim$dta))
  expect_equal(length(sim$dtaL$id), nrow(sim$dta))
  expect_equal(length(sim$dtaL$time), nrow(sim$dta))
})

test_that("simLong returns all four list elements", {
  sim <- boostmtree::simLong(n = 4, ntest = 2, N = 2,
                             model = 0, family = "Continuous")
  expect_named(sim, c("dtaL", "dta", "trn", "f.true"))
  expect_named(sim$dtaL, c("features", "time", "id", "y"))
})
