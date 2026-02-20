test_that("boostmtree fits and predicts small binary data", {
  set.seed(123)
  sim <- boostmtree::simLong(n = 8, ntest = 0, N = 2, rho = 0.2,
                             model = 1, family = "Binary", q = 1)
  y_factor <- factor(ifelse(sim$dtaL$y == 1, "yes", "no"))

  fit <- boostmtree(
    x = sim$dtaL$features,
    tm = sim$dtaL$time,
    id = sim$dtaL$id,
    y = y_factor,
    family = "Binary",
    M = 8,
    nu = 0.15,
    nknots = 3,
    d = 1,
    K = 3,
    verbose = FALSE,
    cv.flag = TRUE
  )

  expect_s3_class(fit, "boostmtree")
  expect_equal(fit$family, "Binary")
  expect_true(fit$Mopt[1] <= fit$M)

  pred <- predict(fit)
  expect_equal(pred$family, "Binary")
  expect_null(pred$rmse)
  expect_true(all(is.finite(unlist(pred$mu))))

  v <- vimp.boostmtree(fit, x.names = colnames(sim$dtaL$features)[1:2])
  expect_equal(length(v), 3)
  expect_true(all(rownames(v[[1]]) %in% colnames(sim$dtaL$features)[1:2]))
})
