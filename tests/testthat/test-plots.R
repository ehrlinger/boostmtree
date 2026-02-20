test_that("plotting helpers execute for binary models", {
  set.seed(2024)
  sim <- boostmtree::simLong(n = 6, ntest = 0, N = 2, rho = 0.25,
                             model = 0, family = "Binary", q = 0)
  fit <- boostmtree(
    x = sim$dtaL$features,
    tm = sim$dtaL$time,
    id = sim$dtaL$id,
    y = sim$dtaL$y,
    family = "Binary",
    M = 6,
    K = 3,
    nknots = 4,
    d = 1,
    nu = 0.1,
    verbose = FALSE
  )

  pp <- boostmtree::partialPlot(fit, xvar.names = colnames(sim$dtaL$features)[1],
                                tm.unq = 1, plot.it = FALSE)
  expect_type(pp$p.obj, "list")

  mp <- boostmtree::marginalPlot(fit, xvar.names = colnames(sim$dtaL$features)[1],
                                 tm.unq = 1, plot.it = FALSE)
  expect_type(mp$p.obj, "list")

  tmp_dir <- tempfile(pattern = "plots")
  dir.create(tmp_dir)
  expect_invisible(plot(fit, path_saveplot = tmp_dir, Verbose = FALSE))

  pred <- predict(fit, x = sim$dtaL$features, tm = sim$dtaL$time,
                  id = sim$dtaL$id, y = sim$dtaL$y)
  v <- boostmtree::vimp.boostmtree(pred, x.names = colnames(sim$dtaL$features)[1:2])
  expect_invisible(boostmtree::vimpPlot(v,
                                        xvar.names = colnames(sim$dtaL$features)[1:2],
                                        path_saveplot = tmp_dir,
                                        Verbose = FALSE))
})
