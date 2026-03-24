## ----------------------------------------------------------------------------
## plot / marginalPlot / partialPlot / vimpPlot tests
## ----------------------------------------------------------------------------

# Shared tiny continuous longitudinal fit
make_cont_fit <- function(seed = 2024) {
  set.seed(seed)
  sim <- boostmtree::simLong(n = 8, ntest = 0, N = 2, rho = 0.25,
                             model = 1, family = "Continuous", q = 0)
  fit <- boostmtree(
    x       = sim$dtaL$features,
    tm      = sim$dtaL$time,
    id      = sim$dtaL$id,
    y       = sim$dtaL$y,
    family  = "Continuous",
    M       = 8,
    K       = 3,
    nknots  = 3,
    d       = 1,
    nu      = 0.1,
    verbose = FALSE,
    cv.flag = TRUE
  )
  list(fit = fit, sim = sim)
}

make_binary_fit_plots <- function(seed = 2024) {
  set.seed(seed)
  sim <- boostmtree::simLong(n = 6, ntest = 0, N = 2, rho = 0.25,
                             model = 0, family = "Binary", q = 0)
  fit <- boostmtree(
    x = sim$dtaL$features, tm = sim$dtaL$time,
    id = sim$dtaL$id, y = sim$dtaL$y,
    family = "Binary", M = 6, K = 3, nknots = 4, d = 1,
    nu = 0.1, verbose = FALSE
  )
  list(fit = fit, sim = sim)
}

# ---------------------------------------------------------------------------
# plot.boostmtree — grow object
# ---------------------------------------------------------------------------

test_that("plot.boostmtree (grow, Continuous) writes a PDF without error", {
  out <- make_cont_fit()
  tmp <- tempfile(fileext = "")
  dir.create(tmp)
  expect_invisible(plot(out$fit, path_saveplot = tmp, Verbose = FALSE))
  pdf_files <- list.files(tmp, pattern = "\\.pdf$")
  expect_length(pdf_files, 1L)
})

test_that("plot.boostmtree (grow, Binary) writes a PDF without error", {
  out <- make_binary_fit_plots()
  tmp <- tempfile(fileext = "")
  dir.create(tmp)
  expect_invisible(plot(out$fit, path_saveplot = tmp, Verbose = FALSE))
  expect_length(list.files(tmp, pattern = "\\.pdf$"), 1L)
})

test_that("plot.boostmtree (predict object) writes a PDF without error", {
  out <- make_cont_fit()
  pred <- predict(out$fit)
  tmp <- tempfile(fileext = "")
  dir.create(tmp)
  expect_invisible(plot(pred, path_saveplot = tmp, Verbose = FALSE))
  expect_length(list.files(tmp, pattern = "\\.pdf$"), 1L)
})

test_that("plot.boostmtree stops on wrong class", {
  expect_error(plot.boostmtree(list()), regexp = "class")
})

# ---------------------------------------------------------------------------
# marginalPlot
# ---------------------------------------------------------------------------

test_that("marginalPlot returns list with p.obj and time for Continuous", {
  out <- make_cont_fit()
  nm <- colnames(out$sim$dtaL$features)[1]
  mp <- marginalPlot(out$fit, xvar.names = nm, tm.unq = 1, plot.it = FALSE)
  expect_type(mp, "list")
  expect_named(mp, c("p.obj", "l.obj", "time"))
  expect_null(mp$l.obj)
  expect_type(mp$p.obj, "list")
})

test_that("marginalPlot writes a PDF when plot.it = TRUE", {
  out <- make_cont_fit()
  tmp <- tempfile(fileext = "")
  dir.create(tmp)
  nm <- colnames(out$sim$dtaL$features)[1]
  mp <- marginalPlot(out$fit, xvar.names = nm, tm.unq = 1,
                     plot.it = TRUE, path_saveplot = tmp, Verbose = FALSE)
  expect_length(list.files(tmp, pattern = "\\.pdf$"), 1L)
  expect_type(mp$l.obj, "list")
})

test_that("marginalPlot stops when xvar.names don't match", {
  out <- make_cont_fit()
  expect_error(
    marginalPlot(out$fit, xvar.names = "nonexistent_var"),
    regexp = "variable names"
  )
})

test_that("marginalPlot stops on wrong object class", {
  expect_error(marginalPlot(list()), regexp = "class")
})

# ---------------------------------------------------------------------------
# partialPlot
# ---------------------------------------------------------------------------

test_that("partialPlot returns p.obj for Continuous", {
  out <- make_cont_fit()
  nm <- colnames(out$sim$dtaL$features)[1]
  pp <- partialPlot(out$fit, xvar.names = nm, tm.unq = 1,
                    plot.it = FALSE, npts = 5)
  expect_type(pp, "list")
  expect_named(pp, c("p.obj", "l.obj", "time"))
  expect_null(pp$l.obj)
  expect_true(is.matrix(pp$p.obj[[nm]]))
})

test_that("partialPlot writes a PDF when plot.it = TRUE", {
  out <- make_cont_fit()
  tmp <- tempfile(fileext = "")
  dir.create(tmp)
  nm <- colnames(out$sim$dtaL$features)[1]
  pp <- partialPlot(out$fit, xvar.names = nm, tm.unq = 1,
                    npts = 5, plot.it = TRUE,
                    path_saveplot = tmp, Verbose = FALSE)
  expect_length(list.files(tmp, pattern = "\\.pdf$"), 1L)
})

test_that("partialPlot stops when xvar.names don't match", {
  out <- make_cont_fit()
  expect_error(
    partialPlot(out$fit, xvar.names = "no_such_var"),
    regexp = "variable names"
  )
})

test_that("partialPlot stops on wrong object class", {
  expect_error(partialPlot(list()), regexp = "class")
})

test_that("partialPlot conditional.xvars length mismatch errors", {
  out <- make_cont_fit()
  nm <- colnames(out$sim$dtaL$features)[1]
  nm2 <- colnames(out$sim$dtaL$features)[2]
  expect_error(
    partialPlot(out$fit, xvar.names = nm,
                conditional.xvars = nm2,
                conditional.values = c(0, 1)),   # length mismatch
    regexp = "same length"
  )
})

# ---------------------------------------------------------------------------
# vimpPlot
# ---------------------------------------------------------------------------

test_that("vimpPlot writes a PDF and returns invisibly", {
  out <- make_binary_fit_plots()
  pred <- predict(out$fit,
                  x  = out$sim$dtaL$features,
                  tm = out$sim$dtaL$time,
                  id = out$sim$dtaL$id,
                  y  = out$sim$dtaL$y)
  v <- vimp.boostmtree(pred, x.names = colnames(out$sim$dtaL$features)[1:2])
  tmp <- tempfile(fileext = "")
  dir.create(tmp)
  expect_invisible(
    vimpPlot(v, xvar.names = colnames(out$sim$dtaL$features)[1:2],
             path_saveplot = tmp, Verbose = FALSE)
  )
  expect_length(list.files(tmp, pattern = "\\.pdf$"), 1L)
})

test_that("vimpPlot stops when vimp is NULL", {
  expect_error(vimpPlot(NULL), regexp = "vimp")
})
