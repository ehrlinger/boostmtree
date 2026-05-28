## ----------------------------------------------------------------------------
## boostmtree core fit/predict/vimp tests
## ----------------------------------------------------------------------------

# Shared small dataset (binary, longitudinal) used across multiple tests.
make_binary_fit <- function(seed = 123) {
  set.seed(seed)
  sim <- boostmtree::simLong(
    n = 8, ntest = 0, N = 2, rho = 0.2,
    model = 1, family = "Binary", q = 1
  )
  y_factor <- factor(ifelse(sim$dtaL$y == 1, "yes", "no"))
  fit <- boostmtree(
    x       = sim$dtaL$features,
    tm      = sim$dtaL$time,
    id      = sim$dtaL$id,
    y       = y_factor,
    family  = "Binary",
    M       = 8,
    nu      = 0.15,
    nknots  = 3,
    d       = 1,
    K       = 3,
    verbose = FALSE,
    cv.flag = TRUE
  )
  list(fit = fit, sim = sim)
}

# ---------------------------------------------------------------------------
# Fit properties
# ---------------------------------------------------------------------------

test_that("boostmtree returns a correctly-classed object (Binary)", {
  out <- make_binary_fit()
  fit <- out$fit
  expect_s3_class(fit, "boostmtree")
  expect_equal(fit$family, "Binary")
  expect_true(fit$Mopt[1] <= fit$M)
  expect_true(fit$Mopt[1] >= 1L)
})

test_that("boostmtree fit contains required slots", {
  out <- make_binary_fit()
  fit <- out$fit
  expected_slots <- c("mu", "x", "time", "id", "family", "M", "Mopt",
                       "n", "K", "nu", "err.rate", "rmse")
  for (s in expected_slots) {
    expect_true(!is.null(fit[[s]]), info = paste("missing slot:", s))
  }
})

test_that("boostmtree Continuous fit returns finite predictions", {
  set.seed(7)
  sim <- boostmtree::simLong(n = 10, ntest = 0, N = 3, model = 0,
                             family = "Continuous", q = 0)
  fit <- boostmtree(
    x = sim$dtaL$features, tm = sim$dtaL$time,
    id = sim$dtaL$id, y = sim$dtaL$y,
    family = "Continuous", M = 10, K = 3, nknots = 3,
    verbose = FALSE
  )
  expect_s3_class(fit, "boostmtree")
  expect_equal(fit$family, "Continuous")
  mu_vals <- unlist(fit$mu)
  expect_true(all(is.finite(mu_vals)))
})

# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

test_that("boostmtree stops on invalid family", {
  sim <- boostmtree::simLong(n = 6, ntest = 0, N = 2,
                             model = 0, family = "Continuous", q = 0)
  expect_error(
    boostmtree(sim$dtaL$features, sim$dtaL$time, sim$dtaL$id,
               sim$dtaL$y, family = "Gaussian", M = 5, verbose = FALSE),
    regexp = "family"
  )
})

test_that("boostmtree stops when multiple families are supplied", {
  sim <- boostmtree::simLong(n = 6, ntest = 0, N = 2,
                             model = 0, family = "Continuous", q = 0)
  expect_error(
    boostmtree(sim$dtaL$features, sim$dtaL$time, sim$dtaL$id,
               sim$dtaL$y, family = c("Continuous", "Binary"), M = 5,
               verbose = FALSE),
    regexp = "famil"
  )
})

# ---------------------------------------------------------------------------
# Prediction
# ---------------------------------------------------------------------------

test_that("predict.boostmtree returns correct class on training data", {
  out <- make_binary_fit()
  pred <- predict(out$fit)
  expect_s3_class(pred, "boostmtree")
  expect_true(inherits(pred, "predict"))
  expect_equal(pred$family, "Binary")
  expect_null(pred$rmse)   # no y supplied → no error metric
  expect_true(all(is.finite(unlist(pred$mu))))
})

test_that("predict.boostmtree computes RMSE when y is supplied", {
  set.seed(42)
  sim <- boostmtree::simLong(n = 10, ntest = 5, N = 3, model = 0,
                             family = "Continuous", q = 0)
  fit <- boostmtree(
    x = sim$dtaL$features[sim$trn, ], tm = sim$dtaL$time[sim$trn],
    id = sim$dtaL$id[sim$trn], y = sim$dtaL$y[sim$trn],
    family = "Continuous", M = 10, K = 3, nknots = 3, verbose = FALSE
  )
  pred <- predict(
    fit,
    x  = sim$dtaL$features[-sim$trn, ],
    tm = sim$dtaL$time[-sim$trn],
    id = sim$dtaL$id[-sim$trn],
    y  = sim$dtaL$y[-sim$trn]
  )
  expect_true(!is.null(pred$rmse))
  expect_true(is.finite(pred$rmse))
  expect_gt(pred$rmse, 0)
})

test_that("predict.boostmtree stops on wrong object class", {
  expect_error(predict.boostmtree(list(x = 1)), regexp = "class")
})

# ---------------------------------------------------------------------------
# VIMP
# ---------------------------------------------------------------------------

test_that("vimp.boostmtree returns a 3-element list for Binary grow", {
  out <- make_binary_fit()
  v <- vimp.boostmtree(out$fit,
                       x.names = colnames(out$sim$dtaL$features)[1:2])
  expect_equal(length(v), 3)
  expect_true(all(rownames(v[[1]]) %in% colnames(out$sim$dtaL$features)[1:2]))
})

test_that("vimp.boostmtree all-variable default works", {
  out <- make_binary_fit()
  v <- vimp.boostmtree(out$fit)
  expect_equal(nrow(v[[1]]), ncol(out$sim$dtaL$features))
})

test_that("joint vimp.boostmtree returns a single row", {
  out <- make_binary_fit()
  v <- vimp.boostmtree(out$fit,
                       x.names = colnames(out$sim$dtaL$features)[1:2],
                       joint   = TRUE)
  expect_equal(nrow(v[[1]]), 1L)
})

# ---------------------------------------------------------------------------
# print method
# ---------------------------------------------------------------------------

test_that("print.boostmtree succeeds for grow and predict objects", {
  out <- make_binary_fit()
  expect_output(print(out$fit), regexp = "boostmtree")
  pred <- predict(out$fit)
  expect_output(print(pred), regexp = "boostmtree")
})

test_that("print.boostmtree stops on wrong class", {
  expect_error(print.boostmtree(list()), regexp = "class")
})

# ---------------------------------------------------------------------------
# boostmtree.news
# ---------------------------------------------------------------------------

test_that("boostmtree.news runs without error", {
  expect_invisible(boostmtree.news())
})

# ---------------------------------------------------------------------------
# Nominal family
# ---------------------------------------------------------------------------

make_nominal_fit <- function(seed = 42) {
  set.seed(seed)
  n_subj <- 8
  n_obs  <- 2
  N      <- n_subj * n_obs
  id     <- rep(seq_len(n_subj), each = n_obs)
  tm     <- rep(c(0, 1), times = n_subj)
  x      <- matrix(rnorm(N * 4), nrow = N,
                   dimnames = list(NULL, paste0("x", 1:4)))
  y      <- factor(sample(c("A", "B", "C"), N, replace = TRUE))
  fit    <- boostmtree(
    x = x, tm = tm, id = id, y = y,
    family = "Nominal", M = 8, K = 3, nknots = 3,
    nu = 0.1, verbose = FALSE
  )
  list(fit = fit, x = x, tm = tm, id = id, y = y)
}

test_that("boostmtree Nominal fit returns correctly-classed object", {
  out <- make_nominal_fit()
  expect_s3_class(out$fit, "boostmtree")
  expect_equal(out$fit$family, "Nominal")
})

test_that("boostmtree Nominal fit contains required slots", {
  out <- make_nominal_fit()
  for (s in c("mu", "x", "time", "id", "family", "M", "n", "K",
              "n.Q", "Q_set")) {
    expect_true(!is.null(out$fit[[s]]), info = paste("missing slot:", s))
  }
})

test_that("boostmtree Nominal fit has n.Q equal to number of response classes minus 1", {
  out <- make_nominal_fit()
  expect_equal(out$fit$n.Q, length(levels(out$y)) - 1L)
})

test_that("predict.boostmtree returns finite mu for Nominal family", {
  out <- make_nominal_fit()
  pred <- predict(out$fit)
  expect_s3_class(pred, "boostmtree")
  mu_vals <- unlist(pred$mu)
  expect_true(all(is.finite(mu_vals)))
})

test_that("print.boostmtree works for Nominal grow object", {
  out <- make_nominal_fit()
  expect_output(print(out$fit), regexp = "Nominal")
})

# ---------------------------------------------------------------------------
# Ordinal family
# ---------------------------------------------------------------------------

make_ordinal_fit <- function(seed = 77) {
  set.seed(seed)
  n_subj <- 8
  n_obs  <- 2
  N      <- n_subj * n_obs
  id     <- rep(seq_len(n_subj), each = n_obs)
  tm     <- rep(c(0, 1), times = n_subj)
  x      <- matrix(rnorm(N * 4), nrow = N,
                   dimnames = list(NULL, paste0("x", 1:4)))
  y      <- factor(sample(c("low", "mid", "high"), N, replace = TRUE),
                   levels = c("low", "mid", "high"), ordered = TRUE)
  fit    <- boostmtree(
    x = x, tm = tm, id = id, y = y,
    family = "Ordinal", M = 8, K = 3, nknots = 3,
    nu = 0.1, verbose = FALSE
  )
  list(fit = fit, x = x, tm = tm, id = id, y = y)
}

test_that("boostmtree Ordinal fit returns correctly-classed object", {
  out <- make_ordinal_fit()
  expect_s3_class(out$fit, "boostmtree")
  expect_equal(out$fit$family, "Ordinal")
})

test_that("boostmtree Ordinal fit contains required slots including Prob_class", {
  out <- make_ordinal_fit()
  for (s in c("mu", "x", "time", "id", "family", "M", "n", "K",
              "n.Q", "Q_set", "Prob_class")) {
    expect_true(!is.null(out$fit[[s]]), info = paste("missing slot:", s))
  }
})

test_that("boostmtree Ordinal fit has n.Q equal to number of levels minus 1", {
  out <- make_ordinal_fit()
  expect_equal(out$fit$n.Q, length(levels(out$y)) - 1L)
})

test_that("predict.boostmtree returns finite mu for Ordinal family", {
  out <- make_ordinal_fit()
  pred <- predict(out$fit)
  expect_s3_class(pred, "boostmtree")
  mu_vals <- unlist(pred$mu)
  expect_true(all(is.finite(mu_vals)))
})

test_that("print.boostmtree works for Ordinal grow object", {
  out <- make_ordinal_fit()
  expect_output(print(out$fit), regexp = "Ordinal")
})

# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

make_valid_args <- function(seed = 99) {
  set.seed(seed)
  sim <- boostmtree::simLong(
    n = 6, ntest = 0, N = 2, rho = 0.2,
    model = 0, family = "Continuous", q = 1
  )
  list(
    x  = sim$dtaL$features,
    tm = sim$dtaL$time,
    id = sim$dtaL$id,
    y  = sim$dtaL$y
  )
}

test_that("boostmtree rejects mismatched y length", {
  a <- make_valid_args()
  expect_error(
    boostmtree(x = a$x, tm = a$tm, id = a$id, y = a$y[-1],
               family = "Continuous", verbose = FALSE),
    regexp = "length of y"
  )
})

test_that("boostmtree rejects mismatched id length", {
  a <- make_valid_args()
  expect_error(
    boostmtree(x = a$x, tm = a$tm, id = a$id[-1], y = a$y,
               family = "Continuous", verbose = FALSE),
    regexp = "length of id"
  )
})

test_that("boostmtree rejects mismatched tm length", {
  a <- make_valid_args()
  expect_error(
    boostmtree(x = a$x, tm = a$tm[-1], id = a$id, y = a$y,
               family = "Continuous", verbose = FALSE),
    regexp = "length of tm"
  )
})

test_that("boostmtree rejects empty x", {
  a <- make_valid_args()
  expect_error(
    boostmtree(x = a$x[0, , drop = FALSE], tm = numeric(0),
               id = integer(0), y = numeric(0),
               family = "Continuous", verbose = FALSE),
    regexp = "at least one observation"
  )
})

test_that("boostmtree rejects non-positive M", {
  a <- make_valid_args()
  expect_error(
    boostmtree(x = a$x, tm = a$tm, id = a$id, y = a$y,
               family = "Continuous", M = 0, verbose = FALSE),
    regexp = "M must be a positive integer"
  )
})

test_that("boostmtree rejects non-positive nu", {
  a <- make_valid_args()
  expect_error(
    boostmtree(x = a$x, tm = a$tm, id = a$id, y = a$y,
               family = "Continuous", nu = 0, verbose = FALSE),
    regexp = "nu must be a positive number"
  )
})

test_that("boostmtree rejects non-positive K", {
  a <- make_valid_args()
  expect_error(
    boostmtree(x = a$x, tm = a$tm, id = a$id, y = a$y,
               family = "Continuous", K = 0, verbose = FALSE),
    regexp = "K must be a positive integer"
  )
})

test_that("boostmtree rejects non-integer M", {
  a <- make_valid_args()
  expect_error(
    boostmtree(x = a$x, tm = a$tm, id = a$id, y = a$y,
               family = "Continuous", M = 1.5, verbose = FALSE),
    regexp = "M must be a positive integer"
  )
})

test_that("boostmtree rejects non-integer K", {
  a <- make_valid_args()
  expect_error(
    boostmtree(x = a$x, tm = a$tm, id = a$id, y = a$y,
               family = "Continuous", K = 2.5, verbose = FALSE),
    regexp = "K must be a positive integer"
  )
})

test_that("boostmtree rejects non-finite M", {
  a <- make_valid_args()
  expect_error(
    boostmtree(x = a$x, tm = a$tm, id = a$id, y = a$y,
               family = "Continuous", M = Inf, verbose = FALSE),
    regexp = "M must be a positive integer"
  )
})
