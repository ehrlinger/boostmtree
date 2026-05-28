test_that("DiagMat handles scalars and vectors", {
  expect_equal(boostmtree:::DiagMat(5), matrix(5, nrow = 1))
  expect_equal(boostmtree:::DiagMat(c(2, 3)), diag(c(2, 3)))
})

test_that("GetMu honors family-specific links", {
  lin <- c(-1, 0, 1)
  expect_equal(boostmtree:::GetMu(lin, "Continuous"), lin)
  expect_equal(
    boostmtree:::GetMu(lin, "Binary"),
    plogis(lin)
  )
  expect_equal(
    boostmtree:::GetMu(lin, "Nominal"),
    exp(lin)
  )
})

test_that("RemoveMiss.Fun drops all-NA rows and columns", {
  mat <- matrix(c(NA, NA, 1, 2, 3, 4), nrow = 3, byrow = TRUE)
  res <- boostmtree:::RemoveMiss.Fun(mat)
  expect_equal(nrow(res$X), 2)
  expect_equal(ncol(res$X), 2)
  expect_equal(res$id.remove, 1)
})

test_that("AppoxMatch locates nearest indices", {
  source <- c(0, 2, 4)
  target <- seq(0, 5, by = 1)
  expect_equal(boostmtree:::AppoxMatch(source, target), c(1, 3, 5))
})

test_that("rho.inv matches analytic inverse correlation", {
  ni <- 5
  rho <- 0.2
  expected <- rho / (1 + (ni - 1) * rho)
  expect_equal(boostmtree:::rho.inv(ni, rho), expected)
})
