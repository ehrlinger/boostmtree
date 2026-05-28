test_that("boostmtree handles nominal responses", {
  set.seed(404)
  x <- data.frame(
    x1 = rep(c(0, 1), each = 3),
    x2 = rep(c(1, 2, 3), 2)
  )
  y <- factor(c("red", "blue", "green", "blue", "green", "red"))

  fit <- boostmtree(
    x = x,
    y = y,
    family = "Nominal",
    M = 5,
    K = 2,
    nu = 0.2,
    verbose = FALSE
  )

  expect_equal(fit$family, "Nominal")
  expect_equal(length(fit$Q_set), length(levels(y)) - 1)

  pred <- predict(fit)
  expect_equal(length(pred$mu), fit$n.Q)
  expect_equal(length(pred$mu[[1]]), fit$n)
})

test_that("boostmtree handles ordinal responses with time trends", {
  x <- data.frame(
    x1 = rep(c(-1, 0, 1), each = 2),
    x2 = rep(c(0, 1), times = 3)
  )
  tm <- rep(c(0, 1), times = 3)
  id <- rep(1:3, each = 2)
  y <- ordered(c(1, 1, 2, 2, 3, 3))

  fit <- boostmtree(
    x = x,
    tm = tm,
    id = id,
    y = y,
    family = "Ordinal",
    M = 6,
    K = 3,
    nknots = 3,
    d = 1,
    nu = 0.15,
    verbose = FALSE
  )

  expect_equal(fit$family, "Ordinal")
  expect_false(is.null(fit$Prob_class))

  pred <- predict(fit, x = x, tm = tm, id = id, y = as.numeric(y))
  expect_equal(length(pred$Prob_class), length(levels(y)))
  expect_true(all(vapply(pred$Prob_class, length, integer(1)) == fit$n))
})
