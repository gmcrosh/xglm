library(glmnet)

test_that("Testing simple GLMnet", {
  x <- model.matrix(I(mpg < 20) ~ -1 + hp + drat + cyl, data = mtcars)
  y <- mtcars$mpg < 20
  family <- binomial_family(link = "logit")
  res <- as.numeric(glmnet_fit(x, y, family, lambda = 0.0001))
  expected <- as.numeric(glmnet(
    x, y, binomial(), lambda = 0.0001, standardize=FALSE, intercept=FALSE
  )$beta@x)
  expect_equivalent(res, expected, tolerance = 1e-2)
})

test_that("Testing simple GLMnet (different lambda)", {
  x <- model.matrix(I(mpg < 20) ~ -1 + hp + drat + cyl, data = mtcars)
  y <- mtcars$mpg < 20
  family <- binomial_family(link = "logit")
  res <- as.numeric(glmnet_fit(x, y, family, lambda = 0.1))
  expected <- as.numeric(glmnet(
    x, y, binomial(), lambda = 0.1, standardize=FALSE, intercept=FALSE
  )$beta)
  expect_equivalent(res, expected, tolerance = 1e-2)
})
