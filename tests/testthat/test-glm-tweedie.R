library(statmod)

test_that("Testing tweedie var 1.5, link 1", {
  X <- model.matrix(hp ~ -1 + drat + cyl, data = mtcars)
  y <- mtcars$hp
  res <- as.numeric(glm_fit(X, y, tweedie_family(1.5, 1)))
  expected <- coef(glm.fit(X, y, family = tweedie(1.5, 1)))
  expect_equivalent(res, expected)
})


test_that("Testing tweedie var 2, link 1", {
  X <- model.matrix(hp ~ -1 + drat + cyl, data = mtcars)
  y <- mtcars$hp
  res <- as.numeric(glm_fit(X, y, tweedie_family(2.0, 1)))
  expected <- coef(glm.fit(X, y, family = tweedie(2.0, 1)))
  expect_equivalent(res, expected, tolerance = 1e-6)
})

test_that("Testing tweedie var 1.5, link 0", {
  X <- model.matrix(hp ~ -1 + drat + cyl, data = mtcars)
  y <- mtcars$hp
  res <- as.numeric(glm_fit(X, y, tweedie_family(1.5, 0)))
  expected <- coef(glm.fit(X, y, family = tweedie(1.5, 0)))
  expect_equivalent(res, expected)
})
