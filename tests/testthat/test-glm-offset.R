library(GLMsData)
data("cervical")

test_that("Testing poisson with offset", {
  x <- model.matrix(Deaths ~ Age + Country, data = cervical)
  y <- cervical$Deaths
  offset <- log(cervical$Wyears)
  family <- poisson_family(link = "log")
  res <- as.numeric(glm_fit(x, y, family, offset = offset))
  expected <- coef(glm.fit(x, y, offset = offset, family = stats::poisson()))
  expect_equivalent(res, expected)
})
