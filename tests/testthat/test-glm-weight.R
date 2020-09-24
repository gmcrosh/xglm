library(GLMsData)
data("germ")

test_that("Testing weighted binomial GLM", {
  x <- model.matrix(Germ ~ Extract + Seeds, data = germ)
  y <- germ$Germ / germ$Total
  weights <- germ$Total
  family <- binomial_family(link = "logit")
  res <- as.numeric(glm_fit(x, y, family, sample_weights = weights))
  expected <- coef(glm.fit(x, y, weights = weights, family = stats::binomial()))
  expect_equivalent(res, expected)
})
