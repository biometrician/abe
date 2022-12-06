
# simulate example data
set.seed(1)
n <- 100
x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
y <- -5 + 5*x1 + 5*x2 + rnorm(n, sd = 5)
dd <- data.frame(y, x1, x2, x3)

# test if abe selects the correct model
test_that("correct model is selected", {
  # use abe on full model
  fit <- lm(y ~ ., data = dd, x = TRUE, y = TRUE)
  abe.fit <- abe(fit, data = dd, include = "x1", active = "x2", tau = 0.05, exp.beta = FALSE,
                 exact = TRUE, criterion = "alpha", alpha = 0.2, type.test = "Chisq", verbose = TRUE)

  # fit expected model (x3 dropped)
  fit_expected <- lm(y ~ x1 + x2, data = dd, x = TRUE, y = TRUE)

  # check if their terms are equal
  expect_equal(abe.fit$terms, fit_expected$terms)
})


