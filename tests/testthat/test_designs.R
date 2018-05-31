test_that("Test that t-design runs through", {

  optimal_design(effect=0.3, alpha=0.05, pow=0.8, t_approx=T)

  expect_equal(1,1)

})


test_that("Test that direct design runs through", {

  optimal_design(effect=0.3, alpha=0.05, pow=0.8, t_approx=F, lagrange=F)

  expect_equal(1,1)

})


test_that("Test that Lagrange design runs through", {

  optimal_design(effect=0.7, alpha=0.05, pow=0.8, effect_null=0.1, sd=2, standardized=F, t_approx=F, lagrange=T)

  expect_equal(1,1)

})
