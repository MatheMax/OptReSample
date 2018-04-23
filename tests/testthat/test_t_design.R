test_that("Test that t-Design runs through", {
  params<-list(
    mualt=0.3,
    mu0=0,
    sigma=1,
    alpha=0.05,
    beta=0.2
  )

  t_design(params)

  expect_equal(1,1)

})
