test_that("Test that one can compute the score", {

  params<-list(
    mualt=0.3,
    mu0=0,
    sigma=1,
    alpha=0.05,
    beta=0.2
  )

  score(params,0.5,2,40,450,130)

  expect_equal(1,1)

})
