# check_rstan_compat() should only stop for the rstan < 2.39 / StanHeaders >= 2.33 combination
test_that("check_rstan_compat() stops for old rstan with new StanHeaders", {
  expect_error(xbwqs:::check_rstan_compat("2.32.7", "2.39.1.1"),
               "Incompatible Stan installation")
})

test_that("check_rstan_compat() passes for compatible version pairs", {
  expect_true(xbwqs:::check_rstan_compat("2.39.0.9000", "2.39.1.1"))
  expect_true(xbwqs:::check_rstan_compat("2.32.7", "2.32.2"))
})

test_that("check_rstan_compat() passes for the installed versions", {
  expect_true(xbwqs:::check_rstan_compat())
})
