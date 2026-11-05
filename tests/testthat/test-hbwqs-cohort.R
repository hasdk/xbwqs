# cohort_index() maps the hbwqs `cohort` argument to consecutive integers 1..J
df <- data.frame(site = c(9, 2, 5, 2, 9), lab = c("SF", "NY", "LA", "NY", "SF"))

test_that("a column name with non-consecutive numeric codes maps to 1..J", {
  res <- xbwqs:::cohort_index(df, "site")
  expect_equal(res$index, c(3L, 1L, 2L, 1L, 3L))
  expect_equal(res$labels, c("2", "5", "9"))
})

test_that("a column of character labels maps to 1..J", {
  res <- xbwqs:::cohort_index(df, "lab")
  expect_equal(res$index, c(3L, 2L, 1L, 2L, 3L))
  expect_equal(res$labels, c("LA", "NY", "SF"))
})

test_that("a factor vector with one entry per row is accepted", {
  res <- xbwqs:::cohort_index(df, factor(df$lab, levels = c("NY", "LA", "SF")))
  expect_equal(res$index, c(3L, 1L, 2L, 1L, 3L))
  expect_equal(res$labels, c("NY", "LA", "SF"))
})

test_that("NULL treats all rows as a single cohort", {
  res <- xbwqs:::cohort_index(df, NULL)
  expect_equal(res$index, rep(1L, 5))
  expect_equal(res$labels, "1")
})

test_that("invalid cohort inputs give informative errors", {
  expect_error(xbwqs:::cohort_index(df, "cohort_id"), "not found in `data`")
  expect_error(xbwqs:::cohort_index(df, c(1, 2)), "one entry per row")
  expect_error(xbwqs:::cohort_index(df, c(1, 2, NA, 1, 2)), "missing values")
})
