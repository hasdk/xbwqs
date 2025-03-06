# Testing DBWQS model fitting with R STAN
set.seed(9987)
C=5
K=3
J=4
N=200

simdt <- gen_dbwqs_simdata(C=C, K=K, J=J, N=N)
obsdf <- simdt$obsdata
obsdf$Cov2 <- as.factor(rep(c('a','b','c','d'), 50))
outcomes <- paste0("Y",1:C)
exposures <- paste0("Ex",1:K)
covars <- paste0("Cov",1:J)
fit_dbwqs <- dbwqs(data=obsdf, out=outcomes, exp=exposures, cov=covars,
                   impute_zeros = TRUE)

test_that("DBWQS returned error_flag=0 for successful model fitting.",{
  expect_equal(fit_dbwqs$error_flag,0)
})
