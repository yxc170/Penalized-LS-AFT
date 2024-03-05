library(microbenchmark)
library(Rcpp)
library(survival)

## #############################################################################
## random data (need to try this on our simulation settings and real data)
## #############################################################################
n <- 100
p <- 10
k <- 4

set.seed(1)
y <- rexp(n * k)
x <- matrix(rnorm(n * k * p), n * k, p)
d <- sample(0:1, n * k, T, c(.2, .8))
w <- rep(rexp(n), each = k)
## w <- rep(1, n * k) 
b0 <- runif(p)
id <- rep(1:n, each = k)
nt <- rep(4, n)


## #############################################################################
## compare point estimates with there is no penalty
## Use args() to see the arguments
## Use str() to see the output structure
## #############################################################################

## seems like we might need to consider higher tolerance for aftgee() and
## lower tolerance for gee()

## aftgeeEst() is used when there is no penalty
paft:::aftgeeEst(log(y), x, d, b0, nt, w, "ind", 1e-4, 500L)$b
paft:::aftgeeEst(log(y), x, d, b0, nt, w, "ex", 1e-4, 500L)$b

## alternatively, paftgeeEst() can be used with a proper pindex or with lambda = 0
pindex1 <- rep(1, p) ## don't penalize; results should match aftgeeEst()
paft:::paftgeeEst(log(y), x, d, b0, nt, pindex1, w, "lasso", "ind", .3, 1e-4, 1e-4, 500L)$b
paft:::paftgeeEst(log(y), x, d, b0, nt, pindex1, w, "lasso", "ex", .3, 1e-4, 1e-4, 500L)$b
paft:::paftgeeEst(log(y), x, d, b0, nt, pindex1, w, "scad", "ind", .3, 1e-4, 1e-4, 500L)$b
paft:::paftgeeEst(log(y), x, d, b0, nt, pindex1, w, "scad", "ex", .3, 1e-4, 1e-4, 500L)$b

pindex2 <- rep(0, p) ## penalize all
paft:::paftgeeEst(log(y), x, d, b0, nt, pindex2, w, "lasso", "ind", 0, 1e-4, 1e-4, 500L)$b
paft:::paftgeeEst(log(y), x, d, b0, nt, pindex2, w, "lasso", "ex", 0, 1e-4, 1e-4, 500L)$b
paft:::paftgeeEst(log(y), x, d, b0, nt, pindex2, w, "scad", "ind", 0, 1e-4, 1e-4, 500L)$b
paft:::paftgeeEst(log(y), x, d, b0, nt, pindex2, w, "scad", "ex", 0, 1e-4, 1e-4, 500L)$b

## compare to exisiting packages
coef(aftgee::aftgee(Surv(y, d) ~ x - 1, B = 0, id = id, weights = w, binit = b0))
paft:::paftgee(Surv(y, d) ~ x - 1, id = id, weights = w)$beta

## approaches a, b, c are faster than d
microbenchmark(a = paft:::aftgeeEst(log(y), x, d, b0, nt, w, "ind", 1e-4, 500L)$b,
               b = paft:::paftgeeEst(log(y), x, d, b0, nt, pindex1, w, "lasso", "ind", .3, 1e-4, 1e-4, 500L)$b,
               c = paft:::paftgeeEst(log(y), x, d, b0, nt, pindex2, w, "lasso", "ind", 0, 1e-4, 1e-4, 500L)$b,
               d = paft:::paftgee(Surv(y, d) ~ x - 1, id = id, weights = w)$beta)



## #############################################################################
## compare point estimates with there is penality
## #############################################################################

## A single lambda

pindex2 <- rep(0, p) ## penalize all 
paft:::paftgeeEst(log(y), x, d, b0, nt, pindex2, w, "lasso", "ind", .3, 1e-4, 1e-4, 500L)$b
paft:::paftgeeEst(log(y), x, d, b0, nt, pindex2, w, "lasso", "ex", .3, 1e-4, 1e-4, 500L)$b
paft:::paftgeeEst(log(y), x, d, b0, nt, pindex2, w, "scad", "ind", .3, 1e-4, 1e-4, 500L)$b
paft:::paftgeeEst(log(y), x, d, b0, nt, pindex2, w, "scad", "ex", .3, 1e-4, 1e-4, 500L)$b
paft:::paftgee(Surv(y, d) ~ x - 1, id = id, weights = w, lambda = .3, control = list(tol = 1e-4, init = b0))$beta

## approaches a is faster than b, but the point estimates are slightly different
microbenchmark(a = paft:::paftgeeEst(log(y), x, d, b0, nt, pindex2, w, "lasso", "ind", .3, 1e-4, 1e-4, 500L)$b, 
               b = paft:::paftgee(Surv(y, d) ~ x - 1, id = id, weights = w, lambda = .3, control = list(tol = 1e-4, init = b0))$beta)

## We see they converged at different steps despite the same tolerance. Please check.
str(paft:::paftgeeEst(log(y), x, d, b0, nt, pindex2, w, "lasso", "ind", .3, 1e-4, 1e-4, 500L))
str(paft:::paftgee(Surv(y, d) ~ x - 1, id = id, weights = w, lambda = .3, control = list(tol = 1e-4, init = b0)))

## A vector of lambda;
## in this case, the user needs to additionally specify either bic, cv, or 1se is used to choose the best lambda

lambda <- 1:20 / 10
pindex2 <- rep(0, p) ## penalize all 

## Different rules paftgeeEst1() uses to prepare for choosing the best lambda
str(paft:::paftgeeCV(log(y), x, d, b0, nt, pindex2, w, "lasso", 5, lambda, 1e-4, 1e-4, 500L))
str(paft:::paftgeeCV(log(y), x, d, b0, nt, pindex2, w, "scad", 5, lambda, 1e-3, 1e-3, 500L))
str(paft:::paftgeeBIC(log(y), x, d, b0, nt, pindex2, w, "scad", lambda, 1e-3, 1e-3, 500L))

## paftgeeEst1() gives the point estimate after selecting the best lambda
## Here are some examples

str(paft:::paftgeeEst1(log(y), x, d, b0, nt, pindex2, w, "lasso", 5, "ind", lambda, "bic", 1e-3, 1e-3, 500L))
str(paft:::paftgeeEst1(log(y), x, d, b0, nt, pindex2, w, "lasso", 5, "ind", lambda, "cv", 1e-3, 1e-3, 500L))
str(paft:::paftgeeEst1(log(y), x, d, b0, nt, pindex2, w, "scad", 5, "ind", lambda, "1se", 1e-3, 1e-3, 500L))
str(paft:::paftgeeEst1(log(y), x, d, b0, nt, pindex2, w, "lasso", 5, "ex", lambda, "bic", 1e-3, 1e-3, 500L))
str(paft:::paftgeeEst1(log(y), x, d, b0, nt, pindex2, w, "lasso", 5, "ex", lambda, "cv", 1e-3, 1e-3, 500L))
str(paft:::paftgeeEst1(log(y), x, d, b0, nt, pindex2, w, "scad", 5, "ex", lambda, "1se", 1e-3, 1e-3, 500L))


## results are different; please check
str(paft:::paftgee(Surv(y, d) ~ x - 1, id = id, weights = w, control = list(tol = 1e-3, init = b0, cv.rule = "bic", lambda.vec = lambda)))

## There are error messages here to clean up
str(paft:::paftgee(Surv(y, d) ~ x - 1, id = id, weights = w, control = list(tol = 1e-3, init = b0, nfold = 5, lambda.vec = lambda)))
str(paft:::paftgee(Surv(y, d) ~ x - 1, id = id, weights = w, control = list(tol = 1e-3, init = b0, cv.rule = "1se", nfold = 5, lambda.vec = lambda)))

