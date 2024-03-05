library(paft)
library(geepack)

## Example from ?geese
gendat <- function() {
    id <- gl(50, 4, 200)
    visit <- rep(1:4, 50)
    x1 <- rbinom(200, 1, 0.6) ## within cluster varying binary covariate
    x2 <- runif(200, 0, 1)   ## within cluster varying continuous covariate
    phi <- 1 + 2 * x1         ## true scale model
    ## the true correlation coefficient rho for an ar(1)
    ## correlation structure is 0.667.
    rhomat <- 0.667 ^ outer(1:4, 1:4, function(x, y) abs(x - y))
    chol.u <- chol(rhomat)
    noise <- as.vector(sapply(1:50, function(x) chol.u %*% rnorm(4)))
    e <- sqrt(phi) * noise
    y <- 1 + 3 * x1 - 2 * x2 + e
    dat <- data.frame(y, id, visit, x1, x2)
    dat
}

set.seed(1); dat <- gendat()

fit1 <- geese(y ~ x1 + x2, id = id, data = dat, corstr = "ind")
fit2 <- geese(y ~ x1 + x2, id = id, data = dat, corstr = "ex")
fit3 <- geese(y ~ x1 + x2, id = id, data = dat, corstr = "ar1")

cbind(fit1$beta, fit2$beta, fit3$beta)

## paft returns zero? 
y <- dat$y
X <- model.matrix(~ x1 + x2, data = dat)
id <- dat$id
weights <- rep(1, length(y))
lambda <- 0
b0 <- rep(0, ncol(X))
pindex <- rep(0, ncol(X))

fit4 <- paft:::pgeeCpp(X, y, weights, id, lambda, b0, pindex, corstr = "independence")
fit5 <- paft:::pgeeCpp(X, y, weights, id, lambda, b0, pindex, corstr = "exchangeable")
fit6 <- paft:::pgeeCpp(X, y, weights, id, lambda, b0, pindex, corstr = "AR1")

cbind(fit4$b1, fit6$b1)
cbind(fit4$b1, fit5$b1, fit6$b1)
