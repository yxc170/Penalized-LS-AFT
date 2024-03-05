library(paft)
library(evd)
library(copula)

simDat <- function(n, cr, marg = c("sn", "sl", "sg")) {
    marg <- match.arg(marg)
    x1 <- rbinom(n, 1, .5)
    x1 <- (x1 - 0.5) / 2
    x2 <- rnorm(n)
    x3 <- rnorm(n)
    ## failure time before censoring
    if (marg == "sn") {
        Y <- exp(2 + 3 * x1 + 1 * x2  + 2.5 * x3 + rnorm(n))
        if (cr == 0) cen <- Inf
        if (cr == .9) cen <- runif(n, 0, 0.4)
        if (cr == .95) cen <- runif(n, 0, 0.15)
    }      
    if (marg == "sl") {
        err <- rlogis(n, location = 0, scale = 0.55)
        Y <- exp(2 + 3 * x1 + 1 * x2  + 2.5 * x3 + err)
        if (cr == 0) cen <- Inf
        if (cr == .9) cen <- runif(n, 0, 0.4)
        if (cr == .95) cen <- runif(n, 0, 0.15)
    }      
    if (marg == "sg") {
        err <- rgumbel(n, loc = -0.44, scale = 0.78)
        Y <- exp(2 + 3 * x1 + 1 * x2  + 2.5 * x3 + err)
        if (cr == 0) cen <- Inf
        if (cr == .9) cen <- runif(n, 0, 0.4)
        if (cr == .95) cen <- runif(n, 0, 0.15)
    }
    dat <- data.frame(Time = pmin(Y, cen), delta = 1 * (Y <= cen),
                      x1 = x1, x2 = x2, x3 = x3, x4 = rbinom(n, 1, .5),
                      x5 = rnorm(n), x6 = rnorm(n),
                      x7 = rexp(n), x8 = rexp(n),
                      x9 = runif(n), x10 = runif(n))
    ## Create case-cohort 
    srs <- sample(1:n, round(n / 9))
    dat <- dat[sort(unique(c(which(dat$delta > 0), srs))),]
    dat$weights <- ifelse(dat$delta > 0, 1, 1 / 9)
    dat$id <- 1:nrow(dat)
    rownames(dat) <- NULL
    return(dat)
}

## set.seed(1)
dat <- simDat(1000, .9, "sg")
head(dat)
## true beta is c(2, 3, 1, 2.5)
fm <- Surv(Time, delta) ~ x1 + x2 + x3

## Testing codes
paftgee(fm, data = dat)
paftgee(fm, data = dat, control = list(init = "lm"))
paftgee(fm, data = dat, control = list(init = "gehan_s"))
paftgee(fm, data = dat, control = list(init = "gehan_ns"))

paftgee(fm, data = dat, weights = weights)
paftgee(fm, data = dat, weights = weights, control = list(init = "lm"))
paftgee(fm, data = dat, weights = weights, control = list(init = "gehan_s"))
paftgee(fm, data = dat, weights = weights, control = list(init = "gehan_ns"))


paftgee(fm, data = dat, weights = weights, id = id)
paftgee(fm, data = dat, weights = weights, id = id, control = list(init = "lm"))
paftgee(fm, data = dat, weights = weights, id = id, control = list(init = "gehan_s"))
paftgee(fm, data = dat, weights = weights, id = id, control = list(init = "gehan_ns"))


paftgee(fm, data = dat, lambda = .5)
paftgee(fm, data = dat, weights = weights, lambda = .5)
paftgee(fm, data = dat, weights = weights, control = list(nfold = 10, lambda.vec = 1:10/10))
paftgee(fm, data = dat, weights = weights,
        control = list(nfold = 10, lambda.vec = 1:10/10, cv.rule = "1se"))
paftgee(fm, data = dat, weights = weights,
        control = list(nfold = 10, lambda.vec = 1:10/10, cv.rule = "bic"))

set.seed(1); paftgee(fm, data = dat, B = 5)
set.seed(1); paftgee(fm, data = dat, B = 5, control = list(parallel = TRUE))

do <- function() {
    dat <- simDat(1000, .9, "sg")
    fm <- Surv(Time, delta) ~ x1 + x2 + x3
    c(paftgee(fm, data = dat)$beta,
      paftgee(fm, data = dat, weights = weights)$beta)
}

foo <- replicate(100, do())

matrix(c(2, 3, 1, 2.5, rowMeans(foo)), 4)
matrix(apply(foo, 1, sd), 4)



cv <- cv.paft(x, y, d, y, id, ctrl$nfold, b0, ctrl$cv.rule, ctrl$lambda.vec, 
              ctrl$pindex, weights = w)

cv1 <- cv.paft(x, y, d, y, id, ctrl$nfold, b0, "1se", ctrl$lambda.vec, 
              ctrl$pindex, weights = w)
