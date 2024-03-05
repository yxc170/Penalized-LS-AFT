library(paft)
library(aftgee)
library(copula)
library(survival)
library(evd)
library(Rcpp)
library(geepack)

simDat <- function(n, tau, cr, marg = "sn", cop = "clayton", full = FALSE) {
  k <- 3
  if (k > 1) {
    if (cop == "clayton")
      invisible(cop <- claytonCopula(iTau(claytonCopula(), tau = tau), dim = k))
    if (marg == "sn") err <- t(qnorm(rCopula(n, cop)))
    if (marg == "sl") {
      err <- t(qlogis(rCopula(n, cop), location = 0, scale = 0.55))
    }
    if (marg == "sg") {
      err <- t(qgumbel(rCopula(n, cop), loc = -0.44, scale = 0.78))
    }
  } else {
    if (marg == "sn") err <- rnorm(n)
    if (marg == "sl") {
      err <- rlogis(n, location = 0, scale = 0.55)
    }
    if (marg == "sg") {
      err <- rgumbel(n, loc = -0.44, scale = 0.78)
    }
  }
  N <- n * k
  ## x1 <- rbinom(N, 1, 0.5)
  x1 <- rep(rbinom(n, 1, 0.5), each = K)
  x1 <- (x1 - 0.5) / 2
  x2 <- rnorm(N)
  x3 <- rnorm(N)
  ## x2 <- rep(rnorm(n), each = K)
  ## x3 <- rep(rnorm(n), each = K)
  Y <- exp(2 + 3 * x1 + 1 * x2 + 2.5 * x3 + c(err))
  if (cr <= 0) csTime <- rep(Inf, n * k)
  if (cr == 0.25 & marg == "sn") csTime <- runif(n * k, 0, 160)
  if (cr == 0.25 & marg == "sl") csTime <- runif(n * k, 0, 160)
  if (cr == 0.25 & marg == "sg") csTime <- runif(n * k, 0, 160)
  if (cr == 0.3 & marg == "sn") csTime <- runif(n * k, 0, 110)
  if (cr == 0.3 & marg == "sl") csTime <- runif(n * k, 0, 110)
  if (cr == 0.3 & marg == "sg") csTime <- runif(n * k, 0, 110)
  if (cr == 0.5 & marg == "sn") csTime <- runif(n * k, 0, 20)
  if (cr == 0.5 & marg == "sl") csTime <- runif(n * k, 0, 20)
  if (cr == 0.5 & marg == "sg") csTime <- runif(n * k, 0, 20)
  if (cr == 0.6 & marg == "sn") csTime <- runif(n * k, 0, 8.5)
  if (cr == 0.6 & marg == "sl") csTime <- runif(n * k, 0, 8.5)
  if (cr == 0.6 & marg == "sg") csTime <- runif(n * k, 0, 8.5)
  if (cr == 0.8 & marg == "sn") csTime <- runif(n * k, 0, 1.4)
  if (cr == 0.8 & marg == "sl") csTime <- runif(n * k, 0, 1.4)
  if (cr == 0.8 & marg == "sg") csTime <- runif(n * k, 0, 1.4)
  if (cr == 0.9 & marg == "sn") csTime <- runif(n * k, 0, 0.4) #0.467
  if (cr == 0.9 & marg == "sg") csTime <- runif(n * k, 0, 0.4) #0.735
  if (cr == 0.9 & marg == "sl") csTime <- runif(n * k, 0, 0.4) #234
  if (cr == 0.95 & marg == "sn") csTime <- runif(n * k, 0, 0.15)
  if (cr == 0.95 & marg == "sl") csTime <- runif(n * k, 0, 0.15)
  if (cr == 0.95 & marg == "sg") csTime <- runif(n * k, 0, 0.15)
  dat <- data.frame(Time = pmin(Y, csTime), delta = 1 * (Y <= csTime),
                    id = rep(1:n, rep(k, n)),
                    x1 = x1, x2 = x2, x3 = x3,
                    x4 = rbinom(n, 1, .5), x5 = rnorm(n), x6 = rnorm(n),
                    x7 = rexp(n), x8 = rexp(n), x9 = runif(n), x10 = runif(n))  
  if (cr <= 0 | full) {
      dat$weights <- 1
      return(dat)
  }
  ## New groups;
  ## tmp <- unlist(lapply(split(dat$delta, dat$id), sum))
  ## groupid <- rep(2, length(tmp))
  ## groupid <- ifelse(tmp == 0, 1, groupid)
  ## groupid <- ifelse(tmp == max(tmp), 3, groupid)
  ## dat$group <- rep(groupid, each = k)
  ## groupid <- lapply(split(names(groupid), groupid), as.numeric)
  ## sizes <- c(round(.2 * length(groupid[[1]])),
  ##            round(.5 * length(groupid[[2]])),
  ##            round(1 * length(groupid[[3]])))
  groupid <- unlist(lapply(split(dat$delta, dat$id), sum)) + 1
  dat$group <- rep(groupid, each = k)
  groupid <- lapply(split(names(groupid), groupid), as.numeric)
  sizes <- c(round(.1 * length(groupid[[1]])),
             round(.2 * length(groupid[[2]])),
             round(.3 * length(groupid[[3]])),
             round(.5 * length(groupid[[4]])))
  dat$weights <- (sizes / unlist(lapply(groupid, length)))[dat$group]
  sampleid <- c(sample(groupid[[1]], sizes[1]),
                sample(groupid[[2]], sizes[2]),
                sample(groupid[[3]], sizes[3]), 
                sample(groupid[[4]], sizes[4]))
  dat2 <- subset(dat, id %in% sort(sampleid))
  dat2$id <- rep(1:sum(sizes), each = k)
  dat2$group <- rownames(dat2) <- NULL
  return(dat2)
}

set.seed(2)
dat <- simDat(n = 200, tau = .6, cr = .5, marg = "sn")
## dat <- simDat(n = 100, tau = .6, cr = .5, marg = "sn", full = TRUE)
dim(dat)
fm <- Surv(Time, delta) ~ x1 + x2 + x3
f1 <- paftgee(fm, data = dat, id = id, weights = 1 / weights,
              B = 100, control = list(init = "gehan_s"))
f2 <- paftgee(fm, data = dat, id = id, weights = 1 / weights, corstr = "ex",
              B = 100, control = list(init = "gehan_s"))
f3 <- paftgee(fm, data = dat, id = id, weights = 1 / weights, corstr = "AR1",
              B = 100, control = list(init = "gehan_s"))

f1$beta
f2$beta
f3$beta

table(f1$bstep > 100)
table(f2$bstep > 100)
table(f3$bstep > 100)

rbind(f1$beta.se, f2$beta.se, f3$beta.se)

rbind(sqrt(diag(var(f1$beta.boot[f1$bstep < 100,]))),
      sqrt(diag(var(f2$beta.boot[f2$bstep < 100,]))),
      sqrt(diag(var(f3$beta.boot[f3$bstep < 100,]))))
      
rbind(apply(f1$beta.boot[f1$bstep < 100,], 2, sd),
      apply(f2$beta.boot[f2$bstep < 100,], 2, sd),
      apply(f3$beta.boot[f3$bstep < 100,], 2, sd))

g1 <- aftgee(fm, data = dat, id = id, weights = 1 / weights, corstr = "ind", B = 100)
g2 <- aftgee(fm, data = dat, id = id, weights = 1 / weights, corstr = "ex", B = 100)
g3 <- aftgee(fm, data = dat, id = id, weights = 1 / weights, corstr = "ar1", B = 100)

coef(g1)
coef(g2)
coef(g3)

rbind(sqrt(diag(vcov(g1))), sqrt(diag(vcov(g2))), sqrt(diag(vcov(g3))))

debug(paft:::est.paftgee)

paftgee(fm, data = dat, id = id, weights = 1 / weights, corstr = "ex", 
        B = 100, control = list(init = "gehan_s"))



tmp1 <- pgeeCpp(X = x, y = Ey, weights = w, id = id, corstr = corstr, 
                lambda = lambda, pindex = pindex, b0 = b0, penalty = penalty, 
                olsb = olsb)
tmp1

tmp2 <- geese.fit(as.matrix(x), as.numeric(Ey), id, corstr = corstr)
with(tmp2, list(beta, alpha))

args(PGEE::mycor_gee2)
foo <- PGEE::mycor_gee2(N = length(unique(id)),
                        nt = as.integer(unlist(lapply(split(id, id), "length"))),
                        y = as.numeric(Ey), X = x, family = gaussian(link = "identity"),
                        beta_new = b0, corstr = "exchangeable", Mv = NULL, 
                        maxclsz = 3, R = NULL, scale.fix = FALSE, scale.value = 1)
foo$Ehat[,,1]
