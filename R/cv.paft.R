#' Function to do cross validation - needs to be reworked
#'
#' @importFrom MASS ginv
#' @noRd
#' 
cv.paft <- function(X, y, dat, id, nb0, scale.fix, scale.value, fold, b0, rule = "1se",
                    lambda.vec, pindex, penalty, olsb, eps = 1e-6, maxiter = 100, tol = 1e-3, weights) {
  N <- length(unique(id))
  nt <- tabulate(id)
  cv.min <- Inf
  bic.min <- Inf
  cv.vec <- NULL
  se.vec <- NULL
  bic.vec <- NULL
  d <- dat$delta
  ## univariate case: select train/test data according to proportion of uncensored/censored data
  index.gr <- lapply(0:1, function(x) dat$id[dat$delta==x])
  set.seed(1)
  cv.index <- lapply(1:2, function(x) sample(rep(1:fold, length = length(index.gr[[x]]))))
  df.temp <- data.frame(groupid=unlist(index.gr),foldindex=unlist(cv.index))
  cv.list <- lapply(1:fold, function(x) df.temp$groupid[which(df.temp$foldindex == x)])
  ## multivariate case: select train/test data according to proportion of groups
  ## groups <- sort(unique(dat$group))
  ## cv.prop <- sapply(groups, function(x) round(table(dat$group)[[x]]/nt[1]/fold) )##only for equal cluster size
  ## index.gr <- lapply(groups, function(x) unique(dat$id[dat$group==x]))
  ## set.seed(48976)
  ## cv.index <- lapply(groups, function(x) sample(rep(1:fold, length = length(index.gr[[x]]))))
  ## df.temp <- data.frame(groupid=unlist(index.gr),foldindex=unlist(cv.index))
  ## cv.list <- lapply(1:fold, function(x) df.temp$groupid[which(df.temp$foldindex == x)])
  for (j in 1:length(lambda.vec))   {
    lam.temp <- lambda.vec[j]
    cv.raw <- cv.w <- NULL
    ## -------------BIC----------------------
    if (rule == "bic") {
      b0.bic <- b0
      b1.bic <- rep(0, nb0)
      H <- NULL
      for (i in 1:100) {
        e <- y - colSums(t(X) * b0.bic)
        eres <- eRes(e, delta = d, z = weights)
        Ey <- d * y + (1 - d) * (eres[[1]] + colSums(t(X) * b0.bic))
        fit.bic <- pgeeCpp(X = X, y = Ey, weights = weights, id = id,
                           corstr = "independence", lambda = lam.temp, pindex = pindex,
                           b0 = b0.bic, penalty = penalty, olsb = olsb)
        b1.bic <- as.numeric(fit.bic$b1)
        ## if (i %% (maxit / 100) == 0) print(paste0("Step:", i, " diff is ", max(abs(b1 - b0))))
        if(max(abs(b1.bic -  b0.bic)) < tol) break
        b0.bic <- b1.bic
        H <-  fit.bic$H
        ## print(b1)
      }
      Nm <- sum(weights)
      N <- sum(weights)
      nV <- ginv(fit.bic$H + N * fit.bic$E)
      error <- y - X %*% b1.bic
      bic.value <-
        log(sum(error * error * weights) / Nm) +
        log(Nm) * sum(diag(nV %*% H)) / Nm
      bic.vec <- c(bic.vec, bic.value)
      if (bic.value < bic.min) {
        lam.bic <- lam.temp
        bic.min <- bic.value
      }
    }
    ## -------------CV & 1se --------------------------
    if ((rule == "1se") || (rule == "cv")) {
      cv.index <- sample(rep(1:fold, length = N))
      cv.list <- lapply(1:fold, function(x) which(cv.index == x))
      for (k in 1:fold) { #k=1
        index.cv <- which(id %in% cv.list[[k]])
        obsy.train <- log(dat$Time[-index.cv])
        y.train <- y[-index.cv] 
        X.train <- X[-index.cv, ]
        id.train <- id[-index.cv]
        d.train <- dat$delta[-index.cv]
        weights.train <- weights[-index.cv]
        b1.train <- rep(0, nb0)
        b.initial <- b0
        for (i in 1:100) {
          e <- obsy.train - X.train %*% b.initial
          eres <- eRes(e, delta = d.train, z = weights.train)
          Ey.train <- d.train * obsy.train + (1 - d.train) * (eres[[1]] + X.train %*% b.initial)
          fit.train <- pgeeCpp(X = X.train, y = Ey.train, weights = weights.train, id = id.train,
                               corstr = "independence", lambda = lam.temp, pindex = pindex, 
                               b0 = b.initial, penalty = penalty, olsb = olsb)
          b1.train <- as.numeric(fit.train$b1)
          if(max(abs(b1.train - b.initial)) < tol) break
          b.initial <- b1.train
        }
        beta.train <- b1.train
        ## ----------new PE---------
        X.test <- X[index.cv, ]
        X.test <- as.matrix(X.test)
        y.test <- log(dat$Time[index.cv])
        d.test <- dat$delta[index.cv]
        weights.test <- weights[index.cv]
        id.test <- id[index.cv]
        b1.test <- rep(0, nb0)
        b.initial <- b0
        for (i in 1:100) {
          e.test <- y.test - X.test %*%b.initial
          eres <- eRes(e.test, delta = d.test, z = weights.test)
          Ey.test <- d.test * y.test + (1 - d.test) * (eres[[1]] + X.test %*% b.initial)
          fit.test <- pgeeCpp(X = X.test, y = Ey.test, weights = weights.test, id = id.test,
                              corstr = "independence", lambda = 0, pindex = pindex,
                              b0 = b.initial, penalty = penalty, olsb = olsb) 
          b1.test <- as.numeric(fit.test$b1)
          if(max(abs(b1.test - b.initial)) < tol) break
          b.initial <- b1.test
        }
        e.test <- y.test - X.test %*% b1.test
        eres <- eRes(e.test, delta = d.test, z = weights.test)
        Ey.test <- d.test * y.test + (1 - d.test) * (eres[[1]] + X.test %*% b1.test)
        error.test <- Ey.test - X.test %*% beta.train #t(X.test)*beta.train
        cv.raw <- c(cv.raw, (error.test) ^2 )
        cv.w <- c(cv.w, weights.test)
      } #k 
      cvm <- weighted.mean(cv.raw, cv.w)
      cvse <- sqrt(weighted.mean((cv.raw - cvm)^2, cv.w ) / (length(y) -1))
      se.vec <- c(se.vec, cvse)
      cv.vec <- c(cv.vec, cvm)
      if (cvm < cv.min) {
        lam.min <- lam.temp
        cv.min <- cvm
        upb <- cv.min + cvse
      }
    }
  } #j
  if (rule == "bic") {
    return(list(lam.bic = lam.bic, bic.min = bic.min))
  }
  df <- data.frame(lambda.vec, pe = cv.vec, se = se.vec)
  df1 <- df[df$lambda.vec>=lam.min,]
  ## lam.1se <- df1$lambda.vec[df1$pe == max(df1$pe[df1$pe < upb])]
  lam.1se <- max(df1$lambda.vec[df1$pe < upb])
  ## --------Plot--------------
  ## p1 <- ggplot(df, aes(x = lambda.vec, y = pe)) + geom_line() +
  ##   geom_errorbar(aes(ymin = pe-se, ymax = pe+se), width = 0.004) +
  ##   geom_point() +
  ##   labs(x = expression(lambda), y = "PE",
  ##        caption = paste("lam.min = ", lam.min, "; lam.1se = ", lam.1se))
  ## p11 <- p1 + geom_vline(xintercept = lam.min, linetype = "dotted") +
  ##   geom_vline(xintercept = lam.1se, linetype = "twodash")
  ## print(p11)
  ## --------------------------
  return(list(lam.opt = lam.min, lam.1se = lam.1se, cv.min = cv.min, cv.vec = cv.vec, se.vec = se.vec))
}
