#' Function to estimate aftgee 
#'
#' @noRd
#' 
est.paftgee <- function(y, x, dat, w, d, id, b0, maxit = 500, tol = 1e-3, fold = 0, pindex = NULL,
                        corstr = c("independence", "exchangeable", "AR1"), rule = c("cv", "1se", "bic"), lam0 = 0,
                        penalty = c("SCAD", "LASSO", "aLASSO")) {
  rule <- match.arg(rule)
  penalty <- match.arg(penalty)
  corstr <- match.arg(corstr)
  nb0 <- length(b0)
  lam <- lam0
  olsb <- b0
  steps <- 0
  if (is.null(pindex)) {
    if (fold > 0) pindex <- c(0, rep(1, ncol(x) - 1))
    else pindex <- rep(0, ncol(x))
  }
  for (i in 1:maxit) {
    steps <- steps + 1
    e <- y - colSums(t(x) * b0)
    eres <- aftgee:::eRes(e, delta = d, z = w)
    Ey <- d * y + (1 - d) * (eres[[1]] + colSums(t(x) * b0))
    if (lam0 != 0 || fold ==0) { #
      foo <- pgeeCpp(X = x, y = Ey, weights = w, id = id,
                     corstr = corstr, lambda = lam, pindex = pindex, b0 = b0,
                     penalty = penalty, olsb = olsb)
    } 
    else if (rule == "cv") {
      if (i == 1 & fold > 0) {
        suppressWarnings(cv <- cv.paft(X = x, y = Ey, dat = dat, id = id, nb0 = nb0,
                                       fold = fold, pindex = pindex, rule = rule,
                                       lambda.vec = 1:30/20, weights = w, b0 = b0,
                                       penalty = penalty, olsb = olsb))
        lam <- cv$lam.opt
        lam.1se <- cv$lam.1se
        
      }
      foo <- pgeeCpp(X = x, y = Ey, weights = w, id = id,
                     corstr = corstr, lambda = lam, pindex = pindex, b0 = b0,
                     penalty = penalty, olsb = olsb)
    }
    else if (rule == "bic") {
      if (i == 1 & fold > 0) {
        suppressWarnings(cv <- cv.paft(X = x, y = Ey, dat = dat, id = id, nb0 = nb0,
                                       fold = fold, pindex = pindex, rule = rule,
                                       lambda.vec = 1:30/20, weights = w, b0 = b0,
                                       penalty = penalty, olsb = olsb))
        lam <- cv$lam.bic
      }
      foo <- pgeeCpp(X = x, y = Ey, weights = w, id = id,
                     corstr = corstr, lambda = lam, pindex = pindex, b0 = b0,
                     penalty = penalty, olsb = olsb)
    }
    
    b1 <- as.numeric(foo$b1)
    ahat <- foo$ahat
    ## if (i %% (maxit / 100) == 0) print(paste0("Step:", i, " diff is ", max(abs(b1 - b0))))
    if(max(abs(b1 -  b0)) < tol) break
    b0 <- b1
    ## print(b1)
  }
  if (fold == 0) 
    return(list(b1 = b1, lam = lam0, steps = steps, ahat = ahat))
  else if (lam0 != 0) 
    return(list(b1 = b1, lam = lam, steps = steps, ahat = ahat))
  else if (fold > 0 & rule == "cv") {
    return(list(b1 = b1, lam = lam, lam.1se = lam.1se, steps = steps, ahat = ahat))
  }
  else if ((fold > 0 & rule == "bic") || (fold > 0 & rule == "1se"))
    return(list(b1 = b1, lam = lam, steps = steps, ahat = ahat))
  
}
