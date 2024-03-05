library(survival)
library(Rcpp)

eRes <- function (e, delta, z = rep(1, length(e))) {
    nobs <- length(e)
    ord <- order(e)
    ei <- e[ord]
    deltai <- delta[ord]
    zi <- z[ord]
    dummy <- 1:nobs
    tmp <- survfit(Surv(ei, deltai) ~ 1, weights = zi)
    Shat <- with(tmp, approx(time, surv, ei))$y
    edif <- c(diff(ei), 0)
    ehat <- rev(cumsum(rev(edif * Shat)))
    ehat <- ehat/Shat + ei
    ehat[is.na(ehat)] <- ei[is.na(ehat)]
    eres <- ehat
    eres[dummy[ord]] <- ehat
    return(eres)
}

n <- 20
set.seed(1)

dat <- data.frame(e = rexp(n), delta = sample(0:1, n, T), z = rexp(n))

## debug(eRes)
eRes(dat$e, dat$delta, dat$z)

e <- dat$e
delta <- dat$delta
z <- dat$z

sourceCpp(code = '
#include <RcppArmadillo.h>
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
arma::vec eResC(arma::vec Time, arma::vec censor, arma::vec wgt) {
  arma::vec T0 = arma::sort(arma::unique(Time));
  int n = T0.n_elem;
  arma::vec d(n, arma::fill::zeros);
  arma::vec r(n, arma::fill::zeros); 
  for (int i = 0; i < n; i++) {
    arma::uvec ind1 = find(Time == T0[i]);
    d[i] = sum(censor.elem(ind1) % wgt.elem(ind1));
    r(span(0, i)) += sum(wgt.elem(ind1));
  }
  arma::vec surv = cumprod(1 - d / r);
  T0.resize(n + 1);
  T0(n) = T0(n - 1);
  arma::vec tmp = diff(T0) % surv;
  tmp.replace(datum::nan, 0);
  arma::vec eres(Time.n_elem, arma::fill::zeros);
  arma::uvec ind = sort_index(Time);
  eres(ind) = (reverse(cumsum(reverse(tmp)))) / surv + T0(span(0, n - 1));
  eres.replace(datum::nan, max(T0));
  return eres;
}
')

drop(eResC(dat$e, dat$delta, dat$z))
