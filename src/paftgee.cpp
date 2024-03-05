#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

//' @noRd
// [[Rcpp::export]]
Rcpp::List pgeeCpp(arma::mat X, arma::vec y, arma::vec weights, arma::vec id, double lambda,
                   arma::vec b0, arma::vec pindex, std::string corstr = "independence", 
                   int maxit = 500, double eps = 1e-6, double tol = 1e-3, bool SHEMonly = false,
                   std::string penalty = "SCAD", arma::vec olsb = 0) {
  Rcpp::List out(7);
  arma::vec w = weights;
  arma::vec u = unique(id);
  int N = u.n_rows;
  double ahat = 0;
  arma::vec nt = zeros<vec>(N); //count unique id
  for (int i = 0; i < N; i++) {
    for (size_t j = 0; j < id.n_elem; j++) {
      if (id(j) == u(i)) nt(i) = nt(i) + 1;
    }
  }
  int maxnt = nt.max();
  if (maxnt == 1) corstr = "independence";
  arma::mat Rhat = zeros<mat>(maxnt, maxnt);
  int p = b0.n_elem;
  arma::mat S = zeros<mat>(p, 1);
  arma::mat E = zeros<mat>(p, p);
  arma::mat H = zeros<mat>(p, p);
  arma::mat M = zeros<mat>(p, p);
  if (corstr == "independence") {
    Rhat = eye<mat>(maxnt, maxnt);
  }
  arma::vec ym = zeros<vec>(y.n_elem);
  arma::vec b1 = zeros<vec>(b0.n_elem);
  int kk = 0;
  for (int i = 0; i < maxit; i++) {
    kk = kk + 1;
    ym = y - X * b0; 
    arma::vec ymd = ym / stddev(X * b0);
    double numer = 0;
    if (corstr == "exchangeable") {
      double sumK = 0;
      arma::vec ym1 = ymd % sqrt(w);
      double phi = sum(square(ym1)) / sum(w);
      for (int j = 0; j < N; j++) {
        double sumById = 0;
        arma::vec ymById = zeros<vec>(nt(j));
        double a = sumK;
        double b = sumK + nt(j) -1;
        ymById = ym1.subvec(a, b);
        double temp = sum(ymById); 
        sumById = (temp * temp - sum(square(ymById)));
        sumK = sumK + nt(j);
        numer = numer + sumById;
      }
      ahat = numer / sum(w) / (maxnt - 1) / phi; //only for equal cluster size
      Rhat = eye<mat>(maxnt, maxnt) + ahat * (1 - eye<mat>(maxnt, maxnt));
    }
    else if (corstr == "AR1") {
      double sumK = 0;
      arma::vec ym1 = ym % sqrt(w);
      double phi = sum(square(ym1)) / sum(w);
      for (int j = 0; j < N; j++) {
        double sumById = 0;
        arma::vec ymById = zeros<vec>(nt(j));
        double a = sumK;
        double b = sumK + nt(j) -1;
        ymById = ym1.subvec(a, b);
        for (int k = 0; k < (nt(j) - 1); k++) {
          sumById = sumById + ymById(k) * ymById(k + 1);
        }
        sumK = sumK + nt(j);
        numer = numer + sumById;
      }
      ahat = numer / sum(w) / phi; // / (maxnt - 1);
      Rhat = eye<mat>(maxnt, maxnt);
      for (int j = 0; j < maxnt; j++) {
        for (int k = 0; k < maxnt; k++) {
          Rhat(j,k) = pow(ahat, abs(j - k));
        }
      }
    }
    arma::mat R = zeros<mat>(nt.n_elem, nt.n_elem);
    arma::mat NT = eye<mat>(nt.n_elem, nt.n_elem);
    R = kron(NT, pinv(Rhat));
    S = X.t() * R * (w % ym);
    arma::mat W = repmat(w, 1, X.n_cols);
    H = (X % W).t() * R * X;
    arma::mat R2 = kron(NT, ones<mat>(maxnt, maxnt)) % (ym * ym.t());
    M = (X % W).t() * R * R2 * R * X;
    arma::vec v1 = zeros<vec>(p);
    arma::vec absb0 = abs(b0);
    v1.elem(find(absb0 > lambda)).ones();
    arma::vec v2 = zeros<vec>(p);
    v2.elem(find(absb0 < (lambda * 3.7))).ones();
    arma::vec qscad = zeros<vec>(p);
    if (penalty == "SCAD") 
      qscad = lambda * (1- v1) + ((lambda * 3.7) - absb0) % v2 / (2.7) % v1;
    if (penalty == "LASSO")     
      qscad = lambda * (1 - v1);    
    if (penalty == "aLASSO") {
      arma::vec wLASSO = 1 / olsb;
      qscad = lambda * wLASSO;
    }
    E = zeros<mat>(p, p);
    arma::vec tmp = qscad / (absb0 + eps);
    for (size_t z = 0; z < pindex.n_elem; z++ ) {
      if (pindex(z) == 0) E(z, z) = 0; 
      else E(z, z) = tmp(z);
    }
    int Nf = sum(w) / maxnt;
    b1 = b0 + pinv(H + Nf * E) * (S - Nf * E * b0);
    //if (SHEMonly == TRUE) {
    out(0) = b1;
    out(1) = S;
    out(2) = H;
    out(3) = E;
    out(4) = M;
    out(5) = i;
    out(6) = ahat;
    // out.names() = Rcpp::CharacterVector::create("b1", "S", "H", "E", "M", "iter", "a");
    // return out;    }
    if (max(abs(b1 - b0)) < tol) break;
    b0 = b1; 
  } 
  // out(0) = b1;
  // out(1) = S;
  // out(2) = H;
  // out(3) = E;
  // out(4) = M;
  // out(5) = i;
  // out(6) = ahat;
  out.names() = Rcpp::CharacterVector::create("b1", "S", "H", "E", "M", "iter", "ahat");
  return out;
}
