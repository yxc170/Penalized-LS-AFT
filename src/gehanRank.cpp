#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// codes from aftgee package
using cmp_par = std::pair<double, arma::uword>;

arma::mat matvec2(arma::mat x, arma::vec y) {
  arma::mat out(x.n_rows, x.n_cols);
  for (size_t i = 0; i < x.n_rows; i++) {
    out.row(i) = x.row(i) % y.t();
  }
  return out;
}

//' @noRd
// [[Rcpp::export(rng = false)]]
arma::vec gehan_ns(const arma::vec& b,
		   const arma::mat& X,
		   const arma::vec& Y,
		   const arma::vec& d,
		   const arma::vec& W) {
  arma::uword const n = Y.n_elem;
  arma::uword const p = b.n_elem;
  arma::vec out(p, arma::fill::zeros);
  if(n < 1) return out;
  arma::vec yexa = Y - X.t() * b;
  yexa.replace(-arma::datum::inf, arma::datum::nan);
  yexa.replace(arma::datum::nan, yexa.min() - 0.01);
  arma::uvec const idx = arma::sort_index(yexa);
  auto cmp = [](cmp_par const &x, cmp_par const &y){
    return x.first <= y.first;
  };
  std::set<cmp_par, decltype(cmp)> indices(cmp);
  double w_sum{};
  arma::vec x_col_sum(p, arma::fill::zeros);
  {
    auto const idx_i = idx[0];
    indices.emplace(yexa[idx_i], idx_i);
    x_col_sum = sum(matvec2(X, W), 1);
    w_sum = sum(W);
    out = d(idx_i) * W(idx_i) * (w_sum * X.col(idx_i) - x_col_sum);
  }
  auto indices_head = indices.begin();
  for(arma::uword i = 1; i < n; ++i) {
    auto const idx_i = idx[i];
    indices.emplace(yexa[idx_i], idx_i);
    if(yexa[idx_i] > indices_head->first) {
      while(yexa[idx_i] > indices_head->first) {
	x_col_sum -= W(indices_head->second) * X.col(indices_head->second);
	w_sum -= W(indices_head->second);
	++indices_head;
      }
    }
    else --indices_head;
    if (w_sum > 0)
      out += d(idx_i) * W(idx_i) * (w_sum * X.col(idx_i) - x_col_sum);
  }
  return out / n / sqrt(n);
}

//' @noRd
// [[Rcpp::export(rng = false)]]
arma::rowvec gehan_s(arma::vec b,
		     arma::mat X,
		     arma::vec y,
		     arma::vec d,
		     arma::vec w) {
  int n = y.size();
  int p = b.size();
  vec xb = X * b;
  double ei = 0;
  double ej = 0;
  double H = 0;
  rowvec out(p, fill::zeros);
  for (int i = 0; i < n; i++) {
    if (d(i) > 0) {
      ei = y(i) - xb(i);
      for (int j = 0; j < n; j++) {
        ej = y(j) - xb(j);
        rowvec xdif = (X.row(i) - X.row(j));
        vec rij = xdif * xdif.t();
        rij(0) = sqrt(rij(0));
        if (rij(0) != 0) {
          H = normcdf(sqrt(n) * (ej - ei) / rij(0));
          out = out + w(i) * w(j) * xdif * H;
        }
      }
    }
  }
  return(out / n / sqrt(n));
}
