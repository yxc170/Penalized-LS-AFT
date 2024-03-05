#include <RcppArmadillo.h>
#include <algorithm>
#include <set>
#include <utility>

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

arma::mat matvec(arma::mat x, arma::vec y) {
  arma::mat out(x.n_rows, x.n_cols);
  for (size_t i = 0; i < x.n_cols; i++) {
    out.col(i) = x.col(i) % y;
  }
  return out;
}

double wmean(arma::vec a, arma::vec b) {
  return sum(a % b) / sum(b);
}

// ############################################################
// functions for gee
// ############################################################

// Functions to calculate parameter in variance-covariance matrix
// Those functions assume equal cluster size and weights at cluster level
double ahatEx(arma::vec a,
							arma::vec nt,
							arma::vec index,
							arma::vec w) {
  vec aa = a % sqrt(w);
  double phi = sum(square(aa)) / sum(w); // this makes scale.fix = FALSE
  int n = nt.n_elem;
  double out = 0;
  for (int i = 0; i < n; i++) {
    arma::vec a2 = a(span(index(i), index(i) + nt(i) - 1));
    out += (sum(a2) * sum(a2) - sum(a2 % a2)) * w(index(i));
  }
  return(out / (nt(0) - 1) / sum(w) / phi);
}

double ahatAR1(arma::vec a,
							 arma::vec nt,
							 arma::vec index,
							 arma::vec w) {
  int n = nt.n_elem;
  vec aa = a % sqrt(w);
  double phi = sum(square(aa)) / sum(w);
  double out = 0;
  for (int i = 0; i < n; i++) {
    arma::vec a2 = aa(span(index(i), index(i) + nt(i) - 1));
    arma::vec w2 = w(span(index(i), index(i) + nt(i) - 1));
    out += sum(a2(span(0, nt(0) - 2)) % a2(span(1, nt(0) - 1))) / sum(w2);
  }
  return(out / (nt(0) - 1) / phi);
}

// Penalty functions
arma::vec qscad(arma::vec b, double lambda) {
  double a = 3.7;
  arma::vec dif = a * lambda - b;
  dif.elem(find(dif < 0)).zeros();
  arma::vec out(b.n_elem, arma::fill::value(lambda));
  // out.fill(lambda);
  arma::uvec ind = find(b > lambda);
  out(ind) = dif(ind) / (a - 1);
  return(out);
}

arma::vec qlasso(arma::vec b, double lambda) {
  arma::vec out(b.n_elem, arma::fill::value(lambda));
  out.elem(find(b > lambda)).zeros();
  return(out);
}

// [[Rcpp::export(rng = false)]]
Rcpp::List pgee(arma::vec y,
								arma::mat X,
								arma::vec b0,
								arma::vec nt,
								arma::vec pindex, // 0 means to penalize
								arma::vec w,
								std::string penalty,
								std::string corstr,
								double lambda,
								double eps,
								double tol,
								int maxit){
  Rcpp::List out(7);
  int N = nt.n_elem;
  int nx = X.n_cols;
  int k = nt(0);
  arma::vec index(N, arma::fill::zeros);
  index(span(1, N - 1)) = cumsum(nt(span(0, N - 2)));
  arma::vec b1 = b0;
  for (int j = 1; j <= maxit; j++) {
    arma::vec E1(nx, arma::fill::zeros);
    if (penalty == "scad")
      E1 = qscad(abs(b0), lambda) / (abs(b0) + eps);
    if (penalty == "lasso")
      E1 = qlasso(abs(b0), lambda) / (abs(b0) + eps);
    E1.elem(find(pindex == 0)).zeros();
    arma::mat E = diagmat(E1);
    arma::mat S(nx, 1, arma::fill::zeros);
    arma::mat H(nx, nx, arma::fill::zeros);
    arma::mat M(nx, nx, arma::fill::zeros);
    arma::vec mu = X * b0;
    arma::vec ym = y - mu;
    arma::mat Rhat(k, k, arma::fill::eye);
    double ahat = 0;
    if (corstr == "ex" && k > 1) {
      ahat = ahatEx(ym / stddev(mu), nt, index, w);
      Rhat = Rhat * (1 - ahat) + ahat;
    }
    if (corstr == "ar1" && k > 1) {
      ahat = ahatAR1(ym / stddev(mu), nt, index, w);
      arma::mat Rhat2(k, k, arma::fill::zeros);
      arma::vec tmp(k - 1, arma::fill::value(ahat));
      tmp = cumprod(tmp);
      for (int i = 0; i < k - 1; i++) {
				Rhat2.submat(i + 1, i, k - 1, i) = tmp(span(0, k - i - 2));
      }
      Rhat = Rhat + Rhat2 + Rhat2.t();
    }
    // Might replace this to matrix operations
    // but I doubt matrix version is faster because it requires more memory
    arma::mat bigD = matvec(X, sqrt(w));
    arma::vec ymw = ym % sqrt(w);
    for (int i = 0; i < N; i++) {
      arma::vec ym2 = ym(span(index(i), index(i) + nt(i) - 1));
      arma::vec ym2w = ymw(span(index(i), index(i) + nt(i) - 1));
      arma::mat bigD2 = bigD.rows(span(index(i), index(i) + nt(i) - 1));
      arma::mat tmp = bigD2.t() * pinv(Rhat);
      S += tmp * ym2w;
      H += tmp * bigD2;
      tmp *= ym2;
      M += tmp * tmp.t();
    }
    // M doesn't match, why?
    int Nf = sum(w) / nt(0);
    b1 = b0 + pinv(H + Nf * E) * (S - Nf * E * b0);
    out(0) = b1;
    out(1) = S;
    out(2) = H;
    out(3) = E;
    out(4) = M;
    out(5) = j;
    out(6) = ahat;
    if(max(abs(b1 - b0)) < tol) break;
    b0 = b1;
  }
  out.names() = Rcpp::CharacterVector::create("b", "S", "H", "E", "M", "iter", "alpha");
  return out;
}

// no penality version
// [[Rcpp::export(rng = false)]]
Rcpp::List gee(arma::vec y,
							 arma::mat X,
							 arma::vec b0,
							 arma::vec nt,
							 arma::vec w,
							 std::string corstr,
							 double tol,
							 int maxit){
  Rcpp::List out(7);
  int N = nt.n_elem;
  int nx = X.n_cols;
  int k = nt(0);
  arma::vec index(N, arma::fill::zeros);
  index(span(1, N - 1)) = cumsum(nt(span(0, N - 2)));
  arma::vec b1 = b0;
  for (int j = 1; j <= maxit; j++) {
    arma::mat S(nx, 1, arma::fill::zeros);
    arma::mat H(nx, nx, arma::fill::zeros);
    arma::mat M(nx, nx, arma::fill::zeros);
    arma::vec mu = X * b0;
    arma::vec ym = y - mu;
    arma::mat Rhat(k, k, arma::fill::eye);
    double ahat = 0;
    if (corstr == "ex" && k > 1) {
      ahat = ahatEx(ym / stddev(mu), nt, index, w);
      Rhat = Rhat * (1 - ahat) + ahat;
    }
    if (corstr == "ar1" && k > 1) {
      ahat = ahatAR1(ym / stddev(mu), nt, index, w);
      arma::vec tmp(k - 1, arma::fill::value(ahat));
      arma::mat Rhat2(k, k, arma::fill::zeros);
      tmp = cumprod(tmp);
      for (int i = 0; i < k - 1; i++) {
				Rhat2.submat(i + 1, i, k - 1, i) = tmp(span(0, k - i - 2));
      }
      Rhat = Rhat + Rhat2 + Rhat2.t();
    }
    arma::mat bigD = matvec(X, sqrt(w));
    arma::vec ymw = ym % sqrt(w);
    for (int i = 0; i < N; i++) {
      arma::vec ym2 = ym(span(index(i), index(i) + nt(i) - 1));
      arma::vec ym2w = ymw(span(index(i), index(i) + nt(i) - 1));
      arma::mat bigD2 = bigD.rows(span(index(i), index(i) + nt(i) - 1));
      arma::mat tmp = bigD2.t() * pinv(Rhat);
      S += tmp * ym2w;
      H += tmp * bigD2;
      tmp *= ym2;
      M += tmp * tmp.t();
    }
    b1 = b0 + pinv(H) * S;
    out(0) = b1;
    out(1) = S;
    out(2) = H;
    out(4) = M;
    out(5) = j;
    out(6) = ahat;
    if(max(abs(b1 - b0)) < tol) break;
    b0 = b1;
  }
  arma::mat E(nx, nx, arma::fill::zeros);
  out(3) = E;
  out.names() = Rcpp::CharacterVector::create("b", "S", "H", "E", "M", "iter", "alpha");
  return out;
}

//Functions for cross-validation according to proportion of groups
Rcpp::List getCVprop(arma::vec index, arma::vec nt, arma::vec w, int nCV) {
  Rcpp::List out(nCV);
  int N = index.n_elem;
  arma::vec uw = unique(w);
  arma::vec grIndex;
  arma::vec cvIndex1;
  arma::vec cv = arma::linspace(1, nCV, nCV);
  arma::vec w1(N, arma::fill::zeros);
  double offset = 0;
  for (int i = 0; i < N; i++) {
    w1(i) = w(offset);
    offset += nt[i];
  }
  for (int i = 0; i < uw.n_elem; i++) {
    grIndex = arma::join_cols(grIndex, index(find(w1 == uw(i))));
  }
  for (int i = 0; i < N / (nCV - 2); i++) {
    cvIndex1 = arma::join_cols(cvIndex1, cv);
  }
  arma::vec cvIndex = cvIndex1(arma::span(0,N-1));
  for (int i = 0; i < nCV; i++) {
    arma::vec g = grIndex(find(cvIndex == (i+1)));
    out(i) = g;
  }
  return(out);
}




// Functions need for cross-validation
Rcpp::List getCV(arma::vec index, int nCV) {
  Rcpp::List out(nCV);
  int N = index.n_elem;
  arma::vec cvIndex = shuffle(index);
  arma::vec nTrain(nCV, arma::fill::value(round(N / nCV)));
  if (sum(nTrain) < N) {
    for (int i = 0; i < N - sum(nTrain); i++) nTrain(i)++;
  }
  if (sum(nTrain) > N) {
    for (int i = 0; i < sum(nTrain) - N; i++) nTrain(i)--;
  }
  double offset = 0;
  for (int i = 0; i < (nCV-1); i++) {
    out(i) = sort(cvIndex(span(0 + offset, nTrain(i) - 1 + offset)));
    offset += nTrain(i);
  }
  out(nCV-1)= sort(cvIndex(span(offset, N - 1)));
  return(out);
}

arma::vec getID(arma::vec nt) {
  arma::vec out(sum(nt));
  int offset = 0;
  int n = nt.n_elem;
  for (int i = 0; i < n; i++) {
    out(span(offset, nt(i) - 1 + offset)).fill(i);
    offset += nt(i);
  }
  return(out);
}

// ############################################################
// functions for aftgee
// ############################################################

// [[Rcpp::export(rng = false)]]
arma::vec eResC(arma::vec const Time,
								arma::vec const censor,
								arma::vec const wgt) {
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


// aftgee algorithm in c, no penalty
// [[Rcpp::export(rng = false)]]
Rcpp::List aftgeeEst(arma::vec y,
										 arma::mat X,
										 arma::vec D,
										 arma::vec b0,
										 arma::vec nt,
										 arma::vec w,
										 std::string corstr,
										 double tol,
										 int maxit) {
  Rcpp::List out(3);
  arma::vec iter_gee(maxit, arma::fill::zeros);
  arma::vec b1 = b0;
  for (int j = 1; j <= maxit; j++) {
    arma::vec e = y - X * b0;
    arma::vec eres = eResC(e, D, w);
    arma::vec Ey = D % y + (1 - D) % (eres + X * b0);
    Rcpp::List fit = gee(Ey, X, b0, nt, w, corstr, tol, maxit);
    arma::vec b1 = fit(0);
    if(max(abs(b1 - b0)) < tol) break;
    b0 = b1;
    out(0) = b1;
    out(1) = j;
    iter_gee(j - 1) = fit(5);
  }
  out(2) = nonzeros(iter_gee);
  out.names() = Rcpp::CharacterVector::create("b", "iter", "iter_gee");
  return out;
}

// aftgee algorithm in c, with penalty, lambda is a scalar
// [[Rcpp::export(rng = false)]]
Rcpp::List paftgeeEst(arma::vec y,
											arma::mat X,
											arma::vec D,
											arma::vec b0,
											arma::vec nt,
											arma::vec pindex,
											arma::vec w,
											std::string penalty,
											std::string corstr,
											double lambda,
											double eps,
											double tol,
											int maxit) {
  Rcpp::List out(5);
  arma::vec iter_gee(maxit, arma::fill::zeros);
  arma::vec b1 = b0;
  for (int j = 1; j <= maxit; j++) {
    arma::vec e = y - X * b0;
    arma::vec eres = eResC(e, D, w);
    arma::vec Ey = D % y + (1 - D) % (eres + X * b0);
    Rcpp::List fit = pgee(Ey, X, b0, nt, pindex, w, penalty, corstr, lambda, eps, tol, maxit);
    arma::vec b1 = fit(0);
    out(0) = b1;
    out(1) = j;
		out(3) = fit(2);
		out(4) = fit(3);
	  iter_gee(j - 1) = fit(5);
    if(max(abs(b1 - b0)) < tol) break;
    b0 = b1;
  }
  out(2) = nonzeros(iter_gee);
	out.names() = Rcpp::CharacterVector::create("b", "iter", "iter_gee", "H", "E");
  return out;
}


// rule = cv
// y0, X0, D0 are from test data
// y1, X1, D1 are from training data
// [[Rcpp::export]]
Rcpp::List paftgeeCV(arma::vec y,
										 arma::mat X,
										 arma::vec D,
										 arma::vec b0,
										 arma::vec nt,
										 arma::vec pindex, // 0 means to penalize
										 arma::vec w,
										 std::string penalty,
										 int nCV,
										 arma::vec lambda,
										 bool CVprop,
										 double eps,
										 double tol,
										 int maxit){
  Rcpp::List out(3);
	int nlambda = lambda.n_elem;
  arma::mat cvRaw(y.n_elem, nlambda, arma::fill::zeros);
  arma::vec id = getID(nt);
  arma::vec uid = unique(id);
  Rcpp::List cvList = getCV(uid, nCV);
  if (CVprop == TRUE) {
    cvList = getCVprop(uid, nt, w, nCV);
    }
  arma::vec wCV;
	for (int i = 0; i < nCV; i++) {
    arma::uvec id0;
    arma::vec id1(id.n_elem, arma::fill::zeros);
    arma::vec idCV = cvList(i);
    //Rcpp::Rcout << "idCV" << size(idCV)<<endl;
    for (int k = 0; k < (int) idCV.n_elem; k++) {
      id0 = arma::join_cols(id0, find(id == idCV(k)));
    }
    id1.elem(id0).ones();
    arma::vec y0 = y(id0); // test data
    arma::mat X0 = X.rows(id0);
    arma::vec D0 = D(id0);
    arma::vec w0 = w(id0);
    wCV = arma::join_cols(wCV, w0);
    arma::vec y1 = y(find(id1 == 0)); // training data
    arma::mat X1 = X.rows(find(id1 == 0));
    arma::vec D1 = D(find(id1 == 0));
    arma::vec w1 = w(find(id1 == 0));
    // assume equal cluster size
    arma::vec nt0 = nt(span(1, idCV.n_elem));
    arma::vec nt1 = nt(span(1, (int) uid.n_elem - idCV.n_elem));
    for (arma::uword j = 0; j < (uword) nlambda; j++) {
      Rcpp::List trainFit = paftgeeEst(y1, X1, D1, b0, nt1, pindex, w1,
																			 penalty, "ind", lambda(j), eps, tol, maxit);
      arma::vec trainBeta = trainFit(0);
      Rcpp::List testFit = paftgeeEst(y0, X0, D0, trainBeta, nt0, pindex, w0,
																			penalty, "ind", 0, eps, tol, maxit);
			arma::vec testBeta = testFit(0);
			arma::vec e0 = y0 - X0 * testBeta;
			arma::vec eres = eResC(e0, D0, w0);
			arma::vec Ey0 = D0 % y0 + (1 - D0) % (eres + X0 * testBeta);
			arma::vec error0 = Ey0 - X0 * trainBeta;
			cvRaw(id0, arma::uvec {j}) = error0 % error0;
    }
  }
  arma::vec cvm(nlambda, arma::fill::zeros);
  arma::vec cvse(nlambda, arma::fill::zeros);
  for (int i = 0; i < nlambda; i++) {
    cvm(i) = wmean(cvRaw.col(i), wCV);
    cvse(i) = sqrt(wmean(pow(cvRaw.col(i) - cvm(i), 2), wCV) / (y.n_elem - 1));
  }
  out(0) = cvm;
  out(1) = cvse;
  out(2) = cvRaw;
  out.names() = Rcpp::CharacterVector::create("cvm", "cvse", "cvRaw");
  return out;
}

// rule = bic
// [[Rcpp::export(rng = false)]]
arma::vec paftgeeBIC(arma::vec y,
										 arma::mat X,
										 arma::vec D,
										 arma::vec b0,
										 arma::vec nt,
										 arma::vec pindex, // 0 means to penalize
										 arma::vec w,
										 std::string penalty,
										 arma::vec lambda,
										 double eps,
										 double tol,
										 int maxit){
	arma::vec bicValues(lambda.n_elem, arma::fill::zeros);
	for (arma::uword j = 0; j < (uword) lambda.n_elem; j++) {
		Rcpp::List fit = paftgeeEst(y, X, D, b0, nt, pindex, w, penalty, "ind", lambda(j), eps, tol, maxit);
		arma::vec b1 = fit(0);
		arma::mat H = fit(3);
		arma::mat E = fit(4);
		double Nm = sum(w);
		double N = sum(w); // why define this twice?
		arma::mat nV = pinv(H + N * E);
		arma::vec error = y - X * b1;
		arma::mat tmp = nV * H;
		bicValues(j) = log(sum(error % error % w) / Nm) + log(Nm) * sum(tmp.diag()) / Nm;
	}
	return bicValues;
}


// aftgee algorithm in c, with penalty,
// lambda is a vector (best lambda is chosen)
// [[Rcpp::export(rng = false)]]
Rcpp::List paftgeeEst1(arma::vec y,
											 arma::mat X,
											 arma::vec D,
											 arma::vec b0,
											 arma::vec nt,
											 arma::vec pindex,
											 arma::vec w,
											 std::string penalty,
											 int nCV,
											 std::string corstr,
											 arma::vec lambda,
											 bool CVprop,
											 std::string rule,
											 double eps,
											 double tol,
											 int maxit) {
  Rcpp::List out(7);
  arma::vec iter_gee(maxit, arma::fill::zeros);
  arma::vec b1 = b0;
	double lambdaBest = 0;
	if (lambda.n_elem == 1) lambdaBest = lambda(0);
	else {
		if (rule == "cv" | rule == "1se") {
			Rcpp::List cvList = paftgeeCV(y, X, D, b0, nt, pindex, w, penalty, nCV, lambda, CVprop, eps, tol, maxit);
			arma::vec cvm = cvList(0);
			arma::vec cvse = cvList(1);
			arma::uword i = cvm.index_min();
			if (rule == "cv") {
				lambdaBest = lambda(i);
			} else {
				double oneSE = cvm(i) + cvse(i);
				arma::vec lambda2 = lambda(span(i, lambda.n_elem - 1));
				arma::vec cvm2 = cvm(span(i, cvm.n_elem - 1));
				if (oneSE > max(cvm2)) lambdaBest = max(lambda2);
				else {
					arma::uvec j = find(cvm2 < oneSE);
					lambdaBest = lambda2(j.n_elem - 1);
				}
			}
			out(5) = cvList;
		}
		if (rule == "bic") {
			arma::vec bics = paftgeeBIC(y, X, D, b0, nt, pindex, w, penalty, lambda, eps, tol, maxit);
			arma::uword i = bics.index_min();
			lambdaBest = lambda(i);
			out(5) = bics;
		}
	}
	out(6) = lambdaBest;
  for (int j = 1; j <= maxit; j++) {
    arma::vec e = y - X * b0;
    arma::vec eres = eResC(e, D, w);
    arma::vec Ey = D % y + (1 - D) % (eres + X * b0);
    Rcpp::List fit = pgee(Ey, X, b0, nt, pindex, w, penalty, corstr, lambdaBest, eps, tol, maxit);
    arma::vec b1 = fit(0);
    out(0) = b1;
    out(1) = j;
		out(3) = fit(2);
		out(4) = fit(3);
	  iter_gee(j - 1) = fit(5);
	  if(max(abs(b1 - b0)) < tol) break;
	  b0 = b1;
  }
  out(2) = nonzeros(iter_gee);
	out.names() =
		Rcpp::CharacterVector::create("b", "iter", "iter_gee", "H", "E", "values", "lambdaBest");
  return out;
}
