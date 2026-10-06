// ecl.cpp - survival-weighted expected credit loss, streamed row by row
//
// For each exposure i and month t = 1..H:
//   ECL_life(i) = sum_t S(t-1) h_t LGD_t EAD_t DF_t,   S(t) = prod_{s<=t} (1 - min(h_s + p_s, 1))
// with the 12-month figure the partial sum at t = hz. Scenario shocks (one
// factor z with correlation rho, multiplier of the hazards, add-on to the
// LGD, multiplier of the EAD) are applied on the fly, so no n x H matrix is
// ever materialised: memory is O(n) whatever the term, where the matrix
// formulation needs several n x H copies (2.9 GB each at n = 1e6, H = 360).
//
// Every input is a scalar, a vector of length n (flat over the months) or a
// column-major n x H matrix.

#include <Rcpp.h>
#include <Rmath.h>
#include <cmath>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

struct Src {
  const double* v = nullptr;
  int mode = -1;   // -1 absent, 0 scalar, 1 vector (n), 2 matrix (n x H)
  std::size_t n = 0;
  inline double at(std::size_t i, std::size_t t) const {
    switch (mode) {
    case 0: return v[0];
    case 1: return v[i];
    default: return v[i + t * n];
    }
  }
};

Src make_src(SEXP x, std::size_t n, std::size_t H, const char* what) {
  Src s;
  s.n = n;
  if (Rf_isNull(x)) return s;
  if (TYPEOF(x) != REALSXP) Rcpp::stop("cpp_ecl_paths(): `%s` must be a double vector or matrix.", what);
  const std::size_t len = static_cast<std::size_t>(Rf_xlength(x));
  s.v = REAL(x);
  if (Rf_isMatrix(x)) {
    if (static_cast<std::size_t>(Rf_nrows(x)) != n || static_cast<std::size_t>(Rf_ncols(x)) != H)
      Rcpp::stop("cpp_ecl_paths(): `%s` must be %d x %d.", what, static_cast<int>(n), static_cast<int>(H));
    s.mode = 2;
  } else if (len == 1) {
    s.mode = 0;
  } else if (len == n) {
    s.mode = 1;
  } else {
    Rcpp::stop("cpp_ecl_paths(): `%s` has length %d, expected 1 or %d.", what, static_cast<int>(len), static_cast<int>(n));
  }
  return s;
}

}  // namespace

//' Expected credit loss paths of one scenario (internal)
//'
//' @param h Marginal monthly default hazards (scalar, length `n` or `n x H`).
//' @param L,E LGD and EAD, same shapes.
//' @param P Prepayment hazard (same shapes) or `NULL`.
//' @param r Annual effective interest rate, scalar or length `n`.
//' @param n,H,hz Rows, months in the lifetime, months in the 12-month figure.
//' @param discount Discount at `r` with monthly compounding.
//' @param stage3 Logical of length `n`: credit-impaired rows carry
//'   `LGD_1 * EAD_1`.
//' @param z,pd_mult,lgd_add,ead_mult Scenario shocks, each `NULL`, a scalar or
//'   length `n`.
//' @param rho Asset correlation of the one-factor shock.
//' @param nthreads OpenMP threads over rows.
//' @return A list with `ecl_12m`, `ecl_life` (scenario) and `pd_12m`,
//'   `pd_life` (cumulative PD of the unshocked hazards).
//' @keywords internal
//' @noRd
// [[Rcpp::export(rng = false)]]
Rcpp::List cpp_ecl_paths(SEXP h, SEXP L, SEXP E, SEXP P, SEXP r, const int n, const int H, const int hz,
                         const bool discount, const Rcpp::LogicalVector& stage3,
                         SEXP z, SEXP pd_mult, SEXP lgd_add, SEXP ead_mult,
                         const double rho, const int nthreads) {
  if (n < 0 || H < 1 || hz < 1 || hz > H) Rcpp::stop("cpp_ecl_paths(): invalid dimensions.");
  const std::size_t nn = static_cast<std::size_t>(n), HH = static_cast<std::size_t>(H);
  const Src sh = make_src(h, nn, HH, "h"), sL = make_src(L, nn, HH, "lgd"), sE = make_src(E, nn, HH, "ead");
  const Src sP = make_src(P, nn, HH, "prepay");
  const Src sr = make_src(r, nn, 1, "eir");
  const Src sz = make_src(z, nn, 1, "z"), sm = make_src(pd_mult, nn, 1, "pd_mult");
  const Src sa = make_src(lgd_add, nn, 1, "lgd_add"), se = make_src(ead_mult, nn, 1, "ead_mult");
  if (sh.mode < 0 || sL.mode < 0 || sE.mode < 0 || sr.mode < 0) Rcpp::stop("cpp_ecl_paths(): `h`, `lgd`, `ead` and `eir` are required.");
  if (sr.mode == 2 || sz.mode == 2 || sm.mode == 2 || sa.mode == 2 || se.mode == 2)
    Rcpp::stop("cpp_ecl_paths(): `eir` and the scenario shocks must be scalars or vectors.");
  if (static_cast<std::size_t>(stage3.size()) != nn) Rcpp::stop("cpp_ecl_paths(): `stage3` must have length n.");
  if (sz.mode >= 0 && !(rho >= 0.0 && rho < 1.0)) Rcpp::stop("cpp_ecl_paths(): `rho` must be in [0, 1).");

  Rcpp::NumericVector ecl12(nn), ecll(nn), pd12(nn), pdl(nn);
  double* o12 = ecl12.begin();
  double* oll = ecll.begin();
  double* p12 = pd12.begin();
  double* pll = pdl.begin();
  const int* s3 = stage3.begin();
  const double srho = std::sqrt(rho), s1rho = std::sqrt(1.0 - rho);
  const long long N = static_cast<long long>(nn);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(std::max(1, nthreads))
#endif
  for (long long ii = 0; ii < N; ++ii) {
    const std::size_t i = static_cast<std::size_t>(ii);
    const double ri = sr.at(i, 0);
    const double lr = std::log1p(ri);
    const double zi = sz.mode >= 0 ? sz.at(i, 0) : 0.0;
    const double mi = sm.mode >= 0 ? sm.at(i, 0) : 1.0;
    const double ai = sa.mode >= 0 ? sa.at(i, 0) : 0.0;
    const double ei = se.mode >= 0 ? se.at(i, 0) : 1.0;
    double S = 1.0, life = 0.0, m12 = 0.0;
    long double lsum = 0.0L, l12 = 0.0L;
    for (std::size_t t = 0; t < HH; ++t) {
      const double hb = sh.at(i, t);
      lsum += std::log1p(-hb);
      if (t + 1 == static_cast<std::size_t>(hz)) l12 = lsum;
      double ht = hb;
      if (sz.mode >= 0) {
        if (ht <= 0.0) ht = 0.0;
        else if (ht >= 1.0) ht = 1.0;
        else ht = R::pnorm((R::qnorm(ht, 0.0, 1.0, 1, 0) - srho * zi) / s1rho, 0.0, 1.0, 1, 0);
      }
      if (sm.mode >= 0) ht = std::min(ht * mi, 1.0);
      double Lt = sL.at(i, t);
      if (sa.mode >= 0) Lt = std::max(Lt + ai, 0.0);
      const double Et = sE.at(i, t) * ei;
      const double df = discount ? std::exp(-static_cast<double>(t + 1) / 12.0 * lr) : 1.0;
      life += S * ht * Lt * Et * df;
      if (t + 1 == static_cast<std::size_t>(hz)) m12 = life;
      const double exit = sP.mode >= 0 ? ht + sP.at(i, t) : ht;
      S *= 1.0 - std::min(exit, 1.0);
    }
    if (s3[i] == TRUE) {
      double L1 = sL.at(i, 0);
      if (sa.mode >= 0) L1 = std::max(L1 + ai, 0.0);
      life = m12 = L1 * sE.at(i, 0) * ei;
    }
    o12[i] = m12;
    oll[i] = life;
    p12[i] = 0.0 - std::expm1(static_cast<double>(l12));
    pll[i] = 0.0 - std::expm1(static_cast<double>(lsum));
  }
  return Rcpp::List::create(Rcpp::_["ecl_12m"] = ecl12, Rcpp::_["ecl_life"] = ecll,
                            Rcpp::_["pd_12m"] = pd12, Rcpp::_["pd_life"] = pdl);
}
