// concordance.cpp - pair counts for Somers' D in O(n log n)
//
// Somers' D of a prediction against a realised continuous outcome (LGD, CCF)
// needs, over all n(n-1)/2 pairs, the concordant minus the discordant count
// and the number of pairs tied on the outcome. Knight's (1966) algorithm
// sorts by the prediction and counts, for each observation, how many earlier
// observations (strictly lower prediction) have a lower / higher outcome with
// a Fenwick tree over the dense ranks of the outcome. Observations that tie
// on the prediction are scored before any of them is inserted, so pairs tied
// on the prediction count neither way.
//
// Knight, W. R. (1966). A computer method for calculating Kendall's tau with
// ungrouped data. Journal of the American Statistical Association, 61(314),
// 436-439.

#include <Rcpp.h>
#include <algorithm>
#include <cstdint>
#include <numeric>
#include <vector>

namespace {

// number of pairs tied within the groups of equal value of x[order]
double tied_pairs(const std::vector<std::size_t>& order, const double* x) {
  const std::size_t n = order.size();
  double t = 0.0;
  std::size_t a = 0;
  while (a < n) {
    std::size_t b = a;
    while (b < n && x[order[b]] == x[order[a]]) ++b;
    const double c = static_cast<double>(b - a);
    t += c * (c - 1.0) / 2.0;
    a = b;
  }
  return t;
}

}  // namespace

//' Concordance counts between a prediction and an outcome (internal)
//'
//' @param p,r Double vectors of equal length, finite (the caller drops
//'   non-finite pairs).
//' @return A named double vector: `cmd` (concordant minus discordant pairs),
//'   `pairs` (`n(n-1)/2`), `ties_p` and `ties_r` (pairs tied on `p` and on
//'   `r`). Counts are exact up to 2^53.
//' @keywords internal
//' @noRd
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector cpp_concordance(const Rcpp::NumericVector& p, const Rcpp::NumericVector& r) {
  const std::size_t n = static_cast<std::size_t>(p.size());
  if (static_cast<std::size_t>(r.size()) != n) Rcpp::stop("cpp_concordance(): `p` and `r` differ in length.");
  const double* pp = p.begin();
  const double* rr = r.begin();
  for (std::size_t i = 0; i < n; ++i) {
    if (!std::isfinite(pp[i]) || !std::isfinite(rr[i])) Rcpp::stop("cpp_concordance(): non-finite value.");
  }

  // dense ranks of r (1..K)
  std::vector<std::size_t> ir(n);
  std::iota(ir.begin(), ir.end(), 0);
  std::sort(ir.begin(), ir.end(), [rr](std::size_t a, std::size_t b) { return rr[a] < rr[b]; });
  std::vector<std::size_t> rk(n);
  std::size_t K = 0;
  for (std::size_t k = 0; k < n; ++k) {
    if (k == 0 || rr[ir[k]] != rr[ir[k - 1]]) ++K;
    rk[ir[k]] = K;
  }

  std::vector<std::size_t> ip(n);
  std::iota(ip.begin(), ip.end(), 0);
  std::sort(ip.begin(), ip.end(), [pp](std::size_t a, std::size_t b) { return pp[a] < pp[b]; });

  std::vector<std::int64_t> fen(K + 1, 0);
  auto add = [&fen, K](std::size_t i) { for (; i <= K; i += i & (~i + 1)) fen[i] += 1; };
  auto below = [&fen](std::size_t i) {   // count of inserted ranks <= i
    std::int64_t s = 0;
    for (; i > 0; i -= i & (~i + 1)) s += fen[i];
    return s;
  };

  std::int64_t conc = 0, disc = 0, inserted = 0;
  std::size_t g = 0;
  while (g < n) {
    std::size_t h = g;
    while (h < n && pp[ip[h]] == pp[ip[g]]) ++h;
    for (std::size_t k = g; k < h; ++k) {
      const std::size_t q = rk[ip[k]];
      conc += below(q - 1);                  // earlier: lower p and lower r
      disc += inserted - below(q);           // earlier: lower p and higher r
    }
    for (std::size_t k = g; k < h; ++k) add(rk[ip[k]]);
    inserted += static_cast<std::int64_t>(h - g);
    g = h;
  }

  const double nn = static_cast<double>(n);
  Rcpp::NumericVector out = Rcpp::NumericVector::create(
    Rcpp::_["cmd"] = static_cast<double>(conc - disc),
    Rcpp::_["pairs"] = nn * (nn - 1.0) / 2.0,
    Rcpp::_["ties_p"] = tied_pairs(ip, pp),
    Rcpp::_["ties_r"] = tied_pairs(ir, rr));
  return out;
}
