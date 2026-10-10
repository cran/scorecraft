// tiers.cpp - optimal contiguous segmentation of ordered blocks into tiers
//
// The blocks are pool-adjacent-violator blocks of score pre-bins, ordered so
// that the event rate does not decrease along the index. The dynamic
// program picks L contiguous segments that maximize the binomial
// log-likelihood, sum_t [e_t log p_t + (n_t - e_t) log(1 - p_t)], or the
// information value, subject to: a minimum volume share per segment; a
// minimum number of events and of non-events per segment; and adjacent
// segments distinct by a one-sided Fisher exact test at alpha / (L - 1).
// The adjacency test needs the previous segment in the state, so the state
// is (number of segments, end of the last segment, start of the last
// segment) and the cost is O(L M^3).

#include <Rcpp.h>
#include <cmath>
#include <limits>
#include <vector>

namespace {

const double NEG_INF = -std::numeric_limits<double>::infinity();

// x log(x / m), with 0 log 0 = 0
inline double xlogr(double x, double m) {
  return (x > 0.0 && m > 0.0) ? x * std::log(x / m) : 0.0;
}

}  // namespace

//' Optimal tier segmentation over ordered blocks (internal)
//'
//' @param e,n Weighted events and weighted volume with a known outcome per
//'   block (the objective).
//' @param vol Weighted volume per block (the share constraint).
//' @param e_raw,n_raw Unweighted events and rows with a known outcome per
//'   block (the event constraints and the Fisher exact test).
//' @param L Number of segments.
//' @param criterion 0 for the binomial log-likelihood, 1 for the
//'   information value.
//' @param min_share,min_events,alpha The constraints.
//' @return A list: `feasible`, `ends` (1-based end block of every segment,
//'   empty when infeasible) and `objective`.
//' @keywords internal
//' @noRd
// [[Rcpp::export(rng = false)]]
Rcpp::List cpp_tier_dp(const Rcpp::NumericVector& e, const Rcpp::NumericVector& n,
                       const Rcpp::NumericVector& vol, const Rcpp::NumericVector& e_raw,
                       const Rcpp::NumericVector& n_raw, const int L, const int criterion,
                       const double min_share, const double min_events, const double alpha) {
  const int M = e.size();
  if (n.size() != M || vol.size() != M || e_raw.size() != M || n_raw.size() != M) {
    Rcpp::stop("cpp_tier_dp(): the block vectors differ in length.");
  }
  Rcpp::List none = Rcpp::List::create(Rcpp::_["feasible"] = false,
                                       Rcpp::_["ends"] = Rcpp::IntegerVector(0),
                                       Rcpp::_["objective"] = NA_REAL);
  if (L < 1 || M < L) return none;

  // prefix sums: segment [a, b] (0-based, inclusive) is P[b + 1] - P[a]
  std::vector<double> E(M + 1, 0.0), N(M + 1, 0.0), V(M + 1, 0.0), ER(M + 1, 0.0), NR(M + 1, 0.0);
  for (int i = 0; i < M; ++i) {
    E[i + 1] = E[i] + e[i];
    N[i + 1] = N[i] + n[i];
    V[i + 1] = V[i] + vol[i];
    ER[i + 1] = ER[i] + e_raw[i];
    NR[i + 1] = NR[i] + n_raw[i];
  }
  const double e_tot = E[M], ne_tot = N[M] - E[M], v_tot = V[M];
  // a relative slack of 1e-12 keeps a share equal to min_share feasible
  // despite rounding in the prefix sums
  const double v_min = (min_share - 1e-12) * v_tot;
  const double thr = L > 1 ? alpha / (L - 1) : alpha;

  // objective of one segment; -Inf when the segment breaks a constraint
  auto seg = [&](int a, int b) -> double {
    const double sv = V[b + 1] - V[a], se = E[b + 1] - E[a], sn = N[b + 1] - N[a];
    const double ser = ER[b + 1] - ER[a], snr = NR[b + 1] - NR[a];
    if (sv < v_min || ser < min_events || snr - ser < min_events) return NEG_INF;
    if (criterion == 0) {
      if (sn <= 0.0) return 0.0;
      return xlogr(se, sn) + xlogr(sn - se, sn);
    }
    // information value: undefined for a segment without events or non-events
    if (!(e_tot > 0.0) || !(ne_tot > 0.0)) return NEG_INF;
    const double pe = se / e_tot, pn = (sn - se) / ne_tot;
    if (!(pe > 0.0) || !(pn > 0.0)) return NEG_INF;
    return (pe - pn) * std::log(pe / pn);
  };

  // one-sided Fisher exact test that segment [a, b] has a higher event rate
  // than the previous segment [p, a - 1]; cached over (p, a, b) when small
  const bool cache = M <= 200;
  std::vector<signed char> memo(cache ? static_cast<std::size_t>(M) * M * M : 0, 0);
  auto distinct = [&](int p, int a, int b) -> bool {
    std::size_t key = 0;
    if (cache) {
      key = (static_cast<std::size_t>(p) * M + a) * M + b;
      if (memo[key] != 0) return memo[key] > 0;
    }
    const double e1 = ER[b + 1] - ER[a], n1 = NR[b + 1] - NR[a];
    const double e2 = ER[a] - ER[p], n2 = NR[a] - NR[p];
    const double pv = (e1 + e2 <= 0.0) ? 1.0 : R::phyper(e1 - 1.0, n1, n2, e1 + e2, false, false);
    const bool ok = pv < thr;
    if (cache) memo[key] = ok ? 1 : -1;
    return ok;
  };

  // best[(l * M + b) * M + a]: best objective of blocks 0..b in l + 1
  // segments, the last one being [a, b]; arg holds the start of the one before
  const std::size_t sz = static_cast<std::size_t>(L) * M * M;
  std::vector<double> best(sz, NEG_INF);
  std::vector<int> arg(sz, -1);
  auto at = [M](int l, int b, int a) { return (static_cast<std::size_t>(l) * M + b) * M + a; };

  for (int b = 0; b < M; ++b) best[at(0, b, 0)] = seg(0, b);
  for (int l = 1; l < L; ++l) {
    for (int b = l; b < M; ++b) {
      for (int a = l; a <= b; ++a) {
        const double s = seg(a, b);
        if (s == NEG_INF) continue;
        double cur = NEG_INF;
        int who = -1;
        for (int p = l - 1; p <= a - 1; ++p) {
          const double prev = best[at(l - 1, a - 1, p)];
          if (prev == NEG_INF) continue;
          if (!distinct(p, a, b)) continue;
          if (prev + s > cur) { cur = prev + s; who = p; }
        }
        if (who >= 0) { best[at(l, b, a)] = cur; arg[at(l, b, a)] = who; }
      }
    }
  }

  // the last segment ends at block M - 1; the first maximum wins
  double opt = NEG_INF;
  int a_opt = -1;
  for (int a = L - 1; a < M; ++a) {
    const double v = best[at(L - 1, M - 1, a)];
    if (v > opt) { opt = v; a_opt = a; }
  }
  if (a_opt < 0) return none;

  Rcpp::IntegerVector ends(L);
  int b = M - 1, a = a_opt;
  for (int l = L - 1; l >= 0; --l) {
    ends[l] = b + 1;
    const int p = arg[at(l, b, a)];
    b = a - 1;
    a = p;
  }
  return Rcpp::List::create(Rcpp::_["feasible"] = true, Rcpp::_["ends"] = ends,
                            Rcpp::_["objective"] = opt);
}
