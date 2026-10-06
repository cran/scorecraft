// correlation.cpp - Pearson / Spearman correlation matrix of many columns
//
// Used by the redundancy pruning of stage 2 on the WOE space. Every column is
// standardised once (ranked first under Spearman, mid-ranks for ties), so the
// correlation matrix is a single cross-product Z'Z computed by BLAS (syrk):
// O(n p log n + n p^2 / BLAS) instead of ranking both columns of every pair,
// O(p^2 n log n).
//
// Columns arrive as a list of double vectors (the columns of a data.table,
// read without copying). Missing values are not accepted: the caller checks
// the stage 1 invariant (the data leave the triage with no NA) and falls back
// to the pairwise routine otherwise.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

// Standardise one column into `z`: centred and scaled to unit Euclidean norm,
// so that the dot product of two standardised columns is their correlation.
// Under `rank`, the values are replaced by their mid-ranks first; the sort is
// done on a copy of the values and each value is mapped back to the average
// rank of its tie group by binary search over the distinct values, which is
// cheap for WOE columns (a handful of distinct values). Returns false when the
// column is constant (zero variance); `z` is then set to zero.
static bool standardise_column(const double* x, const std::size_t n, const bool rank,
                               double* z, std::vector<double>& buf,
                               std::vector<double>& uniq, std::vector<double>& midrank) {
  if (rank) {
    buf.assign(x, x + n);
    std::sort(buf.begin(), buf.end());
    uniq.clear();
    midrank.clear();
    std::size_t a = 0;
    while (a < n) {
      std::size_t b = a;
      while (b < n && buf[b] == buf[a]) ++b;
      uniq.push_back(buf[a]);
      // 1-based ranks a+1 .. b share their average
      midrank.push_back(0.5 * (static_cast<double>(a + 1) + static_cast<double>(b)));
      a = b;
    }
    for (std::size_t i = 0; i < n; ++i) {
      const std::size_t k = static_cast<std::size_t>(
        std::lower_bound(uniq.begin(), uniq.end(), x[i]) - uniq.begin());
      z[i] = midrank[k];
    }
  } else {
    std::copy(x, x + n, z);
  }
  // two-pass centring for numerical stability (no sum-of-squares cancellation)
  long double s = 0.0L;
  for (std::size_t i = 0; i < n; ++i) s += z[i];
  const double m = static_cast<double>(s / static_cast<long double>(n));
  long double ss = 0.0L;
  for (std::size_t i = 0; i < n; ++i) {
    z[i] -= m;
    ss += static_cast<long double>(z[i]) * z[i];
  }
  const double nrm = std::sqrt(static_cast<double>(ss));
  if (!(nrm > 0.0) || !std::isfinite(nrm)) {
    std::fill(z, z + n, 0.0);
    return false;
  }
  for (std::size_t i = 0; i < n; ++i) z[i] /= nrm;
  return true;
}

//' Correlation matrix of a list of numeric columns (internal)
//'
//' @param cols List of `p` double vectors of equal length `n`, no missing
//'   values.
//' @param spearman `TRUE` for Spearman (Pearson on mid-ranks), `FALSE` for
//'   Pearson.
//' @param nthreads Threads used to standardise the columns (OpenMP); the
//'   cross-product runs on the BLAS R is linked to.
//' @return A `p x p` matrix; rows and columns of a constant column are `NA`
//'   (the diagonal included).
//' @keywords internal
//' @noRd
// [[Rcpp::export(rng = false)]]
arma::mat cpp_cor_matrix(const Rcpp::List& cols, const bool spearman, const int nthreads) {
  const std::size_t p = static_cast<std::size_t>(cols.size());
  if (p == 0) return arma::mat(0, 0);
  std::vector<const double*> ptr(p);
  std::size_t n = 0;
  for (std::size_t j = 0; j < p; ++j) {
    SEXP cj = cols[j];
    if (TYPEOF(cj) != REALSXP) Rcpp::stop("cpp_cor_matrix(): column %d is not a double vector.", j + 1);
    const std::size_t nj = static_cast<std::size_t>(Rf_xlength(cj));
    if (j == 0) n = nj;
    else if (nj != n) Rcpp::stop("cpp_cor_matrix(): columns differ in length.");
    ptr[j] = REAL(cj);
    for (std::size_t i = 0; i < nj; ++i) {
      if (!std::isfinite(ptr[j][i])) Rcpp::stop("cpp_cor_matrix(): column %d has a non-finite value.", j + 1);
    }
  }
  if (n < 2) Rcpp::stop("cpp_cor_matrix(): at least two rows are needed.");

  arma::mat Z(n, p);
  std::vector<int> ok(p, 1);
  double* zbase = Z.memptr();
  const long long pp = static_cast<long long>(p);

#ifdef _OPENMP
#pragma omp parallel num_threads(std::max(1, nthreads))
#endif
  {
    std::vector<double> buf, uniq, midrank;   // scratch owned by each thread
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
    for (long long j = 0; j < pp; ++j) {
      ok[j] = standardise_column(ptr[j], n, spearman, zbase + static_cast<std::size_t>(j) * n,
                                 buf, uniq, midrank) ? 1 : 0;
    }
  }

  arma::mat C = Z.t() * Z;   // syrk
  Z.reset();                 // release n x p before building the result
  for (std::size_t j = 0; j < p; ++j) {
    if (!ok[j]) {
      C.row(j).fill(NA_REAL);
      C.col(j).fill(NA_REAL);
    } else {
      C(j, j) = 1.0;
    }
  }
  // clip rounding excursions beyond [-1, 1]
  C.transform([](double v) { return std::isnan(v) ? v : std::max(-1.0, std::min(1.0, v)); });
  return C;
}
