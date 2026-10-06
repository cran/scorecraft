# ============================================================================ #
# cpp.R - thin R wrappers around the compiled kernels (src/*.cpp)
# ============================================================================ #
# The kernels take plain vectors and matrices, never a data.table, and never
# touch the R API from a worker thread. Validation, dispatch and the fall-back
# to the reference R code live here, so every caller sees an R function with
# the contract of the code it replaced.
# ============================================================================ #

#' @useDynLib scorecraft, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' Threads for a compiled kernel
#'
#' `nthread` capped at two under `R CMD check --as-cran`
#' (`_R_CHECK_LIMIT_CORES_`), as `.scr_lapply()` does, and at least one.
#' @keywords internal
#' @noRd
.scr_threads <- function(nthread = 1L) {
  k <- suppressWarnings(as.integer(nthread %||% 1L))
  if (length(k) != 1L || is.na(k) || k < 1L) k <- 1L
  lim <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(lim) && !identical(tolower(lim), "false")) k <- min(k, 2L)
  k
}

#' Correlation matrix of the columns of a list or data.frame
#'
#' Pearson, or Spearman (Pearson on mid-ranks), through one cross-product of
#' the standardised columns (`cpp_cor_matrix()`). The columns are read in
#' place; only integer or logical columns are converted to double.
#' @return A `p x p` matrix with dimnames; constant columns give `NA`.
#' @keywords internal
#' @noRd
.scr_cor_matrix <- function(x, method = c("pearson", "spearman"), nthread = 1L) {
  method <- match.arg(method)
  cols <- lapply(x, function(v) if (is.double(v)) v else as.double(v))
  m <- cpp_cor_matrix(cols, identical(method, "spearman"), .scr_threads(nthread))
  dimnames(m) <- list(names(x), names(x))
  m
}

#' Greedy redundancy pruning on a correlation matrix
#'
#' The algorithm of `OptimalBinningWoE::obwoe_prune()`, reproduced exactly:
#' while a surviving pair has an absolute correlation at or above `cutoff`,
#' drop the worse-ranked member of the strongest one (the first in the
#' pairwise order `(1,2), (1,3), ..., (2,3), ...` on ties); variables absent
#' from `ranking` rank last. Repeated arg-max over the surviving pairs is the
#' same as one sweep of the pairs sorted by decreasing absolute correlation
#' (stable on the pairwise order) that skips pairs with a dead member, so the
#' cost is `O(P log P)` over the `P` pairs above the cutoff instead of one
#' scan of every pair per removal.
#' @param cm Correlation matrix with dimnames (`.scr_cor_matrix()`).
#' @return `list(keep, dropped)`, `dropped` a data.frame of `variable`,
#'   `correlated_with` and `correlation` in the order of removal; `keep` the
#'   survivors in the order of `ranking` (variables outside `ranking` are not
#'   kept, as in `obwoe_prune()`).
#' @keywords internal
#' @noRd
.scr_prune_matrix <- function(cm, ranking, cutoff) {
  vars <- colnames(cm)
  p <- length(vars)
  rank_of <- match(vars, ranking)
  rank_of[is.na(rank_of)] <- length(ranking) + 1L
  ut <- which(upper.tri(cm), arr.ind = TRUE)
  # pairwise order of obcorr(): i ascending, then j ascending
  ut <- ut[order(ut[, "row"], ut[, "col"]), , drop = FALSE]
  r <- cm[ut]
  hit <- which(!is.na(r) & abs(r) >= cutoff)
  alive <- rep(TRUE, p)
  drops <- vector("list", length(hit))
  nd <- 0L
  if (length(hit)) {
    o <- hit[order(-abs(r[hit]), hit)]   # decreasing |r|, stable on the pairwise order
    for (k in o) {
      a <- ut[k, "row"]; b <- ut[k, "col"]
      if (!alive[a] || !alive[b]) next
      loser <- if (rank_of[a] > rank_of[b]) a else b
      winner <- if (loser == a) b else a
      alive[loser] <- FALSE
      nd <- nd + 1L
      drops[[nd]] <- list(variable = vars[loser], correlated_with = vars[winner], correlation = r[k])
    }
  }
  dropped <- if (nd) data.frame(variable = vapply(drops[seq_len(nd)], `[[`, "", "variable"),
                                correlated_with = vapply(drops[seq_len(nd)], `[[`, "", "correlated_with"),
                                correlation = vapply(drops[seq_len(nd)], `[[`, 0, "correlation"),
                                stringsAsFactors = FALSE)
             else data.frame(variable = character(), correlated_with = character(),
                             correlation = numeric(), stringsAsFactors = FALSE)
  list(keep = ranking[ranking %in% vars[alive]], dropped = dropped)
}

#' Somers' D of a prediction against a realised outcome
#'
#' `(C - D) / (n(n-1)/2 - T_r)`: concordant minus discordant pairs over the
#' pairs not tied on the outcome, counted exactly in `O(n log n)`
#' (`cpp_concordance()`, Knight 1966). Pairs with a non-finite member are
#' dropped first.
#' @param p Prediction. @param r Realised outcome.
#' @param const_p What a constant prediction returns: `NA` or `0`.
#' @keywords internal
#' @noRd
.scr_somers <- function(p, r, const_p = NA_real_) {
  ok <- is.finite(p) & is.finite(r)
  p <- as.double(p[ok]); r <- as.double(r[ok])
  if (length(p) < 2L) return(NA_real_)
  cc <- cpp_concordance(p, r)
  den <- cc[["pairs"]] - cc[["ties_r"]]
  if (den <= 0) return(NA_real_)
  if (cc[["pairs"]] - cc[["ties_p"]] <= 0) return(const_p)
  cc[["cmd"]] / den
}

#' @keywords internal
#' @noRd
.onUnload <- function(libpath) {
  library.dynam.unload("scorecraft", libpath)
}
