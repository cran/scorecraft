# ============================================================================ #
# study.R - score-study engine: one pass into a count table, O(K) afterwards
# ============================================================================ #
# A score study never goes back to the rows. The scored samples are
# aggregated once into a table of counts per distinct score cell, keyed by
# sample (and by group when asked); every band, tier, test and light is a
# function of that table, so its cost grows with the number of distinct
# scores K, not with the number of rows. Bands are frozen on the reference
# sample and applied to the other samples with the convention of
# scr_cutoff(): the upper side of a cut is `score >= cut`.
# ============================================================================ #

# Count columns of a study histogram: summed whenever cells are pooled
.study_count_cols <- c("n", "n_y", "e", "w2", "w2_y", "n_raw", "n_y_raw", "e_raw", "v", "ve", "ep")

# -- the count table ---------------------------------------------------------- #

#' Counts per distinct score cell, in one grouped pass
#'
#' Per cell (and per key of `by`): `n` (weighted volume), `n_y` (weighted
#' volume with a known outcome, the rate denominator), `e` (weighted
#' events), `w2` and `w2_y` (sums of squared weights, for the Kish effective
#' size), the unweighted `n_raw`, `n_y_raw` and `e_raw`, and, when given,
#' `v`/`ve` (value and event value) and `ep` (expected events). Rows with a
#' missing or infinite score, or with a zero weight, are not counted; a row
#' with a missing outcome counts in the volume only. `s_lo` and `s_hi` are the
#' smallest and largest score of the cell.
#'
#' With more than `max_cells` distinct scores, the scores are pooled into at
#' most `max_cells` buckets of equal weighted share before the grouped pass
#' (one sort of the scores, so no cell per distinct value is ever built);
#' each bucket edge sits midway between two adjacent distinct scores, and
#' the finite `breaks` are added to the edges so that no cell straddles
#' them. The edges are stored in the attribute `edges` and the bucket of
#' each cell in the column `g`.
#' @keywords internal
#' @noRd
.study_hist <- function(score, y, w = NULL, value = NULL, prob = NULL, by = NULL, max_cells = 1e5,
                        breaks = NULL, fn = "scr_bands") {
  score <- as.double(score)
  y <- .scr_y01(y, fn)
  ok <- is.finite(score)
  wtd <- !is.null(w)
  if (wtd) {
    w <- as.double(w)
    if (anyNA(w) || any(!is.finite(w)) || any(w < 0)) {
      stop(fn, "(): the weights must be finite and non-negative.", call. = FALSE)
    }
    # a row with zero weight does not belong to the population
    ok <- ok & w > 0
  }
  yo <- y[ok]
  k <- !is.na(yo)
  yr <- k & yo == 1L   # FALSE & NA is FALSE: no NA left
  cols <- list(s = score[ok], k = as.integer(k), yr = as.integer(yr))
  for (b in names(by)) cols[[b]] <- by[[b]][ok]
  agg <- alist(n_raw = .N, n_y_raw = sum(k), e_raw = sum(yr))
  if (wtd) {
    wv <- w[ok]
    cols$w <- wv; cols$wk <- wv * k; cols$we <- wv * yr
    cols$wsq <- wv * wv; cols$wsqk <- wv * wv * k
    agg <- c(agg, alist(n = sum(w), n_y = sum(wk), e = sum(we), w2 = sum(wsq), w2_y = sum(wsqk)))
  }
  if (!is.null(value)) {
    vv <- as.double(value)[ok]
    vv[is.na(vv)] <- 0
    if (wtd) vv <- vv * wv
    cols$v <- vv; cols$vy <- vv * yr
    agg <- c(agg, alist(v = sum(v), ve = sum(vy)))
  }
  if (!is.null(prob)) {
    pp <- as.double(prob)[ok]
    if (any(k & (is.na(pp) | pp < 0 | pp > 1))) {
      stop(fn, "(): the expected probability must lie in [0, 1] on every row with a known outcome.", call. = FALSE)
    }
    pp[!k] <- 0
    if (wtd) pp <- pp * wv
    cols$p <- pp
    agg <- c(agg, alist(ep = sum(p)))
  }
  # many distinct scores: pool them into buckets first, so the grouped pass
  # below never builds one cell per distinct value
  edges <- NULL
  if (length(cols$s) > max_cells && data.table::uniqueN(cols$s) > max_cells) {
    edges <- .study_edges(cols$s, if (wtd) wv, max_cells, breaks)
    cols$g <- findInterval(cols$s, edges) + 1L
  }
  dt <- data.table::setDT(cols)
  keys <- c(names(by), "s")
  # one grouped pass, optimized by GForce (sums, min, max and .N only); `j`
  # must be built beforehand, an inline call to as.call() defeats GForce
  if (is.null(edges)) {
    j <- as.call(c(as.name("list"), agg))
    h <- dt[, eval(j), keyby = keys]
    h[, `:=`(s_lo = s, s_hi = s)]
  } else {
    j <- as.call(c(as.name("list"), alist(s_lo = min(s), s_hi = max(s)), agg))
    h <- dt[, eval(j), keyby = c(names(by), "g")]
    h[, s := s_lo]
    data.table::setkeyv(h, keys)
    data.table::setattr(h, "edges", edges)
  }
  for (cn in c("n_raw", "n_y_raw", "e_raw")) data.table::set(h, j = cn, value = as.double(h[[cn]]))
  if (!wtd) h[, `:=`(n = n_raw, n_y = n_y_raw, e = e_raw, w2 = n_raw, w2_y = n_y_raw)]
  data.table::setcolorder(h, c(keys, "s_lo", "s_hi", intersect(.study_count_cols, names(h))))
  data.table::setattr(h, "weighted", wtd)
  h
}

#' Bucket edges of equal weighted share over the distinct scores
#'
#' `u` are the distinct scores in ascending order and `cu` the cumulative
#' weight up to each; an edge sits midway between the two distinct scores
#' where the bucket changes, so a group of tied scores is never split.
#' @keywords internal
#' @noRd
.study_edges_u <- function(u, cu, max_cells, breaks = NULL) {
  b <- pmin(pmax(ceiling(cu / cu[length(cu)] * max_cells), 1), max_cells)
  at <- which(diff(b) != 0)
  edges <- .study_mid(u[at], u[at + 1L])
  br <- as.double(breaks)
  sort(unique(c(edges, br[is.finite(br)])))
}

#' The same edges from the raw scores (and weights), after one sort
#' @keywords internal
#' @noRd
.study_edges <- function(s, w = NULL, max_cells, breaks = NULL) {
  o <- order(s, method = "radix")
  ss <- s[o]
  cw <- if (is.null(w)) as.double(seq_along(ss)) else cumsum(as.double(w[o]))
  last <- c(which(diff(ss) > 0), length(ss))
  .study_edges_u(ss[last], cw[last], max_cells, breaks)
}

#' Pre-aggregated counts as a study histogram
#'
#' Each row is a score cell with `n` rows and `events` events; the counts
#' are taken as frequencies, so they serve both as the weighted and as the
#' unweighted counts (the tests round them). Rows with a missing or infinite
#' score, or a zero count, are dropped; repeated scores are summed.
#' @keywords internal
#' @noRd
.study_hist_counts <- function(score, n, events, value = NULL, value_events = NULL, by = NULL,
                               max_cells = 1e5, breaks = NULL, fn = "scr_bands") {
  score <- as.double(score); n <- as.double(n); events <- as.double(events)
  if (anyNA(n) || anyNA(events) || any(!is.finite(n)) || any(!is.finite(events)) ||
      any(n < 0) || any(events < 0) || any(events > n)) {
    stop(fn, "(): the counts must be finite, non-negative, with events <= n on every row.", call. = FALSE)
  }
  ok <- is.finite(score) & n > 0
  cols <- list(s = score[ok], n = n[ok], e = events[ok])
  for (b in names(by)) cols[[b]] <- by[[b]][ok]
  agg <- alist(n = sum(n), e = sum(e))
  if (!is.null(value)) {
    vv <- as.double(value)[ok]; vv[is.na(vv)] <- 0
    ve <- if (is.null(value_events)) rep(0, length(vv)) else as.double(value_events)[ok]
    ve[is.na(ve)] <- 0
    cols$v <- vv; cols$vy <- ve
    agg <- c(agg, alist(v = sum(v), ve = sum(vy)))
  }
  dt <- data.table::setDT(cols)
  keys <- c(names(by), "s")
  j <- as.call(c(as.name("list"), agg))
  h <- dt[, eval(j), keyby = keys]
  h[, `:=`(n_y = n, w2 = n, w2_y = n, n_raw = n, n_y_raw = n, e_raw = e, s_lo = s, s_hi = s)]
  data.table::setcolorder(h, c(keys, "s_lo", "s_hi", intersect(.study_count_cols, names(h))))
  h <- .study_quantize(h, names(by), max_cells, breaks)
  data.table::setattr(h, "weighted", FALSE)
  h
}

#' Midpoint strictly above `a` and at most `b` (adjacent doubles can round
#' the midpoint onto the lower value)
#' @keywords internal
#' @noRd
.study_mid <- function(a, b) {
  m <- (a + b) / 2
  low <- !(m > a)
  m[low] <- b[low]
  m
}

#' Pool the cells into weighted-quantile buckets when there are too many
#' @keywords internal
#' @noRd
.study_quantize <- function(h, keys, max_cells, breaks = NULL) {
  u <- h[, list(n = sum(n)), keyby = "s"]
  if (nrow(u) <= max_cells) return(h)
  edges <- .study_edges_u(u$s, cumsum(u$n), max_cells, breaks)
  h[, g := findInterval(s, edges) + 1L]
  cnt <- intersect(.study_count_cols, names(h))
  q <- h[, c(list(s = min(s), s_lo = min(s_lo), s_hi = max(s_hi)), lapply(.SD, sum)),
         keyby = c(keys, "g"), .SDcols = cnt]
  data.table::setkeyv(q, c(keys, "s"))
  data.table::setattr(q, "edges", edges)
  q
}

#' Pool the cells of a histogram across samples and groups
#'
#' Quantized cells are pooled by bucket, exact cells by score. The result is
#' keyed by score and keeps the bucket edges.
#' @keywords internal
#' @noRd
.study_collapse <- function(h) {
  cnt <- intersect(.study_count_cols, names(h))
  if (!nrow(h)) {
    # an empty sample: no cell, the same columns
    out <- h[0L, c(if ("g" %in% names(h)) "g", "s", "s_lo", "s_hi", cnt), with = FALSE]
  } else if ("g" %in% names(h)) {
    out <- h[, c(list(s = min(s), s_lo = min(s_lo), s_hi = max(s_hi)), lapply(.SD, sum)),
             keyby = "g", .SDcols = cnt]
  } else {
    out <- h[, c(list(s_lo = min(s_lo), s_hi = max(s_hi)), lapply(.SD, sum)), keyby = "s", .SDcols = cnt]
  }
  data.table::setkeyv(out, "s")
  data.table::setattr(out, "edges", attr(h, "edges"))
  data.table::setattr(out, "weighted", isTRUE(attr(h, "weighted")))
  out
}

#' Rows of a histogram for one sample (and group), pooled to one cell per score
#' @keywords internal
#' @noRd
.study_cells <- function(h, sample, group = NULL) {
  i <- h[["sample"]] %in% sample
  if (!is.null(group)) i <- i & h[["group"]] %in% group
  .study_collapse(h[i])
}

# -- cut points ----------------------------------------------------------------- #

#' Boundaries closest to share targets
#'
#' `w` are cell volumes in the order of accumulation; the boundary `j`
#' (`1..K-1`) leaves the first `j` cells on one side. For every target share
#' the nearest boundary is taken (the lower one on a tie); duplicates are
#' dropped. A target is therefore hit within the share of one cell.
#' @keywords internal
#' @noRd
.study_bounds <- function(w, targets) {
  K <- length(w)
  tot <- sum(w)
  if (K < 2L || !(tot > 0) || !length(targets)) return(integer())
  cb <- cumsum(w)[-K] / tot
  i <- findInterval(targets, cb)
  lo <- pmax(i, 1L); hi <- pmin(i + 1L, K - 1L)
  j <- ifelse(abs(cb[hi] - targets) < abs(cb[lo] - targets), hi, lo)
  sort(unique(as.integer(j)))
}

#' Tie-safe cut points on the reference cells
#'
#' Targets are cumulative weighted shares counted from the event-rich side
#' (`"high"` under higher_is_riskier, `"low"` under higher_is_safer):
#' `k / n_bands` for `"uniform"`, `tail_probs` for `"tail"`. Every cut sits
#' between two adjacent cells, midway between their distinct scores (or on
#' the bucket edge of a quantized table), so no tie is ever split and no cut
#' equals an observed reference score. The band of a score is
#' `findInterval(score, cuts) + 1`, left-closed.
#' @keywords internal
#' @noRd
.study_cuts <- function(hist, n_bands, spacing = c("uniform", "tail"), tail_probs = NULL,
                        event_side = c("high", "low")) {
  spacing <- match.arg(spacing)
  event_side <- match.arg(event_side)
  targets <- if (identical(spacing, "tail")) {
    tp <- as.double(tail_probs %||% c(0.001, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50))
    sort(unique(tp[is.finite(tp) & tp > 0 & tp < 1]))
  } else if (n_bands >= 2L) seq_len(n_bands - 1L) / n_bands else numeric()
  req <- length(targets) + 1L
  K <- nrow(hist)
  o <- if (identical(event_side, "high")) rev(seq_len(K)) else seq_len(K)
  j <- .study_bounds(hist$n[o], targets)
  if (!length(j)) return(list(cuts = numeric(), n_bands_requested = req, n_bands_effective = 1L))
  # lower ascending cell of every boundary
  a <- if (identical(event_side, "high")) K - j else j
  edges <- attr(hist, "edges")
  cuts <- if (!is.null(edges) && "g" %in% names(hist)) edges[hist$g[a]] else
    .study_mid(hist$s_hi[a], hist$s_lo[a + 1L])
  cuts <- sort(unique(cuts))
  list(cuts = cuts, n_bands_requested = req, n_bands_effective = length(cuts) + 1L)
}

#' Event-rich side of the score
#' @keywords internal
#' @noRd
.study_side <- function(direction) if (identical(direction, "higher_is_riskier")) "high" else "low"

#' Short number for an interval label; digits grow until the labels differ
#' @keywords internal
#' @noRd
.study_num <- function(x) {
  for (d in 7:15) {
    s <- vapply(x, format, character(1), digits = d, scientific = FALSE, trim = TRUE)
    if (!anyDuplicated(s[is.finite(x)])) break
  }
  s
}

#' Interval labels `[lo, hi)` in ascending score order
#' @keywords internal
#' @noRd
.study_labels <- function(cuts) {
  b <- .study_num(c(-Inf, cuts, Inf))
  k <- length(b)
  paste0("[", b[-k], ", ", b[-1L], ")")
}

# -- band statistics ------------------------------------------------------------ #

#' Sums of count columns per band index (`1..B`, zero for an empty band)
#' @keywords internal
#' @noRd
.study_sum <- function(h, idx, B, cols) {
  out <- lapply(cols, function(cn) numeric(B))
  names(out) <- cols
  if (length(idx)) {
    m <- rowsum(as.matrix(h[, cols, with = FALSE]), idx, reorder = TRUE)
    at <- as.integer(rownames(m))
    for (cn in cols) out[[cn]][at] <- m[, cn]
  }
  out
}

#' Jeffreys interval of a binomial proportion
#'
#' Beta(x + 1/2, n - x + 1/2) quantiles, with the lower bound set to 0 when
#' `x = 0` and the upper bound to 1 when `x = n` (Brown, Cai and DasGupta,
#' 2001). `n` may be a Kish effective size; `NA` when `n` is not positive.
#' @keywords internal
#' @noRd
.study_jeffreys <- function(x, n, level) {
  a <- (1 - level) / 2
  ok <- is.finite(n) & n > 0 & is.finite(x)
  lo <- hi <- rep(NA_real_, length(n))
  xo <- pmin(pmax(x[ok], 0), n[ok]); no <- n[ok]
  lo[ok] <- ifelse(xo <= 0, 0, stats::qbeta(a, xo + 0.5, no - xo + 0.5))
  hi[ok] <- ifelse(xo >= no, 1, stats::qbeta(1 - a, xo + 0.5, no - xo + 0.5))
  list(lo = lo, hi = hi)
}

#' Kish effective size: (sum w)^2 / sum w^2
#' @keywords internal
#' @noRd
.study_kish <- function(n, w2) ifelse(w2 > 0, n * n / w2, 0)

#' One-sided Fisher exact test that group 1 has the higher event rate
#'
#' `P(X >= e1)` for the events of group 1 under the hypergeometric law with
#' the margins fixed; 1 when the pair has no events. Counts are rounded.
#' @keywords internal
#' @noRd
.study_fisher_gt <- function(e1, n1, e2, n2) {
  e1 <- round(e1); n1 <- round(n1); e2 <- round(e2); n2 <- round(n2)
  p <- suppressWarnings(stats::phyper(e1 - 1, n1, n2, e1 + e2, lower.tail = FALSE))
  p[!is.na(e1 + e2) & (e1 + e2) <= 0] <- 1
  p
}

#' Holm adjustment over the non-missing p-values
#' @keywords internal
#' @noRd
.study_holm <- function(p) {
  ok <- !is.na(p)
  p[ok] <- stats::p.adjust(p[ok], method = "holm")
  p
}

#' Table of one sample over frozen cuts, event-richest band first
#'
#' `ref` (band volumes of the reference and its Kish effective size) adds
#' the band PSI term; its total and the n-adjusted critical value are in the
#' attribute `psi`. The attribute `raw` holds the unweighted events and rows
#' with a known outcome per band, in table order, and `interval` the
#' ascending band index of every row.
#' @keywords internal
#' @noRd
.study_band_table <- function(hist_sample, cuts, ref_shares = NULL, level = 0.95,
                              direction = "higher_is_safer", psi_alpha = 0.05) {
  h <- hist_sample
  B <- length(cuts) + 1L
  idx <- findInterval(h$s, cuts) + 1L
  cnt <- intersect(.study_count_cols, names(h))
  S <- .study_sum(h, idx, B, cnt)
  # event-richest band first: high scores under higher_is_riskier
  high <- identical(direction, "higher_is_riskier")
  o <- if (high) rev(seq_len(B)) else seq_len(B)
  S <- lapply(S, `[`, o)
  n <- S$n; ny <- S$n_y; e <- S$e
  ne <- pmax(ny - e, 0)
  N <- sum(n); NY <- sum(ny); E <- sum(e); NE <- sum(ne)
  rate <- ifelse(ny > 0, e / ny, NA_real_)
  R <- if (NY > 0) E / NY else NA_real_
  # Jeffreys interval on the Kish effective size (the size itself when unweighted)
  neff <- .study_kish(ny, S$w2_y)
  ci <- .study_jeffreys(rate * neff, neff, level)
  cum_ny <- cumsum(ny); cum_e <- cumsum(e)
  cum_rate <- ifelse(cum_ny > 0, cum_e / cum_ny, NA_real_)
  capture <- if (E > 0) cum_e / E else rep(NA_real_, B)
  cum_ne <- if (NE > 0) cumsum(ne) / NE else rep(NA_real_, B)
  bw <- .band_woe(e, ne)
  # odds in the orientation of the scale, 0.5 added to each count (as scr_score_gains())
  odds <- if (high) (e + 0.5) / (ne + 0.5) else (ne + 0.5) / (e + 0.5)
  odds[!(ny > 0)] <- NA_real_
  lab <- .study_labels(cuts)[o]
  tab <- data.table::data.table(
    band = seq_len(B), label = lab, score_lo = c(-Inf, cuts)[o], score_hi = c(cuts, Inf)[o],
    n = n, pct = if (N > 0) n / N else NA_real_, cum_pct = if (N > 0) cumsum(n) / N else NA_real_,
    events = e, rate = rate, rate_lo = ci$lo, rate_hi = ci$hi, cum_rate = cum_rate,
    lift = rate / R, lift_lo = ci$lo / R, lift_hi = ci$hi / R, cum_lift = cum_rate / R,
    capture = capture, cum_nonevent_pct = cum_ne, ks = abs(capture - cum_ne),
    pct_event = bw$pct_event, pct_nonevent = bw$pct_nonevent, woe = bw$log_odds,
    iv = (bw$pct_event - bw$pct_nonevent) * bw$log_odds, odds = odds, log_odds = log(odds))
  ps <- NULL
  tab[, psi := NA_real_]
  if (!is.null(ref_shares)) {
    # shares rescaled to the Kish effective sizes: the smoothing and the
    # n-adjusted critical value then act on effective counts
    nb <- ref_shares$n[o] / sum(ref_shares$n) * ref_shares$n_eff
    nc <- if (N > 0) n / N * .study_kish(N, sum(S$w2)) else n
    ps <- .psi_counts(nb, nc, lab, psi_alpha, c(0.10, 0.25))
    if (!is.null(ps$table)) tab[, psi := ps$table$psi_band]
  }
  # rank order: each band against the previous, event-richer one (unweighted counts)
  er <- S$e_raw; nr <- S$n_y_raw
  p <- c(NA_real_, .study_fisher_gt(er[-1L], nr[-1L], er[-B], nr[-B]))[seq_len(B)]
  tab[, `:=`(p_reversal = p, p_reversal_adj = .study_holm(p))]
  if ("v" %in% cnt) {
    v <- S$v; ve <- S$ve
    tab[, `:=`(value = v, value_events = ve,
               value_capture = if (sum(ve) != 0) cumsum(ve) / sum(ve) else NA_real_,
               value_precision = ifelse(cumsum(v) != 0, cumsum(ve) / cumsum(v), NA_real_))]
  }
  data.table::setattr(tab, "psi", ps)
  data.table::setattr(tab, "raw", list(e = er, n = nr, ep = S$ep, n_y = ny, w2_y = S$w2_y))
  data.table::setattr(tab, "interval", o)
  tab
}

#' Band volumes of the reference, for the PSI of the other samples
#' @keywords internal
#' @noRd
.study_ref_shares <- function(ref_cells, cuts) {
  B <- length(cuts) + 1L
  S <- .study_sum(ref_cells, findInterval(ref_cells$s, cuts) + 1L, B, c("n", "w2"))
  list(n = S$n, n_eff = .study_kish(sum(S$n), sum(S$w2)))
}

# -- discrimination on counts ----------------------------------------------------- #

#' AUC and KS of every column of two count matrices (cells in rows, ascending)
#' @keywords internal
#' @noRd
.study_auc_cols <- function(C1, C0) {
  n1 <- colSums(C1); n0 <- colSums(C0)
  cum0 <- apply(C0, 2L, cumsum)
  cum1 <- apply(C1, 2L, cumsum)
  if (!is.matrix(cum0)) { cum0 <- matrix(cum0, nrow = 1L); cum1 <- matrix(cum1, nrow = 1L) }
  auc <- colSums(C1 * (cum0 - C0 / 2)) / (n1 * n0)
  ks <- apply(abs(sweep(cum1, 2L, n1, "/") - sweep(cum0, 2L, n0, "/")), 2L, max)
  list(auc = auc, ks = ks)
}

#' AUC, KS and Gini from counts, with a stratified bootstrap on the counts
#'
#' `c1` and `c0` are the event and non-event counts per cell in ascending
#' order of the score oriented so that a higher score means more events.
#' Each resample draws the events as Multinomial(`size1`, c1 / sum(c1)) and
#' the non-events as Multinomial(`size0`, c0 / sum(c0)) over the cells: the
#' law of the per-cell counts of a row bootstrap stratified by outcome, at
#' O(K) per resample. `size1` and `size0` are the unweighted class counts
#' (with weights, the cell shares are the weighted ones). The resamples are
#' processed in chunks, so the memory is O(K * chunk). With more than
#' `boot_cells` cells, the resamples run on adjacent cells pooled into
#' `boot_cells` of equal share, and the resampled values are shifted by the
#' difference between the full and the pooled point estimates (exact when
#' K <= `boot_cells`, which covers any points scale). The seed is local to
#' the call, as in scr_metrics(); `keep = TRUE` returns the resampled AUC and
#' KS too.
#' @keywords internal
#' @noRd
.study_auc_boot <- function(c1, c0, n_boot = 0L, level = 0.95, seed = NULL, size1 = sum(c1),
                            size0 = sum(c0), keep = FALSE, boot_cells = 1e4) {
  c1 <- as.double(c1); c0 <- pmax(as.double(c0), 0)
  out <- list(auc = NA_real_, ks = NA_real_, gini = NA_real_, auc_lo = NA_real_, auc_hi = NA_real_,
              ks_lo = NA_real_, ks_hi = NA_real_, gini_lo = NA_real_, gini_hi = NA_real_,
              n_boot = 0L, boot_auc = NULL, boot_ks = NULL)
  if (!length(c1) || !(sum(c1) > 0) || !(sum(c0) > 0)) return(out)
  pt <- .auc_ks_counts(c1, c0)
  out$auc <- pt$auc; out$ks <- pt$ks; out$gini <- pt$gini
  size1 <- round(size1); size0 <- round(size0)
  if (n_boot < 2L || size1 < 1 || size0 < 1) return(out)
  .scr_local_seed(seed)
  shift_auc <- shift_ks <- 0
  if (length(c1) > boot_cells) {
    tot <- c1 + c0
    g <- pmin(pmax(ceiling(cumsum(tot) / sum(tot) * boot_cells), 1), boot_cells)
    c1 <- as.double(rowsum(c1, g, reorder = TRUE)); c0 <- as.double(rowsum(c0, g, reorder = TRUE))
    pp <- .auc_ks_counts(c1, c0)
    shift_auc <- pt$auc - pp$auc; shift_ks <- pt$ks - pp$ks
  }
  K <- length(c1)
  p1 <- c1 / sum(c1); p0 <- c0 / sum(c0)
  chunk <- max(1L, min(as.integer(n_boot), as.integer(floor(1e6 / K))))
  auc <- ks <- numeric(n_boot)
  done <- 0L
  while (done < n_boot) {
    m <- min(chunk, n_boot - done)
    C1 <- stats::rmultinom(m, size1, p1)
    C0 <- stats::rmultinom(m, size0, p0)
    # a single resample goes through the vector kernel: the same values
    # without the column-wise passes over a one-column matrix
    r <- if (m == 1L) .auc_ks_counts(C1, C0) else .study_auc_cols(C1, C0)
    auc[done + seq_len(m)] <- r$auc + shift_auc; ks[done + seq_len(m)] <- r$ks + shift_ks
    done <- done + m
  }
  a <- (1 - level) / 2
  qa <- stats::quantile(auc, c(a, 1 - a), names = FALSE)
  qk <- stats::quantile(ks, c(a, 1 - a), names = FALSE)
  out$auc_lo <- qa[1]; out$auc_hi <- qa[2]
  out$ks_lo <- qk[1]; out$ks_hi <- qk[2]
  out$gini_lo <- 2 * qa[1] - 1; out$gini_hi <- 2 * qa[2] - 1
  out$n_boot <- as.integer(n_boot)
  if (isTRUE(keep)) { out$boot_auc <- auc; out$boot_ks <- ks }
  out
}

#' DeLong standard error of the AUC from counts per cell
#'
#' The structural components of DeLong et al. (1988) are constant within a
#' cell: an event in cell k outranks the non-events below it and half of
#' those tied with it, `V10 = (C0[k-1] + c0[k] / 2) / n0`; a non-event is
#' outranked by the events above it and half of the tied ones,
#' `V01 = (n1 - C1[k] + c1[k] / 2) / n1`. The variances are the
#' count-weighted ones with the `n - 1` denominator, so the result equals
#' the row-level estimator. `size1`/`size0` are the sizes in the
#' denominators (the unweighted class counts under weights).
#' @keywords internal
#' @noRd
.study_delong_counts <- function(c1, c0, size1 = sum(c1), size0 = sum(c0)) {
  c1 <- as.double(c1); c0 <- pmax(as.double(c0), 0)
  n1 <- sum(c1); n0 <- sum(c0)
  if (!(size1 >= 2) || !(size0 >= 2) || !(n1 > 0) || !(n0 > 0)) return(NA_real_)
  cum0 <- cumsum(c0); cum1 <- cumsum(c1)
  v10 <- (cum0 - c0 / 2) / n0
  v01 <- (n1 - cum1 + c1 / 2) / n1
  auc <- sum(c1 * v10) / n1
  var10 <- sum(c1 * (v10 - auc)^2) / n1 * size1 / (size1 - 1)
  var01 <- sum(c0 * (v01 - auc)^2) / n0 * size0 / (size0 - 1)
  sqrt(var10 / size1 + var01 / size0)
}

#' Discrimination of one sample from its cells, oriented by the direction
#' @keywords internal
#' @noRd
.study_discrimination <- function(cells, direction, n_boot, level, keep = FALSE, boot_cells = 1e4) {
  c1 <- cells$e; c0 <- pmax(cells$n_y - cells$e, 0)
  # cells are in ascending score; a higher score means more events only under higher_is_riskier
  if (!identical(direction, "higher_is_riskier")) { c1 <- rev(c1); c0 <- rev(c0) }
  .study_auc_boot(c1, c0, n_boot, level, seed = NULL, size1 = sum(cells$e_raw),
                  size0 = sum(cells$n_y_raw - cells$e_raw), keep = keep, boot_cells = boot_cells)
}

# -- traffic lights ------------------------------------------------------------- #

#' Light from a value and its interval against two thresholds
#'
#' The deviation convention, as for the p-value lights: amber or red only
#' when the interval shows the metric beyond a threshold. With
#' `higher_better`, green when the upper bound reaches `green`, red when it
#' stays below `red`, amber in between; the mirror image on the lower bound
#' when lower is better. Without an interval the value is used for both
#' bounds. A missing value is "grey".
#' @keywords internal
#' @noRd
.rag_light <- function(value, lo = value, hi = value, green, red, higher_better = TRUE) {
  lo <- ifelse(is.na(lo), value, lo); hi <- ifelse(is.na(hi), value, hi)
  out <- if (isTRUE(higher_better)) {
    ifelse(hi >= green, "green", ifelse(hi < red, "red", "amber"))
  } else {
    ifelse(lo <= green, "green", ifelse(lo > red, "red", "amber"))
  }
  out[is.na(value)] <- "grey"
  out
}

#' Worst light: red, then amber, then green; "grey" only when nothing else is lit
#'
#' A "grey" light (no testable result) never rolls up to green, and lights
#' marked `"none"` (reported, not lit) are ignored.
#' @keywords internal
#' @noRd
.rag_worst <- function(lights) {
  l <- lights[!is.na(lights)]
  if (any(l == "red")) "red" else if (any(l == "amber")) "amber" else if (any(l == "green")) "green" else "grey"
}

# -- input normalization -------------------------------------------------------- #

#' Normalize the input of a score study into one count table
#'
#' Dispatches on a scorecard (its scored samples, direction and objective)
#' or on a data.frame (row data or pre-aggregated counts). Returns
#' `hist` (keyed by `sample`, optionally `group`, and `s`), `reference`,
#' `samples` (the reference first), `objective`, `direction`, `target`
#' and `meta` (`weighted`, `has_value`, `has_prob`, `quantized`).
#' @keywords internal
#' @noRd
.study_input <- function(x, ...) {
  if (inherits(x, "scr_scorecard")) return(.study_input_sc(x, ...))
  if (is.data.frame(x)) return(.study_input_df(x, ...))
  stop("a score study needs a scorecard from scr_scorecard() or a data.frame.", call. = FALSE)
}

#' @keywords internal
#' @noRd
.study_input_sc <- function(x, sample = "holdout", reference = "train", by = NULL, max_cells = 1e5,
                            breaks = NULL, fn = "score study") {
  nms <- names(x$samples)
  .study_chr(sample, "sample", fn); .study_chr1(reference, "reference", fn)
  bad <- setdiff(c(reference, sample), nms)
  if (length(bad)) stop(fn, "(): sample(s) ", lst(bad), " not in the scorecard (", lst(nms), ").", call. = FALSE)
  samples <- unique(c(reference, sample))
  s <- lapply(samples, function(nm) x$samples[[nm]])
  score <- unlist(lapply(s, `[[`, "score"), use.names = FALSE)
  y <- unlist(lapply(s, `[[`, "y"), use.names = FALSE)
  lab <- rep(samples, vapply(s, nrow, integer(1)))
  keys <- list(sample = lab)
  if (!is.null(by)) {
    .study_chr1(by, "by", fn)
    miss <- samples[!vapply(s, function(d) by %in% names(d), logical(1))]
    if (length(miss)) stop(fn, "(): column '", by, "' not in the scored sample(s) ", lst(miss), ".", call. = FALSE)
    keys$group <- unlist(lapply(s, function(d) as.character(d[[by]])), use.names = FALSE)
  }
  # expected event probability from the alignment of the scale
  prob <- .score_to_prob(x$alignment, score)
  h <- .study_hist(score, y, prob = prob, by = keys, max_cells = max_cells, breaks = breaks, fn = fn)
  list(hist = h, reference = reference, samples = samples,
       objective = x$config$objective %||% "risk", direction = x$direction, target = x$target,
       group_levels = if (!is.null(by)) .study_num_levels(lapply(s, `[[`, by)),
       sql_table = x$config$sql_table, sql_dialect = x$config$sql_dialect,
       meta = list(weighted = FALSE, has_value = FALSE, has_prob = TRUE, quantized = !is.null(attr(h, "edges"))))
}

#' @keywords internal
#' @noRd
.study_input_df <- function(x, score = "score", y = "y", objective = "risk", direction = NULL, weight = NULL,
                            value = NULL, value_events = NULL, prob = NULL, sample = NULL, reference = NULL,
                            study = NULL, counts = FALSE, n = "n", events = "events", by = NULL,
                            max_cells = 1e5, breaks = NULL, fn = "score study") {
  if (!is.character(objective) || length(objective) != 1L || !objective %in% c("risk", "propensity")) {
    stop(fn, "(): `objective` must be \"risk\" or \"propensity\".", call. = FALSE)
  }
  if (!is.null(direction) && (!is.character(direction) || length(direction) != 1L ||
                              !direction %in% c("higher_is_safer", "higher_is_riskier"))) {
    stop(fn, "(): `direction` must be NULL, \"higher_is_safer\" or \"higher_is_riskier\".", call. = FALSE)
  }
  direction <- resolve_direction(list(objective = objective, direction = direction))
  if (!is.logical(counts) || length(counts) != 1L || is.na(counts)) stop(fn, "(): `counts` must be TRUE or FALSE.", call. = FALSE)
  need <- c(score, if (counts) c(n, events) else y, weight, value, if (counts) value_events, prob, sample, by)
  for (nm in c("score", "y", "weight", "value", "value_events", "prob", "sample", "by", "n", "events")) {
    v <- get(nm)
    if (!is.null(v)) .study_chr1(v, nm, fn)
  }
  miss <- setdiff(need, names(x))
  if (length(miss)) stop(fn, "(): column(s) ", lst(miss), " not in `x`.", call. = FALSE)
  if (counts && !is.null(weight)) stop(fn, "(): `weight` does not apply to pre-aggregated counts.", call. = FALSE)
  for (nm in setdiff(need, c(sample, by, if (!counts) y))) {
    if (!is.numeric(x[[nm]])) stop(fn, "(): column '", nm, "' must be numeric.", call. = FALSE)
  }
  # sample labels: the reference is the first level unless given
  if (is.null(sample)) {
    lab <- rep("all", nrow(x))
    lv <- "all"
    if (!is.null(reference) || !is.null(study)) stop(fn, "(): `reference` and `study` need a `sample` column.", call. = FALSE)
  } else {
    sv <- x[[sample]]
    if (anyNA(sv)) stop(fn, "(): the sample column '", sample, "' has missing values.", call. = FALSE)
    lv <- .study_levels(sv)
    lab <- as.character(sv)
  }
  if (!is.null(reference)) .study_chr1(reference, "reference", fn)
  if (!is.null(study)) .study_chr(study, "study", fn)
  reference <- reference %||% lv[1]
  study <- study %||% setdiff(lv, reference)
  if (!length(study)) study <- reference
  bad <- setdiff(c(reference, study), lv)
  if (length(bad)) stop(fn, "(): sample(s) ", lst(bad), " not in the column '", sample, "'.", call. = FALSE)
  samples <- unique(c(reference, study))
  # rows of the samples in the study (no copy when every row is in)
  keep <- lab %in% samples
  all_in <- all(keep)
  sub <- function(v) if (all_in) v else v[keep]
  keys <- list(sample = sub(lab))
  if (!is.null(by)) keys$group <- sub(as.character(x[[by]]))
  col <- function(nm) if (is.null(nm)) NULL else sub(x[[nm]])
  h <- if (counts) {
    .study_hist_counts(col(score), col(n), col(events), value = col(value), value_events = col(value_events),
                       by = keys, max_cells = max_cells, breaks = breaks, fn = fn)
  } else {
    .study_hist(col(score), col(y), w = col(weight), value = col(value), prob = col(prob), by = keys,
                max_cells = max_cells, breaks = breaks, fn = fn)
  }
  list(hist = h, reference = reference, samples = samples, objective = objective, direction = direction,
       target = if (counts) events else y, group_levels = if (!is.null(by)) .study_num_levels(list(x[[by]])),
       meta = list(weighted = !is.null(weight), has_value = !is.null(value), has_prob = !is.null(prob),
                   quantized = !is.null(attr(h, "edges"))))
}

#' Labels of a sample column, in their reading order
#'
#' The levels of a factor; numbers in numeric order (as text, "10" would
#' sort before "9"); anything else as sorted text.
#' @keywords internal
#' @noRd
.study_levels <- function(v) {
  if (is.factor(v)) return(levels(droplevels(v)))
  if (is.numeric(v)) return(unique(as.character(sort(unique(v)))))
  sort(unique(as.character(v)))
}

#' Labels of a numeric group column in numeric order; NULL when not numeric
#'
#' `vals` is a list of vectors (the column in every sample). A column that
#' is not numeric keeps the order of its labels.
#' @keywords internal
#' @noRd
.study_num_levels <- function(vals) {
  num <- vapply(vals, function(v) is.numeric(v) && !is.factor(v), logical(1))
  if (!length(vals) || !all(num)) return(NULL)
  unique(as.character(sort(unique(unlist(vals, use.names = FALSE)))))
}

#' @keywords internal
#' @noRd
.study_chr1 <- function(x, name, fn) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(fn, "(): `", name, "` must be a single non-empty string.", call. = FALSE)
  }
  invisible(x)
}

#' @keywords internal
#' @noRd
.study_chr <- function(x, name, fn) {
  if (!is.character(x) || !length(x) || anyNA(x) || any(!nzchar(x))) {
    stop(fn, "(): `", name, "` must be a character vector of sample names.", call. = FALSE)
  }
  invisible(x)
}

#' Arguments caught by `...` that no method uses
#' @keywords internal
#' @noRd
.study_dots <- function(dots, fn) {
  if (length(dots)) {
    nm <- names(dots); nm[is.na(nm) | !nzchar(nm)] <- "(unnamed)"
    stop(fn, "(): unused argument(s): ", lst(nm), ".", call. = FALSE)
  }
  invisible(NULL)
}

#' @keywords internal
#' @noRd
.study_whole <- function(x, name, fn, lower = 0) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x != round(x) || x < lower) {
    stop(fn, "(): `", name, "` must be a whole number >= ", lower, ".", call. = FALSE)
  }
  as.integer(x)
}

#' A finite number in a range, with the calling function named in the error
#' @keywords internal
#' @noRd
.study_num1 <- function(x, name, fn, lower = -Inf, upper = Inf, open_lower = FALSE) {
  tryCatch(.scr_num1(x, name, lower = lower, upper = upper, open_lower = open_lower),
           error = function(e) stop(fn, "(): ", conditionMessage(e), call. = FALSE))
}

#' Cells of the bootstrap: a number of at least 2, `Inf` for no pooling
#' @keywords internal
#' @noRd
.study_boot_cells <- function(x, fn) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 2 || (is.finite(x) && x != round(x))) {
    stop(fn, "(): `boot_cells` must be a whole number >= 2, or Inf.", call. = FALSE)
  }
  as.double(x)
}

#' @keywords internal
#' @noRd
.study_level <- function(level, fn) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    stop(fn, "(): `level` must be a number in (0, 1).", call. = FALSE)
  }
  level
}

# -- one score, cuts and shares --------------------------------------------------- #

#' Objective and direction of one score, checked
#'
#' `objective` and `direction` are the arguments as given (`NULL` when not
#' given). A score study passed as `cuts` sets both, and a given argument
#' that disagrees with it is an error; without a study the objective falls
#' back on `default` and the direction follows the objective.
#' @keywords internal
#' @noRd
.study_one_direction <- function(objective, direction, cuts, default, fn) {
  if (!is.null(objective) && (!is.character(objective) || length(objective) != 1L ||
                              !objective %in% c("risk", "propensity"))) {
    stop(fn, "(): `objective` must be \"risk\" or \"propensity\".", call. = FALSE)
  }
  if (!is.null(direction) && (!is.character(direction) || length(direction) != 1L ||
                              !direction %in% c("higher_is_safer", "higher_is_riskier"))) {
    stop(fn, "(): `direction` must be NULL, \"higher_is_safer\" or \"higher_is_riskier\".", call. = FALSE)
  }
  if (inherits(cuts, "scr_study")) {
    for (nm in c("objective", "direction")) {
      given <- get(nm)
      if (!is.null(given) && !identical(given, cuts[[nm]])) {
        stop(fn, "(): `", nm, "` is \"", given, "\" but the study given as `cuts` was fitted with \"",
             cuts[[nm]], "\". Drop the argument, or pass a study of that ", nm, ".", call. = FALSE)
      }
    }
    return(list(objective = cuts$objective, direction = cuts$direction))
  }
  objective <- objective %||% default
  list(objective = objective, direction = resolve_direction(list(objective = objective, direction = direction)))
}

#' Tie-safe cut that selects a share of the volume from the event-rich end
#'
#' The boundary between two distinct scores nearest to the share. A share of
#' 1 selects every row (`-Inf` from the high scores, `Inf` from the low
#' ones); when every score is tied there is no boundary, and the nearer of
#' "no row" and "every row" is taken.
#' @keywords internal
#' @noRd
.study_share_cut <- function(h, shares, side) {
  none <- if (side == "high") Inf else -Inf
  every <- -none
  vapply(shares, function(d) {
    if (d >= 1) return(every)
    cs <- .study_cuts(h, 1L, "tail", d, side)$cuts
    if (length(cs)) cs else if (d > 0.5) every else none
  }, numeric(1))
}

#' Share of the volume on the alert side of each cut, from the score cells
#'
#' The alert side is `score >= cut` from the high scores and `score < cut`
#' from the low ones.
#' @keywords internal
#' @noRd
.study_alert_share <- function(cells, cuts, side) {
  N <- sum(cells$n)
  if (!(N > 0)) return(rep(NA_real_, length(cuts)))
  # volume strictly below each cut: the cells are in ascending score
  below <- c(0, cumsum(cells$n))[findInterval(cuts, cells$s_lo, left.open = TRUE) + 1L]
  if (side == "high") (N - below) / N else below / N
}

# -- production methods ---------------------------------------------------------- #

# Dialects of the production SQL, in the order of OptimalBinningWoE::obwoe_sql()
.sql_dialects <- c("ansi", "postgres", "mysql", "mariadb", "sqlserver", "oracle", "spark", "hive", "databricks",
                   "bigquery", "snowflake", "redshift", "duckdb", "sqlite")

#' @param score For `scr_study`: name of the score column of `newdata`.
#'   `newdata` may also be a numeric vector of scores.
#' @param numbered For a tiers study: `TRUE` (default) returns the tier
#'   labels with their order in front (`"01.very high"`), `FALSE` the plain
#'   labels. Band labels are intervals and never get a prefix.
#' @section Score studies:
#'
#' For a score study ([scr_bands()], [scr_tiers()]), `newdata` is returned
#' (as a copy) with `tier`, the band or tier number, and `tier_label`. The
#' intervals are left-closed: `score >= cut` is the upper side, and a
#' missing score gives a missing tier.
#'
#' The labels of a tiers study carry their order, `"01."` for the tier with
#' the highest event rate (the first row of the tiers table) down to the
#' tier with the lowest, so they sort from the event-richest tier under any
#' objective and direction; `tier` is unchanged and still rises with the
#' event rate. The result joins to the `tier_label` column of the tiers
#' table. `numbered = FALSE` returns the plain labels (its `label` column).
#' @rdname scr_apply
#' @export
scr_apply.scr_study <- function(x, newdata, score = "score", numbered = TRUE, ...) {
  .study_dots(list(...), "scr_apply")
  labs <- .study_code_labels(x, numbered, "scr_apply")
  if (is.numeric(newdata) && is.null(dim(newdata))) {
    newdata <- data.table::data.table(score = newdata)
    score <- "score"
  }
  if (!is.data.frame(newdata)) stop("scr_apply(): `newdata` must be a data.frame or a numeric vector of scores.", call. = FALSE)
  .study_chr1(score, "score", "scr_apply")
  if (!score %in% names(newdata)) stop("scr_apply(): column '", score, "' not in `newdata`.", call. = FALSE)
  if (!is.numeric(newdata[[score]])) stop("scr_apply(): column '", score, "' must be numeric.", call. = FALSE)
  # a copy: the caller's table is never modified by reference
  out <- data.table::copy(data.table::as.data.table(newdata))
  i <- findInterval(as.double(out[[score]]), x$cuts) + 1L
  data.table::set(out, j = "tier", value = x$codes[i])
  data.table::set(out, j = "tier_label", value = labs[i])
  out[]
}

#' Label of every ascending interval, as the production methods assign it
#'
#' Tier labels with their order in front when `numbered`; the plain labels
#' otherwise, and always for a bands study (interval labels). The prefix is
#' rebuilt from the tier labels, so a study fitted before the numbered
#' labels existed gets them too.
#' @keywords internal
#' @noRd
.study_code_labels <- function(x, numbered, fn) {
  if (!is.logical(numbered) || length(numbered) != 1L || is.na(numbered)) {
    stop(fn, "(): `numbered` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!numbered || !inherits(x, "scr_study_tiers")) return(x$code_labels)
  # `labels` are in tier order and `codes` gives the tier of every interval
  .tier_numbered(x$labels)[x$codes]
}

#' @section Score studies:
#'
#' For a score study ([scr_bands()], [scr_tiers()]), the SQL reads the score
#' column of `table` and adds `tier` and `tier_label` with a `CASE` on the
#' frozen cuts (`score >= cut` is the upper side; a `NULL` score gives a
#' `NULL` tier). `table` and `dialect` default to the configuration of the
#' scorecard the study came from, else to `"your_table"` and `"ansi"`. The
#' tiers computed by the SQL match [scr_apply()], by an automated test.
#'
#' The labels of a tiers study carry their order, `'01.very high'` for the
#' tier with the highest event rate down to the tier with the lowest, as in
#' [scr_apply()], so `ORDER BY tier_label` lists the event-richest tier
#' first; `numbered = FALSE` emits the plain labels. Band labels are
#' intervals and never get a prefix.
#' @param score For `scr_study`: name of the score column of `table`.
#' @param numbered For a tiers study: `TRUE` (default) emits the tier
#'   labels with their order in front, `FALSE` the plain labels.
#' @rdname scr_sql
#' @export
scr_sql.scr_study <- function(x, table = NULL, dialect = NULL, file = NULL, score = "score", numbered = TRUE,
                              ...) {
  .study_dots(list(...), "scr_sql")
  labs <- .study_code_labels(x, numbered, "scr_sql")
  table <- table %||% x$sql_table %||% "your_table"
  .study_chr1(table, "table", "scr_sql"); .study_chr1(score, "score", "scr_sql")
  # the dialects of scr_sql() for a scorecard, with the same error
  dialect <- match.arg(dialect %||% x$sql_dialect %||% "ansi", .sql_dialects)
  qs <- paste0("s.", .sql_q(score, dialect))
  cuts <- x$cuts; B <- length(cuts) + 1L
  # left-closed intervals: `score >= cut` is the upper side; a NULL score gives a NULL tier
  whens <- function(vals) paste(c(sprintf("WHEN %s IS NULL THEN NULL", qs),
                                  if (B > 1L) sprintf("WHEN %s < %s THEN %s", qs, .sql_num(cuts), vals[-B]),
                                  sprintf("ELSE %s", vals[B])), collapse = " ")
  kind <- if (inherits(x, "scr_study_tiers")) "tiers" else "bands"
  out <- c("-- =============================================================",
           sprintf("-- scorecraft | score study (%s) | target: %s | %d %s | dialect: %s", kind,
                   .sql_cmt(x$target %||% "?"), B, kind, dialect),
           sprintf("-- Generated on %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
           sprintf("-- Cuts frozen on the reference sample '%s'; %s >= cut is the upper side.",
                   .sql_cmt(x$reference), .sql_cmt(score)),
           "-- =============================================================", "",
           "SELECT", "    s.*,",
           sprintf("    CASE %s END AS tier,", whens(as.character(x$codes))),
           sprintf("    CASE %s END AS tier_label", whens(.sql_str(labs, dialect))),
           sprintf("FROM %s s;", table))
  .sql_out(.sql_lines(out), file)
}

#' @param rag For `scr_study`: an optional [scr_rag()] object whose lights
#'   are written to the same workbook.
#' @rdname scr_export
#' @export
scr_export.scr_study <- function(x, dir, stamp = TRUE, rag = NULL, ...) {
  .study_dots(list(...), "scr_export")
  .need_openxlsx()
  if (!is.null(rag) && !inherits(rag, "scr_rag")) stop("scr_export(): `rag` must come from scr_rag().", call. = FALSE)
  out_dir <- .export_dir(dir, stamp)
  tag <- .file_tag(x$target)
  tiers <- inherits(x, "scr_study_tiers")
  kind <- if (tiers) "tiers" else "bands"
  settings <- list(kind = kind, target = x$target, objective = x$objective, direction = x$direction,
                   reference = x$reference, samples = x$samples, level = x$level)
  settings <- c(settings, if (tiers) list(measure = x$measure, method = x$method, criterion = x$criterion,
                                          n_tiers_requested = x$n_tiers_requested, n_tiers = length(x$cuts) + 1L,
                                          alpha = x$alpha, min_pct = x$min_pct, min_events = x$min_events)
                else list(spacing = x$spacing, n_bands_requested = x$n_bands_requested,
                          n_bands_effective = x$n_bands_effective))
  cuts <- data.frame(cut = seq_along(x$cuts), score = x$cuts)
  if (tiers) cuts$score_raw <- if (length(x$cuts_raw) == length(x$cuts)) x$cuts_raw else NA_real_
  sheets <- list(Summary = x$summary, Table = x$table, Cuts = cuts, Settings = .kv_table(settings))
  names(sheets)[2] <- if (tiers) "Tiers" else "Bands"
  if (tiers) {
    sheets$Ledger <- x$ledger
    if (!is.null(x$stability)) sheets$Stability <- x$stability$cuts
  }
  if (!is.null(rag)) sheets <- c(sheets, .rag_sheets(rag))
  sheets <- lapply(sheets, .study_sheet)
  files <- list(xlsx = .scr_write_xlsx(sheets, file.path(out_dir, sprintf("study_%s_%s.xlsx", kind, tag))))
  for (f in files) msg("  %s", f)
  x$files <- files
  invisible(x)
}

#' Infinite numbers become missing cells (the label column keeps the open ends)
#' @keywords internal
#' @noRd
.study_sheet <- function(d) {
  d <- as.data.frame(d, stringsAsFactors = FALSE)
  for (j in seq_along(d)) if (is.double(d[[j]])) d[[j]][is.infinite(d[[j]])] <- NA_real_
  d
}

# -- shared print helpers --------------------------------------------------------- #

#' @keywords internal
#' @noRd
.study_f <- function(x, fmt, na = "-") ifelse(is.na(x), na, sprintf(fmt, x))

#' A volume with thousands separators; weighted totals can exceed the integer range
#' @keywords internal
#' @noRd
.study_n <- function(x) formatC(round(x), big.mark = ",", format = "f", digits = 0)

# data.table column names used without quotes in this file
utils::globalVariables(c(
  "s", "s_lo", "s_hi", "g", "n", "n_y", "e", "w2", "w2_y", "n_raw", "n_y_raw", "e_raw", "psi",
  "k", "yr", "w", "wk", "we", "wsq", "wsqk", "v", "vy", "p", ".SD", ".N"
))
