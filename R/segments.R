# ============================================================================ #
# segments.R - one score read on many segments
# ============================================================================ #
# One grouped pass counts the rows per segment (and per group of `by`) and
# score cell. Discrimination, the indirect standardization, the stability
# index and the logistic slope of every segment are functions of those
# counts: the slope is a weighted fit on the cells, never on the rows.
# ============================================================================ #

#' One score on many segments
#'
#' Reads one score on the segments of a population (a product, a channel, a
#' region) and says, for each, whether the score ranks and calibrates there
#' as it does on the whole: discrimination with its standard error, the
#' observed events against those expected from the pooled bands, the offset
#' on the log-odds scale, the slope of the score relative to the pooled one
#' and a suggested action.
#'
#' @section Statistics:
#'
#' The rows are counted once per segment and distinct score. Per segment:
#' \itemize{
#'   \item `auc`, `gini` and `ks`, with the DeLong standard error `auc_se`
#'     computed from the counts and the interval `auc +/- z * auc_se` (cut to
#'     \[0, 1\]); with `n_boot > 0`, the interval is the percentile interval
#'     of the count bootstrap of [scr_bands()]. A segment without a positive
#'     standard error (fewer than two events or non-events, or every score
#'     tied) has no interval.
#'   \item `psi`: the population stability index of the score distribution
#'     of the segment against the pooled one, over `n_bands` equal-share
#'     bands frozen on the pooled rows (see [scr_psi()]). The segment is part
#'     of the pooled rows, so the two distributions are not independent:
#'     `psi_critical`, the n-adjusted critical value at `1 - level`, uses
#'     \eqn{(1/n_s - 1/N)\,\chi^2_{B-1}} in place of
#'     \eqn{(1/n_s + 1/N)\,\chi^2_{B-1}}, with \eqn{n_s} and \eqn{N} the
#'     sizes of the segment and of the pooled rows (under weights, the Kish
#'     effective sizes, with the covariance term scaled by the weight share
#'     of the segment). It is `NA` when the segment is the whole pool.
#'   \item `expected`: the events expected by indirect standardization,
#'     \eqn{\sum_b n_{s,b} R_b}, with \eqn{n_{s,b}} the rows of the segment
#'     in the pooled band \eqn{b} and \eqn{R_b} the pooled event rate of the
#'     band; `oe_ratio` = events / `expected`, with the Jeffreys interval of
#'     the observed rate (on the Kish effective size under weights) divided
#'     by the expected rate.
#'   \item `offset` = logit(observed rate) - logit(expected rate): what a
#'     segment intercept would add to the log-odds. `NA` when either rate is
#'     0 or 1.
#'   \item `slope`: the coefficient of the score in a logistic regression of
#'     the outcome within the segment, fitted on the count table with the
#'     counts as weights (equal to the fit on the rows); `slope_ratio` is
#'     `slope` over the pooled slope. `NA` when the segment has a single
#'     score value or its classes are separated by the score. `p_slope` is
#'     the two-sided Wald test of the slope of the segment against the slope
#'     fitted on the rest of its group (the pooled rows without the segment,
#'     an independent sample), \eqn{z = (b_s - b_r) / \sqrt{v_s + v_r}} with
#'     the model-based variances of the two fits (scaled to the Kish
#'     effective size under weights); `p_slope_adj` is its Holm adjustment
#'     across the segments.
#' }
#' Rows with a missing outcome count in the volume and in the PSI only.
#' With more than `max_cells` distinct scores the scores are pooled into
#' cells, and the slope is fitted on the cell midpoints.
#'
#' @section Test of equal discrimination:
#'
#' With \eqn{w_s = 1 / se_s^2}, the inverse-variance weighted mean is
#' \eqn{AUC_w = \sum_s w_s AUC_s / \sum_s w_s}, and
#' \deqn{Q = \sum_s (AUC_s - AUC_w)^2 / se_s^2}
#' is compared with a chi-square law with \eqn{S - 1} degrees of freedom
#' (`test`: `statistic`, `df`, `p_value`). The segments are disjoint, so
#' their AUCs are independent. Per segment, `auc_diff` =
#' \eqn{AUC_s - AUC_w} is tested by a two-sided z test with the variance
#' \eqn{se_s^2 - 1 / \sum_s w_s} (the mean contains the segment); `p_auc` is
#' its p-value and `p_auc_adj` the Holm adjustment across the segments.
#' Segments without a positive standard error (fewer than two events or
#' non-events, or every score tied) are left out of the test. With exactly
#' two segments the two tests are one and the same, so no adjustment is
#' made; the same holds for the slope tests.
#'
#' @section Action:
#'
#' A convention of this package, read in this order:
#' \describe{
#'   \item{`"too few events"`}{Fewer than `min_events` events, or
#'     non-events, in the segment (unweighted counts).}
#'   \item{`"separate model"`}{The score ranks differently: `abs(auc_diff) >
#'     auc_tol` with `p_auc_adj < 1 - level`, or `abs(slope_ratio - 1) >
#'     slope_tol` with `p_slope_adj < 1 - level`. A difference must be both
#'     material and significant: on a small segment, a slope ratio far from
#'     1 is often noise.}
#'   \item{`"offset"`}{The ranking holds but the level does not:
#'     `abs(offset) > offset_tol` and the interval of `oe_ratio` excludes 1.}
#'   \item{`"shared"`}{None of the above: the pooled score serves the
#'     segment.}
#' }
#' The tolerances are starting points, not rules; set them to the policy in
#' force.
#'
#' @section Groups:
#'
#' With `by` (a period, a sample label), the analysis is repeated within
#' each group: the pooled reference of a segment is the whole of its group.
#' The bands are frozen once, on all rows. Groups and segments are listed in
#' the order of their labels; numeric columns in numeric order, and a factor
#' segment column in the order of its levels.
#'
#' @inheritParams scr_bands
#' @param x A `data.frame` with one row per scored case, or an object from
#'   [scr_scorecard()].
#' @param newdata For a scorecard: the rows to score with [scr_apply()],
#'   holding the candidate variables, the target and the segment column.
#' @param segment Name of the segment column. A missing segment is the
#'   segment `"(missing)"`. The segments of a factor column are listed in
#'   the order of its levels.
#' @param by Optional name of a column of periods or groups.
#' @param target For a scorecard: the outcome column of `newdata`; `NULL`
#'   uses the target of the scorecard.
#' @param n_bands Equal-share bands frozen on the pooled rows, for the
#'   expected events and the PSI.
#' @param level Confidence level of the intervals; `1 - level` is the
#'   significance of the tests. For a scorecard, `NULL` uses
#'   `config$study_level` (0.95).
#' @param min_events Fewest events, and non-events, for an action other than
#'   `"too few events"`.
#' @param n_boot Bootstrap resamples of the AUC interval; `0` (default)
#'   keeps the DeLong interval.
#' @param auc_tol,offset_tol,slope_tol Tolerances of the actions: on the
#'   difference of the AUC from the weighted mean, on the absolute offset
#'   and on the distance of the slope ratio from 1.
#'
#' @return An object of class `c("scr_segments", "list")`:
#'   \describe{
#'     \item{`table`}{One row per group and segment: `group`, `segment`,
#'       `n`, `events`, `rate`, `auc`, `auc_se`, `auc_lo`, `auc_hi`, `gini`,
#'       `ks`, `psi`, `psi_critical`, `expected`, `oe_ratio`, `oe_lo`,
#'       `oe_hi`, `offset`, `slope`, `slope_ratio`, `auc_diff`, `p_auc`,
#'       `p_auc_adj`, `p_slope`, `p_slope_adj` and `action`.}
#'     \item{`test`}{One row per group: `group`, `segments` (those in the
#'       test), `auc_w`, `statistic`, `df` and `p_value`.}
#'     \item{`pooled`}{One row per group, the pooled rows: `group`, `n`,
#'       `events`, `rate`, `auc`, `auc_se`, `auc_lo`, `auc_hi`, `gini`, `ks`
#'       and `slope`.}
#'     \item{`cuts`, `segment`, `by`, `level`, `min_events`, `auc_tol`,
#'       `offset_tol`, `slope_tol`, `objective`, `direction`, `target`,
#'       `call`}{The pooled cuts and the settings.}
#'   }
#'
#' @references
#' Brown, L. D., Cai, T. T. and DasGupta, A. (2001). Interval estimation for
#' a binomial proportion. *Statistical Science*, 16(2), 101-133.
#' \doi{10.1214/ss/1009213286}
#'
#' DeLong, E. R., DeLong, D. M. and Clarke-Pearson, D. L. (1988). Comparing
#' the areas under two or more correlated receiver operating characteristic
#' curves: a nonparametric approach. *Biometrics*, 44(3), 837-845.
#'
#' Thomas, L. C., Crook, J. and Edelman, D. (2017). *Credit Scoring and Its
#' Applications*, 2nd edition. SIAM. \doi{10.1137/1.9781611974560}
#'
#' @seealso [scr_mix_shift()] for the change of the event rate between
#'   samples, [scr_rag()] for lights per segment against a reference.
#' @family score-studies
#' @examples
#' local({
#'   set.seed(1)
#'   n <- 6000
#'   seg <- sample(c("app", "store", "web"), n, TRUE, c(0.5, 0.3, 0.2))
#'   x <- rnorm(n)
#'   # the score ranks the same everywhere; "web" defaults more at every score
#'   d <- data.frame(channel = seg, score = round(600 + 50 * x),
#'                   y = rbinom(n, 1, plogis(-2 - x + 0.6 * (seg == "web"))))
#'   sg <- scr_segments(d, segment = "channel")
#'   print(sg)
#'   sg$test
#' })
#' @export
scr_segments <- function(x, ...) UseMethod("scr_segments")

#' @rdname scr_segments
#' @export
scr_segments.data.frame <- function(x, segment, by = NULL, n_bands = 10L, level = 0.95, min_events = 20L,
                                    n_boot = 0L, auc_tol = 0.03, offset_tol = 0.25, slope_tol = 0.25,
                                    score = "score", y = "y", objective = "risk", direction = NULL, weight = NULL,
                                    seed = NULL, max_cells = 1e5, ...) {
  fn <- "scr_segments"
  .study_dots(list(...), fn)
  if (missing(segment)) stop(fn, "(): `segment` is needed: the name of the segment column.", call. = FALSE)
  for (nm in c("segment", "by", "score", "y", "weight")) {
    v <- get(nm)
    if (!is.null(v)) .study_chr1(v, nm, fn)
  }
  rd <- .study_one_direction(objective, direction, NULL, "risk", fn)
  miss <- setdiff(c(segment, by, score, y, weight), names(x))
  if (length(miss)) stop(fn, "(): column(s) ", lst(miss), " not in `x`.", call. = FALSE)
  for (nm in c(score, weight)) {
    if (!is.numeric(x[[nm]])) stop(fn, "(): column '", nm, "' must be numeric.", call. = FALSE)
  }
  .segments_fit(x[[score]], x[[y]], x[[segment]], if (!is.null(by)) x[[by]], if (!is.null(weight)) x[[weight]],
                rd$objective, rd$direction, n_bands, level, min_events, n_boot, seed, auc_tol, offset_tol,
                slope_tol, max_cells, y, segment, by, sys.call())
}

#' @rdname scr_segments
#' @export
scr_segments.scr_scorecard <- function(x, newdata, segment, by = NULL, n_bands = 10L, level = NULL,
                                       min_events = 20L, n_boot = 0L, auc_tol = 0.03, offset_tol = 0.25,
                                       slope_tol = 0.25, target = NULL, seed = NULL, max_cells = 1e5, ...) {
  fn <- "scr_segments"
  .study_dots(list(...), fn)
  if (missing(newdata) || !is.data.frame(newdata)) {
    stop(fn, "(): a scorecard needs `newdata`, the rows to score.", call. = FALSE)
  }
  if (missing(segment)) stop(fn, "(): `segment` is needed: the name of the segment column.", call. = FALSE)
  target <- target %||% x$target
  for (nm in c("segment", "by", "target")) {
    v <- get(nm)
    if (!is.null(v)) .study_chr1(v, nm, fn)
  }
  miss <- setdiff(c(segment, by, target), names(newdata))
  if (length(miss)) stop(fn, "(): column(s) ", lst(miss), " not in `newdata`.", call. = FALSE)
  # the event of the scorecard: a text target by its label, an inverted 0/1 target by 0
  yv <- newdata[[target]]
  lvl <- if (is.factor(yv) || is.character(yv)) x$event$label else if (isTRUE(x$event$inverted)) 0L else NULL
  yy <- .target_as_int(yv, target, lvl)$y
  score <- scr_apply(x, newdata)$score
  .segments_fit(score, yy, newdata[[segment]], if (!is.null(by)) newdata[[by]], NULL,
                x$config$objective %||% "risk", x$direction, n_bands, level %||% x$config$study_level %||% 0.95,
                min_events, n_boot, seed %||% x$config$seed, auc_tol, offset_tol, slope_tol, max_cells, target,
                segment, by, sys.call())
}

#' Slope of the score in a logistic fit on a count table, with its variance
#'
#' A weighted fit on the cells (event share of each cell, its volume with a
#' known outcome as the weight), equal to the fit on the rows. The score is
#' standardized by `mu` and `sdv` for the fit; the slope and its variance
#' are returned per unit of the raw score. The variance is the model-based
#' one, scaled to the Kish effective size under weights (`sum(w^2) /
#' sum(w)`, 1 without weights). `NA` with fewer than two populated cells, a
#' single class, or classes separated by the score (no finite estimate).
#' `cells` needs `s_lo`, `s_hi`, `n_y`, `e` and `w2_y`.
#' @keywords internal
#' @noRd
.seg_slope <- function(cells, mu, sdv) {
  na <- list(slope = NA_real_, var = NA_real_)
  ok <- cells$n_y > 0
  ny <- cells$n_y[ok]; e <- pmin(pmax(cells$e[ok], 0), ny)
  if (length(ny) < 2L || !(sum(e) > 0) || !(sum(ny - e) > 0) || !is.finite(sdv) || !(sdv > 0)) return(na)
  s <- (cells$s_lo[ok] + cells$s_hi[ok]) / 2
  # complete or quasi-complete separation: every event on one side of every non-event
  re <- range(s[e > 0]); rn <- range(s[ny - e > 0])
  if (re[1] >= rn[2] || re[2] <= rn[1]) return(na)
  # iterated to a tight tolerance: the Wald test reads the third digit of the slope
  f <- suppressWarnings(stats::glm.fit(cbind(1, (s - mu) / sdv), e / ny, weights = ny, family = stats::binomial(),
                                       control = list(epsilon = 1e-12, maxit = 50)))
  b <- f$coefficients[2]
  if (!isTRUE(f$converged) || !is.finite(b) || f$rank < 2L) return(na)
  # (X'WX)^-1 from the QR of the fit, in the pivot order of its columns
  V <- chol2inv(f$qr$qr[1:2, 1:2, drop = FALSE])
  j <- match(2L, f$qr$pivot[1:2])
  list(slope = as.double(b) / sdv, var = V[j, j] * sum(cells$w2_y[ok]) / sum(ny) / sdv^2)
}

#' Cells of the rest of a group: the pooled cells less those of one segment
#'
#' Weighted sums leave rounding residue where the segment holds a whole
#' cell; such a cell is emptied.
#' @keywords internal
#' @noRd
.seg_rest <- function(pc, cells) {
  key <- if ("g" %in% names(pc)) "g" else "s"
  i <- match(cells[[key]], pc[[key]])
  ny <- pc$n_y; e <- pc$e; w2 <- pc$w2_y
  ny[i] <- ny[i] - cells$n_y; e[i] <- e[i] - cells$e; w2[i] <- w2[i] - cells$w2_y
  ny[ny < 1e-9 * max(pc$n_y, 0)] <- 0
  list(s_lo = pc$s_lo, s_hi = pc$s_hi, n_y = ny, e = pmin(pmax(e, 0), ny), w2_y = pmax(w2, 0))
}

#' Discrimination of one set of cells: AUC, DeLong error, interval, Gini, KS
#'
#' No interval without a positive standard error (fewer than two events or
#' non-events, or every score tied): there is nothing to rank.
#' @keywords internal
#' @noRd
.seg_disc <- function(cells, direction, level, n_boot) {
  d <- .study_discrimination(cells, direction, n_boot, level)
  c1 <- cells$e; c0 <- pmax(cells$n_y - cells$e, 0)
  if (!identical(direction, "higher_is_riskier")) { c1 <- rev(c1); c0 <- rev(c0) }
  se <- if (nrow(cells)) .study_delong_counts(c1, c0, sum(cells$e_raw), sum(cells$n_y_raw - cells$e_raw)) else NA_real_
  z <- stats::qnorm(1 - (1 - level) / 2)
  lo <- hi <- NA_real_
  if (isTRUE(se > 0)) {
    # the bootstrap interval when asked, the DeLong normal interval otherwise
    boot <- d$n_boot > 0L
    lo <- if (boot) d$auc_lo else max(0, d$auc - z * se)
    hi <- if (boot) d$auc_hi else min(1, d$auc + z * se)
  }
  list(auc = d$auc, se = se, lo = lo, hi = hi, gini = d$gini, ks = d$ks)
}

#' Labels of a segment or group column in their reading order
#'
#' Numbers in numeric order, anything else as sorted text; `extra` (the
#' label of the missing segment) comes last after numbers.
#' @keywords internal
#' @noRd
.seg_levels <- function(v, labels, extra = NULL) {
  num <- .study_num_levels(list(v))
  if (is.null(num)) sort(unique(labels)) else c(num, extra)
}

#' Segment statistics from one count table keyed by segment, group and score
#' @keywords internal
#' @noRd
.segments_fit <- function(score, y, seg, byv, w, objective, direction, n_bands, level, min_events, n_boot, seed,
                          auc_tol, offset_tol, slope_tol, max_cells, target, segment, by, call) {
  fn <- "scr_segments"
  n_bands <- .study_whole(n_bands, "n_bands", fn, lower = 1)
  level <- .study_level(level, fn)
  min_events <- .study_whole(min_events, "min_events", fn)
  n_boot <- .study_whole(n_boot, "n_boot", fn)
  .study_num1(auc_tol, "auc_tol", fn, lower = 0); .study_num1(offset_tol, "offset_tol", fn, lower = 0)
  .study_num1(slope_tol, "slope_tol", fn, lower = 0)
  if (!is.numeric(score)) stop(fn, "(): the score must be numeric.", call. = FALSE)
  sg <- as.character(seg); sg[is.na(sg)] <- "(missing)"
  gv <- if (is.null(byv)) rep("all", length(sg)) else as.character(byv)
  if (anyNA(gv)) stop(fn, "(): the column '", by, "' has missing values.", call. = FALSE)
  # a factor segment in the order of its levels, the missing segment last
  seg_lv <- if (is.factor(seg)) unique(c(.study_levels(seg), "(missing)")) else .seg_levels(seg, sg, "(missing)")
  grp_lv <- if (is.null(byv)) "all" else .seg_levels(byv, gv)
  h <- .study_hist(score, y, w = w, by = list(sample = sg, group = gv), max_cells = max_cells, fn = fn)
  if (!nrow(h)) stop(fn, "(): no row has a score (and a positive weight).", call. = FALSE)
  all_cells <- .study_collapse(h)
  side <- .study_side(direction)
  cuts <- .study_cuts(all_cells, n_bands, "uniform", NULL, side)$cuts
  B <- length(cuts) + 1L
  labs <- .study_labels(cuts)
  # the score is standardized on the pooled cells for the logistic fits
  mid <- (all_cells$s_lo + all_cells$s_hi) / 2
  mu <- sum(all_cells$n * mid) / sum(all_cells$n)
  sdv <- sqrt(sum(all_cells$n * (mid - mu)^2) / sum(all_cells$n))
  band_sums <- function(cells) .study_sum(cells, findInterval(cells$s, cuts) + 1L, B, c("n", "w2", "n_y", "e"))
  alpha <- 1 - level
  .scr_local_seed(seed)
  rows <- list(); tests <- list(); pool <- list()
  for (g in grp_lv[grp_lv %in% h[["group"]]]) {
    # the rows of the group, selected outside `[` (a pooled table has a column `g`)
    in_g <- h[["group"]] == g
    hg <- h[in_g]
    data.table::setattr(hg, "edges", attr(h, "edges"))
    pc <- .study_collapse(hg)
    P <- band_sums(pc)
    Rb <- ifelse(P$n_y > 0, P$e / P$n_y, 0)
    NP <- sum(P$n)
    pn_eff <- P$n / NP * .study_kish(NP, sum(P$w2))
    pd <- .seg_disc(pc, direction, level, n_boot)
    b_pool <- .seg_slope(pc, mu, sdv)$slope
    NYp <- sum(pc$n_y)
    pool[[g]] <- data.table::data.table(
      group = g, n = sum(pc$n), events = sum(pc$e), rate = if (NYp > 0) sum(pc$e) / NYp else NA_real_,
      auc = pd$auc, auc_se = pd$se, auc_lo = pd$lo, auc_hi = pd$hi, gini = pd$gini, ks = pd$ks, slope = b_pool)
    out <- list()
    for (s in seg_lv[seg_lv %in% hg[["sample"]]]) {
      cells <- .study_cells(hg, s)
      S <- band_sums(cells)
      N <- sum(cells$n); NY <- sum(cells$n_y); E <- sum(cells$e)
      d <- .seg_disc(cells, direction, level, n_boot)
      rate <- if (NY > 0) E / NY else NA_real_
      # indirect standardization: the pooled band rates applied to the band volumes of the segment
      expd <- sum(S$n_y * Rb)
      neff <- .study_kish(NY, sum(cells$w2_y))
      ci <- .study_jeffreys(rate * neff, neff, level)
      er <- if (NY > 0) expd / NY else NA_real_
      inside <- function(v) isTRUE(v > 0 && v < 1)
      sn_eff <- .study_kish(N, sum(S$w2))
      ps <- .psi_counts(pn_eff, if (N > 0) S$n / N * sn_eff else S$n, labs, alpha, c(0.10, 0.25))
      # the segment is part of the pooled rows: Var(p_s - p_pool) carries their
      # covariance, 1/n_s - 1/N without weights; no test when it is the whole pool
      m_s <- sn_eff; m_p <- sum(pn_eff); share <- N / NP
      vf <- 1 / m_s + 1 / m_p - 2 * share / m_s
      crit <- if (!is.na(ps$critical) && share < 1 && vf > 0) ps$critical * vf / (1 / m_s + 1 / m_p) else NA_real_
      # the slope of the segment against that of the rest of the group, an independent sample
      sl <- .seg_slope(cells, mu, sdv)
      rs <- .seg_slope(.seg_rest(pc, cells), mu, sdv)
      ev_raw <- sum(cells$e_raw)
      out[[s]] <- data.table::data.table(
        group = g, segment = s, n = N, events = E, rate = rate, auc = d$auc, auc_se = d$se, auc_lo = d$lo,
        auc_hi = d$hi, gini = d$gini, ks = d$ks, psi = ps$psi, psi_critical = crit, expected = expd,
        oe_ratio = if (expd > 0) E / expd else NA_real_, oe_lo = if (isTRUE(er > 0)) ci$lo / er else NA_real_,
        oe_hi = if (isTRUE(er > 0)) ci$hi / er else NA_real_,
        offset = if (inside(rate) && inside(er)) stats::qlogis(rate) - stats::qlogis(er) else NA_real_,
        slope = sl$slope, slope_ratio = if (isTRUE(b_pool != 0)) sl$slope / b_pool else NA_real_,
        p_slope = 2 * stats::pnorm(-abs(sl$slope - rs$slope) / sqrt(sl$var + rs$var)),
        few = ev_raw < min_events || sum(cells$n_y_raw) - ev_raw < min_events)
    }
    tg <- data.table::rbindlist(out)
    # equal AUC across the segments: inverse-variance weighted mean and chi-square
    ok <- is.finite(tg$auc) & is.finite(tg$auc_se) & tg$auc_se > 0
    S_ok <- sum(ok)
    auc_w <- NA_real_; Q <- NA_real_; pq <- NA_real_
    dif <- rep(NA_real_, nrow(tg)); pz <- rep(NA_real_, nrow(tg))
    if (S_ok >= 1L) {
      wt <- 1 / tg$auc_se[ok]^2
      auc_w <- sum(wt * tg$auc[ok]) / sum(wt)
      dif[ok] <- tg$auc[ok] - auc_w
    }
    if (S_ok >= 2L) {
      Q <- sum(wt * (tg$auc[ok] - auc_w)^2)
      pq <- stats::pchisq(Q, S_ok - 1L, lower.tail = FALSE)
      # the mean contains the segment: Var(auc_s - auc_w) = se_s^2 - 1 / sum(w)
      v <- pmax(tg$auc_se[ok]^2 - 1 / sum(wt), 0)
      pz[ok] <- ifelse(v > 0, 2 * stats::pnorm(-abs(dif[ok]) / sqrt(v)), NA_real_)
    }
    # Holm across the segments; two segments give one test twice, left as it is
    p_sl <- tg$p_slope
    data.table::set(tg, j = "p_slope", value = NULL)
    tg[, `:=`(auc_diff = dif, p_auc = pz, p_auc_adj = if (S_ok == 2L) pz else .study_holm(pz), p_slope = p_sl,
              p_slope_adj = if (nrow(tg) == 2L) p_sl else .study_holm(p_sl))]
    # a different ranking needs both a material and a significant difference
    hit <- function(size, p) !is.na(size) & size & !is.na(p) & p < alpha
    sep <- hit(abs(tg$auc_diff) > auc_tol, tg$p_auc_adj) | hit(abs(tg$slope_ratio - 1) > slope_tol, tg$p_slope_adj)
    off <- !is.na(tg$offset) & abs(tg$offset) > offset_tol & !is.na(tg$oe_lo) & (tg$oe_lo > 1 | tg$oe_hi < 1)
    tg[, action := ifelse(few, "too few events", ifelse(sep, "separate model", ifelse(off, "offset", "shared")))]
    tg[, few := NULL]
    rows[[g]] <- tg
    tests[[g]] <- data.table::data.table(group = g, segments = S_ok, auc_w = auc_w, statistic = Q,
                                         df = if (S_ok >= 2L) S_ok - 1L else NA_integer_, p_value = pq)
  }
  structure(list(table = data.table::rbindlist(rows), test = data.table::rbindlist(tests),
                 pooled = data.table::rbindlist(pool), cuts = cuts, segment = segment, by = by, level = level,
                 min_events = min_events, n_boot = n_boot, auc_tol = auc_tol, offset_tol = offset_tol,
                 slope_tol = slope_tol, objective = objective, direction = direction, target = target,
                 weighted = !is.null(w), quantized = !is.null(attr(h, "edges")), call = call),
            class = c("scr_segments", "list"))
}

#' @export
print.scr_segments <- function(x, ...) {
  cat(sprintf("<scr_segments> target \"%s\" | objective %s | %s | segments of '%s'%s\n", x$target, x$objective,
              x$direction, x$segment, if (is.null(x$by)) "" else sprintf(" by '%s'", x$by)))
  cat(sprintf("  tolerances: AUC %s, offset %s, slope %s | fewest events %d | level %s\n", .g3(x$auc_tol),
              .g3(x$offset_tol), .g3(x$slope_tol), x$min_events, paste0(.g3(100 * x$level), "%")))
  for (i in seq_len(nrow(x$pooled))) {
    p <- x$pooled[i]; te <- x$test[i]
    t <- x$table[x$table$group == p$group]
    cat(sprintf("\n%s: n %s | rate %s | AUC %s | equal AUC across %d segments: chi-square %s (df %s), p %s\n",
                if (is.null(x$by)) "Pooled" else sprintf("Group '%s'", p$group), .study_n(p$n), .claims_pct(p$rate),
                .study_f(p$auc, "%.4f"), te$segments, .study_f(te$statistic, "%.2f"),
                if (is.na(te$df)) "-" else as.character(te$df), .study_f(te$p_value, "%.4f")))
    cat(sprintf("  %-14s %9s %8s %-24s %7s %-20s %7s %11s  %s\n", "segment", "n", "rate", "AUC [lo, hi]", "PSI",
                "O/E [lo, hi]", "offset", "slope ratio", "action"))
    for (j in seq_len(nrow(t))) {
      r <- t[j]
      # no interval for a segment left out of the test (no standard error, or every score tied)
      cat(sprintf("  %-14s %9s %8s %-24s %7s %-20s %7s %11s  %s\n", substr(r$segment, 1, 14), .study_n(r$n),
                  .claims_pct(r$rate),
                  if (is.na(r$auc)) "-" else if (is.na(r$auc_lo)) sprintf("%.4f -", r$auc) else
                    sprintf("%.4f [%.3f, %.3f]", r$auc, r$auc_lo, r$auc_hi),
                  .study_f(r$psi, "%.4f"),
                  if (is.na(r$oe_ratio)) "-" else sprintf("%.2f [%.2f, %.2f]", r$oe_ratio, r$oe_lo, r$oe_hi),
                  .study_f(r$offset, "%+.2f"), .study_f(r$slope_ratio, "%.2f"), r$action))
    }
  }
  invisible(x)
}

#' @rdname scr_export
#' @export
scr_export.scr_segments <- function(x, dir, stamp = TRUE, ...) {
  .study_dots(list(...), "scr_export")
  .need_openxlsx()
  out_dir <- .export_dir(dir, stamp)
  tag <- .file_tag(x$target)
  settings <- .kv_table(list(target = x$target, objective = x$objective, direction = x$direction,
                             segment = x$segment, by = x$by %||% "none", level = x$level,
                             min_events = x$min_events, auc_tol = x$auc_tol, offset_tol = x$offset_tol,
                             slope_tol = x$slope_tol, n_boot = x$n_boot))
  cuts <- data.frame(cut = seq_along(x$cuts), score = x$cuts)
  sheets <- lapply(list(Segments = x$table, Test = x$test, Pooled = x$pooled, Cuts = cuts, Settings = settings),
                   .study_sheet)
  files <- list(xlsx = .scr_write_xlsx(sheets, file.path(out_dir, sprintf("segments_%s.xlsx", tag))))
  for (f in files) msg("  %s", f)
  x$files <- files
  invisible(x)
}

# data.table column names used without quotes in this file
utils::globalVariables(c("action", "few", "auc_diff", "p_auc", "p_auc_adj", "p_slope", "p_slope_adj"))
