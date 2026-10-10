# ============================================================================ #
# claims.R - probability statements about the event rate of score groups
# ============================================================================ #
# A claim names a group of the score (a band, a tier or a score range), a
# direction and a rate. It is read on the count table of a score study, so
# its cost is that of the cells, never of the rows.
# ============================================================================ #

#' Probability statements about the event rate of score groups
#'
#' Tests claims such as "the customers scoring 625 or more respond at a rate
#' of at least 60%" or "tier 'low' defaults at no more than 2%" on a study
#' sample, and writes each verdict as one English sentence with the observed
#' rate and its one-sided bound.
#'
#' @section Claims:
#'
#' `claims` is a `data.frame` with one row per claim:
#' \describe{
#'   \item{`op`}{`">="` (the event rate is at least `rate`) or `"<="` (at
#'     most `rate`).}
#'   \item{`rate`}{The claimed event rate, in (0, 1).}
#'   \item{`label`}{A band or tier of the study: its label (the `label`
#'     column of the study table, or the numbered `tier_label` of a tiers
#'     study) or its number (the `band` or `tier` column).}
#'   \item{`score_lo`, `score_hi`}{Instead of a label, a score range
#'     `[score_lo, score_hi)`; a missing end is open. A claim with neither a
#'     label nor a score range is on the whole sample.}
#'   \item{`name`}{Optional: a name for the printed statement.}
#' }
#'
#' @section Test:
#'
#' Each claim is read on the rows of the group with a known outcome in the
#' study sample. With weights, the counts are taken to the Kish effective
#' size \eqn{n = (\sum w)^2 / \sum w^2} and \eqn{x = \hat p\, n} effective
#' events (without weights, the counts themselves). The claim
#' "rate >= r" is the alternative to \eqn{H_0: p \le r}, tested by the exact
#' one-sided binomial p-value \eqn{P(X \ge x)}, \eqn{X \sim
#' \mathrm{Binomial}(n, r)}; \eqn{x} and \eqn{n} are rounded to whole numbers
#' for the binomial only. The refutation is the opposite test,
#' \eqn{P(X \le x)}. For "rate <= r" the two tests swap. With
#' `adjust = "holm"`, the claim p-values and the refutation p-values are each
#' Holm-adjusted across the claims. They are two Holm families, the tests of
#' the claims and the tests of the refutations, each controlled at
#' `1 - level`; a claim and its refutation can never both be significant.
#'
#' The verdict is `"supported"` when the adjusted p-value of the claim is
#' below `1 - level`, `"refuted"` when the adjusted p-value of the
#' refutation is, and `"not proven"` otherwise (also for a group without a
#' row with a known outcome). `bound` is the one-sided Jeffreys bound at
#' `level` in the direction of the claim, on the unrounded effective counts:
#' the lower bound, the Beta(\eqn{x + 1/2, n - x + 1/2}) quantile
#' `1 - level`, for "rate >= r" (0 when \eqn{x = 0}), and the upper bound,
#' its quantile `level`, for "rate <= r" (1 when \eqn{x = n}). The bound
#' describes the uncertainty; the verdict comes from the tests.
#'
#' @section Average and floor:
#'
#' `type = "average"` reads the claim on the event rate of the whole group.
#' `type = "floor"` reads it on the weakest end of the group under a
#' monotone fit of the rate: a claim that holds there holds for the part of
#' the group where it is hardest to meet, not only for the average.
#'
#' The reference rows of the group are cut into `floor_bins` pre-bins of
#' equal share (tie-safe, as the bands of [scr_bands()]), and their event
#' rates are fitted by pool adjacent violators: the rate is made monotone
#' along the score, rising toward its event-rich end. The weakest end is the
#' run of pre-bins with the lowest fitted rate of the group for
#' "rate >= r", with the highest for "rate <= r". Its score range is found
#' on the reference sample and the claim is tested on the rows of the study
#' sample in that range, so the end is not chosen on the rows it is tested
#' on (unless the study sample is the reference). When the fitted rate is
#' flat over the group, the weakest end is the whole group.
#'
#' Two limits follow from the fit. A dip of the rate inside the group is
#' pooled with its neighbors by the monotone fit, so the floor does not see
#' a weak pocket in the middle of the group, only the end the fit calls
#' weakest. And `floor_bins` sets the resolution: the floor speaks for runs
#' of pre-bins, about one part in `floor_bins` of the group or more, never
#' for a single customer (a fit on single score values would spike at its
#' ends and read the floor on a handful of rows).
#'
#' @section Input:
#'
#' An object from [scr_bands()] or [scr_tiers()] is read as it is: labels
#' refer to its bands or tiers, and `sample` picks one of its samples (by
#' default its first study sample). A scorecard or a `data.frame` is first
#' summarized by a percentile study ([scr_bands()], `n_bands` bands frozen
#' on the reference, no bootstrap) whose count table keeps every distinct
#' score up to `max_cells`; the finite ends of the score ranges of `claims`
#' are forced as cell edges, so a score range is exact even when more
#' distinct scores are pooled into cells. A score range read on a study
#' whose cells were pooled must not cut through a cell, or the call is an
#' error.
#'
#' @inheritParams scr_bands
#' @param x An object from [scr_bands()] or [scr_tiers()], an object from
#'   [scr_scorecard()], or a `data.frame` with one row per scored case (or
#'   one row per score value with `counts = TRUE`).
#' @param claims A `data.frame` of claims; see the section Claims.
#' @param level Confidence level of the tests and of the one-sided bounds.
#'   For a scorecard, `NULL` uses `config$study_level` (0.95); for a score
#'   study, `NULL` uses the level of the study.
#' @param adjust `"holm"` (default) or `"none"`: multiplicity adjustment
#'   across the claims.
#' @param type `"average"` (the rate of the whole group) or `"floor"` (the
#'   rate at the weakest end of the group under a monotone fit of the rate);
#'   see the section Average and floor.
#' @param sample For a score study: the sample the claims are read on, `NULL`
#'   for its first study sample. For a scorecard: that sample (`"holdout"`).
#'   For a data.frame: the name of a column with sample labels, as in
#'   [scr_bands()]; the claims are read on `study`.
#' @param study For a data.frame: the label of the sample the claims are read
#'   on; `NULL` takes the first label other than the reference.
#' @param n_bands Bands of the internal percentile study (their labels and
#'   numbers can be used in `claims$label`). For a scorecard, `NULL` uses
#'   `config$study_bands` (20).
#' @param floor_bins Pre-bins of equal share of a group for
#'   `type = "floor"`: the resolution of the floor (default 10; 1 makes the
#'   floor the average).
#'
#' @return An object of class `c("scr_claims", "list")`:
#'   \describe{
#'     \item{`table`}{One row per claim: `name`, `group` (the label or the
#'       score range), `type`, `score_lo` and `score_hi` (the range tested),
#'       `n` (rows with a known outcome; their weighted volume under
#'       weights), `n_eff` (Kish effective size, the `n` of the test; equal
#'       to `n` without weights), `events`, `rate`, `bound` (one-sided
#'       Jeffreys bound), `op`,
#'       `rate_claimed`, `p_value` and `p_adj` (test of the claim),
#'       `p_refute` and `p_refute_adj` (the opposite test), `verdict` and
#'       `statement` (under weights, the sentence quotes the effective n of
#'       the test next to the weighted volume).}
#'     \item{`sample`, `reference`, `level`, `adjust`, `type`, `floor_bins`,
#'       `objective`, `direction`, `target`, `call`}{The settings
#'       (`floor_bins` is used by `type = "floor"` only).}
#'   }
#'
#' @references
#' Brown, L. D., Cai, T. T. and DasGupta, A. (2001). Interval estimation for
#' a binomial proportion. *Statistical Science*, 16(2), 101-133.
#' \doi{10.1214/ss/1009213286}
#'
#' Holm, S. (1979). A simple sequentially rejective multiple test procedure.
#' *Scandinavian Journal of Statistics*, 6(2), 65-70.
#'
#' Kish, L. (1965). *Survey Sampling*. Wiley.
#'
#' @seealso [scr_bands()] and [scr_tiers()] for the groups, [scr_operating()]
#'   to choose a cut under constraints.
#' @family score-studies
#' @examples
#' set.seed(1)
#' x <- rnorm(4000)
#' d <- data.frame(score = round(500 + 50 * x),
#'                 y = rbinom(4000, 1, plogis(-0.5 + 1.2 * x)))
#' claims <- data.frame(op = c(">=", ">=", "<="), rate = c(0.60, 0.80, 0.25),
#'                      score_lo = c(550, 550, NA), score_hi = c(NA, NA, 450),
#'                      name = c("top converts", "top converts well", "bottom is cold"))
#' cl <- scr_claims(d, claims, objective = "propensity")
#' cl
#' cl$table[, c("name", "n", "rate", "bound", "p_adj", "verdict")]
#'
#' # the same claims at the weakest end of each group, not only on its average
#' scr_claims(d, claims, objective = "propensity", type = "floor")$table$verdict
#'
#' # claims on the tiers of a study, by label
#' tr <- scr_tiers(d, objective = "propensity", n_tiers = 3)
#' scr_claims(tr, data.frame(label = "high", op = ">=", rate = 0.5))
#' @export
scr_claims <- function(x, claims, ...) UseMethod("scr_claims")

#' @rdname scr_claims
#' @export
scr_claims.scr_study <- function(x, claims, level = NULL, adjust = c("holm", "none"),
                                 type = c("average", "floor"), sample = NULL, floor_bins = 10L, ...) {
  .study_dots(list(...), "scr_claims")
  level <- .study_level(level %||% x$level %||% 0.95, "scr_claims")
  cl <- .claims_read(claims)
  .claims_fit(x, cl, level, match.arg(adjust), match.arg(type), sample, floor_bins, sys.call())
}

#' @rdname scr_claims
#' @export
scr_claims.scr_scorecard <- function(x, claims, level = NULL, adjust = c("holm", "none"),
                                     type = c("average", "floor"), sample = "holdout", reference = "train",
                                     n_bands = NULL, max_cells = 1e5, floor_bins = 10L, ...) {
  fn <- "scr_claims"
  .study_dots(list(...), fn)
  cfg <- x$config
  level <- .study_level(level %||% cfg$study_level %||% 0.95, fn)
  .study_chr1(sample, "sample", fn)
  cl <- .claims_read(claims)
  # the ends of the score ranges are forced as cell edges: pooled cells never straddle them
  inp <- .study_input(x, sample = sample, reference = reference, max_cells = max_cells,
                      breaks = .claims_edges(cl), fn = fn)
  st <- .bands_fit(inp, n_bands %||% cfg$study_bands %||% 20L, "uniform", NULL, NULL, level, 0L, NULL, NULL)
  .claims_fit(st, cl, level, match.arg(adjust), match.arg(type), sample, floor_bins, sys.call())
}

#' @rdname scr_claims
#' @export
scr_claims.data.frame <- function(x, claims, level = 0.95, adjust = c("holm", "none"),
                                  type = c("average", "floor"), sample = NULL, score = "score", y = "y",
                                  objective = "risk", direction = NULL, weight = NULL, reference = NULL,
                                  study = NULL, counts = FALSE, n = "n", events = "events", n_bands = 20L,
                                  max_cells = 1e5, floor_bins = 10L, ...) {
  fn <- "scr_claims"
  .study_dots(list(...), fn)
  level <- .study_level(level, fn)
  if (!is.null(study)) .study_chr1(study, "study", fn)
  cl <- .claims_read(claims)
  inp <- .study_input(x, score = score, y = y, objective = objective, direction = direction, weight = weight,
                      sample = sample, reference = reference, study = study, counts = counts, n = n,
                      events = events, max_cells = max_cells, breaks = .claims_edges(cl), fn = fn)
  st <- .bands_fit(inp, n_bands, "uniform", NULL, NULL, level, 0L, NULL, NULL)
  .claims_fit(st, cl, level, match.arg(adjust), match.arg(type), study, floor_bins, sys.call())
}

#' Read and check a table of claims
#' @keywords internal
#' @noRd
.claims_read <- function(claims, fn = "scr_claims") {
  if (!is.data.frame(claims) || !nrow(claims)) {
    stop(fn, "(): `claims` must be a data.frame with one row per claim.", call. = FALSE)
  }
  nm <- names(claims)
  miss <- setdiff(c("op", "rate"), nm)
  if (length(miss)) stop(fn, "(): `claims` lacks the column(s) ", lst(miss), ".", call. = FALSE)
  K <- nrow(claims)
  op <- as.character(claims$op)
  if (anyNA(op) || !all(op %in% c(">=", "<="))) stop(fn, "(): `claims$op` must be \">=\" or \"<=\".", call. = FALSE)
  rate <- claims$rate
  if (!is.numeric(rate) || anyNA(rate) || any(rate <= 0 | rate >= 1)) {
    stop(fn, "(): `claims$rate` must be event rates in (0, 1).", call. = FALSE)
  }
  num <- function(cn) {
    if (!cn %in% nm) return(rep(NA_real_, K))
    v <- claims[[cn]]
    if (!all(is.na(v)) && !is.numeric(v)) stop(fn, "(): `claims$", cn, "` must be numeric.", call. = FALSE)
    v <- as.double(v)
    if (any(is.nan(v)) || any(is.infinite(v))) stop(fn, "(): `claims$", cn, "` must be finite or NA (open).", call. = FALSE)
    v
  }
  lo <- num("score_lo"); hi <- num("score_hi")
  if (any(!is.na(lo) & !is.na(hi) & lo >= hi)) stop(fn, "(): a score range needs `score_lo` < `score_hi`.", call. = FALSE)
  # a label is matched as text: a band or tier number becomes "3"
  lab <- as.character(if ("label" %in% nm) claims$label else rep(NA, K))
  if (any(!is.na(lab) & (!is.na(lo) | !is.na(hi)))) {
    stop(fn, "(): a claim takes either a `label` or a score range, not both.", call. = FALSE)
  }
  name <- if ("name" %in% nm) as.character(claims$name) else rep(NA_character_, K)
  data.table::data.table(name = name, op = op, rate = as.double(rate), label = as.character(lab),
                         score_lo = lo, score_hi = hi)
}

#' Finite ends of the score ranges, forced as cell edges of the count table
#' @keywords internal
#' @noRd
.claims_edges <- function(cl) {
  v <- c(cl$score_lo, cl$score_hi)
  v <- sort(unique(v[is.finite(v)]))
  if (length(v)) v else NULL
}

#' Score boundary between adjacent ascending cells `a` and `a + 1`
#'
#' The convention of the cuts of a study: midway between the two distinct
#' scores, or the bucket edge of a pooled table.
#' @keywords internal
#' @noRd
.study_bcut <- function(cells, a = seq_len(max(nrow(cells) - 1L, 0L))) {
  if (!length(a)) return(numeric())
  edges <- attr(cells, "edges")
  if (!is.null(edges) && "g" %in% names(cells)) edges[cells$g[a]] else .study_mid(cells$s_hi[a], cells$s_lo[a + 1L])
}

#' Cells of a score range `[lo, hi)`; an error when a pooled cell straddles an end
#' @keywords internal
#' @noRd
.claims_in <- function(cells, lo, hi, check = FALSE, fn = "scr_claims") {
  lo <- if (is.na(lo)) -Inf else lo
  hi <- if (is.na(hi)) Inf else hi
  if (check && nrow(cells)) {
    cut_in <- function(v) is.finite(v) && any(cells$s_lo < v & cells$s_hi >= v)
    bad <- c(lo, hi)[vapply(c(lo, hi), cut_in, logical(1))]
    if (length(bad)) {
      stop(fn, "(): the score range edge(s) ", lst(format(bad)), " fall inside a pooled cell of the study ",
           "(max_cells). Call scr_claims() on the scores (a scorecard or a data.frame), which keeps the ",
           "edges exact.", call. = FALSE)
    }
  }
  cells$s_lo >= lo & cells$s_lo < hi
}

#' Score range of the weakest end of a group on the reference cells
#'
#' The reference cells of the group are cut into `floor_bins` tie-safe
#' pre-bins of equal share; pool adjacent violators runs on the pre-bins,
#' oriented from the event-poor to the event-rich end of the score; the
#' edge is the run of blocks with the smallest smoothed rate for ">=" and
#' with the largest for "<=" (the first and the last level of the smoothed
#' rate). Returns the score range `[lo, hi)` of that run, `NA` ends when
#' the reference has no cell in the group; a group without a known outcome
#' on the reference is its own edge.
#' @keywords internal
#' @noRd
.claims_floor <- function(rc, inn, lo, hi, op, direction, floor_bins) {
  idx <- which(inn)   # the cells of the group: a contiguous run of the ascending reference cells
  if (!length(idx)) return(c(NA_real_, NA_real_))
  gc <- rc[idx]
  data.table::setattr(gc, "edges", attr(rc, "edges"))
  # tie-safe pre-bins of equal share of the group: an isotonic fit on single
  # score values spikes at its ends, so the edge is read at this resolution
  pc <- .study_cuts(gc, floor_bins, "uniform", NULL, .study_side(direction))$cuts
  P <- length(pc) + 1L
  S <- .study_sum(gc, findInterval(gc$s, pc) + 1L, P, c("n_y", "e"))
  # event-poor end first: low scores under higher_is_riskier, high scores otherwise
  ori <- if (identical(direction, "higher_is_riskier")) seq_len(P) else rev(seq_len(P))
  g <- .study_pav(S$e[ori], S$n_y[ori])
  br <- as.numeric(rowsum(S$e[ori], g, reorder = TRUE)) / as.numeric(rowsum(S$n_y[ori], g, reorder = TRUE))
  # the edge is the level set of the smoothed rate at that end: PAV leaves
  # adjacent blocks with equal rates apart, the smoothed function does not
  B <- length(br)
  blk <- if (!all(is.finite(br))) seq_len(B) else if (identical(op, ">=")) which(br <= br[1] * (1 + 1e-12)) else
    which(br >= br[B] * (1 - 1e-12))
  asc <- sort(ori[g %in% blk])
  i1 <- asc[1]; i2 <- asc[length(asc)]
  c(if (i1 == 1L) lo else pc[i1 - 1L], if (i2 == P) hi else pc[i2])
}

#' Read the claims on one sample of a score study
#' @keywords internal
#' @noRd
.claims_fit <- function(x, cl, level, adjust, type, sample, floor_bins, call) {
  fn <- "scr_claims"
  floor_bins <- .study_whole(floor_bins, "floor_bins", fn, lower = 1)
  smp <- sample %||% setdiff(x$samples, x$reference)[1]
  if (is.na(smp)) smp <- x$reference
  .study_chr1(smp, "sample", fn)
  if (!smp %in% x$samples) stop(fn, "(): sample '", smp, "' is not in the study (", lst(x$samples), ").", call. = FALSE)
  cells <- .study_cells(x$hist, smp)
  rc <- if (identical(type, "floor")) .study_cells(x$hist, x$reference) else NULL
  tb <- x$table[x$table$sample == smp]
  id_col <- if ("tier" %in% names(tb)) "tier" else "band"
  kind <- if (identical(id_col, "tier")) "tier" else "band"
  K <- nrow(cl)
  out <- vector("list", K)
  for (k in seq_len(K)) {
    # the group: a band or tier of the study, or a score range
    if (!is.na(cl$label[k])) {
      # the plain label, the numbered one of a tiers study, or the number
      i <- which(tb$label == cl$label[k])
      if (!length(i) && !is.null(tb$tier_label)) i <- which(tb$tier_label == cl$label[k])
      if (!length(i)) i <- which(as.character(tb[[id_col]]) == cl$label[k])
      if (length(i) != 1L) {
        stop(fn, "(): label '", cl$label[k], "' is ", if (length(i)) "ambiguous" else "not a",
             " ", kind, " of the study (labels: ", lst(tb$label), ").", call. = FALSE)
      }
      lo <- tb$score_lo[i]; hi <- tb$score_hi[i]
      group <- tb$label[i]
      gtext <- sprintf("rows in %s '%s'", kind, tb$label[i])
    } else {
      lo <- cl$score_lo[k]; hi <- cl$score_hi[k]
      group <- .claims_range(lo, hi)
      gtext <- if (identical(group, "all scores")) "all rows" else paste("rows with", group)
      .claims_in(cells, lo, hi, check = TRUE)
    }
    lo <- if (is.na(lo)) -Inf else lo
    hi <- if (is.na(hi)) Inf else hi
    tlo <- lo; thi <- hi
    if (identical(type, "floor")) {
      fr <- .claims_floor(rc, .claims_in(rc, lo, hi), lo, hi, cl$op[k], x$direction, floor_bins)
      tlo <- fr[1]; thi <- fr[2]
      edge <- if (identical(cl$op[k], ">=")) "low-rate" else "high-rate"
      # a flat fitted rate makes the edge the group itself: say so instead of repeating its range
      whole <- !is.na(tlo) && tlo == lo && thi == hi
      gtext <- sprintf("rows at the %s edge (%s) of %s", edge,
                       if (is.na(tlo)) "no reference row" else if (whole) "the whole group" else
                         .claims_range(tlo, thi, both = TRUE), gtext)
    }
    sel <- if (is.na(tlo)) rep(FALSE, nrow(cells)) else .claims_in(cells, tlo, thi)
    ny <- sum(cells$n_y[sel]); e <- sum(cells$e[sel]); w2 <- sum(cells$w2_y[sel])
    out[[k]] <- data.table::data.table(group = group, gtext = gtext, score_lo = tlo, score_hi = thi,
                                       n = ny, n_eff = .study_kish(ny, w2), events = e)
  }
  g <- data.table::rbindlist(out)
  rate <- ifelse(g$n > 0, g$events / g$n, NA_real_)
  te <- .claims_test(rate * g$n_eff, g$n_eff, cl$rate, cl$op, level)
  p_adj <- if (identical(adjust, "holm")) .study_holm(te$p) else te$p
  r_adj <- if (identical(adjust, "holm")) .study_holm(te$p_refute) else te$p_refute
  verdict <- ifelse(!is.na(p_adj) & p_adj < 1 - level, "supported",
                    ifelse(!is.na(r_adj) & r_adj < 1 - level, "refuted", "not proven"))
  claim <- sprintf("rate %s %s%%", cl$op, .g3(100 * cl$rate))
  name <- ifelse(is.na(cl$name), claim, cl$name)
  bound_txt <- sprintf("%.0f%% one-sided %s bound %s", 100 * level, ifelse(cl$op == ">=", "lower", "upper"),
                       .claims_pct(te$bound))
  # under weights the test runs on the Kish effective size: quote it, next to the weighted volume
  n_txt <- if (isTRUE(x$weighted)) sprintf("effective n = %s of a weighted volume of %s", .study_n(g$n_eff),
                                           .study_n(g$n)) else sprintf("n = %s", .study_n(g$n))
  statement <- ifelse(g$n > 0,
    sprintf("On '%s' (%s), %s had an event rate of %s (%s); the claim %s is %s.", smp, n_txt, g$gtext,
            .claims_pct(rate), bound_txt, .claims_quote(cl$name, claim), verdict),
    sprintf("On '%s', %s had no row with a known outcome; the claim %s is %s.", smp, g$gtext,
            .claims_quote(cl$name, claim), verdict))
  tab <- data.table::data.table(
    name = name, group = g$group, type = type, score_lo = g$score_lo, score_hi = g$score_hi, n = g$n,
    n_eff = g$n_eff, events = g$events, rate = rate, bound = te$bound, op = cl$op, rate_claimed = cl$rate,
    p_value = te$p, p_adj = p_adj, p_refute = te$p_refute, p_refute_adj = r_adj, verdict = verdict,
    statement = statement)
  structure(list(table = tab, sample = smp, reference = x$reference, level = level, adjust = adjust, type = type,
                 floor_bins = floor_bins, objective = x$objective, direction = x$direction, target = x$target, call = call),
            class = c("scr_claims", "list"))
}

#' Exact one-sided binomial tests of the claims and their Jeffreys bounds
#'
#' `x` and `n` are the effective events and size; they are rounded for the
#' binomial only. For ">=": p = P(X >= x), refutation P(X <= x); for "<="
#' the two swap. The bound is the one-sided Jeffreys bound at `level`.
#' @keywords internal
#' @noRd
.claims_test <- function(x, n, r, op, level) {
  K <- length(n)
  p <- pr <- bound <- rep(NA_real_, K)
  ok <- is.finite(n) & n > 0 & is.finite(x)
  nr <- round(n); xr <- pmin(pmax(round(x), 0), nr)
  ok <- ok & nr >= 1
  ge <- op == ">="
  up <- stats::pbinom(xr - 1, nr, r, lower.tail = FALSE)   # P(X >= x)
  dn <- stats::pbinom(xr, nr, r)                           # P(X <= x)
  p[ok] <- ifelse(ge, up, dn)[ok]
  pr[ok] <- ifelse(ge, dn, up)[ok]
  # Jeffreys bound on the unrounded effective counts, in the direction of the claim
  okb <- is.finite(n) & n > 0 & is.finite(x)
  xo <- pmin(pmax(x, 0), n)
  lb <- ifelse(xo <= 0, 0, stats::qbeta(1 - level, xo + 0.5, n - xo + 0.5))
  ub <- ifelse(xo >= n, 1, stats::qbeta(level, xo + 0.5, n - xo + 0.5))
  bound[okb] <- ifelse(ge, lb, ub)[okb]
  list(p = p, p_refute = pr, bound = bound)
}

#' Text of a score range: "score >= a", "score < b", "a <= score < b"
#' @keywords internal
#' @noRd
.claims_range <- function(lo, hi, both = FALSE) {
  f <- function(v) format(v, digits = 7, scientific = FALSE, trim = TRUE)
  lo_ok <- !is.na(lo) && is.finite(lo); hi_ok <- !is.na(hi) && is.finite(hi)
  if (lo_ok && hi_ok) sprintf("%s <= score < %s", f(lo), f(hi))
  else if (lo_ok) sprintf("score >= %s", f(lo))
  else if (hi_ok) sprintf("score < %s", f(hi))
  else if (both) "every score" else "all scores"
}

#' A rate in percent: one decimal from 1%, three significant digits below
#' @keywords internal
#' @noRd
.claims_pct <- function(x) {
  ifelse(is.na(x), "-", ifelse(abs(100 * x) >= 1 | x == 0, sprintf("%.1f%%", 100 * x), paste0(.g3(100 * x), "%")))
}

#' The quoted claim, with its name when given
#' @keywords internal
#' @noRd
.claims_quote <- function(name, claim) {
  ifelse(is.na(name), sprintf("'%s'", claim), sprintf("'%s' (%s)", name, claim))
}

#' @export
print.scr_claims <- function(x, ...) {
  cat(sprintf("<scr_claims> target \"%s\" | objective %s | sample '%s' | level %s (one-sided) | adjustment %s | %s\n",
              x$target, x$objective, x$sample, paste0(.g3(100 * x$level), "%"), x$adjust, x$type))
  t <- x$table
  cat(sprintf("  %-24s %-26s %9s %8s %8s %-10s %8s  %s\n", "claim", "group", "n", "rate", "bound", "claimed",
              "p_adj", "verdict"))
  for (i in seq_len(nrow(t))) {
    cat(sprintf("  %-24s %-26s %9s %8s %8s %-10s %8s  %s\n", substr(t$name[i], 1, 24), substr(t$group[i], 1, 26),
                .study_n(t$n[i]), .claims_pct(t$rate[i]), .claims_pct(t$bound[i]),
                paste(t$op[i], .claims_pct(t$rate_claimed[i])), .study_f(t$p_adj[i], "%.4f"), t$verdict[i]))
  }
  cat("\n")
  for (s in t$statement) cat(paste(strwrap(s, width = 100, exdent = 2), collapse = "\n"), "\n", sep = "")
  invisible(x)
}

#' @rdname scr_export
#' @export
scr_export.scr_claims <- function(x, dir, stamp = TRUE, ...) {
  .study_dots(list(...), "scr_export")
  .need_openxlsx()
  out_dir <- .export_dir(dir, stamp)
  tag <- .file_tag(x$target)
  settings <- .kv_table(list(target = x$target, objective = x$objective, direction = x$direction,
                             sample = x$sample, reference = x$reference, level = x$level, adjust = x$adjust,
                             type = x$type))
  sheets <- lapply(list(Claims = x$table, Settings = settings), .study_sheet)
  files <- list(xlsx = .scr_write_xlsx(sheets, file.path(out_dir, sprintf("claims_%s.xlsx", tag))))
  for (f in files) msg("  %s", f)
  x$files <- files
  invisible(x)
}
