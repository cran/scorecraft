# ============================================================================ #
# score-cross.R - two scores on the same rows: cross table, association, overlap
# ============================================================================ #
# The rows are read once: each gets its band and its selection depth under
# both scores, and one grouped pass counts them per (band, band, depth,
# depth). The cross table and every overlap set are sums over that small
# table; only the rank association goes back to the rows, in O(n log n).
# ============================================================================ #

#' Two scores on the same rows
#'
#' Crosses two scores read on the same rows (a credit score and a churn
#' score, a champion and a challenger): the cross table of their bands with
#' the event rate of one or two outcomes, the rank association of the
#' scores, and the overlap of the rows each one selects at a few depths,
#' with the swap-in and swap-out sets.
#'
#' @section Bands:
#'
#' Each score is cut into `n_bands` bands of equal share, tie-safe as in
#' [scr_bands()] (band 1 is the event-richest), unless `cuts_a` or `cuts_b`
#' gives the cuts: a numeric vector, or an object from [scr_bands()] or
#' [scr_tiers()], whose cuts, numbers and labels are then used (the tier
#' labels for tiers). Bands are left-closed, `score >= cut` being the upper
#' side.
#'
#' A study given as cuts also sets the objective and the direction of its
#' score, so the event-rich end of the overlap is the one the study was
#' fitted with. `objective_a`, `direction_a` (or their `_b` counterparts)
#' need not be given then; when given and different from the study, the call
#' is an error.
#'
#' @section Cross table:
#'
#' One row per pair of bands, every pair listed (empty ones with `n = 0`),
#' then the totals of each band of A (`band_b` missing, `label_b =
#' "total"`), of each band of B, and the grand total. Per row: `n`, `pct`
#' (share of all rows) and, for every outcome, `events`, `rate` with its
#' Jeffreys interval `rate_lo`, `rate_hi` (on the Kish effective size under
#' weights) and `lift` (the rate over the overall rate of that outcome).
#' With `y` the outcome columns have no suffix; with `y_a` and `y_b` they
#' end in `_a` and `_b`.
#'
#' @section Association:
#'
#' Spearman's rank correlation (Pearson on mid-ranks, as
#' `cor(method = "spearman")`) and Kendall's tau-b (pairs tied on either
#' score count neither way, as `cor(method = "kendall")`), both on the raw
#' scores and unweighted. The tau-b counts are exact, by Knight's (1966)
#' algorithm in `O(n log n)`. `oriented` multiplies each estimate by the
#' signs of the two directions, so it is positive when the two scores put
#' the same rows at their event-rich ends.
#'
#' @section Overlap:
#'
#' At each depth, each score selects the share `depth` of the rows from its
#' event-rich end (tie-safe: the selection stops at the boundary between two
#' distinct scores nearest to the target, so the share selected can differ
#' from `depth` by the share of one score value; `share_a` and `share_b`
#' report it). `overlap` counts the rows selected by A, by B, by both, by A
#' only and by B only, and the Jaccard index `n_both / (n_a + n_b -
#' n_both)`. With an outcome, `overlap_rates` gives the event rate of each
#' of the five sets with its Jeffreys interval: when B replaces A at the
#' same depth, `"B only"` is the swap-in and `"A only"` the swap-out.
#'
#' Rows with a missing or infinite value of either score, or a zero weight,
#' are left out (`n_dropped`); rows with a missing outcome count in the
#' volume but not in the rates of that outcome.
#'
#' @param x A `data.frame` with both scores on every row.
#' @param score_a,score_b Column names of the two scores.
#' @param y Column name of a 0/1 outcome read under both scores.
#' @param y_a,y_b Instead of `y`: the outcome of score A and of score B
#'   (two different targets on the same rows); either may be given alone.
#' @param objective_a,objective_b `"risk"` or `"propensity"`; `objective_b =
#'   NULL` takes `objective_a`. A study given as cuts sets the objective of
#'   its score.
#' @param direction_a,direction_b `"higher_is_safer"` or
#'   `"higher_is_riskier"`; `NULL` derives each from its objective, or takes
#'   it from the study given as cuts.
#' @param n_bands Bands of each score when its cuts are not given.
#' @param cuts_a,cuts_b Optional cuts of each score: a numeric vector, or an
#'   object from [scr_bands()] or [scr_tiers()], which also sets the
#'   objective and the direction of that score.
#' @param depths Shares of the rows selected by each score for the overlap,
#'   in (0, 1].
#' @param weight Optional column of non-negative case weights.
#' @param level Confidence level of the Jeffreys intervals.
#' @param ... Not used; an unknown argument is an error.
#'
#' @return An object of class `c("scr_score_cross", "list")`:
#'   \describe{
#'     \item{`table`}{The cross table: `band_a`, `label_a`, `band_b`,
#'       `label_b`, `n`, `pct` and the outcome columns (see the section Cross
#'       table).}
#'     \item{`overlap`}{One row per depth: `depth`, `cut_a`, `cut_b`,
#'       `share_a`, `share_b`, `n_a`, `n_b`, `n_both`, `n_a_only`,
#'       `n_b_only` and `jaccard`.}
#'     \item{`overlap_rates`}{With an outcome, one row per depth, outcome and
#'       set (`"A"`, `"B"`, `"both"`, `"A only"`, `"B only"`): `n`, `events`,
#'       `rate`, `rate_lo` and `rate_hi`.}
#'     \item{`association`}{One row per method (`"spearman"`,
#'       `"kendall_tau_b"`): `estimate`, `oriented` and `n`.}
#'     \item{`settings`}{A list: the score and outcome columns, objectives,
#'       directions, cuts, codes and labels of both scores, `depths`,
#'       `level`, `n` (the volume used: the sum of the weights, the number of
#'       rows without weights), `n_rows` (the rows used), `n_dropped` (rows
#'       left out) and `weighted`.}
#'   }
#'
#' @references
#' Brown, L. D., Cai, T. T. and DasGupta, A. (2001). Interval estimation for
#' a binomial proportion. *Statistical Science*, 16(2), 101-133.
#' \doi{10.1214/ss/1009213286}
#'
#' Knight, W. R. (1966). A computer method for calculating Kendall's tau
#' with ungrouped data. *Journal of the American Statistical Association*,
#' 61(314), 436-439. \doi{10.1080/01621459.1966.10480879}
#'
#' @seealso [scr_bands()] and [scr_tiers()] for the cuts of each score.
#' @family score-studies
#' @examples
#' set.seed(1)
#' n <- 4000
#' z <- rnorm(n)
#' d <- data.frame(credit = round(600 + 40 * (-z + rnorm(n, sd = 0.6))),
#'                 churn = round(450 + 30 * (0.4 * z + rnorm(n))),
#'                 default = rbinom(n, 1, plogis(-2 + z)),
#'                 left = rbinom(n, 1, 0.3))
#' cx <- scr_score_cross(d, "credit", "churn", y_a = "default", y_b = "left",
#'                       objective_b = "propensity", n_bands = 4)
#' cx
#' cx$association
#' cx$overlap_rates[cx$overlap_rates$depth == 0.1, ]
#' @export
scr_score_cross <- function(x, ...) UseMethod("scr_score_cross")

#' @rdname scr_score_cross
#' @export
scr_score_cross.data.frame <- function(x, score_a, score_b, y = NULL, y_a = NULL, y_b = NULL, objective_a = "risk",
                                       objective_b = NULL, direction_a = NULL, direction_b = NULL, n_bands = 5L,
                                       cuts_a = NULL, cuts_b = NULL, depths = c(0.05, 0.10, 0.20), weight = NULL,
                                       level = 0.95, ...) {
  fn <- "scr_score_cross"
  .study_dots(list(...), fn)
  .study_chr1(score_a, "score_a", fn); .study_chr1(score_b, "score_b", fn)
  for (nm in c("y", "y_a", "y_b", "weight")) {
    v <- get(nm)
    if (!is.null(v)) .study_chr1(v, nm, fn)
  }
  if (!is.null(y) && (!is.null(y_a) || !is.null(y_b))) stop(fn, "(): give `y`, or `y_a` and `y_b`, not both.", call. = FALSE)
  # a study given as cuts brings its objective and direction; an explicit argument must agree with it
  ra <- .cross_direction(if (!missing(objective_a)) objective_a, direction_a, cuts_a, "risk", "a", fn)
  rb <- .cross_direction(objective_b, direction_b, cuts_b, ra$objective, "b", fn)
  objective_a <- ra$objective; dir_a <- ra$direction
  objective_b <- rb$objective; dir_b <- rb$direction
  miss <- setdiff(c(score_a, score_b, y, y_a, y_b, weight), names(x))
  if (length(miss)) stop(fn, "(): column(s) ", lst(miss), " not in `x`.", call. = FALSE)
  for (nm in c(score_a, score_b, weight)) {
    if (!is.numeric(x[[nm]])) stop(fn, "(): column '", nm, "' must be numeric.", call. = FALSE)
  }
  n_bands <- .study_whole(n_bands, "n_bands", fn, lower = 1)
  if (!is.numeric(depths) || !length(depths) || anyNA(depths) || any(depths <= 0 | depths > 1)) {
    stop(fn, "(): `depths` must be shares in (0, 1].", call. = FALSE)
  }
  depths <- sort(unique(as.double(depths)))
  level <- .study_level(level, fn)
  # outcome columns: no suffix for a common `y`, "_a" and "_b" otherwise
  outs <- if (!is.null(y)) list(y) else c(if (!is.null(y_a)) list(y_a), if (!is.null(y_b)) list(y_b))
  sfx <- if (!is.null(y)) "" else c(if (!is.null(y_a)) "_a", if (!is.null(y_b)) "_b")
  names(outs) <- sfx

  sa <- as.double(x[[score_a]]); sb <- as.double(x[[score_b]])
  ok <- is.finite(sa) & is.finite(sb)
  wtd <- !is.null(weight)
  if (wtd) {
    w <- as.double(x[[weight]])
    if (anyNA(w) || any(!is.finite(w)) || any(w < 0)) stop(fn, "(): the weights must be finite and non-negative.", call. = FALSE)
    # a row with zero weight does not belong to the population
    ok <- ok & w > 0
    w <- w[ok]
  }
  yy <- lapply(outs, function(cn) .scr_y01(x[[cn]], fn)[ok])
  n_drop <- sum(!ok)
  sa <- sa[ok]; sb <- sb[ok]
  if (!length(sa)) stop(fn, "(): no row has both scores (and a positive weight).", call. = FALSE)
  side_a <- .study_side(dir_a); side_b <- .study_side(dir_b)
  # one count table per score: the cuts of the bands and of the depths
  ha <- .study_hist(sa, integer(length(sa)), w = if (wtd) w, fn = fn)
  hb <- .study_hist(sb, integer(length(sb)), w = if (wtd) w, fn = fn)
  ba <- .cross_bands(cuts_a, ha, n_bands, side_a, "cuts_a", fn)
  bb <- .cross_bands(cuts_b, hb, n_bands, side_b, "cuts_b", fn)
  da <- .cross_depth_cuts(ha, depths, side_a); db <- .cross_depth_cuts(hb, depths, side_b)

  # one grouped pass: band and depth level of both scores
  D <- length(depths)
  cols <- list(ba = ba$codes[findInterval(sa, ba$cuts) + 1L], bb = bb$codes[findInterval(sb, bb$cuts) + 1L],
               la = .cross_level(sa, da, side_a), lb = .cross_level(sb, db, side_b))
  agg <- if (wtd) alist(n = sum(w)) else alist(n = .N)
  if (wtd) cols$w <- w
  for (k in seq_along(outs)) {
    s <- sfx[k]; yk <- yy[[k]]
    kn <- !is.na(yk); ev <- kn & yk == 1L
    ww <- if (wtd) w else 1
    cols[[paste0("k", s)]] <- ww * kn; cols[[paste0("v", s)]] <- ww * ev; cols[[paste0("q", s)]] <- ww * ww * kn
    agg[[paste0("n_y", s)]] <- call("sum", as.name(paste0("k", s)))
    agg[[paste0("e", s)]] <- call("sum", as.name(paste0("v", s)))
    agg[[paste0("w2", s)]] <- call("sum", as.name(paste0("q", s)))
  }
  dt <- data.table::setDT(cols)
  j <- as.call(c(as.name("list"), agg))
  G <- dt[, eval(j), keyby = c("ba", "bb", "la", "lb")]
  cnt <- setdiff(names(G), c("ba", "bb", "la", "lb"))
  for (cn in cnt) data.table::set(G, j = cn, value = as.double(G[[cn]]))
  N <- sum(G$n)
  # overall event rate of each outcome, by position (the suffix of `y` is empty)
  tot <- vapply(sfx, function(s) {
    NY <- sum(G[[paste0("n_y", s)]])
    if (NY > 0) sum(G[[paste0("e", s)]]) / NY else NA_real_
  }, numeric(1), USE.NAMES = FALSE)

  tab <- .cross_table(G, cnt, ba, bb, sfx, tot, N, level)
  ov <- .cross_overlap(G, cnt, depths, da, db, sfx, unlist(outs), N, level)
  sg <- (if (side_a == "high") 1 else -1) * (if (side_b == "high") 1 else -1)
  sp <- if (length(sa) >= 2L) suppressWarnings(stats::cor(data.table::frank(sa), data.table::frank(sb))) else NA_real_
  kt <- .scr_kendall_tau_b(sa, sb)
  assoc <- data.table::data.table(method = c("spearman", "kendall_tau_b"), estimate = c(sp, kt),
                                  oriented = sg * c(sp, kt), n = length(sa))
  settings <- list(score_a = score_a, score_b = score_b, y = y, y_a = y_a, y_b = y_b, objective_a = objective_a,
                   objective_b = objective_b, direction_a = dir_a, direction_b = dir_b, cuts_a = ba$cuts,
                   cuts_b = bb$cuts, codes_a = ba$codes, codes_b = bb$codes, labels_a = ba$labels,
                   labels_b = bb$labels, bands_a = ba$source, bands_b = bb$source, n_bands = n_bands,
                   depths = depths, level = level, n = N, n_rows = length(sa), n_dropped = n_drop, weighted = wtd)
  structure(list(table = tab, overlap = ov$overlap, overlap_rates = ov$rates, association = assoc,
                 settings = settings, call = sys.call()),
            class = c("scr_score_cross", "list"))
}

#' Objective and direction of one score, checked
#'
#' `objective` and `direction` are the arguments as given (`NULL` when not
#' given). A score study passed as `cuts` sets both, and a given argument
#' that disagrees with it is an error: the bands of the study and the
#' event-rich end of the overlap must read the score the same way. Without
#' a study, the objective falls back on `default` and the direction follows
#' the objective.
#' @keywords internal
#' @noRd
.cross_direction <- function(objective, direction, cuts, default, which, fn) {
  if (!is.null(objective) && (!is.character(objective) || length(objective) != 1L ||
                              !objective %in% c("risk", "propensity"))) {
    stop(fn, "(): `objective_", which, "` must be \"risk\" or \"propensity\".", call. = FALSE)
  }
  if (!is.null(direction) && (!is.character(direction) || length(direction) != 1L ||
                              !direction %in% c("higher_is_safer", "higher_is_riskier"))) {
    stop(fn, "(): `direction_", which, "` must be NULL, \"higher_is_safer\" or \"higher_is_riskier\".", call. = FALSE)
  }
  if (inherits(cuts, "scr_study")) {
    for (nm in c("objective", "direction")) {
      given <- get(nm)
      if (!is.null(given) && !identical(given, cuts[[nm]])) {
        stop(fn, "(): `", nm, "_", which, "` is \"", given, "\" but the study given as `cuts_", which,
             "` was fitted with \"", cuts[[nm]], "\". Drop the argument, or pass a study of that ", nm, ".",
             call. = FALSE)
      }
    }
    return(list(objective = cuts$objective, direction = cuts$direction))
  }
  objective <- objective %||% default
  list(objective = objective, direction = resolve_direction(list(objective = objective, direction = direction)))
}

#' Cuts, band numbers and labels of one score
#'
#' From a score study (its cuts, codes and labels), from numeric cuts, or
#' tie-safe equal shares of the count table. Band 1 is the event-richest
#' unless a study says otherwise (tiers number the lowest rate first).
#' @keywords internal
#' @noRd
.cross_bands <- function(cuts, h, n_bands, side, name, fn) {
  if (inherits(cuts, "scr_study")) {
    return(list(cuts = cuts$cuts, codes = cuts$codes, labels = cuts$code_labels,
                source = if (inherits(cuts, "scr_study_tiers")) "tiers" else "bands"))
  }
  if (!is.null(cuts)) {
    if (!is.numeric(cuts) || anyNA(cuts)) stop(fn, "(): `", name, "` must be numeric cuts or a score study.", call. = FALSE)
    cc <- sort(unique(as.double(cuts[is.finite(cuts)])))
    src <- "cuts"
  } else {
    cc <- .study_cuts(h, n_bands, "uniform", NULL, side)$cuts
    src <- "equal shares"
  }
  B <- length(cc) + 1L
  list(cuts = cc, codes = if (side == "high") rev(seq_len(B)) else seq_len(B), labels = .study_labels(cc),
       source = src)
}

#' Tie-safe cut of each depth from the event-rich end; NA when no boundary exists
#' @keywords internal
#' @noRd
.cross_depth_cuts <- function(h, depths, side) {
  vapply(depths, function(d) {
    if (d >= 1) return(if (side == "high") -Inf else Inf)
    cs <- .study_cuts(h, 1L, "tail", d, side)$cuts
    if (length(cs)) cs else NA_real_
  }, numeric(1))
}

#' Smallest depth index at which each row is selected (D + 1 when never)
#'
#' The selections are nested (a deeper depth selects a superset), so the
#' count of depths that select a row is a binary search on the sorted cuts.
#' @keywords internal
#' @noRd
.cross_level <- function(s, cuts, side) {
  D <- length(cuts)
  # a depth without a boundary selects nothing
  cc <- cuts; cc[is.na(cc)] <- if (side == "high") Inf else -Inf
  hit <- if (side == "high") findInterval(s, sort(cc)) else D - findInterval(s, sort(cc))
  as.integer(D + 1L - hit)
}

#' Event rate, Jeffreys interval and lift of count columns with suffix `s`
#' @keywords internal
#' @noRd
.cross_rates <- function(d, s, overall, level, lift = TRUE) {
  ny <- d[[paste0("n_y", s)]]; e <- d[[paste0("e", s)]]
  rate <- ifelse(ny > 0, e / ny, NA_real_)
  neff <- .study_kish(ny, d[[paste0("w2", s)]])
  ci <- .study_jeffreys(rate * neff, neff, level)
  out <- list(e, rate, ci$lo, ci$hi)
  names(out) <- paste0(c("events", "rate", "rate_lo", "rate_hi"), s)
  # no lift without an overall rate (no event, or no known outcome)
  if (lift) out[[paste0("lift", s)]] <- if (isTRUE(overall > 0)) rate / overall else rep(NA_real_, length(rate))
  out
}

#' The cross table with every pair of bands and the totals
#' @keywords internal
#' @noRd
.cross_table <- function(G, cnt, ba, bb, sfx, tot, N, level) {
  # every band and pair of bands, empty ones included (built outside `[`, where ba and bb are columns)
  ga <- data.table::data.table(ba = sort(ba$codes)); gb <- data.table::data.table(bb = sort(bb$codes))
  grid <- data.table::CJ(ba = ga$ba, bb = gb$bb)
  cell <- G[, lapply(.SD, sum), keyby = c("ba", "bb"), .SDcols = cnt][grid, on = c("ba", "bb")]
  for (cn in cnt) data.table::set(cell, i = which(is.na(cell[[cn]])), j = cn, value = 0)
  row_a <- G[, lapply(.SD, sum), keyby = "ba", .SDcols = cnt][ga, on = "ba"]
  row_b <- G[, lapply(.SD, sum), keyby = "bb", .SDcols = cnt][gb, on = "bb"]
  for (d in list(row_a, row_b)) for (cn in cnt) data.table::set(d, i = which(is.na(d[[cn]])), j = cn, value = 0)
  all_ <- G[, lapply(.SD, sum), .SDcols = cnt]
  tb <- data.table::rbindlist(list(cell, row_a[, bb := NA_integer_], row_b[, ba := NA_integer_],
                                   all_[, `:=`(ba = NA_integer_, bb = NA_integer_)]), use.names = TRUE)
  lab <- function(code, b) ifelse(is.na(code), "total", b$labels[match(code, b$codes)])
  out <- data.table::data.table(band_a = as.integer(tb$ba), label_a = lab(tb$ba, ba), band_b = as.integer(tb$bb),
                                label_b = lab(tb$bb, bb), n = tb$n, pct = if (N > 0) tb$n / N else NA_real_)
  for (k in seq_along(sfx)) out <- cbind(out, data.table::as.data.table(.cross_rates(tb, sfx[k], tot[k], level)))
  out[]
}

#' Overlap counts and the rates of the five sets at every depth
#' @keywords internal
#' @noRd
.cross_overlap <- function(G, cnt, depths, da, db, sfx, ynames, N, level) {
  rows <- list(); rates <- list()
  sets <- c("A", "B", "both", "A only", "B only")
  M <- as.matrix(G[, cnt, with = FALSE])
  for (k in seq_along(depths)) {
    A <- G$la <= k; B <- G$lb <= k
    m <- list(A, B, A & B, A & !B, B & !A)
    # sums of the count columns over each set (an empty set sums to zero)
    S <- lapply(m, function(i) data.table::as.data.table(as.list(colSums(M[i, , drop = FALSE]))))
    nn <- vapply(S, function(d) d$n, numeric(1))
    rows[[k]] <- data.table::data.table(depth = depths[k], cut_a = da[k], cut_b = db[k], share_a = nn[1] / N,
                                        share_b = nn[2] / N, n_a = nn[1], n_b = nn[2], n_both = nn[3],
                                        n_a_only = nn[4], n_b_only = nn[5],
                                        jaccard = if (nn[1] + nn[2] - nn[3] > 0) nn[3] / (nn[1] + nn[2] - nn[3]) else NA_real_)
    Sd <- data.table::rbindlist(S)
    for (o in seq_along(sfx)) {
      s <- sfx[o]
      r <- .cross_rates(Sd, s, NA_real_, level, lift = FALSE)
      rates[[length(rates) + 1L]] <- data.table::data.table(depth = depths[k], outcome = ynames[o], set = sets, n = nn,
                                                            events = r[[1]], rate = r[[2]], rate_lo = r[[3]],
                                                            rate_hi = r[[4]])
    }
  }
  list(overlap = data.table::rbindlist(rows), rates = if (length(rates)) data.table::rbindlist(rates) else NULL)
}

#' @export
print.scr_score_cross <- function(x, ...) {
  st <- x$settings
  cat(sprintf("<scr_score_cross> A \"%s\" (%s, %s) | B \"%s\" (%s, %s)\n", st$score_a, st$objective_a, st$direction_a,
              st$score_b, st$objective_b, st$direction_b))
  outs <- c(st$y, st$y_a, st$y_b)
  cat(sprintf("  %s rows%s | outcome%s: %s | bands: A %s, B %s\n", .study_n(st$n_rows),
              if (st$n_dropped > 0) sprintf(" (%s left out)", .study_n(st$n_dropped)) else "",
              if (length(outs) > 1L) "s" else "", if (length(outs)) paste(outs, collapse = ", ") else "none",
              st$bands_a, st$bands_b))
  a <- x$association
  cat(sprintf("  association (raw scores): Spearman %s, Kendall tau-b %s | oriented to the event-rich ends: %s, %s\n",
              .study_f(a$estimate[1], "%.3f"), .study_f(a$estimate[2], "%.3f"), .study_f(a$oriented[1], "%.3f"),
              .study_f(a$oriented[2], "%.3f")))
  t <- x$table
  ca <- sort(st$codes_a); cb <- sort(st$codes_b)
  mat <- function(col, fmt) {
    cat(sprintf("  %-24s%s %9s\n", "A \\ B", paste(sprintf(" %9s", paste0("B", cb)), collapse = ""), "total"))
    for (i in c(ca, NA)) {
      r <- t[if (is.na(i)) is.na(t$band_a) else t$band_a %in% i]
      v <- r[[col]][match(c(cb, NA), r$band_b)]
      lab <- if (is.na(i)) "total" else sprintf("A%d %s", i, st$labels_a[match(i, st$codes_a)])
      cat(sprintf("  %-24s%s\n", substr(lab, 1, 24), paste(sprintf(" %9s", fmt(v)), collapse = "")))
    }
  }
  cat("\nShare of rows\n")
  mat("pct", function(v) .study_f(100 * v, "%.1f%%"))
  sfx <- if (!is.null(st$y)) "" else c(if (!is.null(st$y_a)) "_a", if (!is.null(st$y_b)) "_b")
  for (k in seq_along(sfx)) {
    cat(sprintf("\nEvent rate of \"%s\"\n", outs[k]))
    mat(paste0("rate", sfx[k]), .claims_pct)
  }
  cat(sprintf("  B bands: %s\n", paste(sprintf("B%d %s", cb, st$labels_b[match(cb, st$codes_b)]), collapse = ", ")))
  o <- x$overlap
  cat("\nOverlap (each score selects from its event-rich end)\n")
  cat(sprintf("  %6s %8s %8s %9s %9s %9s %9s %9s %8s\n", "depth", "share_a", "share_b", "n_a", "n_b", "both",
              "A only", "B only", "Jaccard"))
  for (i in seq_len(nrow(o))) {
    cat(sprintf("  %5.1f%% %7.1f%% %7.1f%% %9s %9s %9s %9s %9s %8s\n", 100 * o$depth[i], 100 * o$share_a[i],
                100 * o$share_b[i], .study_n(o$n_a[i]), .study_n(o$n_b[i]), .study_n(o$n_both[i]),
                .study_n(o$n_a_only[i]), .study_n(o$n_b_only[i]), .study_f(o$jaccard[i], "%.3f")))
  }
  r <- x$overlap_rates
  if (!is.null(r) && nrow(r)) {
    cat("\nEvent rate of each set (B only = swap-in, A only = swap-out)\n")
    cat(sprintf("  %6s %-12s %9s %9s %9s %9s %9s\n", "depth", "outcome", "A", "B", "both", "A only", "B only"))
    for (g in split(r, list(r$depth, r$outcome), drop = TRUE)) {
      cat(sprintf("  %5.1f%% %-12s%s\n", 100 * g$depth[1], substr(g$outcome[1], 1, 12),
                  paste(sprintf(" %9s", .claims_pct(g$rate[match(c("A", "B", "both", "A only", "B only"), g$set)])),
                        collapse = "")))
    }
  }
  invisible(x)
}

#' @rdname scr_export
#' @export
scr_export.scr_score_cross <- function(x, dir, stamp = TRUE, ...) {
  .study_dots(list(...), "scr_export")
  .need_openxlsx()
  out_dir <- .export_dir(dir, stamp)
  st <- x$settings
  tag <- paste(.file_tag(st$score_a), .file_tag(st$score_b), sep = "_")
  settings <- .kv_table(st[c("score_a", "score_b", "y", "y_a", "y_b", "objective_a", "objective_b", "direction_a",
                             "direction_b", "cuts_a", "cuts_b", "labels_a", "labels_b", "depths", "level", "n",
                             "n_rows", "n_dropped", "weighted")])
  sheets <- list(Cross_Table = x$table, Overlap = x$overlap, Overlap_Rates = x$overlap_rates,
                 Association = x$association, Settings = settings)
  sheets <- lapply(sheets[!vapply(sheets, is.null, logical(1))], .study_sheet)
  files <- list(xlsx = .scr_write_xlsx(sheets, file.path(out_dir, sprintf("score_cross_%s.xlsx", tag))))
  for (f in files) msg("  %s", f)
  x$files <- files
  invisible(x)
}

# data.table column names used without quotes in this file
utils::globalVariables(c("ba", "bb", "la", "lb"))
