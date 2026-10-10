# ============================================================================ #
# mix-shift.R - why the event rate changed: the mix of the bands or their rates
# ============================================================================ #
# Both samples are read on bands frozen on the base, from the count table of
# a score study: the decomposition is a function of the band counts, so its
# cost is that of the cells, never of the rows.
# ============================================================================ #

#' Mix and rate effects of a change in the event rate
#'
#' Splits the change of the event rate between a base and a comparison
#' sample into the part due to the mix (the population moved across the
#' score bands) and the part due to the rates (the bands themselves have a
#' different event rate), on bands frozen on the base. With `by`, every
#' period or segment is compared with the base.
#'
#' @section Decomposition:
#'
#' With \eqn{p_k} the share of band \eqn{k} among the rows with a known
#' outcome and \eqn{r_k} its event rate, on the base (\eqn{b}) and on the
#' comparison (\eqn{c}):
#' \deqn{mix_k = (p_{c,k} - p_{b,k}) (r_{b,k} + r_{c,k}) / 2, \qquad
#'       rate_k = (r_{c,k} - r_{b,k}) (p_{b,k} + p_{c,k}) / 2.}
#' The weights are the midpoints of the two samples, so the effects carry no
#' interaction term and add up exactly:
#' \eqn{\sum_k mix_k + \sum_k rate_k = R_c - R_b}, the change of the overall
#' rate. A band empty in one sample has no rate there; it takes the rate of
#' the other sample, so its rate effect is 0 and its whole contribution is
#' a mix effect (the columns `rate_base` and `rate_cmp` keep the missing
#' rate as `NA`). `total` is `mix_effect + rate_effect`.
#'
#' The summary adds, per comparison, `share_mix = mix_total / delta` (`NA`
#' when the rate did not change; it can fall outside \[0, 1\] when the two
#' effects have opposite signs) and the PSI of the band shares against the
#' base with its n-adjusted critical value at `1 - level` (see [scr_psi()]).
#'
#' @section Tests:
#'
#' `p_rate` is the two-sided p-value of the change of the rate of the band,
#' on the unweighted counts: Fisher's exact test when the smallest expected
#' count of the two-by-two table (sample by outcome) is below 5, and the
#' pooled two-proportion z test otherwise. `p_rate_adj` is Holm-adjusted
#' across the bands of the comparison. A band empty in either sample has no
#' test.
#'
#' @section Input:
#'
#' \describe{
#'   \item{A score study}{An object from [scr_bands()] or [scr_tiers()] with
#'     at least two samples: its cuts and labels are used, `base` is one of
#'     its samples (default its reference) and `compare` the others.}
#'   \item{A scorecard}{`base` and `compare` are scored samples (`"train"`
#'     against `"holdout"`). With `by`, a column of the scored samples such
#'     as `"date"`, the rows of both samples are pooled and grouped by that
#'     column; `base` is then one of its values (default the first) and
#'     `compare` the others.}
#'   \item{A data.frame}{`by` names the column that tells the groups apart (a
#'     period, a segment or a sample label); `base` is one of its values
#'     (default the first level) and `compare` the others.}
#' }
#' The values of `by` are read in the order of the levels of a factor, in
#' numeric order for a numeric column and as sorted text otherwise (dates
#' included); the first is the default base. Rows with a missing `by` value
#' are left out. Rows with a missing outcome are left out of the shares and
#' of the rates, so that the effects add up to the change of the rate.
#'
#' @inheritParams scr_bands
#' @param x An object from [scr_bands()], [scr_tiers()] or [scr_scorecard()],
#'   or a `data.frame` with one row per scored case (or one row per score
#'   value with `counts = TRUE`).
#' @param base The base: a sample or a value of `by` (see the section
#'   Input). `NULL` takes the reference of a study, `"train"` for a
#'   scorecard, and the first value of `by` otherwise.
#' @param compare The comparisons: samples or values of `by`. `NULL` takes
#'   every one but the base.
#' @param by Name of the column of periods or segments. Needed for a
#'   data.frame; for a scorecard, a column of its scored samples.
#' @param n_bands Equal-share bands frozen on the base. `NULL` uses
#'   `config$study_bands` (20) for a scorecard and 10 for a data.frame.
#' @param breaks Explicit ascending cut points; overrides `n_bands`.
#' @param level Confidence level: `1 - level` is the significance of the PSI
#'   critical value and of the count of changed bands. `NULL` uses the level
#'   of the study or `config$study_level` (0.95).
#' @param counts `TRUE` when `x` is pre-aggregated: one row per group and
#'   score value with the columns `score`, `n` and `events`.
#'
#' @return An object of class `c("scr_mix_shift", "list")`:
#'   \describe{
#'     \item{`table`}{One row per comparison and band, event-richest band
#'       first: `group`, `band`, `label`, `n_base`, `n_cmp` (rows with a
#'       known outcome), `pct_base`, `pct_cmp`, `rate_base`, `rate_cmp`,
#'       `mix_effect`, `rate_effect`, `total`, `p_rate` and `p_rate_adj`.}
#'     \item{`summary`}{One row per comparison: `group`, `n_base`, `n_cmp`,
#'       `rate_base`, `rate_cmp`, `delta`, `mix_total`, `rate_total`,
#'       `share_mix`, `psi`, `psi_critical` and `bands_changed` (bands with
#'       `p_rate_adj` below `1 - level`).}
#'     \item{`cuts`, `base`, `groups`, `by`, `level`, `objective`,
#'       `direction`, `target`, `call`}{The cuts and the settings.}
#'   }
#'
#' @references
#' Holm, S. (1979). A simple sequentially rejective multiple test procedure.
#' *Scandinavian Journal of Statistics*, 6(2), 65-70.
#'
#' Yurdakul, B. and Naranjo, J. (2020). Statistical properties of the
#' population stability index. *Journal of Risk Model Validation*, 14(4),
#' 89-100.
#'
#' @seealso [scr_bands()] for the bands, [scr_rag()] for the lights of a
#'   sample against its reference, [scr_segments()] for one score read on
#'   many segments.
#' @family score-studies
#' @examples
#' local({
#'   set.seed(1)
#'   n <- 6000
#'   period <- rep(c("2025", "2026"), each = n / 2)
#'   # 2026 has riskier applicants (mix) and a higher rate at every score (rate)
#'   x <- rnorm(n, mean = ifelse(period == "2026", -0.3, 0))
#'   d <- data.frame(period = period, score = round(600 + 50 * x),
#'                   y = rbinom(n, 1, plogis(-2 - x + 0.3 * (period == "2026"))))
#'   ms <- scr_mix_shift(d, by = "period", n_bands = 5)
#'   print(ms)
#'   ms$table[, c("band", "label", "pct_base", "pct_cmp", "mix_effect", "rate_effect")]
#' })
#' @export
scr_mix_shift <- function(x, ...) UseMethod("scr_mix_shift")

#' @rdname scr_mix_shift
#' @export
scr_mix_shift.scr_study <- function(x, base = NULL, compare = NULL, level = NULL, ...) {
  fn <- "scr_mix_shift"
  .study_dots(list(...), fn)
  base <- base %||% x$reference
  .study_chr1(base, "base", fn)
  if (!is.null(compare)) .study_chr(compare, "compare", fn)
  compare <- compare %||% setdiff(x$samples, base)
  bad <- setdiff(c(base, compare), x$samples)
  if (length(bad)) stop(fn, "(): sample(s) ", lst(bad), " not in the study (", lst(x$samples), ").", call. = FALSE)
  .mix_fit(x$hist, base, setdiff(compare, base), x$cuts, x$codes, x$code_labels, x$objective, x$direction,
           x$target, level %||% x$level %||% 0.95, NULL, sys.call())
}

#' @rdname scr_mix_shift
#' @export
scr_mix_shift.scr_scorecard <- function(x, base = NULL, compare = NULL, by = NULL, n_bands = NULL, breaks = NULL,
                                        level = NULL, max_cells = 1e5, ...) {
  fn <- "scr_mix_shift"
  .study_dots(list(...), fn)
  cfg <- x$config
  nms <- names(x$samples)
  if (!is.null(base)) .study_chr1(base, "base", fn)
  if (!is.null(compare)) .study_chr(compare, "compare", fn)
  if (is.null(by)) {
    base <- base %||% "train"
    inp <- .study_input(x, sample = compare %||% setdiff(nms, base), reference = base, max_cells = max_cells,
                        breaks = breaks, fn = fn)
    h <- inp$hist
    groups <- setdiff(inp$samples, base)
  } else {
    # the rows of both scored samples, grouped by the column `by`
    inp <- .study_input(x, sample = nms, reference = nms[1], by = by, max_cells = max_cells, breaks = breaks,
                        fn = fn)
    h <- .mix_regroup(inp$hist)
    # a numeric column in numeric order, any other by its labels
    lv <- sort(unique(h[["sample"]]))
    if (!is.null(inp$group_levels)) lv <- inp$group_levels[inp$group_levels %in% lv]
    base <- base %||% lv[1]
    compare <- compare %||% setdiff(lv, base)
    bad <- setdiff(c(base, compare), lv)
    if (length(bad)) stop(fn, "(): value(s) ", lst(bad), " not in the column '", by, "' (", lst(lv), ").", call. = FALSE)
    groups <- setdiff(compare, base)
  }
  bd <- .mix_cuts(h, base, n_bands %||% cfg$study_bands %||% 20L, breaks, inp$direction, fn)
  .mix_fit(h, base, groups, bd$cuts, bd$codes, bd$labels, inp$objective, inp$direction, inp$target,
           level %||% cfg$study_level %||% 0.95, by, sys.call())
}

#' @rdname scr_mix_shift
#' @export
scr_mix_shift.data.frame <- function(x, base = NULL, compare = NULL, by = NULL, n_bands = NULL, breaks = NULL,
                                     level = 0.95, score = "score", y = "y", objective = "risk", direction = NULL,
                                     weight = NULL, counts = FALSE, n = "n", events = "events", max_cells = 1e5,
                                     ...) {
  fn <- "scr_mix_shift"
  .study_dots(list(...), fn)
  if (is.null(by)) {
    stop(fn, "(): a data.frame needs `by`, the column that tells the base from the comparisons.", call. = FALSE)
  }
  .study_chr1(by, "by", fn)
  if (!by %in% names(x)) stop(fn, "(): column(s) ", by, " not in `x`.", call. = FALSE)
  # rows without a period or segment belong to no group. The columns read are
  # filtered one by one into a plain data.frame: `[` on a data.table would
  # evaluate the filter among its columns
  keep <- !is.na(x[[by]])
  if (!all(keep)) {
    cols <- unlist(Filter(function(v) is.character(v) && length(v) == 1L, list(score, y, weight, n, events, by)))
    cols <- intersect(unique(cols), names(x))
    x <- structure(lapply(stats::setNames(cols, cols), function(cn) x[[cn]][keep]), class = "data.frame",
                   row.names = .set_row_names(sum(keep)))
  }
  if (!nrow(x)) stop(fn, "(): `x` has no row with a value of '", by, "'.", call. = FALSE)
  if (!is.null(base)) base <- as.character(base)
  if (!is.null(compare)) compare <- as.character(compare)
  inp <- .study_input(x, score = score, y = y, objective = objective, direction = direction, weight = weight,
                      sample = by, reference = base, study = compare, counts = counts, n = n, events = events,
                      max_cells = max_cells, breaks = breaks, fn = fn)
  bd <- .mix_cuts(inp$hist, inp$reference, n_bands %||% 10L, breaks, inp$direction, fn)
  .mix_fit(inp$hist, inp$reference, setdiff(inp$samples, inp$reference), bd$cuts, bd$codes, bd$labels,
           inp$objective, inp$direction, inp$target, level, by, sys.call())
}

#' Count table of a grouped scorecard input, keyed by the group instead of the sample
#'
#' The cells of the samples are pooled within each group; rows without a
#' group are dropped.
#' @keywords internal
#' @noRd
.mix_regroup <- function(h) {
  key <- if ("g" %in% names(h)) "g" else "s"
  cnt <- intersect(.study_count_cols, names(h))
  hh <- h[!is.na(h[["group"]])]
  out <- hh[, c(list(s_lo = min(s_lo), s_hi = max(s_hi)), lapply(.SD, sum)), keyby = c("group", key), .SDcols = cnt]
  if (identical(key, "g")) data.table::set(out, j = "s", value = out[["s_lo"]])
  data.table::setnames(out, "group", "sample")
  data.table::setkeyv(out, c("sample", "s"))
  data.table::setattr(out, "edges", attr(h, "edges"))
  data.table::setattr(out, "weighted", isTRUE(attr(h, "weighted")))
  out
}

#' Cuts, band numbers and labels frozen on the base
#' @keywords internal
#' @noRd
.mix_cuts <- function(h, base, n_bands, breaks, direction, fn) {
  n_bands <- .study_whole(n_bands, "n_bands", fn, lower = 1)
  side <- .study_side(direction)
  if (!is.null(breaks)) {
    if (!is.numeric(breaks) || anyNA(breaks)) stop(fn, "(): `breaks` must be numeric without NA.", call. = FALSE)
    cuts <- sort(unique(as.double(breaks[is.finite(breaks)])))
  } else {
    ref <- .study_cells(h, base)
    if (!nrow(ref)) stop(fn, "(): the base '", base, "' has no scored row.", call. = FALSE)
    cuts <- .study_cuts(ref, n_bands, "uniform", NULL, side)$cuts
  }
  B <- length(cuts) + 1L
  list(cuts = cuts, codes = if (side == "high") rev(seq_len(B)) else seq_len(B), labels = .study_labels(cuts))
}

#' Two-sided test of a change of rate between two samples, per band
#'
#' Fisher's exact test (the sum of the hypergeometric probabilities not
#' above the observed one) when the smallest expected count is below 5, the
#' pooled two-proportion z test otherwise. Counts are rounded; `NA` when a
#' sample is empty.
#' @keywords internal
#' @noRd
.mix_rate_test <- function(e1, n1, e2, n2) {
  e1 <- round(e1); n1 <- round(n1); e2 <- round(e2); n2 <- round(n2)
  p <- rep(NA_real_, length(n1))
  ok <- n1 > 0 & n2 > 0
  N <- n1 + n2; k <- e1 + e2
  small <- ok & pmin(n1, n2) * pmin(k, N - k) / N < 5
  big <- ok & !small
  pp <- k / N
  z <- (e2 / n2 - e1 / n1) / sqrt(pp * (1 - pp) * (1 / n1 + 1 / n2))
  p[big] <- 2 * stats::pnorm(-abs(z[big]))
  for (i in which(small)) {
    lo <- max(0, k[i] - n2[i])
    d <- stats::dhyper(lo:min(k[i], n1[i]), n1[i], n2[i], k[i])
    # the relative tolerance of fisher.test(), so that ties in probability count
    p[i] <- min(1, sum(d[d <= d[e1[i] - lo + 1] * (1 + 1e-7)]))
  }
  p
}

#' Mix and rate effects of every comparison against the base
#' @keywords internal
#' @noRd
.mix_fit <- function(h, base, groups, cuts, codes, labels, objective, direction, target, level, by, call) {
  fn <- "scr_mix_shift"
  level <- .study_level(level, fn)
  if (!length(groups)) stop(fn, "(): nothing to compare with the base '", base, "'.", call. = FALSE)
  B <- length(cuts) + 1L
  # event-richest band first: the high scores under higher_is_riskier
  o <- if (identical(direction, "higher_is_riskier")) rev(seq_len(B)) else seq_len(B)
  cols <- c("n_y", "e", "w2_y", "n_y_raw", "e_raw")
  sums <- function(nm) {
    cells <- .study_cells(h, nm)
    lapply(.study_sum(cells, findInterval(cells$s, cuts) + 1L, B, cols), `[`, o)
  }
  b <- sums(base)
  NYb <- sum(b$n_y)
  if (!(NYb > 0)) stop(fn, "(): the base '", base, "' has no row with a known outcome.", call. = FALSE)
  pb <- b$n_y / NYb
  rb <- ifelse(b$n_y > 0, b$e / b$n_y, NA_real_)
  Rb <- sum(b$e) / NYb
  tabs <- list(); summ <- list()
  for (g in groups) {
    cc <- sums(g)
    NYc <- sum(cc$n_y)
    has <- NYc > 0
    pc <- if (has) cc$n_y / NYc else rep(NA_real_, B)
    rc <- ifelse(cc$n_y > 0, cc$e / cc$n_y, NA_real_)
    # a band empty in one sample takes the rate of the other (rate effect 0);
    # empty in both, it has no share and no effect
    rb2 <- ifelse(is.na(rb), rc, rb); rc2 <- ifelse(is.na(rc), rb, rc)
    rb2[is.na(rb2)] <- 0; rc2[is.na(rc2)] <- 0
    mix <- (pc - pb) * (rb2 + rc2) / 2
    rte <- (rc2 - rb2) * (pb + pc) / 2
    pr <- .mix_rate_test(b$e_raw, b$n_y_raw, cc$e_raw, cc$n_y_raw)
    pa <- .study_holm(pr)
    tabs[[g]] <- data.table::data.table(
      group = g, band = codes[o], label = labels[o], n_base = b$n_y, n_cmp = cc$n_y, pct_base = pb, pct_cmp = pc,
      rate_base = rb, rate_cmp = rc, mix_effect = mix, rate_effect = rte, total = mix + rte, p_rate = pr,
      p_rate_adj = pa)
    Rc <- if (has) sum(cc$e) / NYc else NA_real_
    delta <- Rc - Rb
    # PSI of the band shares, on the Kish effective sizes
    ps <- if (has) .psi_counts(pb * .study_kish(NYb, sum(b$w2_y)), pc * .study_kish(NYc, sum(cc$w2_y)), labels[o],
                               1 - level, c(0.10, 0.25)) else list(psi = NA_real_, critical = NA_real_)
    summ[[g]] <- data.table::data.table(
      group = g, n_base = NYb, n_cmp = NYc, rate_base = Rb, rate_cmp = Rc, delta = delta, mix_total = sum(mix),
      rate_total = sum(rte), share_mix = if (isTRUE(delta != 0)) sum(mix) / delta else NA_real_, psi = ps$psi,
      psi_critical = ps$critical, bands_changed = sum(pa < 1 - level, na.rm = TRUE))
  }
  structure(list(table = data.table::rbindlist(tabs), summary = data.table::rbindlist(summ), cuts = cuts,
                 base = base, groups = groups, by = by, level = level, objective = objective,
                 direction = direction, target = target, weighted = isTRUE(attr(h, "weighted")), call = call),
            class = c("scr_mix_shift", "list"))
}

#' A signed difference of rates in percentage points
#' @keywords internal
#' @noRd
.mix_pp <- function(x) ifelse(is.na(x), "-", sprintf("%+.2f pp", 100 * x))

#' @export
print.scr_mix_shift <- function(x, ...) {
  cat(sprintf("<scr_mix_shift> target \"%s\" | objective %s | %s\n", x$target, x$objective, x$direction))
  cat(sprintf("  base '%s'%s | %d bands frozen on the base | %d comparison%s\n", x$base,
              if (is.null(x$by)) "" else sprintf(" of '%s'", x$by), length(x$cuts) + 1L, length(x$groups),
              if (length(x$groups) == 1L) "" else "s"))
  s <- x$summary
  for (i in seq_len(nrow(s))) {
    cat(sprintf("  %-12s rate %s -> %s (%s): mix %s, rate %s | PSI %s (critical %s)\n", substr(s$group[i], 1, 12),
                .claims_pct(s$rate_base[i]), .claims_pct(s$rate_cmp[i]), .mix_pp(s$delta[i]),
                .mix_pp(s$mix_total[i]), .mix_pp(s$rate_total[i]), .study_f(s$psi[i], "%.4f"),
                .study_f(s$psi_critical[i], "%.4f")))
  }
  # the bands that moved the rate most: every comparison when few, else the largest change
  show <- if (nrow(s) <= 3L) s$group else s$group[which.max(abs(s$delta))]
  for (g in show) {
    t <- x$table[x$table$group == g]
    t <- t[utils::head(order(-abs(t$total)), 5L)]
    cat(sprintf("\nLargest band effects on '%s' (p_adj: Holm-adjusted test of the band rate)\n", g))
    cat(sprintf("  %4s %-24s %16s %18s %10s %10s %10s %7s\n", "band", "score", "share", "rate", "mix", "rate",
                "total", "p_adj"))
    for (i in seq_len(nrow(t))) {
      cat(sprintf("  %4d %-24s %16s %18s %10s %10s %10s %7s\n", t$band[i], substr(t$label[i], 1, 24),
                  sprintf("%s -> %s", .study_f(100 * t$pct_base[i], "%.1f%%"), .study_f(100 * t$pct_cmp[i], "%.1f%%")),
                  sprintf("%s -> %s", .claims_pct(t$rate_base[i]), .claims_pct(t$rate_cmp[i])),
                  .mix_pp(t$mix_effect[i]), .mix_pp(t$rate_effect[i]), .mix_pp(t$total[i]),
                  .study_f(t$p_rate_adj[i], "%.3f")))
    }
  }
  invisible(x)
}

#' @rdname scr_export
#' @export
scr_export.scr_mix_shift <- function(x, dir, stamp = TRUE, ...) {
  .study_dots(list(...), "scr_export")
  .need_openxlsx()
  out_dir <- .export_dir(dir, stamp)
  tag <- .file_tag(x$target)
  settings <- .kv_table(list(target = x$target, objective = x$objective, direction = x$direction, base = x$base,
                             groups = x$groups, by = x$by %||% "sample", level = x$level))
  cuts <- data.frame(cut = seq_along(x$cuts), score = x$cuts)
  sheets <- lapply(list(Summary = x$summary, Bands = x$table, Cuts = cuts, Settings = settings), .study_sheet)
  files <- list(xlsx = .scr_write_xlsx(sheets, file.path(out_dir, sprintf("mix_shift_%s.xlsx", tag))))
  for (f in files) msg("  %s", f)
  x$files <- files
  invisible(x)
}
