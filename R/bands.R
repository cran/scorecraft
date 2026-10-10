# ============================================================================ #
# bands.R - percentile study: bands frozen on the reference, read on a study
# ============================================================================ #

#' Percentile study of a score
#'
#' Cuts the score into bands of equal share (or into tail percentiles)
#' **frozen on the reference sample**, and reads every band on the reference
#' and on the study samples: volume, event rate with a Jeffreys interval,
#' lift, capture (recall), the cumulative non-event share (false positive
#' rate), KS, the band WOE and IV, odds, the PSI term against the reference
#' and a one-sided Fisher exact test of rank order against the previous band.
#' The summary adds the AUC, Gini and KS of every sample with a bootstrap
#' interval.
#'
#' @section Method:
#'
#' The scored rows are aggregated once into a table of counts per distinct
#' score (per sample); every statistic is then computed from that table, so
#' the cost is one pass over the rows plus work proportional to the number
#' of distinct scores. With more than `max_cells` distinct scores, the
#' scores are first pooled into `max_cells` cells of equal weighted share.
#'
#' The cut targets are cumulative shares counted from the event-rich side
#' of the score (the high scores under `higher_is_riskier`, the low ones
#' under `higher_is_safer`): `k / n_bands` with `spacing = "uniform"`, or
#' the shares in `tail_probs` with `spacing = "tail"` (default 0.1%, 0.5%,
#' 1%, 2%, 5%, 10%, 20% and 50%). Each cut is placed midway between two
#' adjacent distinct reference scores, at the boundary nearest to its
#' target, so a group of tied scores is never split and a target is hit
#' within the share of one score value. Targets that land on the same
#' boundary give one cut: `n_bands_effective` can be smaller than
#' `n_bands_requested`, and both are reported. A band is left-closed,
#' `[lo, hi)`: `score >= cut` is the upper side, the convention of
#' [scr_cutoff()]. `breaks` given explicitly are used as they are; the bands
#' of [scr_score_gains()] (`breaks = sc$breaks`) are right-closed, so a
#' score equal to a break falls one band higher here.
#'
#' The band table lists the event-richest band first (`band = 1`). With
#' \eqn{e_b} events and \eqn{m_b} non-events in band \eqn{b}, totals \eqn{E}
#' and \eqn{M}, and overall rate \eqn{R}:
#' \itemize{
#'   \item `rate` = \eqn{e_b / (e_b + m_b)}, with the Jeffreys interval
#'     `rate_lo`, `rate_hi`: the Beta(\eqn{e_b + 1/2, m_b + 1/2}) quantiles,
#'     0 and 1 at the edges (Brown, Cai and DasGupta, 2001). Under weights,
#'     the counts are scaled to the Kish effective size
#'     \eqn{(\sum w)^2 / \sum w^2}.
#'   \item `lift` = `rate` / \eqn{R} (its interval divides the rate bounds by
#'     \eqn{R}); `cum_rate` and `cum_lift` accumulate from the first band.
#'   \item `capture` = \eqn{\sum_{j \le b} e_j / E} (recall),
#'     `cum_nonevent_pct` = \eqn{\sum_{j \le b} m_j / M} (false positive
#'     rate) and `ks` their absolute difference.
#'   \item `pct_event`, `pct_nonevent` and `woe` = \eqn{\ln(e_b / E) -
#'     \ln(m_b / M)} as in [scr_strategy()] (0.5 is added to every band only
#'     when a band lacks events or non-events); `iv` =
#'     (`pct_event` - `pct_nonevent`) * `woe`.
#'   \item `odds` and `log_odds` in the orientation of the scale
#'     (non-events per event under `higher_is_safer`, events per non-event
#'     under `higher_is_riskier`, 0.5 added to each count), as in
#'     [scr_score_gains()].
#'   \item `psi`: the band term of the PSI against the reference shares
#'     (see [scr_psi()]); `NA` on the reference itself.
#'   \item `p_reversal`: one-sided Fisher exact test that the band has a
#'     **higher** event rate than the previous, event-richer band (a
#'     reversal of the rank order); `p_reversal_adj` is Holm-adjusted over
#'     the bands of the sample. The tests use the unweighted counts.
#' }
#'
#' Rows with a missing or infinite score, or a zero weight, are not counted;
#' rows with a missing outcome count in the volume (`n`, `pct`, the PSI) but
#' not in the rates.
#'
#' The bootstrap of the AUC draws the event and non-event counts of every
#' score value from multinomial laws with the observed shares (the law of a
#' row bootstrap stratified by outcome), with the unweighted class counts as
#' sizes; a given `seed` is local to the call, while `seed = NULL` draws
#' from, and advances, the user's random stream. Up to `boot_cells` cells of
#' the count table the bootstrap is exact. Above it, the resamples run on
#' `boot_cells` cells of equal share (adjacent cells pooled) and are shifted
#' to the point estimate of the full table, so the interval is an
#' approximation. The point estimates use every cell of the count table:
#' every distinct score, unless `max_cells` pooled the scores into cells.
#' `boot_cells = Inf` keeps the bootstrap exact, at a cost per resample
#' proportional to the number of cells.
#'
#' @param x An object from [scr_scorecard()], or a `data.frame` with one row
#'   per scored case (or one row per score value with `counts = TRUE`).
#' @param ... Passed on to the methods; an unknown argument is an error.
#' @param n_bands Number of equal-share bands. For a scorecard, `NULL` uses
#'   `config$study_bands` (20).
#' @param spacing `"uniform"` (equal shares) or `"tail"` (the shares of
#'   `tail_probs`, counted from the event-rich side).
#' @param tail_probs Cumulative shares from the event-rich side for
#'   `spacing = "tail"`. `NULL` uses 0.001, 0.005, 0.01, 0.02, 0.05, 0.10,
#'   0.20 and 0.50.
#' @param sample For a scorecard: the study sample(s), `"holdout"` (default)
#'   and/or `"train"`. For a data.frame: the name of a column with sample
#'   labels, or `NULL` (all rows are one sample, reference and study at once).
#' @param reference For a scorecard: the sample the bands are frozen on
#'   (`"train"`). For a data.frame: the label of the reference sample;
#'   `NULL` takes the first level of the `sample` column (the levels of a
#'   factor in their order, numbers in numeric order, text sorted).
#' @param breaks Explicit ascending cut points; overrides `n_bands` and
#'   `spacing`. Infinite values are ignored.
#' @param level Confidence level of the intervals. For a scorecard, `NULL`
#'   uses `config$study_level` (0.95).
#' @param n_boot Bootstrap resamples of the AUC interval (`0` skips it). For
#'   a scorecard, `NULL` uses `config$n_boot`.
#' @param seed Seed of the bootstrap. A number is local to the call (the
#'   user's random stream is restored on exit); `NULL` draws from the user's
#'   stream and advances it. For a scorecard, `NULL` uses `config$seed`.
#' @param max_cells Largest number of distinct score values kept exactly.
#' @param boot_cells Largest number of score cells resampled exactly by the
#'   bootstrap (default 10,000; `Inf` for no pooling). See the section
#'   Method.
#' @param score,y Column names of the score and of the 0/1 outcome (`NA`
#'   allowed).
#' @param objective `"risk"` (the event is the bad case) or `"propensity"`
#'   (the event is the good case).
#' @param direction `"higher_is_safer"` or `"higher_is_riskier"`; `NULL`
#'   derives it from `objective`.
#' @param weight Optional column of non-negative case weights.
#' @param value Optional column of a value per case (an amount, a balance):
#'   adds the value captured per band. With `counts = TRUE`, the value per
#'   score cell.
#' @param study Labels of the study samples; `NULL` takes every level other
#'   than the reference.
#' @param counts `TRUE` when `x` is pre-aggregated: one row per score value
#'   with the columns `score`, `n` and `events` (and, optionally, `value`
#'   and `value_events`).
#' @param n,events Column names of the counts when `counts = TRUE`.
#' @param value_events With `counts = TRUE`: the column of the value of the
#'   events per score cell.
#'
#' @return An object of class `c("scr_study_bands", "scr_study", "list")`:
#'   \describe{
#'     \item{`table`}{One row per sample and band, event-richest band first:
#'       `sample`, `band`, `label`, `score_lo`, `score_hi`, `n`, `pct`,
#'       `cum_pct`, `events`, `rate`, `rate_lo`, `rate_hi`, `cum_rate`,
#'       `lift`, `lift_lo`, `lift_hi`, `cum_lift`, `capture`,
#'       `cum_nonevent_pct`, `ks`, `pct_event`, `pct_nonevent`, `woe`, `iv`,
#'       `odds`, `log_odds`, `psi`, `p_reversal`, `p_reversal_adj` and, with
#'       `value`, `value`, `value_events`, `value_capture` (cumulative share
#'       of the event value) and `value_precision` (cumulative event value
#'       over cumulative value).}
#'     \item{`summary`}{One row per sample: `sample`, `n`, `events`, `rate`,
#'       `auc`, `auc_lo`, `auc_hi`, `gini`, `gini_lo`, `gini_hi`, `ks`, `iv`,
#'       `psi` (against the reference), `n_bands_requested`,
#'       `n_bands_effective` and `reversals` (bands with `p_reversal_adj`
#'       below 0.05).}
#'     \item{`cuts`}{The ascending cut points.}
#'     \item{`codes`, `code_labels`}{Band number and label of every interval
#'       in ascending score order, used by [scr_apply()] and [scr_sql()].}
#'     \item{`objective`, `direction`, `level`, `spacing`, `reference`,
#'       `samples`, `target`, `call`}{The settings.}
#'     \item{`n_bands_requested`, `n_bands_effective`, `quantized`, `hist`}{
#'       The band counts, whether the scores were pooled, and the count
#'       table the study was computed from.}
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
#' Siddiqi, N. (2006). *Credit Risk Scorecards: Developing and Implementing
#' Intelligent Credit Scoring*. Wiley.
#'
#' Yurdakul, B. and Naranjo, J. (2020). Statistical properties of the
#' population stability index. *Journal of Risk Model Validation*, 14(4),
#' 89-100.
#'
#' @seealso [scr_tiers()] for a small number of policy tiers, [scr_rag()]
#'   for traffic lights, [scr_apply()] and [scr_sql()] to assign the bands in
#'   production.
#' @family score-studies
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   use_lightgbm = FALSE, xgb_rounds = 40, n_boot = 20)
#' res <- scr_select(scr_demo, "default", config = cfg, drop = c("id", "churn"),
#'                   date_col = "ref_date")
#' sc <- scr_scorecard(res)
#' b <- scr_bands(sc, n_bands = 10)
#' b
#' b$table[sample == "holdout", .(band, label, n, rate, lift, capture, ks)]
#'
#' # tail percentiles, from a data.frame
#' d <- data.frame(score = sc$samples$holdout$score, y = sc$samples$holdout$y)
#' scr_bands(d, spacing = "tail", n_boot = 0)$table[, .(band, label, pct, rate, capture)]
#' @export
scr_bands <- function(x, ...) UseMethod("scr_bands")

#' @rdname scr_bands
#' @export
scr_bands.scr_scorecard <- function(x, n_bands = NULL, spacing = c("uniform", "tail"), tail_probs = NULL,
                                    sample = "holdout", reference = "train", breaks = NULL, level = NULL,
                                    n_boot = NULL, seed = NULL, max_cells = 1e5, boot_cells = 1e4, ...) {
  .study_dots(list(...), "scr_bands")
  cfg <- x$config
  inp <- .study_input(x, sample = sample, reference = reference, max_cells = max_cells, breaks = breaks,
                      fn = "scr_bands")
  .bands_fit(inp, n_bands %||% cfg$study_bands %||% 20L, match.arg(spacing), tail_probs, breaks,
             level %||% cfg$study_level %||% 0.95, n_boot %||% cfg$n_boot, seed %||% cfg$seed, sys.call(),
             boot_cells)
}

#' @rdname scr_bands
#' @export
scr_bands.data.frame <- function(x, score = "score", y = "y", objective = "risk", direction = NULL,
                                 weight = NULL, value = NULL, sample = NULL, reference = NULL, study = NULL,
                                 counts = FALSE, n = "n", events = "events", n_bands = 20L,
                                 spacing = c("uniform", "tail"), tail_probs = NULL, breaks = NULL,
                                 level = 0.95, n_boot = 200L, seed = NULL, max_cells = 1e5,
                                 value_events = NULL, boot_cells = 1e4, ...) {
  .study_dots(list(...), "scr_bands")
  inp <- .study_input(x, score = score, y = y, objective = objective, direction = direction, weight = weight,
                      value = value, value_events = value_events, sample = sample, reference = reference,
                      study = study, counts = counts, n = n, events = events, max_cells = max_cells,
                      breaks = breaks, fn = "scr_bands")
  .bands_fit(inp, n_bands, match.arg(spacing), tail_probs, breaks, level, n_boot, seed, sys.call(), boot_cells)
}

#' Bands frozen on the reference, read on every sample
#' @keywords internal
#' @noRd
.bands_fit <- function(inp, n_bands, spacing, tail_probs, breaks, level, n_boot, seed, call, boot_cells = 1e4) {
  fn <- "scr_bands"
  n_bands <- .study_whole(n_bands, "n_bands", fn, lower = 1)
  n_boot <- .study_whole(n_boot, "n_boot", fn)
  boot_cells <- .study_boot_cells(boot_cells, fn)
  level <- .study_level(level, fn)
  h <- inp$hist
  ref_cells <- .study_cells(h, inp$reference)
  if (!nrow(ref_cells)) stop(fn, "(): the reference sample '", inp$reference, "' has no scored row.", call. = FALSE)
  side <- .study_side(inp$direction)
  if (!is.null(breaks)) {
    if (!is.numeric(breaks) || anyNA(breaks)) stop(fn, "(): `breaks` must be numeric without NA.", call. = FALSE)
    cuts <- sort(unique(as.double(breaks[is.finite(breaks)])))
    cs <- list(cuts = cuts, n_bands_requested = length(cuts) + 1L, n_bands_effective = length(cuts) + 1L)
  } else {
    cs <- .study_cuts(ref_cells, n_bands, spacing, tail_probs, side)
  }
  cuts <- cs$cuts
  B <- length(cuts) + 1L
  ref_sh <- .study_ref_shares(ref_cells, cuts)
  # band number of every ascending interval: 1 is the event-richest band
  codes <- if (side == "high") rev(seq_len(B)) else seq_len(B)
  .scr_local_seed(seed)
  tabs <- list(); summ <- list()
  for (nm in inp$samples) {
    cells <- .study_cells(h, nm)
    tb <- .study_band_table(cells, cuts, if (!identical(nm, inp$reference)) ref_sh, level, inp$direction)
    ps <- attr(tb, "psi")
    d <- .study_discrimination(cells, inp$direction, n_boot, level, boot_cells = boot_cells)
    tabs[[nm]] <- data.table::data.table(sample = nm, tb)
    NY <- sum(cells$n_y)
    summ[[nm]] <- data.table::data.table(
      sample = nm, n = sum(cells$n), events = sum(cells$e), rate = if (NY > 0) sum(cells$e) / NY else NA_real_,
      auc = d$auc, auc_lo = d$auc_lo, auc_hi = d$auc_hi, gini = d$gini, gini_lo = d$gini_lo, gini_hi = d$gini_hi,
      ks = d$ks, iv = if (all(is.na(tb$iv))) NA_real_ else sum(tb$iv, na.rm = TRUE),
      psi = if (is.null(ps)) NA_real_ else ps$psi,
      n_bands_requested = cs$n_bands_requested, n_bands_effective = cs$n_bands_effective,
      reversals = sum(tb$p_reversal_adj < 0.05, na.rm = TRUE))
  }
  structure(list(
    table = data.table::rbindlist(tabs), summary = data.table::rbindlist(summ), cuts = cuts,
    codes = codes, code_labels = .study_labels(cuts), objective = inp$objective, direction = inp$direction,
    level = level, spacing = if (is.null(breaks)) spacing else "breaks", reference = inp$reference,
    samples = inp$samples, target = inp$target, n_bands_requested = cs$n_bands_requested,
    n_bands_effective = cs$n_bands_effective, n_boot = n_boot, boot_cells = boot_cells,
    quantized = inp$meta$quantized,
    weighted = inp$meta$weighted, sql_table = inp$sql_table, sql_dialect = inp$sql_dialect,
    hist = h, call = call),
    class = c("scr_study_bands", "scr_study", "list"))
}

#' @export
print.scr_study_bands <- function(x, ...) {
  cat(sprintf("<scr_study_bands> target \"%s\" | objective %s | %s\n", x$target, x$objective, x$direction))
  study <- setdiff(x$samples, x$reference)
  if (!length(study)) study <- x$reference
  cat(sprintf("  bands frozen on '%s', read on %s | %d requested, %d effective (%s)%s\n", x$reference,
              paste0("'", study, "'", collapse = ", "), x$n_bands_requested, x$n_bands_effective, x$spacing,
              if (isTRUE(x$quantized)) " | scores pooled into cells" else ""))
  .study_print_summary(x$summary, x$level)
  for (nm in study) {
    t <- x$table[x$table$sample == nm]
    cat(sprintf("\nBands on '%s' (event-richest first; %.0f%% Jeffreys interval of the rate)\n", nm, 100 * x$level))
    cat(sprintf("  %4s %-24s %7s %-26s %7s %8s %6s %9s\n", "band", "score", "pct", "rate [lo, hi]", "lift",
                "capture", "KS", "p_rev_adj"))
    for (i in seq_len(nrow(t))) {
      # the rate with its interval; p_rev_adj is the Holm-adjusted p-value of the reversal test
      ri <- if (is.na(t$rate[i])) "-" else if (is.na(t$rate_lo[i])) sprintf("%.2f%%", 100 * t$rate[i]) else
        sprintf("%.2f%% [%.2f%%, %.2f%%]", 100 * t$rate[i], 100 * t$rate_lo[i], 100 * t$rate_hi[i])
      cat(sprintf("  %4d %-24s %7s %-26s %7s %8s %6s %9s\n", t$band[i], substr(t$label[i], 1, 24),
                  .study_f(100 * t$pct[i], "%.1f%%"), ri, .study_f(t$lift[i], "%.2f"),
                  .study_f(100 * t$capture[i], "%.1f%%"), .study_f(t$ks[i], "%.3f"),
                  .study_f(t$p_reversal_adj[i], "%.3f")))
    }
  }
  invisible(x)
}

#' Summary lines shared by the band and tier printers
#' @keywords internal
#' @noRd
.study_print_summary <- function(s, level) {
  has_auc <- "auc" %in% names(s)
  cat(sprintf("  %-10s %9s %9s %8s %s\n", "sample", "n", "events", "rate",
              if (has_auc) sprintf("AUC [%.0f%% CI]           Gini     KS     IV     PSI  reversals", 100 * level)
              else "    IV     PSI  tiers  monotone  distinct"))
  for (i in seq_len(nrow(s))) {
    r <- s[i]
    head <- sprintf("  %-10s %9s %9s %8s ", substr(r$sample, 1, 10), .study_n(r$n), .study_n(r$events),
                    .study_f(100 * r$rate, "%.2f%%"))
    tail <- if (has_auc) {
      sprintf("%s %-14s %6s %6s %6s %7s %9d", .study_f(r$auc, "%.4f", "  -   "),
              if (is.na(r$auc_lo)) "" else sprintf("[%.3f, %.3f]", r$auc_lo, r$auc_hi),
              .study_f(r$gini, "%.3f"), .study_f(r$ks, "%.3f"), .study_f(r$iv, "%.3f"), .study_f(r$psi, "%.4f"),
              r$reversals)
    } else {
      sprintf("%6s %7s %6d %9s %9s", .study_f(r$iv, "%.3f"), .study_f(r$psi, "%.4f"), r$n_tiers,
              if (is.na(r$monotone)) "-" else if (r$monotone) "yes" else "no",
              if (is.na(r$all_distinct)) "-" else if (r$all_distinct) "yes" else "no")
    }
    cat(head, tail, "\n", sep = "")
  }
}

# data.table column names used without quotes in this file
utils::globalVariables(c("band", "label", "rate", "lift", "capture", "ks"))
