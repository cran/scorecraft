# ============================================================================ #
# cutoff.R - Stage 6: cut-off sweep, strategy table, reject inference
# ============================================================================ #

#' @keywords internal
#' @noRd
check_scorecard <- function(x, fn) {
  if (!inherits(x, "scr_scorecard")) stop(sprintf("%s() expects an object from scr_scorecard().", fn), call. = FALSE)
  invisible(TRUE)
}

#' Stage 6: cut-off sweep with frozen cuts
#'
#' For each candidate cut, what happens in each sample: the fraction of the
#' population on the safe side (approval), the event rate on both sides, the
#' events avoided (share of events falling on the risky side), the
#' non-events lost and the KS at the cut. The candidate cuts are quantiles
#' of the score **on train**, applied frozen to the hold-out: both samples
#' answer on the same numbers, and the comparison between them measures the
#' stability of the decision, not a sample difference.
#'
#' The "safe side" is the high-score side under `higher_is_safer` (credit)
#' and the low-score side under `higher_is_riskier` (fraud, propensity).
#'
#' @param x An object from [scr_scorecard()].
#' @param n_cuts Number of candidate cuts. `NULL` uses `config$cutoff_n`.
#' @param cuts Explicit vector of cuts; overrides `n_cuts`.
#'
#' @return An `scr_cutoff` object with `table` (one row per sample and cut)
#'   and `direction`.
#'
#' @family stages
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   xgb_rounds = 60, n_boot = 20)
#' res <- scr_select(scr_demo, "default", config = cfg, drop = "id",
#'                   date_col = "ref_date")
#' sc <- scr_scorecard(res)
#' ct <- scr_cutoff(sc, n_cuts = 10)
#' ct
#' st <- scr_strategy(sc, revenue_good = 1080, loss_bad = 4500)
#' st
#' rj <- scr_reject(sc)
#' rj
#' @export
scr_cutoff <- function(x, n_cuts = NULL, cuts = NULL) {
  check_scorecard(x, "scr_cutoff")
  n_cuts <- n_cuts %||% x$config$cutoff_n
  tr <- x$samples$train$score
  if (is.null(cuts)) {
    probs <- seq(0, 1, length.out = n_cuts + 2L)[-c(1L, n_cuts + 2L)]
    cuts <- unique(round(stats::quantile(tr, probs = probs, names = FALSE), 1))
  }
  dir <- x$direction
  safer <- identical(dir, "higher_is_safer")
  tb <- data.table::rbindlist(lapply(names(x$samples), function(nm) {
    s <- x$samples[[nm]]
    # one sort per sample, then every cut is a binary search on the sorted
    # scores and a lookup in the cumulative events: O(n log n + cuts log n)
    # instead of a pass over the sample per cut
    o <- order(s$score); sc <- s$score[o]; y <- as.integer(s$y[o])
    n <- length(y); e <- sum(y); ne <- n - e
    cum_e <- c(0L, cumsum(y))
    n_lo <- findInterval(cuts, sc, left.open = TRUE)   # rows with score < cut
    e_lo <- cum_e[n_lo + 1L]
    # safe side: score >= cut under higher_is_safer, score < cut otherwise
    n_safe <- if (safer) n - n_lo else n_lo
    e_safe <- if (safer) e - e_lo else e_lo
    n_risk <- n - n_safe; e_risk <- e - e_safe
    data.table::data.table(
      sample = nm, cut = cuts, n_safe = n_safe, pct_safe = n_safe / n,
      event_rate_safe = data.table::fifelse(n_safe > 0L, e_safe / n_safe, NA_real_),
      event_rate_risky = data.table::fifelse(n_risk > 0L, e_risk / n_risk, NA_real_),
      events_avoided_pct = e_risk / max(1L, e),
      nonevents_lost_pct = (n_risk - e_risk) / max(1L, ne),
      ks_at_cut = abs(e_risk / max(1L, e) - (n_risk - e_risk) / max(1L, ne)))
  }))
  structure(list(table = tb[], cuts = cuts, direction = dir, target = x$target), class = c("scr_cutoff", "list"))
}

#' @export
print.scr_cutoff <- function(x, ...) {
  cat(sprintf("<scr_cutoff> target \"%s\" | %d cuts frozen on train | safe side: %s score\n",
              x$target, length(x$cuts), if (x$direction == "higher_is_safer") "high" else "low"))
  h <- x$table[sample == "holdout"]
  if (!nrow(h)) h <- x$table[sample == x$table$sample[1]]
  cat(sprintf("  %8s %9s %10s %10s %10s %8s\n", "cut", "%safe", "ev.safe", "ev.risky", "ev.avoid", "KS"))
  for (i in seq_len(nrow(h))) cat(sprintf("  %8.1f %8.1f%% %9.2f%% %9.2f%% %9.1f%% %8.3f\n", h$cut[i], 100 * h$pct_safe[i],
                                          100 * h$event_rate_safe[i], 100 * h$event_rate_risky[i], 100 * h$events_avoided_pct[i], h$ks_at_cut[i]))
  invisible(x)
}

#' Stage 6: strategy table per band, with marginal expected profit
#'
#' Score bands (by default the deciles frozen on train) with volume, event
#' rate, the event and non-event distributions, decision and the expected
#' result per account. The good case is the non-event under
#' `objective = "risk"` (credit, fraud) and the event under
#' `"propensity"`; the bad case is the other one. With \eqn{p} the rate of
#' the bad case in the band (the event rate under risk, one minus it under
#' propensity):
#' \deqn{EP = (1 - p)\,\mathrm{revenue\_good} - p\,\mathrm{loss\_bad},}
#' which makes visible the band that is profitable **at the margin** even
#' with a high rate of the bad case. `EP = 0` at the break-even rate of the
#' bad case, `revenue_good / (revenue_good + loss_bad)`. The object stores
#' it as an event rate (`breakeven`): the same value under risk, and
#' `loss_bad / (revenue_good + loss_bad)` under propensity, where a band is
#' targeted at or above it.
#'
#' The table runs from the band richest in the good case to the poorest:
#' the safest band first under risk, the most likely first under
#' propensity.
#'
#' @section Event and non-event distributions:
#'
#' With \eqn{e_k} events and \eqn{m_k} non-events in band \eqn{k}, and
#' \eqn{E} and \eqn{M} their totals over the sample:
#' \deqn{\mathrm{pct\_event}_k = e_k / E, \qquad
#'       \mathrm{pct\_nonevent}_k = m_k / M,}
#' \deqn{\mathrm{odds\_event}_k = \mathrm{pct\_event}_k / \mathrm{pct\_nonevent}_k,
#'       \qquad \mathrm{log\_odds}_k = \ln \mathrm{odds\_event}_k.}
#' `log_odds` is the WOE of the band, event-oriented like the WOE of the
#' variables: `log_odds > 0` if and only if the band event rate is above the
#' overall event rate, that is, the lift of the band is above 1 (exact when
#' every band has both classes; under the smoothing below, a band at the
#' overall rate can fall on either side). When a band has no events or no
#' non-events, 0.5 is added to the counts of every band for `odds_event`
#' and `log_odds`; the shares stay exact. With a single class in the
#' sample, the shares of the missing class and every ratio are `NA`. This
#' `log_odds` is the `woe` column of [scr_score_gains()], not its
#' `log_odds`, which is the log of the band odds in the orientation of the
#' scale.
#'
#' @section Decision rules:
#'
#' `rule = "breakeven"` (default) gives the good label (`"approve"` under
#' risk, `"target"` under propensity) to a band whose rate of the bad case
#' is at or below break-even, `"review"` to one up to 25% above it, and the
#' bad label (`"decline"` or `"skip"`) to the rest.
#'
#' `rule = "crossing"` cuts where the event and non-event distributions are
#' furthest apart. With
#' \deqn{D_k = \left|\sum_{j \le k} \mathrm{pct\_event}_j -
#'       \sum_{j \le k} \mathrm{pct\_nonevent}_j\right|}
#' over the first \eqn{k} rows of the table, the first maximum of \eqn{D_k}
#' over the boundaries between rows is the KS of the table; the rows up to
#' it get the good label and the rest the bad label, with no review band.
#' When `log_odds` is monotone along the table this is where it changes
#' sign, the band event rate crossing the overall rate; when it is not, the
#' cut still gives a contiguous set of bands. The boundary is always
#' computed and stored in `crossing`. It is undefined with fewer than two
#' bands or a single class in the sample, and `rule = "crossing"` is then an
#' error. Scores outside `breaks` form a last row with a missing `band`,
#' which gets no decision (`NA`) under the crossing rule; the shares, and
#' hence `ks`, stay relative to the whole sample, that row included.
#'
#' `decisions`, when given, overrides either rule.
#'
#' @param x An object from [scr_scorecard()].
#' @param breaks Band cut points. `NULL` uses the deciles frozen on train.
#' @param decisions Vector of decisions, one per band (from the first row
#'   of the table to the last). `NULL` derives them from `rule`; when
#'   given, it overrides `rule`.
#' @param revenue_good Expected revenue per account of the good case (the
#'   non-event under risk, the event under propensity; default `1`).
#' @param loss_bad Expected loss per account of the bad case (default `1`;
#'   with both defaults the break-even is 50%). `revenue_good` and
#'   `loss_bad` cannot both be 0.
#' @param sample `"holdout"` (default) or `"train"`.
#' @param rule `"breakeven"` (default) or `"crossing"`; see the section
#'   Decision rules.
#'
#' @return An `scr_strategy` object with
#'   \describe{
#'     \item{`table`}{One row per band: `id`, `band`, `min_score`,
#'       `max_score`, `n`, `pct`, `events`, `event_rate`, `pct_event`,
#'       `pct_nonevent`, `odds_event`, `log_odds`, `decision`,
#'       `ep_per_account`, `band_profit`, `cum_pct`, `cum_event_rate` and
#'       `cum_profit`.}
#'     \item{`breakeven`}{The break-even event rate.}
#'     \item{`crossing`}{A list: `cut`, the score where the upper side of
#'       the crossing starts (`score >= cut`, the convention of
#'       [scr_cutoff()]), frozen on the training scores like the bands:
#'       midway between the largest training score at or below the band
#'       edge of the crossing and the smallest training score above it.
#'       `score >= cut` then reproduces the split of the bands on train and
#'       on any score seen in training; a score of another sample strictly
#'       between those two training scores can fall on the other side. When
#'       `breaks` is a number of intervals (whose edges come from `sample`),
#'       or no training score lies on one side of the edge, the cut is the
#'       midpoint between the bands on `sample`; `ks`, the distance
#'       \eqn{D_k} at it; `after_band`, the last band on the good side;
#'       `single_crossing`, whether `log_odds` changes sign exactly once
#'       along the table. All `NA` when undefined.}
#'     \item{`objective`, `rule`}{The objective of the scorecard and the
#'       rule used.}
#'     \item{`revenue_good`, `loss_bad`, `sample`, `direction`, `target`}{
#'       The parameters and the scorecard's direction and target.}
#'   }
#'
#' @family stages
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   xgb_rounds = 60, n_boot = 20)
#' res <- scr_select(scr_demo, "default", config = cfg, drop = "id",
#'                   date_col = "ref_date")
#' sc <- scr_scorecard(res)
#' scr_strategy(sc, revenue_good = 1080, loss_bad = 4500)
#' # approve down to where the event and non-event distributions cross
#' st <- scr_strategy(sc, rule = "crossing")
#' st$crossing
#' st$table[, .(band, event_rate, log_odds, decision)]
#' @export
scr_strategy <- function(x, breaks = NULL, decisions = NULL, revenue_good = 1, loss_bad = 1,
                         sample = "holdout", rule = c("breakeven", "crossing")) {
  check_scorecard(x, "scr_strategy")
  rule <- match.arg(rule)
  .scr_num1(revenue_good, "revenue_good", lower = 0); .scr_num1(loss_bad, "loss_bad", lower = 0)
  if (revenue_good + loss_bad <= 0) stop("`revenue_good` and `loss_bad` cannot both be 0: the break-even is undefined.", call. = FALSE)
  # the good case is the non-event under risk and the event under propensity
  prop <- identical(x$config$objective, "propensity")
  breaks <- breaks %||% x$breaks
  s <- x$samples[[sample]]
  if (is.null(s)) stop("sample '", sample, "' does not exist.", call. = FALSE)
  # aggregated by the integer band index; the labels are attached to the table
  b <- .score_band(s$score, breaks)
  d <- data.table::data.table(band = b$idx, y = s$y, score = s$score)[
    , .(n = .N, events = sum(y), event_rate = mean(y), min_score = min(score), max_score = max(score)), by = band]
  # band richest in the good case first: descending score when the good case sits at the high end
  d <- if (xor(identical(x$direction, "higher_is_safer"), prop)) d[order(-band)] else d[order(band)]
  idx <- d$band   # band index, kept for the score boundary of the crossing
  d[, `:=`(id = seq_len(.N), pct = n / sum(n), band = b$labels[band])]
  bw <- .band_woe(d$events, d$n - d$events)
  d[, `:=`(pct_event = bw$pct_event, pct_nonevent = bw$pct_nonevent, odds_event = bw$odds_event, log_odds = bw$log_odds)]
  # EP and break-even on the rate of the bad case; `breakeven` is stored as an event rate
  p_bad <- if (prop) 1 - d$event_rate else d$event_rate
  be_bad <- revenue_good / (revenue_good + loss_bad)
  breakeven <- if (prop) loss_bad / (revenue_good + loss_bad) else be_bad
  d[, ep_per_account := (1 - p_bad) * revenue_good - p_bad * loss_bad]
  d[, band_profit := n * ep_per_account]
  cr <- .strategy_crossing(d$pct_event, d$pct_nonevent, d$log_odds, idx, d$band, d$min_score, d$max_score,
                           breaks, x$samples$train$score)
  lab <- if (prop) c("target", "review", "skip") else c("approve", "review", "decline")
  if (is.null(decisions)) {
    if (rule == "crossing") {
      if (is.na(cr$k)) stop("scr_strategy(): rule = \"crossing\" needs at least two bands and both classes in the sample.", call. = FALSE)
      # scores outside the breaks (a last row without band) get no decision
      d[, decision := data.table::fifelse(is.na(idx), NA_character_,
                       data.table::fifelse(seq_len(.N) <= cr$k, lab[1], lab[3]))]
    } else {
      d[, decision := data.table::fifelse(p_bad <= be_bad, lab[1],
                       data.table::fifelse(p_bad <= 1.25 * be_bad, lab[2], lab[3]))]
    }
  } else {
    if (length(decisions) != nrow(d)) stop("`decisions` needs one decision per band (", nrow(d), ").", call. = FALSE)
    d[, decision := as.character(decisions)]
  }
  d[, `:=`(cum_pct = cumsum(pct), cum_event_rate = cumsum(events) / cumsum(n), cum_profit = cumsum(band_profit))]
  data.table::setcolorder(d, c("id", "band", "min_score", "max_score", "n", "pct", "events", "event_rate",
                               "pct_event", "pct_nonevent", "odds_event", "log_odds",
                               "decision", "ep_per_account", "band_profit", "cum_pct", "cum_event_rate", "cum_profit"))
  structure(list(table = d[], breakeven = breakeven, revenue_good = revenue_good, loss_bad = loss_bad,
                 sample = sample, direction = x$direction, target = x$target,
                 objective = if (prop) "propensity" else "risk", rule = rule, crossing = cr$crossing),
            class = c("scr_strategy", "list"))
}

#' Boundary where the cumulative event and non-event shares of the strategy
#' table are furthest apart (the KS of the table), in the row order given
#'
#' `idx` is the band index of every row (`NA` for scores outside the
#' breaks). The cut is frozen on the reference (training) scores `ref`, like
#' the bands: the midpoint between the largest reference score at or below
#' the band edge of the boundary (the upper edge of the lower band) and the
#' smallest one above it, so that `score >= cut` (the convention of
#' scr_cutoff()) reproduces the right-closed split on the reference and on
#' any score seen in it. When `breaks` is a number of intervals (the edges
#' come from the evaluated sample), or no reference score lies on one side
#' of the edge, the midpoint is taken on the evaluated sample: the largest
#' score of the lower band and the smallest of the upper one. `k` is the
#' last row on the good side, `NA` when the crossing is undefined.
#' @keywords internal
#' @noRd
.strategy_crossing <- function(pct_event, pct_nonevent, log_odds, idx, band, min_score, max_score,
                               breaks = NULL, ref = NULL) {
  none <- list(crossing = list(cut = NA_real_, ks = NA_real_, after_band = NA_character_, single_crossing = NA),
               k = NA_integer_)
  v <- which(!is.na(idx))
  pe <- pct_event[v]; pn <- pct_nonevent[v]
  # fewer than two bands, or a single class (shares NA): no boundary to find
  if (length(v) < 2L || anyNA(pe) || anyNA(pn)) return(none)
  D <- abs(cumsum(pe) - cumsum(pn))[-length(v)]
  k <- which.max(D)
  # rows at the edge of a right-closed band belong to the lower band: cut
  # strictly between the two bands, by band index (either table order)
  lo <- if (idx[v[k]] < idx[v[k + 1L]]) v[k] else v[k + 1L]
  hi <- if (lo == v[k]) v[k + 1L] else v[k]
  a <- max_score[lo]; b <- min_score[hi]
  # frozen like the bands: the reference scores on either side of the band
  # edge (cut() sorts the breaks; band i is (edge_i, edge_i+1])
  if (length(breaks) > 1L && length(ref)) {
    edge <- sort(as.double(breaks))[idx[lo] + 1L]
    ref <- ref[!is.na(ref)]
    if (any(ref <= edge) && any(ref > edge)) { a <- max(ref[ref <= edge]); b <- min(ref[ref > edge]) }
  }
  cut_at <- (a + b) / 2
  # adjacent doubles: the midpoint can round onto the lower score
  if (!(cut_at > a)) cut_at <- b
  sg <- sign(log_odds[v]); sg <- sg[!is.na(sg) & sg != 0]
  list(crossing = list(cut = cut_at, ks = D[k], after_band = band[v[k]],
                       single_crossing = sum(diff(sg) != 0) == 1L),
       k = v[k])
}

#' @export
print.scr_strategy <- function(x, ...) {
  prop <- identical(x$objective, "propensity")
  cat(sprintf("<scr_strategy> target \"%s\" | objective %s | rule %s | sample %s\n",
              x$target, x$objective %||% "risk", x$rule %||% "breakeven", x$sample))
  cat(sprintf("  break-even event rate: %s%.2f%% (revenue %s, loss %s)\n", if (prop) "target at or above " else "",
              100 * x$breakeven, format(x$revenue_good), format(x$loss_bad)))
  d <- x$table
  # an object saved before log_odds existed prints without that column
  lo <- if (is.null(d$log_odds)) rep("", nrow(d)) else sprintf("%8.3f ", d$log_odds)
  cat(sprintf("  %-24s %6s %8s %s%-9s %10s %12s\n", "band", "vol%", "event", if (is.null(d$log_odds)) "" else "log_odds ",
              "decision", "EP/acct", "profit"))
  for (i in seq_len(nrow(d))) cat(sprintf("  %-24s %5.1f%% %7.2f%% %s%-9s %10.2f %12.0f\n", substr(d$band[i], 1, 24), 100 * d$pct[i],
                                          100 * d$event_rate[i], lo[i], d$decision[i], d$ep_per_account[i], d$band_profit[i]))
  cr <- x$crossing
  if (!is.null(cr)) {
    if (is.na(cr$ks)) {
      cat("  event and non-event distributions: crossing undefined (fewer than two bands or a single class)\n")
    } else {
      cat(sprintf("  event and non-event distributions cross %s (KS %.3f)%s\n",
                  if (is.na(cr$cut)) paste("after band", cr$after_band) else sprintf("at score %.1f", cr$cut), cr$ks,
                  if (isFALSE(cr$single_crossing)) "; log_odds does not change sign exactly once" else ""))
    }
  }
  invisible(x)
}

#' Stage 6: honest reject inference through a sensitivity band
#'
#' Does not ship parceling as the default behavior: instead
#' of inventing a single multiplier and reweighting, it declares the
#' **population scope** of the scorecard, measures the **coverage per band**
#' (where an observed outcome exists, and in what volume) and presents a
#' **sensitivity band**: the event rate each band would have if the
#' population without an outcome were 2, 4 or 8 times worse than the
#' observed one, with the effect on the total. The analyst reads the band;
#' no single number is fabricated.
#'
#' @param x An object from [scr_scorecard()].
#' @param population Optional: a table of the full population (accepted and
#'   rejected, without outcome), scored by [scr_apply()]. `NULL` restricts
#'   the scope to the population with an outcome.
#' @param accepted Optional: a logical vector, of the length of
#'   `population`, marking the rows with an observed outcome. `NULL` treats
#'   the whole `population` as without an outcome beyond the development sample.
#' @param multipliers Sensitivity band. `NULL` uses the configuration.
#' @param sample Reference sample of the observed outcomes.
#'
#' @return An `scr_reject` object with `scope`, `coverage` (per band) and
#'   `sensitivity` (per band and multiplier, plus the `TOTAL` row).
#'
#' @family stages
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   xgb_rounds = 60, n_boot = 20)
#' res <- scr_select(scr_demo, "default", config = cfg, drop = "id",
#'                   date_col = "ref_date")
#' sc <- scr_scorecard(res)
#' scr_reject(sc)
#' # with a through-the-door population: rows with an outcome are the hold-out
#' acc <- seq_len(nrow(scr_demo)) %in% res$split$holdout_idx
#' scr_reject(sc, population = scr_demo, accepted = acc)
#' @export
scr_reject <- function(x, population = NULL, accepted = NULL, multipliers = NULL, sample = "holdout") {
  check_scorecard(x, "scr_reject")
  multipliers <- multipliers %||% x$config$reject_multipliers
  if (!is.numeric(multipliers) || !length(multipliers) || any(!is.finite(multipliers)) || any(multipliers <= 0)) {
    stop("scr_reject(): `multipliers` must be positive numbers.", call. = FALSE)
  }
  s <- x$samples[[sample]]
  if (is.null(s)) stop("sample '", sample, "' does not exist.", call. = FALSE)
  breaks <- x$breaks
  # bands as integer indices; the labels are attached to the coverage table
  b_dev <- .score_band(s$score, breaks)
  dev <- data.table::data.table(band = b_dev$idx, y = s$y)[, .(n_dev = .N, events_dev = sum(y), rate_dev = mean(y)), by = band]

  n_pop <- NA_integer_; n_unk <- 0L; pop_tb <- NULL
  if (!is.null(population)) {
    sp <- scr_apply(x, population)$score
    acc <- if (is.null(accepted)) rep(FALSE, length(sp)) else as.logical(accepted)
    if (length(acc) != length(sp)) stop("`accepted` must have the length of `population`.", call. = FALSE)
    if (anyNA(acc)) stop("`accepted` must be TRUE or FALSE on every row (no NA).", call. = FALSE)
    pop_tb <- data.table::data.table(band = .score_band(sp, breaks)$idx, acc = acc)[, .(n_pop = .N, n_unknown = sum(!acc)), by = band]
    n_pop <- length(sp); n_unk <- sum(!acc)
  }
  lv <- b_dev$labels
  cov <- data.table::data.table(band = seq_along(lv))
  cov <- merge(cov, dev, by = "band", all.x = TRUE)
  if (!is.null(pop_tb)) cov <- merge(cov, pop_tb, by = "band", all.x = TRUE) else cov[, `:=`(n_pop = NA_integer_, n_unknown = 0L)]
  for (cn in c("n_dev", "events_dev", "n_unknown")) cov[is.na(get(cn)), (cn) := 0L]
  cov[, coverage := if (all(is.na(n_pop))) NA_real_ else n_dev / pmax(1L, n_pop)]
  cov[, coverage_flag := data.table::fifelse(n_dev == 0L, "no_outcome",
                          data.table::fifelse(events_dev < 30L, "few_events", "ok"))]
  cov <- if (identical(x$direction, "higher_is_safer")) cov[order(-band)] else cov[order(band)]
  cov[, band := lv[band]]

  sens <- data.table::rbindlist(lapply(multipliers, function(m) {
    r <- data.table::copy(cov)
    r[, multiplier := m]
    r[, rate_unknown := pmin(1, rate_dev * m)]
    r[, events_implied := events_dev + n_unknown * data.table::fifelse(is.na(rate_unknown), 0, rate_unknown)]
    r[, rate_implied := events_implied / pmax(1L, n_dev + n_unknown)]
    tot <- data.table::data.table(band = "TOTAL", n_dev = sum(r$n_dev), events_dev = sum(r$events_dev),
                                  rate_dev = sum(r$events_dev) / max(1L, sum(r$n_dev)), n_pop = sum(r$n_pop),
                                  n_unknown = sum(r$n_unknown), coverage = NA_real_, coverage_flag = "",
                                  multiplier = m, rate_unknown = NA_real_, events_implied = sum(r$events_implied),
                                  rate_implied = sum(r$events_implied) / max(1L, sum(r$n_dev + r$n_unknown)))
    data.table::rbindlist(list(r, tot), use.names = TRUE, fill = TRUE)
  }))
  scope <- list(
    n_with_outcome = nrow(s), n_population = n_pop, n_without_outcome = n_unk,
    share_with_outcome = if (is.na(n_pop)) NA_real_ else nrow(s) / n_pop,
    statement = if (is.null(population))
      "The scorecard describes the population WITH an observed outcome. No extrapolation to rejects was made; the sensitivity band shows the effect of declared assumptions, not an inferred number."
    else sprintf("The full population has %s rows, of which %s (%.1f%%) have an observed outcome. The rest enter only the sensitivity band, under declared multipliers.",
                 n_fmt(n_pop), n_fmt(n_pop - n_unk), 100 * (n_pop - n_unk) / n_pop))
  structure(list(scope = scope, coverage = cov[, .(band, n_dev, events_dev, rate_dev, n_pop, n_unknown, coverage, coverage_flag)],
                 sensitivity = sens[, .(multiplier, band, n_dev, events_dev, rate_dev, n_unknown, rate_unknown, events_implied, rate_implied)],
                 multipliers = multipliers, target = x$target), class = c("scr_reject", "list"))
}

#' @export
print.scr_reject <- function(x, ...) {
  cat(sprintf("<scr_reject> target \"%s\" | multipliers %s\n", x$target, paste0(x$multipliers, "x", collapse = ", ")))
  cat("  ", x$scope$statement, "\n", sep = "")
  tot <- x$sensitivity[band == "TOTAL"]
  cat(sprintf("  observed event rate: %.2f%%\n", 100 * tot$rate_dev[1]))
  for (i in seq_len(nrow(tot))) cat(sprintf("  implied rate if the population without outcome is %gx worse: %.2f%%\n", tot$multiplier[i], 100 * tot$rate_implied[i]))
  cf <- x$coverage[coverage_flag != "ok"]
  if (nrow(cf)) cat(sprintf("  bands with weak coverage: %s\n", lst(paste0(cf$band, " (", cf$coverage_flag, ")"))))
  invisible(x)
}
