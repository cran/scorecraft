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
#' rate, decision and the expected result per account:
#' \deqn{EP = (1 - p)\,\mathrm{revenue\_good} - p\,\mathrm{loss\_bad},}
#' which makes visible the band that is profitable **at the margin** even
#' with a high event rate. The break-even event rate, where `EP = 0`, is
#' `revenue_good / (revenue_good + loss_bad)`.
#'
#' The automatic decision approves a band whose rate is below break-even,
#' sends to review a band up to 25% above it and declines the rest; pass
#' `decisions` to fix the policy.
#'
#' @param x An object from [scr_scorecard()].
#' @param breaks Band cut points. `NULL` uses the deciles frozen on train.
#' @param decisions Vector of decisions, one per band (from the safest to
#'   the riskiest). `NULL` derives them from break-even.
#' @param revenue_good Expected revenue per account without the event
#'   (default `1`).
#' @param loss_bad Expected loss per account with the event (default `1`;
#'   with both defaults the break-even event rate is 50%).
#' @param sample `"holdout"` (default) or `"train"`.
#'
#' @return An `scr_strategy` object with `table`, `breakeven` and the parameters.
#'
#' @family stages
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   xgb_rounds = 60, n_boot = 20)
#' res <- scr_select(scr_demo, "default", config = cfg, drop = "id",
#'                   date_col = "ref_date")
#' sc <- scr_scorecard(res)
#' scr_strategy(sc, revenue_good = 1080, loss_bad = 4500)
#' @export
scr_strategy <- function(x, breaks = NULL, decisions = NULL, revenue_good = 1, loss_bad = 1,
                         sample = "holdout") {
  check_scorecard(x, "scr_strategy")
  .scr_num1(revenue_good, "revenue_good", lower = 0); .scr_num1(loss_bad, "loss_bad", lower = 0)
  breaks <- breaks %||% x$breaks
  s <- x$samples[[sample]]
  if (is.null(s)) stop("sample '", sample, "' does not exist.", call. = FALSE)
  band <- cut(s$score, breaks = breaks, include.lowest = TRUE)
  d <- data.table::data.table(band = band, y = s$y, score = s$score)[
    , .(n = .N, events = sum(y), event_rate = mean(y), min_score = min(score), max_score = max(score)), by = band]
  d <- if (identical(x$direction, "higher_is_safer")) d[order(-as.integer(band))] else d[order(band)]
  d[, `:=`(id = seq_len(.N), pct = n / sum(n), band = as.character(band))]
  breakeven <- revenue_good / (revenue_good + loss_bad)
  d[, ep_per_account := (1 - event_rate) * revenue_good - event_rate * loss_bad]
  d[, band_profit := n * ep_per_account]
  if (is.null(decisions)) {
    d[, decision := data.table::fifelse(event_rate <= breakeven, "approve",
                     data.table::fifelse(event_rate <= 1.25 * breakeven, "review", "decline"))]
  } else {
    if (length(decisions) != nrow(d)) stop("`decisions` needs one decision per band (", nrow(d), ").", call. = FALSE)
    d[, decision := as.character(decisions)]
  }
  d[, `:=`(cum_pct = cumsum(pct), cum_event_rate = cumsum(events) / cumsum(n), cum_profit = cumsum(band_profit))]
  data.table::setcolorder(d, c("id", "band", "min_score", "max_score", "n", "pct", "events", "event_rate",
                               "decision", "ep_per_account", "band_profit", "cum_pct", "cum_event_rate", "cum_profit"))
  structure(list(table = d[], breakeven = breakeven, revenue_good = revenue_good, loss_bad = loss_bad,
                 sample = sample, direction = x$direction, target = x$target), class = c("scr_strategy", "list"))
}

#' @export
print.scr_strategy <- function(x, ...) {
  cat(sprintf("<scr_strategy> target \"%s\" | sample %s | break-even event rate: %.2f%% (revenue %s, loss %s)\n",
              x$target, x$sample, 100 * x$breakeven, format(x$revenue_good), format(x$loss_bad)))
  d <- x$table
  cat(sprintf("  %-24s %6s %8s %-9s %10s %12s\n", "band", "vol%", "event", "decision", "EP/acct", "profit"))
  for (i in seq_len(nrow(d))) cat(sprintf("  %-24s %5.1f%% %7.2f%% %-9s %10.2f %12.0f\n", substr(d$band[i], 1, 24), 100 * d$pct[i],
                                          100 * d$event_rate[i], d$decision[i], d$ep_per_account[i], d$band_profit[i]))
  invisible(x)
}

#' Stage 6: honest reject inference through a sensitivity band
#'
#' Does not ship parcelling as the default behaviour: instead
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
  band_dev <- cut(s$score, breaks = breaks, include.lowest = TRUE)
  dev <- data.table::data.table(band = band_dev, y = s$y)[, .(n_dev = .N, events_dev = sum(y), rate_dev = mean(y)), by = band]

  n_pop <- NA_integer_; n_unk <- 0L; pop_tb <- NULL
  if (!is.null(population)) {
    sp <- scr_apply(x, population)$score
    acc <- if (is.null(accepted)) rep(FALSE, length(sp)) else as.logical(accepted)
    if (length(acc) != length(sp)) stop("`accepted` must have the length of `population`.", call. = FALSE)
    if (anyNA(acc)) stop("`accepted` must be TRUE or FALSE on every row (no NA).", call. = FALSE)
    band_pop <- cut(sp, breaks = breaks, include.lowest = TRUE)
    pop_tb <- data.table::data.table(band = band_pop, acc = acc)[, .(n_pop = .N, n_unknown = sum(!acc)), by = band]
    n_pop <- length(sp); n_unk <- sum(!acc)
  }
  lv <- levels(band_dev)
  cov <- data.table::data.table(band = factor(lv, levels = lv))
  cov <- merge(cov, dev, by = "band", all.x = TRUE)
  if (!is.null(pop_tb)) cov <- merge(cov, pop_tb, by = "band", all.x = TRUE) else cov[, `:=`(n_pop = NA_integer_, n_unknown = 0L)]
  for (cn in c("n_dev", "events_dev", "n_unknown")) cov[is.na(get(cn)), (cn) := 0L]
  cov[, coverage := if (all(is.na(n_pop))) NA_real_ else n_dev / pmax(1L, n_pop)]
  cov[, coverage_flag := data.table::fifelse(n_dev == 0L, "no_outcome",
                          data.table::fifelse(events_dev < 30L, "few_events", "ok"))]
  cov <- if (identical(x$direction, "higher_is_safer")) cov[order(-as.integer(band))] else cov[order(band)]
  cov[, band := as.character(band)]

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
