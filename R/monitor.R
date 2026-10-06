# ============================================================================ #
# monitor.R - PSI/CSI over time (exported, never runs by itself)
# ============================================================================ #

#' Monitoring plan read by scr_monitor()
#'
#' A small `item`/`value` table with the thresholds and the frozen score
#' bands of a scorecard. It is created by [scr_scorecard()] from the
#' configuration, written to the `Monitoring_Plan` sheet of the strategy
#' workbook by [scr_export()], and **read back** by [scr_monitor()]: change a
#' threshold in the sheet, pass the file (or the edited table) as `plan`, and
#' the flags follow the plan, not the configuration.
#'
#' @param x An object from [scr_scorecard()], or a configuration from
#'   [scr_config()] plus `breaks`.
#' @param breaks Frozen score bands, when `x` is a configuration.
#'
#' @return A `data.frame` of class `scr_monitoring_plan` with the items
#'   `psi_score_fixed_moderate`, `psi_score_fixed_action`, `psi_adjusted_alpha`,
#'   `csi_variable_fixed_moderate`, `csi_variable_fixed_action`,
#'   `score_bands`, `min_events_per_period` and `threshold_source`.
#'
#' @family production
#' @examples
#' plan <- scr_monitoring_plan(scr_config(), breaks = c(-Inf, 500, 550, 600, Inf))
#' plan
#' @export
scr_monitoring_plan <- function(x, breaks = NULL) {
  if (inherits(x, "scr_scorecard")) { cfg <- x$config; breaks <- x$breaks }
  else if (inherits(x, "scr_config")) { cfg <- x; if (is.null(breaks)) stop("`breaks` is needed with a configuration.", call. = FALSE) }
  else stop("scr_monitoring_plan() expects an scr_scorecard or an scr_config.", call. = FALSE)
  d <- data.frame(
    item = c("psi_score_fixed_moderate", "psi_score_fixed_action", "psi_adjusted_alpha",
             "csi_variable_fixed_moderate", "csi_variable_fixed_action", "score_bands",
             "min_events_per_period", "threshold_source"),
    value = c("0.10", "0.25", as.character(cfg$psi_alpha), "0.10", "0.25",
              paste(round(breaks[is.finite(breaks)], 2), collapse = " | "), "100",
              "0.10/0.25: market convention, no published authority; adjusted: Yurdakul & Naranjo (2020)"),
    stringsAsFactors = FALSE)
  class(d) <- c("scr_monitoring_plan", "data.frame")
  d
}

#' Read a monitoring plan from a table or from the strategy workbook
#' @keywords internal
#' @noRd
.read_plan <- function(plan) {
  if (is.character(plan) && length(plan) == 1L && file.exists(plan)) {
    .need_openxlsx()
    sheets <- openxlsx::getSheetNames(plan)
    if (!"Monitoring_Plan" %in% sheets) stop("no 'Monitoring_Plan' sheet in ", plan, call. = FALSE)
    plan <- openxlsx::read.xlsx(plan, sheet = "Monitoring_Plan")
  }
  plan <- as.data.frame(plan, stringsAsFactors = FALSE)
  if (!all(c("item", "value") %in% names(plan))) stop("a monitoring plan needs `item` and `value` columns.", call. = FALSE)
  get <- function(nm, default) {
    v <- plan$value[match(nm, plan$item)]
    if (is.na(v)) return(default)
    out <- suppressWarnings(as.numeric(v))
    if (is.na(out)) stop("monitoring plan: item '", nm, "' is not numeric (", v, ").", call. = FALSE)
    out
  }
  th <- list(psi = c(get("psi_score_fixed_moderate", 0.10), get("psi_score_fixed_action", 0.25)),
             csi = c(get("csi_variable_fixed_moderate", 0.10), get("csi_variable_fixed_action", 0.25)),
             alpha = get("psi_adjusted_alpha", 0.05), min_events = get("min_events_per_period", 100))
  if (any(diff(th$psi) <= 0) || any(diff(th$csi) <= 0)) stop("monitoring plan: the action threshold must exceed the moderate one.", call. = FALSE)
  list(table = plan, thresholds = th)
}

#' CSI of one variable against its training bin shares, with both thresholds
#' @keywords internal
#' @noRd
.csi_row <- function(pt, cmp, alpha, thresholds = c(0.10, 0.25)) {
  n_new <- sum(cmp); n_tr <- sum(pt$count_train); k <- nrow(pt)
  csi_v <- if (n_new > 0L) {
    sm <- if (any(pt$count_train == 0L) || any(cmp == 0L)) 0.5 else 0
    pb <- (pt$count_train + sm) / (n_tr + sm * k); pc <- (cmp + sm) / (n_new + sm * k)
    sum((pb - pc) * log(pb / pc))
  } else NA_real_
  crit <- (1 / n_tr + 1 / max(1L, n_new)) * stats::qchisq(1 - alpha, df = max(1L, k - 1L))
  list(n = n_new, csi = csi_v,
       flag_fixed = if (is.na(csi_v)) NA_character_ else if (csi_v < thresholds[1]) "stable" else if (csi_v < thresholds[2]) "moderate" else "shift",
       critical = crit,
       flag_adjusted = if (is.na(csi_v)) NA_character_ else if (csi_v < crit) "stable" else "shift",
       points_shift = if (n_new > 0L) .points_shift(pt$count_train / n_tr, cmp / n_new, pt$points) else NA_real_)
}

#' @keywords internal
#' @noRd
.csi_dt <- function(period, variable, r) {
  data.table::data.table(period = period, variable = variable, n = r$n, csi = r$csi, flag_fixed = r$flag_fixed,
                         critical = r$critical, flag_adjusted = r$flag_adjusted, points_shift = r$points_shift)
}

#' Monitor the scorecard on new data
#'
#' Recomputes, per period of `date_col` (or for the whole data), the score
#' PSI against train with frozen bands, the CSI of every variable with
#' frozen bins plus the signed points shift and, when the target is
#' present, the performance by vintage (event rate, AUC/KS/Gini with CI).
#' Always reports both thresholds (fixed and n-adjusted). Schedules
#' nothing: the analyst calls it when needed.
#'
#' @param x An object from [scr_scorecard()].
#' @param newdata New table with the source columns.
#' @param date_col Period column. `NULL` treats `newdata` as a single period.
#'   Rows with a missing date form a period of their own (`NA`, last).
#' @param target Target column in `newdata`, for the performance by vintage.
#'   `NULL` skips it.
#' @param alpha Level of the adjusted threshold. `NULL` (default) takes it
#'   from the plan.
#' @param n_boot CI resamples per vintage. `NULL` uses the configuration.
#' @param plan The monitoring contract: `NULL` (default) uses the plan stored
#'   in the scorecard ([scr_monitoring_plan()]); otherwise an `item`/`value`
#'   table, or the path of a strategy workbook written by [scr_export()],
#'   whose `Monitoring_Plan` sheet is read. The fixed thresholds of the PSI
#'   and CSI flags, the alpha of the adjusted threshold and
#'   `min_events_per_period` come from it.
#'
#' @return An `scr_monitor` object with `psi` (score, per period), `csi`
#'   (per variable and period), `vintage` (or `NULL`; `status` says
#'   `"insufficient"` when a period has fewer events than the plan requires)
#'   and `plan` (the contract actually used).
#'
#' @family production
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   xgb_rounds = 60, n_boot = 20)
#' res <- scr_select(scr_demo, "default", config = cfg, drop = "id",
#'                   date_col = "ref_date")
#' sc <- scr_scorecard(res)
#' mo <- scr_monitor(sc, scr_demo, date_col = "ref_date", target = "default")
#' mo
#' mo$psi
#' head(mo$csi)
#' @export
scr_monitor <- function(x, newdata, date_col = NULL, target = NULL, alpha = NULL, n_boot = NULL, plan = NULL) {
  check_scorecard(x, "scr_monitor")
  n_boot <- n_boot %||% x$config$n_boot
  pl <- .read_plan(plan %||% x$monitoring_plan %||% scr_monitoring_plan(x))
  th <- pl$thresholds
  alpha <- alpha %||% th$alpha
  .scr_num1(alpha, "alpha", lower = 0, upper = 1, open_lower = TRUE)
  dt <- data.table::as.data.table(newdata)
  if (!is.null(date_col) && !date_col %in% names(dt)) stop("`date_col` does not exist in newdata.", call. = FALSE)
  if (!is.null(target) && !target %in% names(dt)) stop("`target` does not exist in newdata.", call. = FALSE)
  base <- .scr_preprocess(dt, x$ledger, x$features, x$config$special_values)
  w <- apply_woe(x$fit, base, x$features, "both")
  link <- .glm_link(x$coef, w, x$features)
  score <- x$alignment$a + x$alignment$b * link
  period <- if (is.null(date_col)) rep("all", nrow(dt)) else as.character(dt[[date_col]])
  # rows without a date form their own period (NA, listed last) instead of
  # silently leaving every table
  periods <- sort(unique(period), na.last = TRUE)
  tr_score <- x$samples$train$score
  # rows of every period, grouped once; the training bands, the score band
  # of every new row and the bin index of every variable are also computed
  # once, and each period only re-tabulates its rows
  rows <- split(seq_along(period), factor(match(period, periods), levels = seq_along(periods)))
  gb <- cut(tr_score, breaks = x$breaks, include.lowest = TRUE)
  lv <- levels(gb)
  nb <- tabulate(as.integer(gb), nbins = length(lv))
  band_new <- as.integer(cut(score, breaks = x$breaks, include.lowest = TRUE))

  psi <- data.table::rbindlist(lapply(seq_along(periods), function(k) {
    i <- rows[[k]]
    r <- .psi_counts(nb, tabulate(band_new[i], nbins = length(lv)), lv, alpha, th$psi)
    data.table::data.table(period = periods[k], n = length(i), mean_score = mean(score[i]), psi = r$psi,
                           flag_fixed = r$flag_fixed, critical = r$critical, flag_adjusted = r$flag_adjusted)
  }))
  pts_f <- lapply(stats::setNames(x$features, x$features), function(f) x$points[variable == f])
  # the base distribution is the training bin share stored in the points
  # table; the adjusted critical value uses the real training size
  bin_idx <- lapply(x$features, function(f) match(w[[paste0(f, "_bin")]], pts_f[[f]]$bin))
  names(bin_idx) <- x$features
  csi <- data.table::rbindlist(lapply(seq_along(periods), function(k) {
    i <- rows[[k]]
    data.table::rbindlist(lapply(x$features, function(f) {
      pt <- pts_f[[f]]
      .csi_dt(periods[k], f, .csi_row(pt, tabulate(bin_idx[[f]][i], nbins = nrow(pt)), alpha, th$csi))
    }))
  }))
  vintage <- NULL
  if (!is.null(target)) {
    y <- .target_as_int(dt[[target]], target, if (isTRUE(x$event$inverted)) 0L else NULL)$y
    hie <- identical(x$direction, "higher_is_riskier")
    vintage <- data.table::rbindlist(lapply(seq_along(periods), function(k) {
      i <- rows[[k]]
      m <- scr_metrics(score[i], y[i], higher_is_event = hie, ci = TRUE, n_boot = n_boot,
                       level = x$config$ci_level, seed = x$config$seed, nthread = x$config$nthread)
      data.table::data.table(period = periods[k], n = length(i), events = sum(y[i]), event_rate = mean(y[i]),
                             mean_score = mean(score[i]), auc = m$auc, auc_lo = m$auc_lo, auc_hi = m$auc_hi,
                             ks = m$ks, ks_lo = m$ks_lo, ks_hi = m$ks_hi, gini = m$gini,
                             status = if (sum(y[i]) < th$min_events) "insufficient" else "ok")
    }))
  }
  structure(list(psi = psi[], csi = csi[], vintage = vintage, plan = pl$table, thresholds = th,
                 periods = periods, target = x$target),
            class = c("scr_monitor", "list"))
}

#' Platform-stable short number (sprintf("%.2f", 0.005) differs between C libraries)
#' @keywords internal
#' @noRd
.g3 <- function(v) formatC(v, digits = 3, format = "g", flag = "#") |> sub(pattern = "\\.?0+$", replacement = "") |> sub(pattern = "\\.$", replacement = "")

#' @export
print.scr_monitor <- function(x, ...) {
  th <- x$thresholds
  cat(sprintf("<scr_monitor> target \"%s\" | %d period(s) | plan: PSI %s/%s, CSI %s/%s, alpha %s, min events %g\n",
              x$target, length(x$periods), .g3(th$psi[1]), .g3(th$psi[2]), .g3(th$csi[1]), .g3(th$csi[2]), .g3(th$alpha), th$min_events))
  cat(sprintf("  %-12s %8s %10s %8s %-9s %9s %-8s\n", "period", "n", "score", "PSI", "fixed", "critical", "adj."))
  p <- x$psi
  for (i in seq_len(nrow(p))) cat(sprintf("  %-12s %8s %10.1f %8.4f %-9s %9.4f %-8s\n", p$period[i], n_fmt(p$n[i]), p$mean_score[i],
                                          p$psi[i], p$flag_fixed[i], p$critical[i], p$flag_adjusted[i]))
  top <- x$csi[order(-abs(points_shift))][seq_len(min(5L, nrow(x$csi)))]
  if (nrow(top)) {
    cat("  largest points shifts (variable @ period):\n")
    for (i in seq_len(nrow(top))) cat(sprintf("    %-28s %-12s CSI %.4f  shift %+.2f pts\n", top$variable[i], top$period[i], top$csi[i], top$points_shift[i]))
  }
  if (!is.null(x$vintage)) {
    v <- x$vintage
    cat("  performance by vintage:\n")
    for (i in seq_len(nrow(v))) cat(sprintf("    %-12s n %-7s event %6.2f%%  AUC %.4f [%.4f, %.4f]  KS %.4f%s\n", v$period[i], n_fmt(v$n[i]),
                                            100 * v$event_rate[i], v$auc[i], v$auc_lo[i], v$auc_hi[i], v$ks[i],
                                            if (identical(v$status[i], "insufficient")) "  (insufficient events)" else ""))
  }
  invisible(x)
}
