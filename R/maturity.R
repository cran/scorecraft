# ============================================================================ #
# maturity.R - events over time by score band, with censoring
# ============================================================================ #
# The units are counted once per band and time (events and censorings); the
# Kaplan-Meier curve of every band is a cumulative product over its distinct
# times, and every horizon is read from the curve by a rolling join. The
# discrimination at each horizon comes from one count table per score cell.
# ============================================================================ #

#' Maturity of the event by score band
#'
#' Follows the units of each score band over time and estimates, at each
#' horizon, the cumulative share that has had the event, allowing for units
#' whose follow-up ends before the horizon (censoring). It shows how fast
#' each band matures, where the curve flattens (the performance window) and
#' how the discrimination of the score changes with the horizon.
#'
#' @section Data:
#'
#' One row per unit (a loan, a customer), with the score at the origin,
#' `time` and `event`:
#' \describe{
#'   \item{`time`}{A non-negative number of periods (days, months) from the
#'     origin to the event, or to the end of the follow-up.}
#'   \item{`event`}{1 when the event was observed at `time`, 0 when the
#'     follow-up ended there without it (censored).}
#' }
#' From a monthly panel, take per unit the months from its opening date to
#' its first month in default (`event = 1`), or to its last observed month
#' when it never defaulted (`event = 0`); the example does it for
#' [scr_demo_panel].
#'
#' Rows with a missing or infinite score, a missing time or event, or a zero
#' weight are left out (`n_dropped`).
#'
#' @section Incidence:
#'
#' Per band and for all units, the Kaplan-Meier estimate of the survival is
#' \eqn{S(t) = \prod_{t_j \le t} (1 - d_j / n_j)}, with \eqn{d_j} the events
#' at the distinct time \eqn{t_j} and \eqn{n_j} the units still followed
#' just before it (units censored at \eqn{t_j} are at risk at \eqn{t_j}).
#' `incidence` is \eqn{1 - S(h)} at the horizon \eqn{h}, and `se` its
#' standard error by Greenwood's formula,
#' \eqn{S(h) \sqrt{\sum_{t_j \le h} d_j / (n_j (n_j - d_j))}}. The interval
#' `lo`, `hi` is the complementary log-log interval of \eqn{S(h)},
#' \eqn{S^{\exp(\pm z \sigma)}} with
#' \eqn{\sigma = \sqrt{\sum d_j / (n_j (n_j - d_j))} / |\ln S|}, which stays
#' inside \[0, 1\]. It is undefined when no event has happened by the
#' horizon (`lo = 0`, `hi = NA`) and when every unit has had the event
#' (`NA`).
#'
#' Under weights, \eqn{d_j} and \eqn{n_j} are weighted sums and each term of
#' Greenwood's sum is computed on the Kish effective size of the risk set,
#' \eqn{d_j / ((n_j - d_j)\, n^{eff}_j)} with
#' \eqn{n^{eff}_j = n_j^2 / \sum w^2}; equal weights give the unweighted
#' result.
#'
#' Per band and horizon the table also counts `events` (events up to the
#' horizon), `censored` (units censored before it) and `at_risk` (the other
#' units: still followed at the horizon without the event), which add up to
#' `n`. `pct_of_final` is the incidence over the incidence at the largest
#' horizon: the share of the final events already seen.
#'
#' Past the last follow-up time of a band nothing more is observed: the
#' curve is carried flat, and a horizon beyond it repeats the last estimate
#' with `at_risk = 0`. Read such a row as "no information", not as a
#' plateau of the event rate.
#'
#' @section Discrimination:
#'
#' At each horizon, the units with a complete window are those with the
#' event by the horizon and those followed for at least the horizon without
#' it; units censored earlier are left out. `auc` and `gini` are those of
#' the score for "event by the horizon" among them, from the counts per
#' score value. The complete window ignores the censored units, so it is
#' unbiased only when censoring does not depend on the score.
#'
#' @inheritParams scr_bands
#' @param x A `data.frame` with one row per unit.
#' @param score,time,event Column names of the score, of the time to the
#'   event or to the end of the follow-up, and of the 0/1 event indicator.
#' @param horizons Times at which the incidence is read, in the unit of
#'   `time` (positive numbers).
#' @param n_bands Equal-share bands of the score (tie-safe, the event-richest
#'   first) when `cuts` is not given.
#' @param cuts Optional ascending cuts of the score, or an object from
#'   [scr_bands()] or [scr_tiers()], whose cuts, numbers and labels are then
#'   used, with its objective and direction.
#' @param level Confidence level of the intervals.
#' @param weight Optional column of non-negative case weights.
#'
#' @return An object of class `c("scr_maturity", "list")`:
#'   \describe{
#'     \item{`table`}{One row per band and horizon, event-richest band
#'       first, then all units (`band = NA`, `label = "all"`): `band`,
#'       `label`, `n`, `horizon`, `at_risk`, `events`, `censored`,
#'       `incidence`, `se`, `lo`, `hi` and `pct_of_final`.}
#'     \item{`discrimination`}{One row per horizon: `horizon`, `n` (units
#'       with a complete window), `events`, `auc` and `gini`.}
#'     \item{`cuts`, `codes`, `labels`}{The ascending cuts, and the band
#'       number and label of every interval in ascending score order.}
#'     \item{`horizons`, `level`, `objective`, `direction`, `score`, `time`,
#'       `event`, `n`, `n_rows`, `n_dropped`, `weighted`, `call`}{The
#'       settings: `n` is the volume used (the sum of the weights), `n_rows`
#'       the rows used and `n_dropped` those left out.}
#'   }
#'
#' @references
#' Greenwood, M. (1926). The natural duration of cancer. *Reports on Public
#' Health and Medical Subjects*, 33, 1-26. HMSO.
#'
#' Kalbfleisch, J. D. and Prentice, R. L. (2002). *The Statistical Analysis
#' of Failure Time Data*, 2nd edition. Wiley. \doi{10.1002/9781118032985}
#'
#' Kaplan, E. L. and Meier, P. (1958). Nonparametric estimation from
#' incomplete observations. *Journal of the American Statistical
#' Association*, 53(282), 457-481. \doi{10.1080/01621459.1958.10501452}
#'
#' @seealso [scr_bands()] and [scr_tiers()] for the bands, [scr_default()]
#'   and [scr_default_rate()] for the default flag and the cohort rates of a
#'   monthly panel.
#' @family score-studies
#' @examples
#' # time and event from a monthly panel: months from the opening date to
#' # the first default (90 days past due), or to the last observed month
#' p <- scr_demo_panel[order(scr_demo_panel$id, scr_demo_panel$ref_date), ]
#' month <- 12 * as.integer(format(p$ref_date, "%Y")) + as.integer(format(p$ref_date, "%m"))
#' first <- !duplicated(p$id)
#' bad <- which(p$dpd >= 90)
#' bad <- bad[!duplicated(p$id[bad])]
#' u <- data.frame(id = p$id[first], score = p$score[first], open = month[first])
#' u$event <- as.integer(u$id %in% p$id[bad])
#' u$end <- as.vector(tapply(month, p$id, max)[as.character(u$id)])
#' u$end[u$event == 1] <- month[bad][match(u$id[u$event == 1], p$id[bad])]
#' u$time <- u$end - u$open
#'
#' mt <- scr_maturity(u, horizons = c(6, 12, 24, 35), n_bands = 4)
#' mt
#' mt$discrimination
#' @export
scr_maturity <- function(x, ...) UseMethod("scr_maturity")

#' @rdname scr_maturity
#' @export
scr_maturity.data.frame <- function(x, score = "score", time = "time", event = "event", horizons,
                                    objective = "risk", direction = NULL, n_bands = 5L, cuts = NULL, level = 0.95,
                                    weight = NULL, max_cells = 1e5, ...) {
  fn <- "scr_maturity"
  .study_dots(list(...), fn)
  for (nm in c("score", "time", "event", "weight")) {
    v <- get(nm)
    if (!is.null(v)) .study_chr1(v, nm, fn)
  }
  if (missing(horizons) || !is.numeric(horizons) || !length(horizons) || any(!is.finite(horizons)) ||
      any(horizons <= 0)) {
    stop(fn, "(): `horizons` must be positive numbers, in the unit of `time`.", call. = FALSE)
  }
  hz <- sort(unique(as.double(horizons)))
  rd <- .study_one_direction(if (!missing(objective)) objective, direction, cuts, "risk", fn)
  miss <- setdiff(c(score, time, event, weight), names(x))
  if (length(miss)) stop(fn, "(): column(s) ", lst(miss), " not in `x`.", call. = FALSE)
  for (nm in c(score, time, weight)) {
    if (!is.numeric(x[[nm]])) stop(fn, "(): column '", nm, "' must be numeric.", call. = FALSE)
  }
  n_bands <- .study_whole(n_bands, "n_bands", fn, lower = 1)
  level <- .study_level(level, fn)
  sc <- as.double(x[[score]]); tm <- as.double(x[[time]])
  ev <- .scr_y01(x[[event]], fn)
  ok <- is.finite(sc) & !is.na(tm) & !is.na(ev)
  wtd <- !is.null(weight)
  w <- NULL
  if (wtd) {
    w <- as.double(x[[weight]])
    if (anyNA(w) || any(!is.finite(w)) || any(w < 0)) stop(fn, "(): the weights must be finite and non-negative.", call. = FALSE)
    # a row with zero weight does not belong to the population
    ok <- ok & w > 0
    w <- w[ok]
  }
  n_drop <- sum(!ok)
  sc <- sc[ok]; tm <- tm[ok]; ev <- ev[ok]
  if (!length(sc)) stop(fn, "(): no row has a score, a time and an event (and a positive weight).", call. = FALSE)
  if (any(!is.finite(tm)) || any(tm < 0)) stop(fn, "(): `time` must be finite and non-negative.", call. = FALSE)
  side <- .study_side(rd$direction)

  # one count table per score cell and horizon level: an event row is an
  # event from the first horizon at or after its time, a censored row has a
  # complete window up to the last horizon at or before its time
  lev <- ifelse(ev == 1L, findInterval(tm, hz, left.open = TRUE), findInterval(tm, hz))
  h <- .study_hist(sc, ev, w = w, by = list(sample = lev), max_cells = max_cells, fn = fn)
  cells <- .study_collapse(h)
  bd <- .cross_bands(cuts, cells, n_bands, side, "cuts", fn)
  B <- length(bd$cuts) + 1L

  # one grouped pass: events and censorings per band and time
  cols <- list(b = findInterval(sc, bd$cuts) + 1L, t = tm, d = if (wtd) w * ev else ev)
  agg <- if (wtd) { cols$w <- w; cols$q <- w * w; alist(d = sum(d), m = sum(w), q = sum(q)) } else alist(d = sum(d), m = .N)
  dt <- data.table::setDT(cols)
  j <- as.call(c(as.name("list"), agg))
  a <- dt[, eval(j), keyby = c("b", "t")]
  for (cn in setdiff(names(a), "b")) data.table::set(a, j = cn, value = as.double(a[[cn]]))
  if (!wtd) a[, q := m]
  # band 0 is every unit
  a <- data.table::rbindlist(list(a[, list(d = sum(d), m = sum(m), q = sum(q)), keyby = "t"][, b := 0L], a),
                             use.names = TRUE)
  data.table::setkeyv(a, c("b", "t"))
  km <- .maturity_km(a, c(0L, seq_len(B)), hz, level)

  # rows in table order: the event-richest band first, then every unit
  o <- if (side == "high") rev(seq_len(B)) else seq_len(B)
  ord <- unlist(lapply(c(o, 0L), function(b) which(km$b == b)))
  km <- km[ord]
  tab <- data.table::data.table(
    band = ifelse(km$b == 0L, NA_integer_, bd$codes[pmax(km$b, 1L)]),
    label = ifelse(km$b == 0L, "all", bd$labels[pmax(km$b, 1L)]), n = km$n, horizon = km$t, at_risk = km$at_risk,
    events = km$events, censored = km$censored, incidence = km$incidence, se = km$se, lo = km$lo, hi = km$hi,
    pct_of_final = km$pct_of_final)

  disc <- .maturity_discrimination(h, cells, hz, rd$direction)
  structure(list(table = tab, discrimination = disc, cuts = bd$cuts, codes = bd$codes, labels = bd$labels,
                 bands = bd$source, horizons = hz, level = level, objective = rd$objective,
                 direction = rd$direction, score = score, time = time, event = event,
                 n = if (wtd) sum(w) else as.double(length(sc)), n_rows = length(sc), n_dropped = n_drop,
                 weighted = wtd, quantized = !is.null(attr(h, "edges")), call = sys.call()),
            class = c("scr_maturity", "list"))
}

#' Kaplan-Meier incidence of every band at every horizon
#'
#' `a` is keyed by band `b` and time `t`, with the events `d`, the volume
#' `m` (events and censorings) and the squared weights `q` of each time. The
#' curve of a band is a cumulative product over its times; a horizon is read
#' at the last time at or before it. Greenwood's sum is taken on the Kish
#' effective size of the risk set, and the interval is the complementary
#' log-log one.
#' @keywords internal
#' @noRd
.maturity_km <- function(a, bands, hz, level) {
  a[, cens := m - d]
  # the risk set of a time: every unit with that time or a later one
  a[, `:=`(risk = rev(cumsum(rev(m))), q_risk = rev(cumsum(rev(q)))), by = "b"]
  a[, rem := pmax(risk - d, 0)]
  a[, `:=`(surv = cumprod(rem / risk),
           gw = cumsum(ifelse(d > 0, d / (rem * risk * risk / q_risk), 0)),
           cum_d = cumsum(d), cum_c = cumsum(cens), cum_c_before = cumsum(cens) - cens, t0 = t), by = "b"]
  tot <- a[, list(n = sum(m)), keyby = "b"]
  qry <- data.table::CJ(b = bands, t = hz)
  # the last time at or before each horizon; none when the horizon comes first
  r <- a[qry, on = c("b", "t"), roll = TRUE]
  n <- tot$n[match(r$b, tot$b)]
  n[is.na(n)] <- 0
  none <- is.na(r$t0)
  surv <- ifelse(none, 1, r$surv)
  gw <- ifelse(none, 0, r$gw)
  events <- ifelse(none, 0, r$cum_d)
  # censored strictly before the horizon: a unit censored at the horizon was followed that long
  cens <- ifelse(none, 0, ifelse(r$t0 == r$t, r$cum_c_before, r$cum_c))
  empty <- !(n > 0)
  surv[empty] <- NA_real_
  z <- stats::qnorm(1 - (1 - level) / 2)
  inc <- 1 - surv
  se <- ifelse(surv > 0, surv * sqrt(gw), NA_real_)
  # complementary log-log interval of the survival, turned into the incidence
  mid <- !is.na(surv) & surv > 0 & surv < 1 & is.finite(gw)
  sg <- rep(NA_real_, length(surv))
  sg[mid] <- sqrt(gw[mid]) / abs(log(surv[mid]))
  lo <- hi <- rep(NA_real_, length(surv))
  lo[mid] <- 1 - surv[mid]^exp(-z * sg[mid])
  hi[mid] <- 1 - surv[mid]^exp(z * sg[mid])
  # no event yet: the incidence is 0 and the interval has no upper end
  lo[!is.na(surv) & surv >= 1] <- 0
  out <- data.table::data.table(b = r$b, t = r$t, n = n, at_risk = pmax(n - events - cens, 0), events = events,
                                censored = cens, incidence = inc, se = se, lo = lo, hi = hi)
  out[empty, `:=`(at_risk = 0, events = 0, censored = 0)]
  # share of the incidence at the largest horizon
  out[, pct_of_final := if (isTRUE(incidence[.N] > 0)) incidence / incidence[.N] else NA_real_, by = "b"]
  out[]
}

#' AUC of the score for "event by the horizon" among the complete windows
#'
#' `h` is keyed by the horizon level and the score cell: an event row of
#' level `l` is an event at the horizons after `l` and a non-event with a
#' complete window at the first `l`; a censored row of level `l` has a
#' complete window at the first `l` horizons only.
#' @keywords internal
#' @noRd
.maturity_discrimination <- function(h, cells, hz, direction) {
  H <- length(hz); K <- nrow(cells)
  key <- if ("g" %in% names(cells)) "g" else "s"
  i <- match(h[[key]], cells[[key]])
  E <- N <- matrix(0, K, H + 1L)
  at <- cbind(i, as.integer(h[["sample"]]) + 1L)
  E[at] <- h$e; N[at] <- h$n_y
  flip <- !identical(direction, "higher_is_riskier")
  out <- lapply(seq_len(H), function(j) {
    # events by horizon j: event rows of level below j; complete non-events: every row of level j or more
    c1 <- rowSums(E[, seq_len(j), drop = FALSE])
    c0 <- rowSums(N[, (j + 1L):(H + 1L), drop = FALSE])
    if (flip) { c1 <- rev(c1); c0 <- rev(c0) }
    ok <- sum(c1) > 0 && sum(c0) > 0
    auc <- if (ok) .auc_ks_counts(c1, c0)$auc else NA_real_
    data.table::data.table(horizon = hz[j], n = sum(c1) + sum(c0), events = sum(c1), auc = auc, gini = 2 * auc - 1)
  })
  data.table::rbindlist(out)
}

#' A horizon as short text, each on its own
#' @keywords internal
#' @noRd
.maturity_h <- function(h) vapply(h, format, character(1), digits = 7, trim = TRUE)

#' @export
print.scr_maturity <- function(x, ...) {
  cat(sprintf("<scr_maturity> event \"%s\" over \"%s\" | objective %s | %s\n", x$event, x$time, x$objective,
              x$direction))
  cat(sprintf("  %s units%s | %d bands (%s) | horizons %s | %s%% complementary log-log intervals\n",
              .study_n(x$n), if (x$n_dropped > 0) sprintf(" (%s rows left out)", .study_n(x$n_dropped)) else "",
              length(x$cuts) + 1L, x$bands, paste(.maturity_h(x$horizons), collapse = ", "), .g3(100 * x$level)))
  t <- x$table
  hz <- x$horizons
  first <- t[t$horizon == hz[1]]
  cat("\nCumulative incidence (Kaplan-Meier), event-richest band first\n")
  cat(sprintf("  %4s %-24s %9s%s\n", "band", "score", "n", paste(sprintf(" %9s", paste0("t=", .maturity_h(hz))), collapse = "")))
  for (i in seq_len(nrow(first))) {
    r <- t[if (is.na(first$band[i])) is.na(t$band) else t$band %in% first$band[i]]
    cat(sprintf("  %4s %-24s %9s%s\n", if (is.na(first$band[i])) "" else as.character(first$band[i]),
                substr(first$label[i], 1, 24), .study_n(first$n[i]),
                paste(sprintf(" %9s", .claims_pct(r$incidence[match(hz, r$horizon)])), collapse = "")))
  }
  al <- t[is.na(t$band)]
  cat(sprintf("  %4s %-24s %9s%s\n", "", "share of the final (all)", "",
              paste(sprintf(" %9s", .study_f(100 * al$pct_of_final[match(hz, al$horizon)], "%.0f%%")), collapse = "")))
  d <- x$discrimination
  cat("\nDiscrimination at each horizon (units with a complete window)\n")
  cat(sprintf("  %9s %9s %9s %8s %8s\n", "horizon", "n", "events", "AUC", "Gini"))
  for (i in seq_len(nrow(d))) {
    cat(sprintf("  %9s %9s %9s %8s %8s\n", .maturity_h(d$horizon[i]), .study_n(d$n[i]), .study_n(d$events[i]),
                .study_f(d$auc[i], "%.4f"), .study_f(d$gini[i], "%.4f")))
  }
  invisible(x)
}

#' @rdname scr_export
#' @export
scr_export.scr_maturity <- function(x, dir, stamp = TRUE, ...) {
  .study_dots(list(...), "scr_export")
  .need_openxlsx()
  out_dir <- .export_dir(dir, stamp)
  tag <- .file_tag(x$event)
  settings <- .kv_table(list(score = x$score, time = x$time, event = x$event, objective = x$objective,
                             direction = x$direction, horizons = x$horizons, level = x$level, n = x$n,
                             n_rows = x$n_rows, n_dropped = x$n_dropped, weighted = x$weighted))
  cuts <- data.frame(cut = seq_along(x$cuts), score = x$cuts)
  sheets <- lapply(list(Incidence = x$table, Discrimination = x$discrimination, Cuts = cuts, Settings = settings),
                   .study_sheet)
  files <- list(xlsx = .scr_write_xlsx(sheets, file.path(out_dir, sprintf("maturity_%s.xlsx", tag))))
  for (f in files) msg("  %s", f)
  x$files <- files
  invisible(x)
}

# data.table column names used without quotes in this file
utils::globalVariables(c("b", "d", "m", "q", "cens", "risk", "q_risk", "rem", "surv", "gw", "cum_d", "cum_c",
                         "cum_c_before", "t0", "incidence", "pct_of_final", "at_risk", "censored"))
