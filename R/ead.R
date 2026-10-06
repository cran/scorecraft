# ============================================================================ #
# ead.R - exposure at default: realised CCF data set, pools, downturn,
#         application, production SQL, validation and export
# ============================================================================ #
# Notation: R the reference date, D the default
# date, E and L the drawn amount and the limit. The realised measures are
#   ULF / CCF = (E_D - E_R) / (L_R - E_R)      undrawn part material at R,
#   LF        = E_D / L_R                       region of instability, over
#                                               limit or nothing undrawn,
#   EADF      = E_D / E_R                       no usable limit,
# with the identity LF = u + CCF (1 - u), u = E_R / L_R. The realised EAD
# is never capped at the limit. Pools are cells of driver bins fitted with
# the continuous binner of R/irb-binning.R, so that obwoe_apply() and
# obwoe_sql() reproduce the bin index in R and in every SQL dialect. The
# applied CCF is the long-run average plus the estimation-error margin,
# raised by the downturn value when one is supplied, and floored at
# `ccf_floor_fraction * ccf_sa_ccf`; the predicted EAD is never below the
# drawn amount.
# ============================================================================ #

# -- 1. Reference data set -------------------------------------------------- #

#' Build the realised-CCF reference data set from facility snapshots
#'
#' One row per default event (or per event and reference date under the
#' variable-horizon comparison), with the facility as it stood at the
#' reference date and the realised exposure at default (EAD) at the default
#' date, from which the realised credit conversion factor (CCF) follows. The
#' reference date follows `config$ccf_horizon`: `"fixed"` takes the
#' snapshot `ccf_horizon_months` before the default month (the nearest
#' earlier snapshot when that month is missing; the first snapshot for a
#' facility younger than the horizon, flagged `FAST_DEFAULT`); `"cohort"`
#' takes the start of the calendar cohort window in which the default
#' falls; `"variable"` takes every snapshot in the horizon before the
#' default, for comparison only.
#'
#' The realised measure per row follows `config$ccf_measure`: under
#' `"auto"` the undrawn-limit factor (CCF) when the utilisation at the
#' reference date is below `ccf_u_star` and the limit factor (LF) at or
#' above it; rows with nothing undrawn or over the limit at the reference
#' date are always routed to the limit factor (`ZERO_UNDRAWN`,
#' `OVER_LIMIT_AT_REF`), never dropped. The raw realised value is kept in
#' `ccf_raw`; `ccf` carries the value after the optional floor
#' (`ccf_floor_realised`) and cap (`ccf_cap_realised`), both logged in the
#' funnel (`NEGATIVE_CCF_FLOORED`, `CCF_ABOVE_ONE`). The realised EAD is
#' the drawn amount at the default date, uncapped; with
#' `post_default_drawings_in = "ccf"` it is the maximum drawn amount over
#' the default event when `defaulted` is given.
#'
#' @param snapshots A `data.frame` or `data.table` with one row per
#'   facility and month.
#' @param facility_id,date_col,limit,drawn Column names of the facility
#'   identifier, the snapshot month, the limit and the drawn amount.
#' @param obligor_id Column name of the obligor, optional. With
#'   `config$default_level = "obligor"` a default of any facility of the
#'   obligor is a default of all its facilities observed at that date.
#' @param default_date Either the name of a column of `snapshots` holding
#'   the default date of the facility (`NA` when it never defaults) or a
#'   `data.frame` with the facility identifier column (same name as
#'   `facility_id`) and a `default_date` column, one row per event.
#' @param defaulted Column name of a 0/1 default flag per snapshot,
#'   alternative to `default_date`: every run of ones opens an event at
#'   its first month.
#' @param drivers Column names measured at the reference date and carried
#'   into the data set (candidate drivers of the pools). `utilisation_ref`,
#'   `limit_ref`, `drawn_ref` and `horizon_months` are always available.
#' @param config A [scr_config()]; keys `ccf_*`,
#'   `post_default_drawings_in`, `default_level`.
#' @param keep_rows If `TRUE`, keeps every candidate event with its
#'   exclusion rule in `rows`.
#'
#' @return An object of class `scr_ead_data`: `rds` (one row per event:
#'   `event_id`, `facility_id`, `obligor_id`, `ref_date`, `default_date`,
#'   `cohort`, `horizon_months`, `fast_default`, `limit_ref`, `drawn_ref`,
#'   `undrawn_ref`, `utilisation_ref`, `limit_default`, `limit_change`,
#'   `ead_realised`, `measure`, `ccf_raw`, `ccf`, `rule`, drivers),
#'   `funnel` (`rule`, `action`, `n`, `share`), `summary` (simple and
#'   exposure-weighted averages by cohort and measure, with a total row),
#'   `lra` (long-run averages and shares), `meta`, `ledger`, `config` and
#'   `rows` (with `keep_rows = TRUE`).
#'
#' @references
#' Basel Committee on Banking Supervision (2023). *The Basel Framework*,
#' CRE32 and CRE36. Moral, G. (2006). EAD estimates for facilities with
#' explicit limits. In Engelmann, B. and Rauhmeier, R. (eds), *The Basel
#' II Risk Parameters*. Springer.
#'
#' @family irb-ead
#' @examples
#' cfg <- scr_config(verbose = FALSE)
#' ed <- scr_ead_data(scr_demo_ead, facility_id = "facility_id", obligor_id = "obligor_id",
#'                    date_col = "ref_date", limit = "limit", drawn = "drawn",
#'                    defaulted = "defaulted",
#'                    drivers = c("product", "months_on_book", "dpd"), config = cfg)
#' ed
#' ed$funnel
#' head(ed$rds[, c("facility_id", "ref_date", "default_date", "utilisation_ref", "measure", "ccf")])
#' @export
scr_ead_data <- function(snapshots, facility_id, obligor_id = NULL, date_col, limit, drawn,
                         default_date = NULL, defaulted = NULL, drivers = NULL,
                         config = scr_config(), keep_rows = FALSE) {
  check_config(config, "scr_ead_data")
  cfg <- config
  old <- scr_verbose(isTRUE(cfg$verbose)); on.exit(scr_verbose(old), add = TRUE)
  t0 <- Sys.time()
  dt <- data.table::as.data.table(snapshots)
  need <- c(facility_id, obligor_id, date_col, limit, drawn, defaulted, drivers,
            if (is.character(default_date)) default_date)
  miss <- setdiff(need, names(dt))
  if (length(miss)) stop("scr_ead_data(): column(s) not found: ", lst(miss), call. = FALSE)
  if (is.null(default_date) && is.null(defaulted)) stop("scr_ead_data(): give `default_date` or `defaulted`.", call. = FALSE)
  reserved <- c("facility_id", "obligor_id", "ref_date", "default_date", "limit_ref", "drawn_ref", "undrawn_ref",
                "utilisation_ref", "ead_realised", "ccf", "ccf_raw", "measure", "rule", "event_id", "cohort",
                "horizon_months", "fast_default", "limit_default", "limit_change")
  clash <- intersect(drivers, reserved)
  if (length(clash)) stop("scr_ead_data(): driver name(s) clash with reserved columns: ", lst(clash), call. = FALSE)
  # names of the internal panel columns: a driver may only reuse `limit` or
  # `drawn` when it is that very column
  clash <- c(intersect(drivers, c("fid", "oid", "date", "def", ".ddate")),
             if ("limit" %in% drivers && !identical(limit, "limit")) "limit",
             if ("drawn" %in% drivers && !identical(drawn, "drawn")) "drawn")
  if (length(clash)) stop("scr_ead_data(): driver name(s) clash with internal columns: ", lst(clash), call. = FALSE)
  H <- as.integer(cfg$ccf_horizon_months)
  if (is.na(H) || H < 1L) stop("scr_ead_data(): `ccf_horizon_months` must be a positive integer.", call. = FALSE)

  msg_stage(11, sprintf("EAD reference data set (%s horizon, %d months, measure %s)", cfg$ccf_horizon, H, cfg$ccf_measure))
  p <- data.table::data.table(fid = as.character(dt[[facility_id]]),
                              oid = if (is.null(obligor_id)) as.character(dt[[facility_id]]) else as.character(dt[[obligor_id]]),
                              date = as.Date(dt[[date_col]]),
                              limit = as.double(dt[[limit]]), drawn = as.double(dt[[drawn]]))
  if (anyNA(p$fid) || anyNA(p$date)) stop("scr_ead_data(): the facility identifier and the date cannot be missing.", call. = FALSE)
  p[, date := .month_floor(date)]
  # checked after the month floor: two snapshots of one facility in the same
  # month would otherwise both survive and match() would pick one silently
  if (anyDuplicated(p, by = c("fid", "date"))) stop("scr_ead_data(): duplicated (facility, month) rows.", call. = FALSE)
  for (d in drivers) data.table::set(p, j = d, value = dt[[d]])
  # the default-date column travels with its row through the sort below
  if (is.character(default_date)) data.table::set(p, j = ".ddate", value = .month_floor(as.Date(dt[[default_date]])))
  p[, def := if (is.null(defaulted)) NA_integer_ else as.integer(isTRUE_vec(dt[[defaulted]]))]
  data.table::setorder(p, fid, date)

  # -- default events ----------------------------------------------------- #
  events <- .ead_events(p, dt, facility_id, default_date, defaulted, obligor_level = !is.null(obligor_id) && identical(cfg$default_level, "obligor"))
  if (!nrow(events)) stop("scr_ead_data(): no default event found.", call. = FALSE)
  msg("  %s facilities x %s months | %s default events on %s facilities", n_fmt(data.table::uniqueN(p$fid)),
      n_fmt(data.table::uniqueN(p$date)), n_fmt(nrow(events)), n_fmt(data.table::uniqueN(events$fid)))

  # -- reference rows ----------------------------------------------------- #
  all_dates <- sort(unique(p$date))
  cohort_starts <- all_dates[seq(1L, length(all_dates), by = H)]
  # rolling joins over all events at once (no per-event closure or subset)
  rows <- .ead_reference_table(events, p, cfg$ccf_horizon, H, cohort_starts, cfg$post_default_drawings_in, drivers)
  if (is.null(rows) || !nrow(rows)) stop("scr_ead_data(): no reference row could be built.", call. = FALSE)

  # -- admission rules and measures --------------------------------------- #
  rows[, rule := "OK"]
  rows[is.na(limit_ref) | is.na(drawn_ref) | is.na(ead_realised) | limit_ref < 0 | drawn_ref < 0 | ead_realised < 0, rule := "DATA_ERROR"]
  rows[rule == "OK" & horizon_months < 1L, rule := "FAST_DEFAULT_EXCLUDED"]
  rows[rule == "OK" & limit_ref <= 0, rule := "NOT_IN_SCOPE"]
  rows[, undrawn_ref := limit_ref - drawn_ref]
  rows[, utilisation_ref := data.table::fifelse(limit_ref > 0, drawn_ref / limit_ref, NA_real_)]
  rows[, limit_change := data.table::fifelse(limit_ref > 0, limit_default / limit_ref - 1, NA_real_)]
  rows[, fast_default := horizon_months < H]
  keep <- rows$rule == "OK"
  meas <- .ead_measure(rows$utilisation_ref, rows$undrawn_ref, rows$drawn_ref, cfg$ccf_measure, cfg$ccf_u_star)
  rows[, measure := meas]
  rows[keep & drawn_ref > limit_ref, rule := "OVER_LIMIT_AT_REF"]
  rows[keep & drawn_ref == limit_ref, rule := "ZERO_UNDRAWN"]
  rows[, ccf_raw := data.table::fcase(
    measure == "ulf", (ead_realised - drawn_ref) / undrawn_ref,
    measure == "lf", ead_realised / limit_ref,
    measure == "eadf", ead_realised / drawn_ref,
    default = NA_real_)]
  rows[, ccf := ccf_raw]
  fl <- cfg$ccf_floor_realised; cp <- cfg$ccf_cap_realised
  neg <- keep & rows$measure == "ulf" & !is.na(fl) & rows$ccf_raw < fl
  if (any(neg)) rows[neg, `:=`(ccf = fl, rule = data.table::fifelse(rule == "OK", "NEGATIVE_CCF_FLOORED", rule))]
  above <- keep & rows$measure == "ulf" & rows$ccf_raw > 1
  if (any(above)) rows[above & rule == "OK", rule := "CCF_ABOVE_ONE"]
  capped <- keep & rows$measure == "ulf" & !is.na(cp) & rows$ccf_raw > cp
  if (any(capped)) rows[capped, ccf := cp]
  rows[rule == "OK" & fast_default, rule := "FAST_DEFAULT"]
  rows[, cohort := format(ref_date, "%Y")]
  excluded <- c("DATA_ERROR", "FAST_DEFAULT_EXCLUDED", "NOT_IN_SCOPE")
  rows[, action := data.table::fcase(
    rule %in% excluded, "excluded",
    rule %in% c("OVER_LIMIT_AT_REF", "ZERO_UNDRAWN"), "routed_to_lf",
    rule == "NEGATIVE_CCF_FLOORED", "floored",
    rule == "CCF_ABOVE_ONE", if (!is.na(cp)) "capped" else "kept",
    rule == "FAST_DEFAULT", "kept",
    default = "kept")]
  rds <- data.table::copy(rows[!(rule %in% excluded)])   # explicit copy: rows is kept for the funnel
  rds[, action := NULL]
  data.table::setcolorder(rds, c("event_id", "facility_id", "obligor_id", "ref_date", "default_date", "cohort", "horizon_months",
                                 "fast_default", "limit_ref", "drawn_ref", "undrawn_ref", "utilisation_ref", "limit_default",
                                 "limit_change", "ead_realised", "measure", "ccf_raw", "ccf", "rule"))
  data.table::setorder(rds, ref_date, facility_id)

  # -- funnel ------------------------------------------------------------- #
  rule_order <- c("NOT_IN_SCOPE", "DATA_ERROR", "FAST_DEFAULT_EXCLUDED", "FAST_DEFAULT", "OVER_LIMIT_AT_REF", "ZERO_UNDRAWN",
                  "NEGATIVE_CCF_FLOORED", "CCF_ABOVE_ONE", "OK")
  fun <- rows[, list(n = .N), by = c("rule", "action")]
  # shares of the rows built (before exclusion); CCF_ABOVE_ONE counted on every row above one, floored counted where floored
  fun[, share := n / nrow(rows)]
  fun[, ord := match(rule, rule_order)]
  data.table::setorder(fun, ord); fun[, ord := NULL]
  n_above <- sum(above); n_neg <- sum(neg)
  for (i in seq_len(nrow(fun))) msg("  %-24s %-13s %6d (%s)", fun$rule[i], fun$action[i], fun$n[i], fmt_pct(fun$share[i]))
  msg("  reference data set: %s rows | %s above one | %s floored | %s in the LF measure", n_fmt(nrow(rds)), n_fmt(n_above),
      n_fmt(n_neg), n_fmt(sum(rds$measure == "lf")))

  # -- summary by cohort and measure -------------------------------------- #
  summary <- .ead_summary(rds)
  main <- if (identical(cfg$ccf_measure, "auto")) "ulf" else cfg$ccf_measure
  mr <- rds[measure == main]
  yrs <- if (nrow(rds)) as.numeric(difftime(max(rds$ref_date), min(rds$ref_date), units = "days")) / 365.25 else NA_real_
  lra <- list(measure = main, n = nrow(mr), simple = if (nrow(mr)) mean(mr$ccf) else NA_real_,
              exposure_weighted = .ead_ew(mr), n_total = nrow(rds),
              share_lf = mean(rds$measure == "lf"), share_above_one = if (nrow(mr)) mean(mr$ccf_raw > 1) else NA_real_,
              share_negative = if (nrow(mr)) mean(mr$ccf_raw < 0) else NA_real_, share_fast = mean(rds$fast_default),
              n_cohorts = data.table::uniqueN(rds$cohort), years = yrs)
  ledger <- data.table::data.table(
    step = "reference_data", action = "build",
    detail = sprintf("horizon %s (%d months); measure %s with u* = %s; floor %s; cap %s; post-default drawings in %s; default level %s",
                     cfg$ccf_horizon, H, cfg$ccf_measure, format(cfg$ccf_u_star), format(fl), format(cp),
                     cfg$post_default_drawings_in, if (!is.null(obligor_id)) cfg$default_level else "facility"),
    reason = "configuration", date = format(Sys.Date()))
  meta <- list(horizon = cfg$ccf_horizon, horizon_months = H, measure_rule = cfg$ccf_measure, main_measure = main,
               u_star = cfg$ccf_u_star, floor_realised = fl, cap_realised = cp,
               post_default_drawings_in = cfg$post_default_drawings_in,
               n_snapshots = nrow(p), n_facilities = data.table::uniqueN(p$fid), n_events = nrow(events), n_rows = nrow(rds),
               date_min = min(p$date), date_max = max(p$date), drivers = drivers %||% character(),
               cols = list(facility_id = facility_id, obligor_id = obligor_id, date_col = date_col, limit = limit, drawn = drawn),
               seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  msg("  done (%.2fs)", meta$seconds)
  structure(list(rds = rds, funnel = fun[, list(rule, action, n, share)], summary = summary, lra = lra, meta = meta,
                 ledger = ledger, config = cfg, rows = if (isTRUE(keep_rows)) rows else NULL),
            class = c("scr_ead_data", "list"))
}

#' First day of the month of a Date
#' @keywords internal
#' @noRd
.month_floor <- function(d) {
  # format/parse only the distinct dates: a panel has few months and many rows
  d <- as.Date(d)
  u <- unique(d)
  as.Date(format(u, "%Y-%m-01"))[match(as.double(d), as.double(u))]
}

# .months_between() is defined once, in R/lgd.R (day-aware); on the
# month-floored dates of this file it is the plain difference in months.
# A second definition here was silently masked by the one of R/lgd.R
# (collated later).

#' Default events: one row per (facility, default date)
#'
#' With a 0/1 flag, every run of ones opens an event at its first month; at
#' obligor level the flag is first raised to every facility of the obligor
#' observed in a month where any of its facilities is flagged, so that a
#' facility already in default is not given a second event. With a default
#' date, an obligor-level default is propagated to the facilities of the
#' obligor observed at that month.
#' @keywords internal
#' @noRd
.ead_events <- function(p, dt, facility_id, default_date, defaulted, obligor_level) {
  if (!is.null(defaulted)) {
    q <- p[, list(fid, oid, date, def)]
    if (isTRUE(obligor_level)) q[, def := max(def), by = c("oid", "date")]
    ev <- q[, {
      d <- def; s <- which(d == 1L & c(0L, d[-length(d)]) != 1L)
      list(default_date = date[s])
    }, by = "fid"]
  } else if (is.character(default_date)) {
    # p is sorted by (fid, date): read the default date from its own column,
    # never from `dt`, whose rows are in the caller's order
    ev <- unique(data.table::data.table(fid = p$fid, default_date = p$.ddate))
    ev <- ev[!is.na(default_date)]
  } else if (is.data.frame(default_date)) {
    dd <- data.table::as.data.table(default_date)
    if (!all(c(facility_id, "default_date") %in% names(dd))) {
      stop("scr_ead_data(): `default_date` table needs columns '", facility_id, "' and 'default_date'.", call. = FALSE)
    }
    ev <- unique(data.table::data.table(fid = as.character(dd[[facility_id]]), default_date = .month_floor(as.Date(dd$default_date))))
    ev <- ev[!is.na(default_date)]
  } else stop("scr_ead_data(): `default_date` must be a column name or a data.frame.", call. = FALSE)
  if (is.null(defaulted) && isTRUE(obligor_level) && nrow(ev)) {
    ob <- unique(p[, list(fid, oid)])
    ev <- merge(ev, ob, by = "fid")
    ev_o <- unique(ev[, list(oid, default_date)])
    ev <- merge(ev_o, p[, list(fid, oid, date)], by.x = c("oid", "default_date"), by.y = c("oid", "date"))[, list(fid, default_date)]
    ev <- unique(ev)
  }
  data.table::setorder(ev, fid, default_date)
  ev[, event_id := paste0(fid, "#", seq_len(.N)), by = "fid"]
  ev[]
}

#' Reference rows of every default event under the configured horizon
#'
#' Vectorised over the events with rolling joins on the (facility, month)
#' panel `p`, sorted by `fid` and `date`:
#' * the default row `i_d` is the first snapshot at or after the default
#'   month (next observation carried backward); an event without one is
#'   dropped;
#' * `"fixed"`: the last snapshot at or before `D - H` months (last
#'   observation carried forward), else the first snapshot of the facility;
#' * `"cohort"`: the last snapshot at or before the latest cohort start
#'   strictly before `D`, else the first snapshot;
#' * `"variable"`: every snapshot before `D` at most `H` months before it,
#'   else the last snapshot before `D`;
#' * a facility with no snapshot before `D` takes the default row itself
#'   (`horizon_months = 0`, excluded later as `FAST_DEFAULT_EXCLUDED`).
#' With `post_default = "ccf"` and a default flag, the realised EAD is the
#' maximum drawn amount from the default row to the end of the run of
#' flagged months that follows it.
#' @keywords internal
#' @noRd
.ead_reference_table <- function(events, p, horizon, H, cohort_starts, post_default, drivers) {
  q <- p[fid %in% events$fid]
  q[, .r := seq_len(.N)]
  q[, `:=`(.first = .r[1L], .last = .r[.N]), by = "fid"]
  ev <- data.table::data.table(event_id = events$event_id, fid = events$fid, D = events$default_date)
  ev[, i_d := q[ev, .r, on = c(fid = "fid", date = "D"), roll = -Inf]]
  ev <- ev[!is.na(i_d)]
  if (!nrow(ev)) return(NULL)
  first <- q$.first[ev$i_d]
  has_before <- ev$i_d > first
  # realised EAD: the drawn amount at the default row, or the running maximum
  # over the flagged months that follow it
  ead <- q$drawn[ev$i_d]
  if (identical(post_default, "ccf") && !all(is.na(q$def))) {
    q[, .run := data.table::rleid(fid, def)]
    q[, .smax := rev(cummax(rev(data.table::fifelse(is.na(drawn), -Inf, drawn)))), by = ".run"]
    nxt <- ev$i_d + 1L
    ok <- nxt <= q$.last[ev$i_d] & q$def[pmin(nxt, nrow(q))] %in% 1L
    m <- pmax(data.table::fifelse(is.na(ead), -Inf, ead), data.table::fifelse(ok, q$.smax[pmin(nxt, nrow(q))], -Inf))
    ead <- data.table::fifelse(m == -Inf, NA_real_, m)
  }
  locf <- function(target) q[data.table::data.table(fid = ev$fid, t = target), .r, on = c(fid = "fid", date = "t"), roll = Inf]
  if (identical(horizon, "variable")) {
    nb <- data.table::fifelse(has_before, ev$i_d - first, 0L)
    e <- rep(seq_len(nrow(ev)), nb)
    r <- rep(first, nb) + sequence(nb) - 1L
    ok <- .months_between(q$date[r], ev$D[e]) <= H
    pairs <- data.table::data.table(e = e[ok], r = r[ok])
    none <- setdiff(seq_len(nrow(ev)), pairs$e)
    fb <- data.table::data.table(e = none, r = data.table::fifelse(has_before[none], ev$i_d[none] - 1L, ev$i_d[none]))
    pairs <- data.table::rbindlist(list(pairs, fb))
    data.table::setorder(pairs, e, r)
    e <- pairs$e; refs <- pairs$r
  } else {
    target <- if (identical(horizon, "fixed")) .add_months(ev$D, -H) else {
      k <- findInterval(as.double(ev$D) - 0.5, as.double(cohort_starts))
      cohort_starts[replace(k, k == 0L, NA_integer_)]
    }
    refs <- locf(target)
    refs <- data.table::fifelse(is.na(refs), first, refs)
    refs <- data.table::fifelse(has_before, refs, ev$i_d)
    e <- seq_len(nrow(ev))
  }
  out <- data.table::data.table(event_id = ev$event_id[e], facility_id = ev$fid[e], obligor_id = q$oid[first[e]],
                                ref_date = q$date[refs], default_date = ev$D[e],
                                horizon_months = data.table::fifelse(has_before[e], as.integer(.months_between(q$date[refs], ev$D[e])), 0L),
                                limit_ref = q$limit[refs], drawn_ref = q$drawn[refs], limit_default = q$limit[ev$i_d[e]],
                                ead_realised = ead[e])
  for (d in drivers) data.table::set(out, j = d, value = q[[d]][refs])
  out
}

#' Measure of every row under the measure rule
#' @keywords internal
#' @noRd
.ead_measure <- function(u, undrawn, drawn, rule, u_star) {
  n <- length(u)
  base <- switch(rule, auto = "ulf", ulf = "ulf", lf = "lf", eadf = "eadf")
  m <- rep(base, n)
  if (identical(rule, "auto")) m[!is.na(u) & u >= u_star] <- "lf"
  if (base == "ulf") m[!is.na(undrawn) & undrawn <= 0] <- "lf"
  if (base == "eadf") m[!is.na(drawn) & drawn <= 0] <- "lf"
  m
}

#' Exposure-weighted realised value of a set of rows, by their measure
#' @keywords internal
#' @noRd
.ead_ew <- function(r, m = NULL) {
  if (!nrow(r)) return(NA_real_)
  m <- m %||% r$measure[1]
  # weighted mean of the (floored, capped) realised value, weights = the
  # denominator of the measure; equals the ratio of sums when nothing is floored
  w <- switch(m, ulf = r$undrawn_ref, lf = r$limit_ref, eadf = r$drawn_ref)
  ok <- is.finite(w) & w > 0 & is.finite(r$ccf)
  if (!any(ok)) NA_real_ else sum(w[ok] * r$ccf[ok]) / sum(w[ok])
}

#' Simple and exposure-weighted averages by cohort and measure, plus totals
#' @keywords internal
#' @noRd
.ead_summary <- function(rds) {
  if (!nrow(rds)) return(data.table::data.table(cohort = character(), measure = character(), n = integer(),
                                                ccf_simple = numeric(), ccf_ew = numeric(), ead_realised = numeric(),
                                                share_above_one = numeric()))
  f <- function(r, m) list(n = nrow(r), ccf_simple = mean(r$ccf), ccf_ew = .ead_ew(r, m), ead_realised = sum(r$ead_realised),
                           share_above_one = mean(r$ccf_raw > 1))
  by_c <- rds[, f(.SD, measure), by = c("cohort", "measure")]
  tot <- rds[, f(.SD, measure), by = "measure"][, cohort := "ALL"]
  out <- rbind(by_c, tot, use.names = TRUE)
  data.table::setcolorder(out, c("cohort", "measure"))
  data.table::setorder(out, cohort, measure)
  out[]
}

#' @export
print.scr_ead_data <- function(x, ...) {
  m <- x$meta; l <- x$lra
  cat(sprintf("<scr_ead_data> %s rows from %s default events | %s facilities x %s snapshots\n",
              n_fmt(m$n_rows), n_fmt(m$n_events), n_fmt(m$n_facilities), n_fmt(m$n_snapshots)))
  cat(sprintf("  horizon: %s (%d months) | measure: %s (u* = %s) | reference dates %s to %s (%.1f years)\n",
              m$horizon, m$horizon_months, m$measure_rule, format(m$u_star), format(min(x$rds$ref_date)),
              format(max(x$rds$ref_date)), l$years))
  cat(sprintf("  %s: LRA simple %.4f | exposure-weighted %.4f | n = %s\n", toupper(l$measure), l$simple, l$exposure_weighted, n_fmt(l$n)))
  cat(sprintf("  LF rows %s | above one %s | negative (raw) %s | fast defaults %s\n", fmt_pct(l$share_lf), fmt_pct(l$share_above_one),
              fmt_pct(l$share_negative), fmt_pct(l$share_fast)))
  cat("  funnel:\n")
  for (i in seq_len(nrow(x$funnel))) cat(sprintf("    %-24s %-13s %6d (%s)\n", x$funnel$rule[i], x$funnel$action[i],
                                                 x$funnel$n[i], fmt_pct(x$funnel$share[i])))
  if (is.finite(l$years) && l$years < 5) cat("  note: fewer than five years of reference dates; the average is not a long-run one yet\n")
  invisible(x)
}

# -- 2. Pools ----------------------------------------------------------------- #

#' Estimate CCF pools from the reference data set
#'
#' Splits the reference data set by reference date (the most recent
#' cohorts form the hold-out), bins every driver against the realised CCF
#' with the continuous binner ([scr_bin_continuous()]) on the training
#' rows, revalidates the frozen bins on the hold-out and admits a driver
#' when it passes the named rules `TOO_FEW_DEFAULTS`, `NO_SEPARATION`,
#' `NOT_MONOTONIC` and `UNSTABLE_HOLDOUT`. The cells of the cross of the
#' admitted drivers are ordered by their predicted CCF and merged, adjacent
#' cells first, down to `config$ccf_n_pools` pools with at least
#' `config$ccf_min_defaults` defaults each. Rows in the limit-factor
#' measure form their own pool `LF`.
#'
#' Per pool the estimate is the long-run (default-weighted) average of the
#' realised values on the training rows, `lra`; `moc_est` is the one-sided
#' normal estimation-error margin at `config$ccf_moc_alpha`;
#' `ccf_dt` is the downturn value (equal to `lra` until
#' [scr_ead_downturn()] is run); `ccf_final = max(lra, ccf_dt) + moc_est`;
#' `ccf_floor = params$ccf_floor_fraction * config$ccf_sa_ccf`; and
#' `ccf_applied = max(ccf_final, ccf_floor)`. For the `LF` pool the floor
#' depends on the utilisation and is applied per row by [scr_apply()].
#'
#' @param x An [scr_ead_data()] object.
#' @param drivers Column names of the candidate drivers (columns of
#'   `x$rds`).
#' @param config A [scr_config()].
#' @param holdout Hold-out share, by whole reference dates.
#' @param params An [scr_irb_params()] object; `NULL` uses the preset of
#'   `config$framework`.
#'
#' @return An object of class `scr_ead`: `pools` (the pool table),
#'   `cells` (every cell of the cross with its pool), `bins` (the
#'   `obwoe`-shaped fit of the admitted drivers), `bins_all` (the fit of
#'   every driver), `drivers` (admission table), `holdout` (frozen bins on
#'   the hold-out), `rds` (the reference rows with `sample` and `pool`),
#'   `metrics` (per sample: `rmse`, `mae`, `gauc` with a bootstrap interval,
#'   `spearman`, `ead_rmse`, `ead_mae`, `adequacy`, `cear`), `split`,
#'   `funnel`, `data_summary`, `downturn`, `ledger`, `model_card`,
#'   `params`, `config`, `meta`.
#'   Also `survivors` (the admitted drivers) and `lra` (the long-run averages
#'   of the data set); `metrics` also carries `n`, `n_main`, `gauc_se`,
#'   `somers_d` and `share_floor_binding`.
#'
#' @family irb-ead
#' @examples
#' cfg <- scr_config(verbose = FALSE, n_boot = 20, nthread = 1)
#' ed <- scr_ead_data(scr_demo_ead, facility_id = "facility_id", date_col = "ref_date",
#'                    limit = "limit", drawn = "drawn", defaulted = "defaulted",
#'                    drivers = c("product", "months_on_book", "dpd"), config = cfg)
#' m <- scr_ead(ed, drivers = c("utilisation_ref", "product", "months_on_book"), config = cfg)
#' m
#' m$pools
#' m$drivers
#' @export
scr_ead <- function(x, drivers, config = scr_config(), holdout = 0.3, params = NULL) {
  if (!inherits(x, "scr_ead_data")) stop("scr_ead(): `x` must come from scr_ead_data().", call. = FALSE)
  check_config(config, "scr_ead")
  cfg <- config
  old <- scr_verbose(isTRUE(cfg$verbose)); on.exit(scr_verbose(old), add = TRUE)
  params <- .check_params(params %||% scr_irb_params(cfg$framework), "scr_ead")
  .scr_num1(holdout, "holdout", lower = 0, upper = 1, open_lower = TRUE)
  rds <- data.table::copy(x$rds)
  miss <- setdiff(drivers, names(rds))
  if (length(miss)) stop("scr_ead(): driver(s) not in the reference data set: ", lst(miss), call. = FALSE)
  if (!length(drivers)) stop("scr_ead(): give at least one driver.", call. = FALSE)
  t0 <- Sys.time()
  main <- x$meta$main_measure
  min_def <- as.integer(cfg$ccf_min_defaults)
  msg_stage(11, sprintf("CCF pools (%d driver(s), target %d pools, >= %d defaults per pool)", length(drivers), cfg$ccf_n_pools, min_def))

  # -- split by reference date ------------------------------------------ #
  sp <- split_train_holdout(rds, "ccf", "ref_date", holdout, cfg$seed)
  rds[, sample := "train"]; rds[sp$holdout_idx, sample := "holdout"]
  tr <- rds[sample == "train"]; ho <- rds[sample == "holdout"]
  tr_m <- tr[measure == main]; ho_m <- ho[measure == main]
  if (nrow(tr_m) < 2L * min_def) {
    stop(sprintf("scr_ead(): %d training rows in the %s measure, fewer than twice ccf_min_defaults (%d).", nrow(tr_m), main, min_def), call. = FALSE)
  }
  msg("  split by reference date: train %s (%s in %s) | hold-out %s (%s in %s)%s", n_fmt(nrow(tr)), n_fmt(nrow(tr_m)), main,
      n_fmt(nrow(ho)), n_fmt(nrow(ho_m)), main, if (!is.na(sp$cutoff)) sprintf(" from %s", sp$cutoff) else "")

  # -- driver bins on train, revalidated frozen on hold-out --------------- #
  for (d in drivers) if (is.numeric(tr_m[[d]]) && anyNA(rds[[d]])) {
    stop("scr_ead(): driver '", d, "' has missing values; impute before pooling.", call. = FALSE)
  }
  fit_all <- time_it("driver bins (train)", .scr_bin_continuous(tr_m, "ccf", drivers, min_bins = 2L, max_bins = 6L,
                                                                 min_share = 0.05, min_n = min_def, monotone = "auto",
                                                                 scale = "mean", nthread = cfg$nthread))
  hold <- if (nrow(ho_m)) .cbins_holdout(fit_all, tr_m, ho_m, "ccf", alpha = cfg$psi_alpha) else NULL
  adm <- .ead_admission(fit_all, hold, tr_m, min_def, cfg$psi_alpha)
  for (i in seq_len(nrow(adm))) msg("  %-24s %d bins | eta2 %.3f | %s", adm$feature[i], adm$n_bins[i], adm$eta2[i],
                                    if (adm$admitted[i]) "admitted" else adm$reason[i])
  survivors <- adm[admitted == TRUE, feature]
  ledger <- data.table::copy(x$ledger)
  if (!length(survivors)) {
    msg("  no driver admitted: a single pool carries the long-run average")
    ledger <- rbind(ledger, data.table::data.table(step = "pools", action = "single_pool", detail = "no driver passed the admission rules",
                                                   reason = "NO_DRIVER_ADMITTED", date = format(Sys.Date())))
  }
  fit <- fit_all
  fit$results <- fit_all$results[survivors]
  fit$summary <- fit_all$summary[fit_all$summary$feature %in% survivors, , drop = FALSE]
  fit$n_features <- length(survivors)

  # -- cells and pools ---------------------------------------------------- #
  pl <- .ead_pools(fit, tr_m, survivors, as.integer(cfg$ccf_n_pools), min_def)
  cells <- pl$cells
  # every reference row gets its pool (main measure through the cells, the rest to LF)
  rds[, pool := .ead_pool_of(rds, fit, survivors, cells, main, x$meta)]
  tr <- rds[sample == "train"]

  pools <- .ead_pool_table(tr, main, cfg, params)
  msg("  %d pool(s)%s: applied CCF from %.4f to %.4f", sum(pools$pool != "LF"), if ("LF" %in% pools$pool) " + LF" else "",
      min(pools$ccf_applied), max(pools$ccf_applied))
  model <- structure(list(pools = pools, cells = cells, bins = fit, bins_all = fit_all, drivers = adm, holdout = hold,
                          survivors = survivors, meta = c(x$meta, list(u_star = cfg$ccf_u_star, sa_ccf = cfg$ccf_sa_ccf,
                                                                       floor_fraction = params$ccf_floor_fraction,
                                                                       floor = params$ccf_floor_fraction * cfg$ccf_sa_ccf,
                                                                       split_cutoff = sp$cutoff, n_train = nrow(tr), n_holdout = nrow(ho))),
                          config = cfg, params = params), class = c("scr_ead", "list"))
  metrics <- time_it("metrics (train and hold-out, bootstrap CI)", .ead_metrics(model, rds))
  model$rds <- rds
  model$metrics <- metrics
  model$split <- list(method = sp$method, cutoff = sp$cutoff, ratio = holdout, n_train = nrow(tr), n_holdout = nrow(ho))
  model$funnel <- x$funnel; model$data_summary <- x$summary; model$lra <- x$lra
  model$downturn <- NULL
  ledger <- rbind(ledger, data.table::data.table(
    step = "pools", action = "fit",
    detail = sprintf("drivers admitted: %s; %d pools; MoC alpha %s; floor %s x %s = %s%s",
                     if (length(survivors)) paste(survivors, collapse = ", ") else "(none)", sum(pools$pool != "LF"),
                     format(cfg$ccf_moc_alpha), format(params$ccf_floor_fraction), format(cfg$ccf_sa_ccf),
                     format(params$ccf_floor_fraction * cfg$ccf_sa_ccf), if (isTRUE(params$modified)) "; params_modified = TRUE" else ""),
    reason = "configuration", date = format(Sys.Date())))
  model$ledger <- ledger
  model$model_card <- .ead_model_card(model)
  msg("  done (%.2fs)", as.numeric(difftime(Sys.time(), t0, units = "secs")))
  model
}

#' Admission rules per driver
#' @keywords internal
#' @noRd
.ead_admission <- function(fit_all, hold, tr_m, min_def, alpha) {
  s <- data.table::as.data.table(fit_all$summary)
  idx <- .cbins_apply_idx(fit_all, tr_m)
  y <- tr_m$ccf
  rows <- lapply(s$feature, function(f) {
    e <- fit_all$results[[f]]
    k <- length(e$bin)
    reasons <- character()
    if (nrow(tr_m) < 2L * min_def || any(e$count < min_def)) reasons <- c(reasons, "TOO_FEW_DEFAULTS")
    p_f <- NA_real_
    if (k >= 2L) {
      g <- factor(idx[[f]], levels = seq_len(k))
      ok <- !is.na(g)
      if (data.table::uniqueN(g[ok]) >= 2L) {
        a <- tryCatch(stats::anova(stats::lm(y[ok] ~ g[ok])), error = function(err) NULL)
        if (!is.null(a)) p_f <- a[["Pr(>F)"]][1]
      }
    }
    if (k < 2L || is.na(p_f) || p_f > alpha) reasons <- c(reasons, "NO_SEPARATION")
    h <- if (!is.null(hold)) hold$summary[feature == f] else NULL
    if (!is.null(h) && nrow(h)) {
      hr <- strsplit(h$holdout_reason, ";", fixed = TRUE)[[1]]
      if ("UNSTABLE_HOLDOUT" %in% hr) reasons <- c(reasons, "NOT_MONOTONIC")
      if (any(c("PSI_ACTION", "UNBINNED_HOLDOUT") %in% hr)) reasons <- c(reasons, "UNSTABLE_HOLDOUT")
    }
    data.table::data.table(feature = f, type = e$type, n_bins = k, eta2 = sum(e$iv), direction = e$direction,
                           p_anova = p_f, eta2_holdout = if (!is.null(h) && nrow(h)) h$eta2_holdout else NA_real_,
                           psi = if (!is.null(h) && nrow(h)) h$psi else NA_real_,
                           psi_flag = if (!is.null(h) && nrow(h)) h$psi_flag else NA_character_,
                           admitted = !length(reasons), reason = if (length(reasons)) paste(reasons, collapse = ";") else "OK")
  })
  data.table::rbindlist(rows)
}

#' Cells of the cross of the admitted drivers, merged into pools by predicted CCF
#' @keywords internal
#' @noRd
.ead_pools <- function(fit, tr_m, survivors, n_pools, min_def) {
  y <- tr_m$ccf
  if (!length(survivors)) {
    cells <- data.table::data.table(cell = 1L, n = length(y), mean = mean(y), pred = mean(y), pool = "P1", seen = TRUE)
    return(list(cells = cells))
  }
  ks <- vapply(survivors, function(f) length(fit$results[[f]]$bin), integer(1))
  if (prod(ks) > 5000) stop("scr_ead(): the cross of the admitted drivers has more than 5,000 cells; drop a driver.", call. = FALSE)
  grid <- do.call(data.table::CJ, c(lapply(ks, seq_len), list(sorted = TRUE)))
  data.table::setnames(grid, survivors)
  grid[, cell := seq_len(.N)]
  idx <- .cbins_apply_idx(fit, tr_m, survivors)
  key <- .ead_cell_key(idx, survivors, ks)
  st <- data.table::data.table(cell = key, y = y)[!is.na(cell), list(n = .N, s = sum(y)), by = "cell"]
  grid <- merge(grid, st, by = "cell", all.x = TRUE)
  grid[is.na(n), `:=`(n = 0L, s = 0)]
  grid[, mean := data.table::fifelse(n > 0L, s / pmax(n, 1L), NA_real_)]
  # additive prediction from the driver bin means, used for small and unseen cells
  add <- Reduce(`+`, lapply(survivors, function(f) fit$results[[f]]$mean[grid[[f]]])) / length(survivors)
  grid[, pred := data.table::fifelse(n >= 10L, mean, add)]
  seen <- grid[n > 0L][order(pred, cell)]
  mg <- .cbin_merge(seen$n, seen$s, max_bins = n_pools, min_n = min_def, min_share = 0)
  seen[, grp := mg$group]
  gm <- seen[, list(m = sum(s) / sum(n)), by = "grp"][order(m)]
  gm[, pool := sprintf("P%d", seq_len(.N))]
  seen <- merge(seen, gm[, list(grp, pool)], by = "grp")
  grid <- merge(grid, seen[, list(cell, pool)], by = "cell", all.x = TRUE)
  # unseen cells: the pool whose mean is nearest to the additive prediction
  un <- is.na(grid$pool)
  if (any(un)) grid$pool[un] <- gm$pool[vapply(grid$pred[un], function(v) which.min(abs(gm$m - v)), integer(1))]
  grid[, seen := n > 0L]
  grid[, s := NULL]
  data.table::setcolorder(grid, c("cell", survivors, "n", "mean", "pred", "pool", "seen"))
  list(cells = grid[order(cell)])
}

#' Integer key of the cell of every row from the driver bin indices (NA when any is NA)
#' @keywords internal
#' @noRd
.ead_cell_key <- function(idx, survivors, ks) {
  # CJ(sorted = TRUE) enumerates with the LAST column varying fastest
  key <- rep(1L, nrow(idx)); mult <- 1L
  for (j in rev(seq_along(survivors))) {
    v <- idx[[survivors[j]]]
    key <- key + (v - 1L) * mult
    mult <- mult * ks[j]
  }
  key
}

#' Pool of every row of a reference data set (or of scored data, see scr_apply)
#' @keywords internal
#' @noRd
.ead_pool_of <- function(d, fit, survivors, cells, main, meta) {
  n <- nrow(d)
  pool <- rep(NA_character_, n)
  # the LF pool collects the limit-factor rows of a model whose main measure
  # is another one; when the main measure is the limit factor itself its rows
  # go through the cells like any other (as the SQL does: LF condition 1 = 0)
  to_lf <- d$measure == "lf" & !identical(main, "lf")
  if (!length(survivors)) {
    pool[!to_lf] <- "P1"
  } else {
    ks <- vapply(survivors, function(f) length(fit$results[[f]]$bin), integer(1))
    idx <- .cbins_apply_idx(fit, d, survivors)
    key <- .ead_cell_key(idx, survivors, ks)
    pool <- cells$pool[match(key, cells$cell)]
    # unbinned rows (missing driver, unseen category): the highest pool, the
    # most conservative one since pools are labelled in increasing mean
    pool[is.na(pool)] <- .ead_top_pool(cells)
  }
  pool[to_lf] <- "LF"
  pool
}

#' Highest pool label (pools are labelled in increasing mean)
#' @keywords internal
#' @noRd
.ead_top_pool <- function(cells) {
  lv <- unique(cells$pool)
  lv[which.max(as.integer(sub("^P", "", lv)))]
}

#' The pool table from the training rows
#' @keywords internal
#' @noRd
.ead_pool_table <- function(tr, main, cfg, params) {
  z <- stats::qnorm(1 - cfg$ccf_moc_alpha)
  floor_v <- params$ccf_floor_fraction * cfg$ccf_sa_ccf
  f <- function(r, ms) {
    n <- nrow(r); m <- if (n) mean(r$ccf) else NA_real_
    se <- if (n >= 2L) stats::sd(r$ccf) / sqrt(n) else NA_real_
    list(n = n, lra = m, lra_ew = .ead_ew(r, ms), se = se, ccf_min = if (n) min(r$ccf) else NA_real_,
         ccf_max = if (n) max(r$ccf) else NA_real_, share_above_one = if (n) mean(r$ccf > 1) else NA_real_)
  }
  pt <- tr[, f(.SD, measure), by = c("pool", "measure")]
  lv <- unique(tr$pool)
  if (!"LF" %in% pt$pool && (identical(cfg$ccf_measure, "auto") || main != "lf")) {
    # no LF row in training: a conservative fully-drawn limit factor
    pt <- rbind(pt, data.table::data.table(pool = "LF", measure = "lf", n = 0L, lra = 1, lra_ew = NA_real_, se = NA_real_,
                                           ccf_min = NA_real_, ccf_max = NA_real_, share_above_one = NA_real_))
  }
  pt[, moc_est := data.table::fifelse(is.na(se), 0, z * se)]
  pt[, ccf_dt := lra]
  pt[, downturn := "none"]
  pt[, ccf_final := pmax(lra, ccf_dt) + moc_est]
  pt[, ccf_floor := data.table::fifelse(pool == "LF", NA_real_, floor_v)]
  pt[, ccf_applied := data.table::fifelse(pool == "LF", ccf_final, pmax(ccf_final, ccf_floor))]
  pt[, floor_binding := !is.na(ccf_floor) & ccf_floor > ccf_final]
  pt[, ord := data.table::fifelse(pool == "LF", Inf, suppressWarnings(as.numeric(sub("^P", "", pool))))]
  data.table::setorder(pt, ord); pt[, ord := NULL]
  data.table::setcolorder(pt, c("pool", "measure", "n", "lra", "lra_ew", "se", "moc_est", "ccf_dt", "downturn", "ccf_final",
                                "ccf_floor", "ccf_applied", "floor_binding"))
  pt[]
}

#' Predicted CCF of pools looked up by label
#' @keywords internal
#' @noRd
.ead_lookup <- function(pools, pool, what = "ccf_applied") pools[[what]][match(pool, pools$pool)]

#' Predicted EAD per row from a pool assignment and the drawn/limit amounts
#' @keywords internal
#' @noRd
.ead_predict_rows <- function(pool, measure, drawn, limit, pools, floor_v) {
  ccf <- .ead_lookup(pools, pool, "ccf_applied")
  undrawn <- pmax(limit - drawn, 0)
  # vector branches through fifelse(): fcase() takes a vector `default` only
  # from data.table 1.15.0, and DESCRIPTION allows 1.14
  ead_model <- data.table::fifelse(measure == "lf", ccf * limit,
                                   data.table::fifelse(measure == "eadf", ccf * drawn, drawn + ccf * undrawn))
  ead_floor <- drawn + floor_v * undrawn
  ead <- pmax(drawn, ead_model, ead_floor)
  list(ccf_applied = ccf, undrawn = undrawn, ead_model = ead_model, ead_floor = ead_floor, ead_predicted = ead,
       ead_floor_binding = ead_floor >= ead_model & undrawn > 0)
}

#' Metrics per sample: CCF error, gAUC with bootstrap CI, EAD error, adequacy, CEAR
#' @keywords internal
#' @noRd
.ead_metrics <- function(model, rds) {
  cfg <- model$config; pools <- model$pools; main <- model$meta$main_measure
  rows <- lapply(c("train", "holdout"), function(sm) {
    r <- rds[sample == sm]
    if (!nrow(r)) return(NULL)
    pr <- .ead_predict_rows(r$pool, r$measure, r$drawn_ref, r$limit_ref, pools, model$meta$floor)
    rm_ <- r$measure == main
    pred <- .ead_lookup(pools, r$pool, "lra")
    g <- .ead_gauc(pred[rm_], r$ccf[rm_], n_boot = cfg$n_boot, level = cfg$ci_level, seed = cfg$seed, nthread = cfg$nthread)
    data.table::data.table(
      sample = sm, n = nrow(r), n_main = sum(rm_),
      rmse = sqrt(mean((r$ccf[rm_] - pred[rm_])^2)), mae = mean(abs(r$ccf[rm_] - pred[rm_])),
      gauc = g$gauc, gauc_lo = g$lo, gauc_hi = g$hi, gauc_se = g$se, somers_d = g$d,
      spearman = if (sum(rm_) > 2L) suppressWarnings(stats::cor(pred[rm_], r$ccf[rm_], method = "spearman")) else NA_real_,
      ead_rmse = sqrt(mean((r$ead_realised - pr$ead_predicted)^2)), ead_mae = mean(abs(r$ead_realised - pr$ead_predicted)),
      adequacy = sum(r$ead_realised) / sum(pr$ead_predicted),
      cear = .ead_cear(.ead_lookup(pools, r$pool, "ccf_applied"), pmax(r$ead_realised - r$drawn_ref, 0)),
      share_floor_binding = mean(pr$ead_floor_binding))
  })
  data.table::rbindlist(rows)
}

#' Somers' D between a prediction and a realised value, over pairs with distinct realised values
#'
#' Concordant when the prediction orders the pair as the realised values do;
#' pairs tied on the prediction count as neither. gAUC = (D + 1) / 2. Pairs
#' are counted exactly in `O(n log n)` by the compiled kernel
#' (`.scr_somers()`), whatever the number of distinct predictions; a
#' constant prediction gives 0.
#' @keywords internal
#' @noRd
.somers_d <- function(pred, y) .scr_somers(pred, y, const_p = 0)

#' gAUC with a percentile bootstrap interval (seeds drawn in the parent)
#' @keywords internal
#' @noRd
.ead_gauc <- function(pred, y, n_boot = 200L, level = 0.95, seed = NULL, nthread = 1L) {
  d <- .somers_d(pred, y)
  out <- list(d = d, gauc = (d + 1) / 2, lo = NA_real_, hi = NA_real_, se = NA_real_)
  n <- length(y)
  if (is.na(d) || n_boot < 2L || n < 3L) return(out)
  .scr_local_seed(seed)
  seeds <- sample.int(.Machine$integer.max, n_boot)
  .scr_rng_guard()   # the replicates re-seed; the user's stream moves by the n_boot draws above only
  reps <- .scr_lapply(seeds, function(sd) {
    set.seed(sd)
    j <- sample.int(n, n, replace = TRUE)
    .somers_d(pred[j], y[j])
  }, nthread = nthread)
  b <- (unlist(reps) + 1) / 2
  a <- (1 - level) / 2
  q <- stats::quantile(b, c(a, 1 - a), na.rm = TRUE, names = FALSE)
  out$lo <- q[1]; out$hi <- q[2]; out$se <- stats::sd(b, na.rm = TRUE)
  out
}

#' Cumulative EAD accuracy ratio: share of realised additional drawing captured by the predicted ranking
#' @keywords internal
#' @noRd
.ead_cear <- function(pred, add) {
  ok <- is.finite(pred) & is.finite(add); pred <- pred[ok]; add <- add[ok]
  n <- length(add); tot <- sum(add)
  if (n < 2L || tot <= 0) return(NA_real_)
  area <- function(p) {
    # tie groups on p share their total linearly: the curve is evaluated at the group boundaries
    g <- data.table::data.table(p = p, a = add)[, list(n = .N, a = sum(a)), by = "p"][order(-p)]
    xs <- c(0, cumsum(g$n) / n); ys <- c(0, cumsum(g$a) / tot)
    sum(diff(xs) * (ys[-1] + ys[-length(ys)]) / 2)
  }
  a_model <- area(pred); a_perfect <- area(add)
  if (a_perfect <= 0.5) return(NA_real_)
  (a_model - 0.5) / (a_perfect - 0.5)
}

#' @keywords internal
#' @noRd
.ead_model_card <- function(m) {
  mt <- m$metrics; ho <- mt[sample == "holdout"]
  list(
    package = sprintf("scorecraft %s", as.character(utils::packageVersion("scorecraft"))),
    fitted_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    horizon = m$meta$horizon, horizon_months = m$meta$horizon_months, measure_rule = m$meta$measure_rule,
    main_measure = m$meta$main_measure, u_star = m$meta$u_star, floor_realised = m$meta$floor_realised,
    cap_realised = m$meta$cap_realised, post_default_drawings_in = m$meta$post_default_drawings_in,
    framework = m$params$framework, params_modified = isTRUE(m$params$modified),
    sa_ccf = m$meta$sa_ccf, floor_fraction = m$meta$floor_fraction, ccf_floor = m$meta$floor,
    moc_alpha = m$config$ccf_moc_alpha, downturn = if (is.null(m$downturn)) "none" else m$downturn$method,
    n_events = m$meta$n_events, n_rows = m$meta$n_rows, n_train = m$meta$n_train, n_holdout = m$meta$n_holdout,
    split_cutoff = m$meta$split_cutoff, years_of_data = m$lra$years,
    drivers_candidate = paste(m$drivers$feature, collapse = ", "),
    drivers_admitted = if (length(m$survivors)) paste(m$survivors, collapse = ", ") else "(none)",
    n_pools = sum(m$pools$pool != "LF"), lf_pool = "LF" %in% m$pools$pool,
    lra_simple = m$lra$simple, lra_exposure_weighted = m$lra$exposure_weighted,
    ccf_applied_min = min(m$pools$ccf_applied), ccf_applied_max = max(m$pools$ccf_applied),
    gauc_train = mt[sample == "train", gauc], gauc_holdout = if (nrow(ho)) ho$gauc else NA_real_,
    rmse_holdout = if (nrow(ho)) ho$rmse else NA_real_, adequacy_holdout = if (nrow(ho)) ho$adequacy else NA_real_,
    references = "Basel Framework CRE32.29, CRE32.36, CRE36.92-36.97; EBA/GL/2017/16 section 4.4 (margin of conservatism)"
  )
}

#' @export
print.scr_ead <- function(x, ...) {
  m <- x$meta; p <- x$pools
  cat(sprintf("<scr_ead> %d pool(s)%s from %s reference rows | %s horizon (%d months) | measure %s\n",
              sum(p$pool != "LF"), if ("LF" %in% p$pool) " + LF" else "", n_fmt(nrow(x$rds)), m$horizon, m$horizon_months, m$measure_rule))
  cat(sprintf("  split by reference date: train %s | hold-out %s%s | drivers admitted: %s\n", n_fmt(m$n_train), n_fmt(m$n_holdout),
              if (!is.na(m$split_cutoff)) sprintf(" (from %s)", m$split_cutoff) else "",
              if (length(x$survivors)) paste(x$survivors, collapse = ", ") else "(none)"))
  cat(sprintf("  floor %.4f (= %s x SA-CCF %s) | MoC alpha %s | downturn %s\n", m$floor, format(m$floor_fraction), format(m$sa_ccf),
              format(x$config$ccf_moc_alpha), if (is.null(x$downturn)) "none" else x$downturn$method))
  cat(sprintf("  %-5s %-5s %6s %8s %8s %8s %8s %8s %8s %8s\n", "pool", "meas", "n", "lra", "lra_ew", "moc", "ccf_dt", "final", "floor", "applied"))
  for (i in seq_len(nrow(p))) {
    cat(sprintf("  %-5s %-5s %6d %8.4f %8s %8.4f %8.4f %8.4f %8s %8.4f%s\n", p$pool[i], p$measure[i], p$n[i], p$lra[i],
                if (is.na(p$lra_ew[i])) "-" else sprintf("%.4f", p$lra_ew[i]), p$moc_est[i], p$ccf_dt[i], p$ccf_final[i],
                if (is.na(p$ccf_floor[i])) "row" else sprintf("%.4f", p$ccf_floor[i]), p$ccf_applied[i],
                if (isTRUE(p$floor_binding[i])) " *" else ""))
  }
  mt <- x$metrics
  for (i in seq_len(nrow(mt))) {
    cat(sprintf("  %-8s n %5d | RMSE %.4f | MAE %.4f | gAUC %.4f [%.4f, %.4f] | EAD adequacy %.4f | CEAR %s\n", mt$sample[i], mt$n[i],
                mt$rmse[i], mt$mae[i], mt$gauc[i], mt$gauc_lo[i], mt$gauc_hi[i], mt$adequacy[i],
                if (is.na(mt$cear[i])) "-" else sprintf("%.4f", mt$cear[i])))
  }
  if (any(p$floor_binding)) cat("  * the standardised floor binds\n")
  invisible(x)
}

# -- 3. Downturn ---------------------------------------------------------------- #

#' Downturn CCF per pool
#'
#' Quantifies the downturn component of the CCF from user-supplied downturn
#' periods. `"type1"` (observed impact) takes, per pool, the default-weighted
#' average of the realised values of the training events whose default date
#' falls in the periods (the hold-out stays independent) and sets
#' `ccf_dt = max(lra, observed)`; `"type3"` (long-run average plus add-on)
#' sets `ccf_dt = lra + add_on`; `"none"` resets
#' `ccf_dt = lra`. The pool table is recomputed (`ccf_final`, `ccf_applied`)
#' and the ledger records the periods, the method and the reason.
#'
#' @param x An [scr_ead()] object.
#' @param periods A `data.frame` with `start` and `end` dates of the
#'   downturn periods (needed for `"type1"`).
#' @param method `"type1"`, `"type3"` or `"none"`; `NULL` uses
#'   `config$ccf_downturn`.
#' @param add_on Add-on of the `"type3"` method, in CCF units.
#' @param reason Text justifying the periods and the method; mandatory.
#'
#' @return The `scr_ead` object with `downturn` (a list with `method`,
#'   `periods`, `add_on` and the per-pool `table`: `pool`, `lra`,
#'   `n_downturn`, `dt_observed`, `dt_type3`, `ccf_dt`), the updated
#'   `pools` and a new ledger row.
#'   The table also carries `ccf_final` and `ccf_applied`; the object
#'   `n_rows_in_periods` and `reason`.
#'
#' @family irb-ead
#' @examples
#' cfg <- scr_config(verbose = FALSE, n_boot = 20, nthread = 1)
#' ed <- scr_ead_data(scr_demo_ead, facility_id = "facility_id", date_col = "ref_date",
#'                    limit = "limit", drawn = "drawn", defaulted = "defaulted",
#'                    drivers = c("product", "months_on_book"), config = cfg)
#' m <- scr_ead(ed, drivers = c("utilisation_ref", "product"), config = cfg)
#' m2 <- scr_ead_downturn(m, periods = data.frame(start = as.Date("2024-01-01"),
#'                                                end = as.Date("2024-12-01")),
#'                        reason = "2024 chosen as the stress year of the demo panel")
#' m2$downturn$table
#' @export
scr_ead_downturn <- function(x, periods = NULL, method = NULL, add_on = 0.15, reason = NULL) {
  if (!inherits(x, "scr_ead")) stop("scr_ead_downturn(): `x` must come from scr_ead().", call. = FALSE)
  method <- method %||% x$config$ccf_downturn
  method <- match.arg(method, c("type1", "type3", "none"))
  if (is.null(reason) || !nzchar(trimws(reason))) stop("scr_ead_downturn(): a `reason` is mandatory (it goes in the ledger).", call. = FALSE)
  .scr_num1(add_on, "add_on", lower = 0)
  cfg <- x$config
  old <- scr_verbose(isTRUE(cfg$verbose)); on.exit(scr_verbose(old), add = TRUE)
  rds <- x$rds
  in_dt <- rep(FALSE, nrow(rds))
  if (!is.null(periods)) {
    pd <- data.table::as.data.table(periods)
    if (!all(c("start", "end") %in% names(pd))) stop("scr_ead_downturn(): `periods` needs `start` and `end` columns.", call. = FALSE)
    pd[, `:=`(start = as.Date(start), end = as.Date(end))]
    # shared with the LGD downturn: refuses missing dates and end < start
    in_dt <- .lgd_in_periods(rds$default_date, pd)
  } else if (identical(method, "type1")) {
    stop("scr_ead_downturn(): `periods` is needed for the observed-impact method (type1).", call. = FALSE)
  }
  # training rows only, like the long-run averages: the hold-out stays independent
  in_dt <- in_dt & rds$sample == "train"
  obs <- rds[in_dt, list(n_downturn = .N, dt_observed = mean(ccf)), by = "pool"]
  p <- data.table::copy(x$pools)
  p <- merge(p, obs, by = "pool", all.x = TRUE, sort = FALSE)
  p[is.na(n_downturn), n_downturn := 0L]
  p[, dt_type3 := lra + add_on]
  p[, ccf_dt := switch(method,
    type1 = data.table::fifelse(is.na(dt_observed), lra, pmax(lra, dt_observed)),
    type3 = dt_type3,
    none = lra)]
  p[, downturn := method]
  p[, ccf_final := pmax(lra, ccf_dt) + moc_est]
  p[, ccf_applied := data.table::fifelse(pool == "LF", ccf_final, pmax(ccf_final, ccf_floor))]
  p[, floor_binding := !is.na(ccf_floor) & ccf_floor > ccf_final]
  tab <- p[, list(pool, lra, n_downturn, dt_observed, dt_type3, ccf_dt, ccf_final, ccf_applied)]
  p[, c("n_downturn", "dt_observed", "dt_type3") := NULL]
  x$pools <- p
  x$downturn <- list(method = method, periods = if (is.null(periods)) NULL else pd, add_on = add_on, table = tab,
                     n_rows_in_periods = sum(in_dt), reason = reason)
  x$metrics <- .ead_metrics(x, rds)
  x$ledger <- rbind(x$ledger, data.table::data.table(
    step = "downturn", action = method,
    detail = sprintf("periods: %s; %d training reference rows in the periods; add-on %s; applied CCF now %.4f to %.4f",
                     if (is.null(periods)) "(none)" else paste(sprintf("%s to %s", format(pd$start), format(pd$end)), collapse = ", "),
                     sum(in_dt), format(add_on), min(p$ccf_applied), max(p$ccf_applied)),
    reason = reason, date = format(Sys.Date())))
  x$model_card <- .ead_model_card(x)
  msg("  downturn %s: %d training rows in the periods; applied CCF %.4f to %.4f", method, sum(in_dt), min(p$ccf_applied), max(p$ccf_applied))
  x
}

# -- 4. Application ----------------------------------------------------------- #

#' @rdname scr_apply
#' @export
scr_apply.scr_ead <- function(x, newdata, what = c("all", "ead", "pool"), ...) {
  what <- match.arg(what)
  dt <- data.table::as.data.table(newdata)
  cols <- x$meta$cols
  need <- c(cols$limit, cols$drawn, setdiff(x$survivors, "utilisation_ref"))
  miss <- setdiff(need, names(dt))
  if (length(miss)) stop("newdata lacks the column(s): ", lst(miss), call. = FALSE)
  d <- data.table::data.table(limit = as.double(dt[[cols$limit]]), drawn = as.double(dt[[cols$drawn]]))
  d[, utilisation := data.table::fifelse(limit > 0, drawn / limit, NA_real_)]
  d[, undrawn := pmax(limit - drawn, 0)]
  for (f in setdiff(x$survivors, "utilisation_ref")) data.table::set(d, j = f, value = dt[[f]])
  if ("utilisation_ref" %in% x$survivors) d[, utilisation_ref := utilisation]
  d[, measure := .ead_measure(utilisation, limit - drawn, drawn, x$meta$measure_rule, x$meta$u_star)]
  d[, pool := .ead_pool_of(d, x$bins, x$survivors, x$cells, x$meta$main_measure, x$meta)]
  pr <- .ead_predict_rows(d$pool, d$measure, d$drawn, d$limit, x$pools, x$meta$floor)
  out <- data.table::data.table(pool = d$pool, measure = d$measure, utilisation = d$utilisation, undrawn = d$undrawn,
                                ccf_applied = pr$ccf_applied, ead_model = pr$ead_model, ead_floor = pr$ead_floor,
                                ead_predicted = pr$ead_predicted, ead_floor_binding = pr$ead_floor_binding)
  keep <- switch(what,
    pool = c("pool", "measure"),
    ead  = c("pool", "measure", "ccf_applied", "ead_predicted", "ead_floor_binding"),
    all  = names(out))
  out[, keep, with = FALSE][]
}

# -- 5. Production SQL -------------------------------------------------------- #

#' @rdname scr_sql
#' @export
scr_sql.scr_ead <- function(x, table = NULL, dialect = NULL, file = NULL, ...) {
  .scr_local_verbose(x)
  cfg <- x$config
  if (!is.null(table)) cfg$sql_table <- table
  if (!is.null(dialect)) cfg$sql_dialect <- dialect
  dl <- cfg$sql_dialect
  dialects <- c("ansi", "databricks", "spark", "hive", "mysql", "mariadb", "sqlserver", "bigquery", "postgres",
                "oracle", "snowflake", "redshift", "duckdb", "sqlite")
  if (!dl %in% dialects) stop("unknown SQL dialect '", dl, "'.", call. = FALSE)
  cols <- x$meta$cols
  lim <- .sql_ident(cols$limit, dl); drw <- .sql_ident(cols$drawn, dl)
  keep <- cfg$sql_keep_columns
  surv <- x$survivors
  raw_drivers <- setdiff(surv, "utilisation_ref")
  greatest <- function(...) .sql_greatest(c(...), dl)
  # block 1: utilisation and undrawn
  base <- c(if (length(keep)) sprintf("    %s,", keep),
            sprintf("    %s AS limit_amt,", lim), sprintf("    %s AS drawn_amt,", drw),
            sprintf("    CASE WHEN %s > 0 THEN %s / %s ELSE NULL END AS utilisation,", lim, drw, lim),
            sprintf("    CASE WHEN %s > %s THEN %s - %s ELSE 0 END AS undrawn%s", lim, drw, lim, drw, if (length(raw_drivers)) "," else ""),
            if (length(raw_drivers)) paste0("    ", vapply(raw_drivers, .sql_ident, character(1), dl), collapse = ",\n"))
  # block 2: driver bin index from the frozen cut points
  fit <- x$bins
  if ("utilisation_ref" %in% surv) {
    # the driver lives on the derived column: rename the feature for the engine
    fit$results$utilisation_ref$feature <- "utilisation"
    names(fit$results)[names(fit$results) == "utilisation_ref"] <- "utilisation"
    fit$summary$feature[fit$summary$feature == "utilisation_ref"] <- "utilisation"
  }
  feats_sql <- sub("^utilisation_ref$", "utilisation", surv)
  idx_exprs <- if (length(surv)) {
    ix <- OptimalBinningWoE::obwoe_sql(obj = fit, table = "base_ead", features = feats_sql, output = "index",
                                       style = "select", dialect = dl, digits = NULL, comment = FALSE,
                                       bin_separator = cfg$bin_separator)
    .sql_select_exprs(ix)
  } else character()
  # block 3: pool from the cells, LF branch from the measure rule
  lf_cond <- .sql_lf_condition(x$meta$measure_rule, x$meta$u_star)
  pool_case <- if (length(surv)) {
    cells <- x$cells
    whens <- vapply(unique(cells$pool), function(pl) {
      cc <- cells[pool == pl]
      conds <- vapply(seq_len(nrow(cc)), function(i) {
        paste(sprintf("%s_idx = %d", feats_sql, unlist(cc[i, surv, with = FALSE])), collapse = " AND ")
      }, character(1))
      sprintf("WHEN %s THEN %s", paste0("(", conds, ")", collapse = " OR "), .sql_str(pl))
    }, character(1))
    sprintf("CASE WHEN %s THEN 'LF' %s ELSE %s END", lf_cond, paste(whens, collapse = " "), .sql_str(.ead_top_pool(cells)))
  } else sprintf("CASE WHEN %s THEN 'LF' ELSE 'P1' END", lf_cond)
  p <- x$pools
  ccf_case <- sprintf("CASE pool %s ELSE NULL END", paste(sprintf("WHEN %s THEN %s", .sql_str(p$pool), .sql_num(p$ccf_applied)), collapse = " "))
  floor_v <- .sql_num(x$meta$floor)
  main <- x$meta$main_measure
  model_expr <- switch(main,
    ulf = "drawn_amt + ccf_applied * undrawn",
    lf = "ccf_applied * limit_amt",
    eadf = "ccf_applied * drawn_amt")
  ead_expr <- sprintf("CASE WHEN pool = 'LF' THEN %s ELSE %s END",
                      greatest("drawn_amt", "ccf_applied * limit_amt", sprintf("drawn_amt + %s * undrawn", floor_v)),
                      greatest("drawn_amt", model_expr, sprintf("drawn_amt + %s * undrawn", floor_v)))
  final <- c("SELECT",
             if (length(keep)) sprintf("    %s,", keep),
             "    limit_amt, drawn_amt, utilisation, undrawn, pool, ccf_applied,",
             sprintf("    %s AS ead_predicted", ead_expr),
             "FROM (", "  SELECT", "    *,",
             sprintf("    %s AS ccf_applied", ccf_case),
             "  FROM pool_ead", ") ead;")
  c("-- =============================================================",
    sprintf("-- scorecraft | EAD/CCF pools | %d pool(s)%s | %s horizon (%d months) | dialect: %s", sum(p$pool != "LF"),
            if ("LF" %in% p$pool) " + LF" else "", x$meta$horizon, x$meta$horizon_months, dl),
    sprintf("-- Generated on %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf("-- measure rule %s | u* = %s | floor %s (= %s x SA-CCF %s)", x$meta$measure_rule, format(x$meta$u_star), floor_v,
            format(x$meta$floor_fraction), format(x$meta$sa_ccf)),
    "-- Block 1 (CTE base_ead): utilisation and undrawn amount.",
    "-- Block 2 (CTE bins_ead): driver bin index from the frozen cut points (OptimalBinningWoE::obwoe_sql()).",
    "-- Block 3 (CTE pool_ead): pool from the cells; LF branch where the limit factor applies.",
    "-- Block 4: applied CCF and predicted EAD = GREATEST(drawn, model, drawn + floor x undrawn).",
    "-- =============================================================", "",
    "WITH base_ead AS (", "  SELECT", base, sprintf("  FROM %s", cfg$sql_table), "),",
    "bins_ead AS (", "  SELECT", if (length(idx_exprs)) c("    *,", paste0("    ", idx_exprs)) else "    *", "  FROM base_ead", "),",
    "pool_ead AS (", "  SELECT", "    *,", sprintf("    %s AS pool", pool_case), "  FROM bins_ead", ")", "", final) |>
    .sql_lines() |> .sql_out(file)
}

#' SQL condition of the limit-factor branch under the measure rule
#' @keywords internal
#' @noRd
.sql_lf_condition <- function(rule, u_star) {
  switch(rule,
    auto = sprintf("undrawn <= 0 OR utilisation >= %s", .sql_num(u_star)),
    ulf  = "undrawn <= 0",
    eadf = "drawn_amt <= 0",
    lf   = "1 = 0")
}

#' GREATEST() of several expressions in every dialect
#' @keywords internal
#' @noRd
.sql_greatest <- function(exprs, dialect) {
  if (identical(dialect, "sqlite")) return(sprintf("MAX(%s)", paste(exprs, collapse = ", ")))
  if (identical(dialect, "sqlserver")) {
    return(sprintf("(SELECT MAX(v) FROM (VALUES %s) AS g(v))", paste(sprintf("(%s)", exprs), collapse = ", ")))
  }
  sprintf("GREATEST(%s)", paste(exprs, collapse = ", "))
}

#' Quote an identifier when it is a reserved word or not a plain name
#' @keywords internal
#' @noRd
.sql_ident <- function(x, dialect) {
  reserved <- c("limit", "default", "order", "group", "select", "from", "where", "date", "user", "table", "index",
                "key", "value", "values", "desc", "asc", "count", "offset", "case", "when", "end", "level", "rank",
                "row", "rows", "start", "type", "position", "partition", "current", "primary", "check", "column")
  plain <- grepl("^[A-Za-z_][A-Za-z0-9_]*$", x)
  need <- !plain | tolower(x) %in% reserved
  q <- if (dialect %in% c("mysql", "mariadb", "spark", "hive", "databricks", "bigquery")) c("`", "`")
       else if (identical(dialect, "sqlserver")) c("[", "]") else c("\"", "\"")
  ifelse(need, paste0(q[1], x, q[2]), x)
}

# -- 6. Validation --------------------------------------------------------- #

#' Validate CCF pools: calibration, discrimination, back-testing and stability
#'
#' Per pool and in total, compares realised and predicted values on the
#' validation rows (the hold-out of the model by default): simple and
#' exposure-weighted averages, the one-sided t-test of realised above
#' predicted (under-estimation) with its p-value, the EAD adequacy ratio
#' (sum of realised EAD over sum of predicted EAD) and traffic lights
#' (red at or below `lights[1]`, amber at or below `lights[2]`, green above;
#' adequacy green at or below `adequacy_lights[1]`, amber up to
#' `adequacy_lights[2]`, red above). Adds the
#' discrimination block (gAUC with a bootstrap interval against the
#' development value, Spearman correlation, cumulative EAD accuracy
#' ratio), the back-test by cohort and the stability of the pool
#' distribution and of the driver bins ([scr_psi()], fixed and
#' sample-size-adjusted thresholds). The numeric limits of the lights are
#' a convention of the package, stated as such in the output.
#'
#' @param x An [scr_ead()] object.
#' @param newdata `NULL` (the hold-out rows of `x`), an [scr_ead_data()]
#'   object or its `rds` table.
#' @param lights Two increasing p-value thresholds: red at or below the
#'   first, amber at or below the second.
#' @param adequacy_lights Two increasing adequacy-ratio thresholds.
#'
#' @return An object of class `scr_ead_validation`: `calibration`,
#'   `discrimination`, `backtest`, `stability`, `summary` (test,
#'   statistic, p, light), `n`, `source`.
#'
#' @family irb-ead
#' @examples
#' cfg <- scr_config(verbose = FALSE, n_boot = 20, nthread = 1)
#' ed <- scr_ead_data(scr_demo_ead, facility_id = "facility_id", date_col = "ref_date",
#'                    limit = "limit", drawn = "drawn", defaulted = "defaulted",
#'                    drivers = c("product", "months_on_book"), config = cfg)
#' m <- scr_ead(ed, drivers = c("utilisation_ref", "product"), config = cfg)
#' v <- scr_ead_validate(m)
#' v
#' v$calibration
#' @export
scr_ead_validate <- function(x, newdata = NULL, lights = c(0.01, 0.05), adequacy_lights = c(1.00, 1.05)) {
  if (!inherits(x, "scr_ead")) stop("scr_ead_validate(): `x` must come from scr_ead().", call. = FALSE)
  if (length(lights) != 2L || any(!is.finite(lights)) || lights[1] >= lights[2]) stop("`lights` must be two increasing p-values.", call. = FALSE)
  cfg <- x$config; main <- x$meta$main_measure; pools <- x$pools
  if (is.null(newdata)) {
    v <- x$rds[sample == "holdout"]; source <- "holdout"
    # copy: `:=` below must not add columns to the model's own table
    if (!nrow(v)) { v <- data.table::copy(x$rds); source <- "train (no hold-out)" }
  } else {
    v <- if (inherits(newdata, "scr_ead_data")) data.table::copy(newdata$rds) else data.table::as.data.table(newdata)
    need <- c("ccf", "measure", "ead_realised", "limit_ref", "drawn_ref", "undrawn_ref", "ref_date", "cohort", setdiff(x$survivors, names(v)))
    miss <- setdiff(need, names(v))
    if (length(miss)) stop("scr_ead_validate(): newdata lacks the column(s): ", lst(miss), call. = FALSE)
    v[, pool := .ead_pool_of(v, x$bins, x$survivors, x$cells, main, x$meta)]
    source <- "newdata"
  }
  light_p <- function(p) data.table::fcase(is.na(p), NA_character_, p <= lights[1], "red", p <= lights[2], "amber", default = "green")
  light_a <- function(a) data.table::fcase(is.na(a), NA_character_, a <= adequacy_lights[1], "green", a <= adequacy_lights[2], "amber", default = "red")
  pr <- .ead_predict_rows(v$pool, v$measure, v$drawn_ref, v$limit_ref, pools, x$meta$floor)
  v[, `:=`(predicted = .ead_lookup(pools, pool, "ccf_applied"), ead_predicted = pr$ead_predicted)]
  calib_row <- function(r, label) {
    rm_ <- if (label == "TOTAL") r$measure == main else rep(TRUE, nrow(r))
    n <- sum(rm_)
    real <- if (n) mean(r$ccf[rm_]) else NA_real_
    pred <- if (n) mean(r$predicted[rm_]) else NA_real_
    se <- if (n >= 2L) stats::sd(r$ccf[rm_]) / sqrt(n) else NA_real_
    t <- if (!is.na(se) && se > 0) (real - pred) / se else NA_real_
    p <- if (!is.na(t)) stats::pt(t, df = n - 1L, lower.tail = FALSE) else NA_real_
    ew_real <- if (n) .ead_ew(r[rm_]) else NA_real_
    ew_pred <- if (n) {
      rr <- r[rm_]; m <- rr$measure[1]
      num <- if (m == "ulf") sum(rr$ead_predicted - rr$drawn_ref) else sum(rr$ead_predicted)
      den <- switch(m, ulf = sum(rr$undrawn_ref), lf = sum(rr$limit_ref), eadf = sum(rr$drawn_ref))
      if (den > 0) num / den else NA_real_
    } else NA_real_
    adequacy <- sum(r$ead_realised) / sum(r$ead_predicted)
    data.table::data.table(pool = label, n = nrow(r), n_main = n, realised = real, predicted = pred, realised_ew = ew_real,
                           predicted_ew = ew_pred, se = se, t = t, p = p, light_p = light_p(p),
                           ead_realised = sum(r$ead_realised), ead_predicted = sum(r$ead_predicted),
                           adequacy = adequacy, light_adequacy = light_a(adequacy))
  }
  by_pool <- data.table::rbindlist(lapply(pools$pool, function(pl) calib_row(v[pool == pl], pl)))
  calibration <- rbind(by_pool, calib_row(v, "TOTAL"))
  # discrimination against development
  rm_ <- v$measure == main
  dev <- x$metrics[sample == "train"]
  g <- .ead_gauc(.ead_lookup(pools, v$pool[rm_], "lra"), v$ccf[rm_], n_boot = cfg$n_boot, level = cfg$ci_level, seed = cfg$seed, nthread = cfg$nthread)
  z <- if (nrow(dev) && is.finite(g$se) && is.finite(dev$gauc_se) && (g$se^2 + dev$gauc_se^2) > 0) (dev$gauc - g$gauc) / sqrt(g$se^2 + dev$gauc_se^2) else NA_real_
  discrimination <- data.table::data.table(
    n = sum(rm_), gauc = g$gauc, gauc_lo = g$lo, gauc_hi = g$hi, gauc_dev = if (nrow(dev)) dev$gauc else NA_real_,
    z_vs_dev = z, p_vs_dev = if (is.na(z)) NA_real_ else stats::pnorm(z, lower.tail = FALSE),
    spearman = if (sum(rm_) > 2L) suppressWarnings(stats::cor(v$predicted[rm_], v$ccf[rm_], method = "spearman")) else NA_real_,
    cear = .ead_cear(v$predicted, pmax(v$ead_realised - v$drawn_ref, 0)))
  discrimination[, light := light_p(p_vs_dev)]
  # back-test by cohort
  backtest <- data.table::rbindlist(lapply(sort(unique(v$cohort)), function(ch) {
    r <- calib_row(v[cohort == ch], "TOTAL"); r[, pool := NULL]; r[, cohort := ch]
    data.table::setcolorder(r, "cohort"); r
  }))
  # stability: pools and driver bins against the training distribution
  tr <- x$rds[sample == "train"]
  ps <- scr_psi(tr$pool, v$pool, levels = pools$pool, alpha = cfg$psi_alpha)
  stab <- data.table::data.table(item = "pool", n_base = ps$n_base, n_compare = ps$n_compare, psi = ps$psi,
                                 flag_fixed = ps$flag_fixed, critical = ps$critical, flag_adjusted = ps$flag_adjusted)
  if (length(x$survivors)) {
    tr_m <- tr[measure == main]; v_m <- v[rm_]
    i_tr <- .cbins_apply_idx(x$bins, tr_m, x$survivors); i_v <- .cbins_apply_idx(x$bins, v_m, x$survivors)
    stab <- rbind(stab, data.table::rbindlist(lapply(x$survivors, function(f) {
      k <- length(x$bins$results[[f]]$bin)
      q <- scr_psi(as.character(i_tr[[f]]), as.character(i_v[[f]]), levels = as.character(seq_len(k)), alpha = cfg$psi_alpha)
      data.table::data.table(item = f, n_base = q$n_base, n_compare = q$n_compare, psi = q$psi, flag_fixed = q$flag_fixed,
                             critical = q$critical, flag_adjusted = q$flag_adjusted)
    })))
  }
  stab[, light := data.table::fcase(is.na(flag_fixed), NA_character_, flag_fixed == "stable", "green", flag_fixed == "moderate", "amber", default = "red")]
  tot <- calibration[pool == "TOTAL"]
  summary <- data.table::data.table(
    test = c("calibration_t_total", "ead_adequacy_total", "gauc_vs_development", "pool_psi"),
    statistic = c(tot$t, tot$adequacy, z, ps$psi),
    p = c(tot$p, NA_real_, discrimination$p_vs_dev, NA_real_),
    light = c(tot$light_p, tot$light_adequacy, discrimination$light, stab$light[1]),
    convention = c(sprintf("one-sided t; red <= %s, amber <= %s (package convention)", lights[1], lights[2]),
                   sprintf("sum realised / sum predicted EAD; amber above %s, red above %s (package convention)", adequacy_lights[1], adequacy_lights[2]),
                   "z of development minus current gAUC over the bootstrap standard errors",
                   "PSI 0.10/0.25 fixed; adjusted: Yurdakul & Naranjo (2020)"))
  structure(list(calibration = calibration, discrimination = discrimination, backtest = backtest, stability = stab,
                 summary = summary, n = nrow(v), source = source, lights = lights, adequacy_lights = adequacy_lights),
            class = c("scr_ead_validation", "list"))
}

#' @export
print.scr_ead_validation <- function(x, ...) {
  cat(sprintf("<scr_ead_validation> %s rows (%s)\n", n_fmt(x$n), x$source))
  cat(sprintf("  %-6s %6s %9s %9s %8s %8s %6s %9s %6s\n", "pool", "n", "realised", "predicted", "t", "p", "light", "adequacy", "light"))
  c <- x$calibration
  for (i in seq_len(nrow(c))) {
    cat(sprintf("  %-6s %6d %9s %9s %8s %8s %-6s %9.4f %-6s\n", c$pool[i], c$n[i],
                if (is.na(c$realised[i])) "-" else sprintf("%.4f", c$realised[i]),
                if (is.na(c$predicted[i])) "-" else sprintf("%.4f", c$predicted[i]),
                if (is.na(c$t[i])) "-" else sprintf("%.3f", c$t[i]), if (is.na(c$p[i])) "-" else sprintf("%.4f", c$p[i]),
                c$light_p[i] %||% "-", c$adequacy[i], c$light_adequacy[i]))
  }
  d <- x$discrimination
  cat(sprintf("  gAUC %.4f [%.4f, %.4f] vs development %.4f (p %s) | Spearman %s | CEAR %s\n", d$gauc, d$gauc_lo, d$gauc_hi, d$gauc_dev,
              if (is.na(d$p_vs_dev)) "-" else sprintf("%.4f", d$p_vs_dev), if (is.na(d$spearman)) "-" else sprintf("%.4f", d$spearman),
              if (is.na(d$cear)) "-" else sprintf("%.4f", d$cear)))
  s <- x$stability
  cat(sprintf("  stability: %s\n", paste(sprintf("%s PSI %.4f (%s)", s$item, s$psi, s$flag_fixed), collapse = " | ")))
  cat(sprintf("  lights: %s\n", paste(sprintf("%s %s", x$summary$test, x$summary$light), collapse = " | ")))
  invisible(x)
}

# -- 7. Export ---------------------------------------------------------------- #

#' @rdname scr_export
#' @export
scr_export.scr_ead <- function(x, dir, stamp = TRUE, validation = NULL, tag = "ccf", ...) {
  .scr_local_verbose(x)
  .need_openxlsx()
  out_dir <- .export_dir(dir, stamp)
  val <- validation %||% scr_ead_validate(x)
  bins_tab <- if (!is.null(x$holdout)) data.table::copy(x$holdout$bins) else
    data.table::rbindlist(lapply(names(x$bins_all$results), function(f) {
      e <- x$bins_all$results[[f]]
      data.table::data.table(feature = f, bin = e$bin, n_train = e$count, mean_train = e$mean, share_train = e$count / sum(e$count))
    }))
  bins_tab[, admitted := feature %in% x$survivors]
  dt_tab <- if (is.null(x$downturn)) x$pools[, list(pool, lra, moc_est, ccf_dt, downturn, ccf_final, ccf_floor, ccf_applied)] else
    cbind(x$downturn$table, method = x$downturn$method, add_on = x$downturn$add_on,
          periods = if (is.null(x$downturn$periods)) "" else paste(sprintf("%s to %s", format(x$downturn$periods$start), format(x$downturn$periods$end)), collapse = "; "))
  sheets <- list(
    "RDS_Funnel"      = x$funnel,
    "RDS_Summary"     = x$data_summary,
    "Driver_Bins"     = bins_tab,
    "Driver_Admission" = x$drivers,
    "Pools"           = x$pools,
    "Cells"           = x$cells,
    "Downturn_MoC"    = dt_tab,
    "Holdout"         = x$metrics,
    "Calibration"     = val$calibration,
    "Discrimination"  = val$discrimination,
    "Backtest_Cohort" = val$backtest,
    "Stability"       = val$stability,
    "Validation_Summary" = val$summary,
    "Model_Card"      = .kv_table(x$model_card),
    "Decision_Ledger" = x$ledger)
  files <- list(xlsx = .scr_write_xlsx(sheets, file.path(out_dir, sprintf("ead_%s.xlsx", .file_tag(tag)))),
                sql = file.path(out_dir, sprintf("sql_ead_%s.sql", .file_tag(tag))))
  writeLines(scr_sql(x), files$sql)
  for (f in files) msg("  %s", f)
  x$files <- files
  invisible(x)
}

# NSE column names used in data.table expressions of this file
utils::globalVariables(c(
  ".first",
  ".last",
  ".r",
  ".run",
  ".smax",
  "i_d",
  "a",
  "admitted",
  "ccf",
  "ccf_applied",
  "ccf_dt",
  "ccf_final",
  "ccf_floor",
  "ccf_raw",
  "cell",
  "def",
  "downturn",
  "drawn",
  "drawn_ref",
  "ead_realised",
  "end",
  "event_id",
  "fast_default",
  "fid",
  "flag_fixed",
  "floor_binding",
  "grp",
  "horizon_months",
  "limit",
  "limit_change",
  "limit_default",
  "limit_ref",
  "m",
  "measure",
  "moc_est",
  "oid",
  "ord",
  "p_vs_dev",
  "pred",
  "ref_date",
  "rule",
  "start",
  "undrawn",
  "undrawn_ref",
  "utilisation",
  "utilisation_ref"
))
