# ============================================================================ #
# lgd.R - Stage 10: loss given default under the IRB approach
# ============================================================================ #
# Five steps, one object carried through: scr_workout() turns default events
# and cash flows into the reference data set (realised LGD per default);
# scr_lgd() fits the two-stage model (cure x severity) and derives the pools;
# scr_lgd_downturn() and scr_lgd_floor() add the downturn and the input
# floor to the pools; scr_elbe() derives the in-default grid. scr_apply(),
# scr_sql(), scr_lgd_validate() and scr_export() close the production and
# validation contracts. Every judgement lands in the ledger.
# ============================================================================ #

# -- small helpers ---------------------------------------------------------- #

#' Whole months elapsed from one date to another (day-aware, vectorised)
#' @keywords internal
#' @noRd
.months_between <- function(from, to) {
  f <- as.POSIXlt(as.Date(from)); t <- as.POSIXlt(as.Date(to))
  m <- 12L * (t$year - f$year) + (t$mon - f$mon)
  as.integer(m - as.integer(t$mday < f$mday))
}

#' Reference rate in force at each date (last observation at or before it)
#' @keywords internal
#' @noRd
.lgd_rate_at <- function(rates, dates, flat) {
  if (is.null(rates)) {
    if (!is.numeric(flat) || length(flat) != 1L || !is.finite(flat)) {
      stop("scr_workout(): give a `rates` table (date, rate) or set `lgd_discount_rate` in the configuration.", call. = FALSE)
    }
    return(rep(as.double(flat), length(dates)))
  }
  r <- data.table::as.data.table(rates)
  if (!all(c("date", "rate") %in% names(r))) stop("scr_workout(): `rates` needs the columns `date` and `rate`.", call. = FALSE)
  r <- r[, list(date = as.Date(date), rate = as.double(rate))][order(date)]
  if (!nrow(r) || anyNA(r$date) || any(!is.finite(r$rate))) {
    stop("scr_workout(): `rates` needs at least one row and no missing date or rate.", call. = FALSE)
  }
  idx <- findInterval(as.numeric(as.Date(dates)), as.numeric(r$date))
  idx[idx == 0L] <- 1L
  r$rate[idx]
}

#' Root default of every event: defaults of one facility closer than the window are one event
#' @keywords internal
#' @noRd
.lgd_merge_map <- function(d, window) {
  d <- d[order(facility_id, default_date)]
  # only facilities with more than one default need the sequential scan; the
  # rest are their own root (keeps the grouped R loop off millions of groups)
  multi <- d$facility_id %in% d$facility_id[duplicated(d$facility_id)]
  single <- d[!multi, list(default_id, root_id = default_id)]
  if (!any(multi)) return(single)
  d <- d[multi]
  out <- d[, {
    root <- default_id[1]; root_close <- close_date[1]; roots <- character(.N); roots[1] <- root
    if (.N > 1L) for (i in 2:.N) {
      gap <- if (is.na(root_close)) 0L else .months_between(root_close, default_date[i])
      if (is.na(root_close) || gap < window) {
        roots[i] <- root; root_close <- close_date[i]
      } else {
        root <- default_id[i]; root_close <- close_date[i]; roots[i] <- root
      }
    }
    list(default_id = default_id, root_id = roots)
  }, by = facility_id]
  data.table::rbindlist(list(single, out[, list(default_id, root_id)]))
}

#' Cumulative discounted recovery rate by product and month, from closed defaults
#' @keywords internal
#' @noRd
.lgd_recovery_profile <- function(cf, rds, t_max) {
  closed <- rds[status == "closed", list(root_id, product, ead)]
  months <- 0:t_max
  build <- function(ids, label) {
    e <- sum(closed$ead[closed$root_id %in% ids])
    x <- cf[type == "recovery" & root_id %in% ids, list(pv = sum(pv)), by = list(m = pmin(month, t_max))]
    v <- numeric(length(months)); v[match(x$m, months)] <- x$pv
    data.table::data.table(product = label, month = months, n_closed = length(ids),
                           cum_recovery = if (e > 0) cumsum(v) / e else numeric(length(months)))
  }
  rows <- lapply(unique(closed$product), function(p) build(closed$root_id[closed$product == p], p))
  data.table::rbindlist(c(rows, list(build(closed$root_id, "all"))))
}

#' Allocate indirect costs over the events by the configured key
#' @keywords internal
#' @noRd
.lgd_allocate_costs <- function(rds, indirect, key) {
  if (is.null(indirect)) return(rep(0, nrow(rds)))
  keyv <- switch(key, ead = rds$ead, count = rep(1, nrow(rds)), duration = pmax(1, rds$months_in_default))
  out <- numeric(nrow(rds))
  if (is.numeric(indirect) && length(indirect) == 1L) {
    if (indirect > 0 && sum(keyv) > 0) out <- indirect * keyv / sum(keyv)
    return(out)
  }
  tb <- data.table::as.data.table(indirect)
  if (!all(c("product", "amount") %in% names(tb))) {
    stop("scr_workout(): `indirect_costs` must be a single total or a table with `product` and `amount`.", call. = FALSE)
  }
  for (i in seq_len(nrow(tb))) {
    w <- rds$product == tb$product[i]
    if (any(w) && sum(keyv[w]) > 0) out[w] <- tb$amount[i] * keyv[w] / sum(keyv[w])
  }
  out
}

#' Ledger row in the house format
#' @keywords internal
#' @noRd
.lgd_ledger_row <- function(action, detail, reason = "") {
  data.table::data.table(action = action, detail = detail, reason = reason, date = format(Sys.Date()))
}

# ============================================================================ #
# scr_workout
# ============================================================================ #

#' Workout LGD: the reference data set from default events and cash flows
#'
#' Builds the reference data set (RDS) of realised loss given default, one
#' row per default event, from a table of default events and the long table
#' of their post-default cash flows. Every cash flow is discounted to the
#' default date at the reference rate in force at that date plus
#' `lgd_discount_add_on`, with monthly compounding over whole months:
#' \deqn{\mathrm{PV} = \frac{A}{(1 + r/12)^{t}}}
#' where `t` is the number of whole months between the default date and the
#' cash-flow date. The realised LGD is the economic loss
#' \deqn{\mathrm{LGD} = \frac{E - \mathrm{PV}(R) + \mathrm{PV}(C) + \mathrm{PV}(D) + C^{\mathrm{ind}}}{E}}
#' with `E` the exposure at default, `R` recoveries, `C` direct costs, `D`
#' drawings after default and `C^ind` the indirect costs allocated by
#' `lgd_cost_allocation`.
#'
#' @section Rules:
#'
#' * **Cures.** An event with `status == "cured"` returns to performing: the
#'   balance outstanding at the cure date (`ead` plus the drawings after
#'   default, net of the cash recovered) enters as an artificial recovery
#'   on the cure date, so the cure carries its costs and the discount
#'   effect, never a zero loss by decree.
#' * **Multiple defaults.** Two defaults of one facility separated by fewer
#'   than `lgd_cure_window` months (from the close of the first to the start
#'   of the second), or a new default while the first is still open, are
#'   one event: the first default date and exposure are kept, the cash
#'   flows of both spells are pooled and the status of the last spell rules.
#' * **Incomplete workouts.** An open event younger than `lgd_t_max` months
#'   receives the expected further recovery read from the recovery profile
#'   of the closed defaults of the same product (cumulative discounted
#'   recovery rate by month in default); an open event at or beyond
#'   `lgd_t_max` is treated as closed with no further recovery.
#' * **Bounds.** With `lgd_floor_at_zero` the realised LGD used in the
#'   averages is floored at zero and with `lgd_cap_at_one` capped at one;
#'   `lgd_raw` always keeps the unbounded value.
#'
#' The long-run average is reported default-weighted (the arithmetic mean
#' over events) and exposure-weighted, overall, by product and by calendar
#' year of default.
#'
#' @param defaults A `data.frame` or `data.table` with one row per default
#'   event: `default_id`, `facility_id`, `default_date`, `ead`, `product`,
#'   `status` (`"closed"`, `"cured"` or `"open"`), optionally `close_date`,
#'   plus any driver columns, which are carried into the RDS.
#' @param cashflows Long table: `default_id`, `date`, `amount`, `type`
#'   (`"recovery"`, `"direct_cost"` or `"drawing"`).
#' @param rates Optional table `(date, rate)` of the annual reference rate;
#'   the rate in force at the default date is used. `NULL` uses the flat
#'   `lgd_discount_rate` of the configuration.
#' @param config A [scr_config()]; keys `lgd_*`.
#' @param indirect_costs Total indirect workout cost to allocate: a single
#'   number, or a table `(product, amount)` allocated within product.
#' @param obs_date Observation date; `NULL` uses the latest date seen.
#' @param keep_rows Keep the cash-flow table with its present values.
#'
#' @return An object of class `scr_workout`: `rds` (one row per default
#'   event: identifiers, `default_date`, `ead`, `product`, drivers,
#'   `status`, `months_in_default`, `discount_rate`, `pv_recovery`,
#'   `pv_cost`, `pv_drawing`, `cost_indirect`, `recovery_extrapolated`,
#'   `lgd_raw`, `lgd_real`, `is_cure`, `is_incomplete`, `merged_n`),
#'   `recovery_profile` (product x month: `cum_recovery`), `extrapolation`
#'   (per open event), `funnel` (rule, n, action), `summary` (`n`,
#'   `cure_rate`, `lra_default_weighted`, `lra_exposure_weighted`,
#'   `share_incomplete`, `discount_rate_mean`, `by_product`, `by_year`),
#'   `ledger`, `config`, `obs_date`, and `cashflows` with `keep_rows`.
#'   `rds` also carries `year`, `close_date`, `recovery_nominal`,
#'   `recovery_artificial`, `cost_nominal`, `drawing_nominal`,
#'   `closed_at_t_max`, `absorbed`, `n_cashflows` and `last_month`; `funnel`
#'   has `kept`; `summary` has `n_cure`, `lra_raw`, `ead_total` and `years`.
#'
#' @family irb-lgd
#' @examples
#' cfg <- scr_config(verbose = FALSE)
#' wo <- scr_workout(scr_demo_lgd, scr_demo_lgd_cashflows, rates = scr_demo_rates, config = cfg)
#' wo
#' wo$funnel
#' head(wo$rds[, c("default_id", "product", "status", "lgd_raw", "lgd_real", "is_cure")])
#' @export
scr_workout <- function(defaults, cashflows, rates = NULL, config = scr_config(), indirect_costs = 0,
                        obs_date = NULL, keep_rows = FALSE) {
  check_config(config, "scr_workout")
  cfg <- config
  old <- scr_verbose(isTRUE(cfg$verbose)); on.exit(scr_verbose(old), add = TRUE)
  t0 <- Sys.time()
  d <- data.table::as.data.table(defaults)
  cf <- data.table::as.data.table(cashflows)
  need_d <- c("default_id", "facility_id", "default_date", "ead", "product", "status")
  miss <- setdiff(need_d, names(d))
  if (length(miss)) stop("scr_workout(): `defaults` lacks column(s): ", lst(miss), call. = FALSE)
  need_cf <- c("default_id", "date", "amount", "type")
  miss <- setdiff(need_cf, names(cf))
  if (length(miss)) stop("scr_workout(): `cashflows` lacks column(s): ", lst(miss), call. = FALSE)
  if (!"close_date" %in% names(d)) d[, close_date := as.Date(NA)]
  d[, `:=`(default_id = as.character(default_id), facility_id = as.character(facility_id),
           default_date = as.Date(default_date), ead = as.double(ead), product = as.character(product),
           status = as.character(status), close_date = as.Date(close_date))]
  if (anyDuplicated(d$default_id)) stop("scr_workout(): duplicated `default_id`.", call. = FALSE)
  bad <- setdiff(unique(d$status), c("closed", "cured", "open"))
  if (length(bad)) stop("scr_workout(): `status` must be closed, cured or open (got ", lst(bad), ").", call. = FALSE)
  cf[, `:=`(default_id = as.character(default_id), date = as.Date(date), amount = as.double(amount), type = as.character(type))]
  bad <- setdiff(unique(cf$type), c("recovery", "direct_cost", "drawing"))
  if (length(bad)) stop("scr_workout(): cash-flow `type` must be recovery, direct_cost or drawing (got ", lst(bad), ").", call. = FALSE)
  if (anyNA(cf$amount) || anyNA(cf$date)) stop("scr_workout(): cash flows with a missing amount or date.", call. = FALSE)
  obs_date <- as.Date(obs_date %||% max(c(d$default_date, d$close_date, cf$date), na.rm = TRUE))
  d[!is.na(close_date) & close_date > obs_date, `:=`(close_date = as.Date(NA), status = "open")]

  msg_stage(10, "workout LGD (reference data set)")
  msg("  %s default events | %s cash-flow rows | observation date %s", n_fmt(nrow(d)), n_fmt(nrow(cf)), format(obs_date))
  funnel <- list()

  # -- scope --------------------------------------------------------------- #
  in_scope <- !is.na(d$ead) & d$ead > 0 & !is.na(d$default_date) & d$default_date <= obs_date
  funnel$scope <- data.table::data.table(rule = "NOT_IN_SCOPE", n = sum(!in_scope),
                                         action = "excluded: EAD <= 0, missing default date or default after the observation date")
  d <- d[in_scope]
  if (!nrow(d)) stop("scr_workout(): no default event in scope.", call. = FALSE)

  # -- multiple defaults of one facility ----------------------------------- #
  map <- .lgd_merge_map(d, cfg$lgd_cure_window)
  d <- merge(d, map, by = "default_id", sort = FALSE)
  last <- d[order(default_date), list(status_last = status[.N], close_last = close_date[.N], merged_n = .N,
                                      absorbed = paste(default_id[-1], collapse = ";")), by = root_id]
  rds <- merge(d[default_id == root_id], last, by = "root_id", sort = FALSE)
  rds[, `:=`(status = status_last, close_date = close_last)][, c("status_last", "close_last") := NULL]
  n_merged <- sum(d$default_id != d$root_id)
  funnel$merge <- data.table::data.table(rule = "MULTIPLE_DEFAULT_MERGED", n = n_merged,
                                         action = sprintf("merged into the earlier default of the facility (window %d months)", cfg$lgd_cure_window))
  cf[, root_id := map$root_id[match(default_id, map$default_id)]]
  orphan <- is.na(cf$root_id)
  funnel$orphan <- data.table::data.table(rule = "CASHFLOW_WITHOUT_DEFAULT", n = sum(orphan),
                                          action = "cash-flow rows dropped: no default event in scope")
  cf <- cf[!orphan]

  # -- discounting --------------------------------------------------------- #
  rds[, rate_base := .lgd_rate_at(rates, default_date, cfg$lgd_discount_rate)]
  rds[, discount_rate := rate_base + cfg$lgd_discount_add_on]
  cf <- merge(cf, rds[, list(root_id, default_date, discount_rate)], by = "root_id", sort = FALSE)
  cf[, month := pmax(0L, .months_between(default_date, date))]
  cf[, pv := amount / (1 + discount_rate / 12)^month]
  # grouped sums on indicator-weighted columns: plain sum()/max() by group are
  # GForce-optimised, per-group subsetting (amount[type == ...]) is not
  sums <- c("recovery_nominal", "pv_recovery_cash", "cost_nominal", "pv_cost", "drawing_nominal", "pv_drawing")
  w <- cf[, list(root_id, month,
                 recovery_nominal = amount * (type == "recovery"), pv_recovery_cash = pv * (type == "recovery"),
                 cost_nominal = amount * (type == "direct_cost"), pv_cost = pv * (type == "direct_cost"),
                 drawing_nominal = amount * (type == "drawing"), pv_drawing = pv * (type == "drawing"))]
  agg <- w[, c(lapply(.SD, sum), list(n_cashflows = .N)), by = root_id, .SDcols = sums]
  agg <- merge(agg, w[, list(last_month = max(month)), by = root_id], by = "root_id", sort = FALSE)
  rm(w)
  rds <- merge(rds, agg, by = "root_id", all.x = TRUE, sort = FALSE)
  for (nm in c("recovery_nominal", "pv_recovery_cash", "cost_nominal", "pv_cost", "drawing_nominal", "pv_drawing", "last_month")) {
    data.table::set(rds, which(is.na(rds[[nm]])), nm, 0)
  }
  data.table::set(rds, which(is.na(rds$n_cashflows)), "n_cashflows", 0L)
  rds[, months_in_default := pmax(0L, data.table::fifelse(status == "open" | is.na(close_date),
                                                           .months_between(default_date, obs_date),
                                                           .months_between(default_date, close_date)))]
  rds[, is_cure := status == "cured"]
  # outstanding at the cure date: exposure at default plus drawings after
  # default less the cash recovered (the drawings are in the loss, so they
  # must also be in the balance returned to performing)
  rds[, recovery_artificial := data.table::fifelse(is_cure, pmax(0, ead + drawing_nominal - recovery_nominal), 0)]
  rds[, pv_artificial := recovery_artificial / (1 + discount_rate / 12)^months_in_default]
  rds[, pv_recovery := pv_recovery_cash + pv_artificial]

  # -- incomplete workouts ------------------------------------------------- #
  beyond <- rds$status == "open" & rds$months_in_default >= cfg$lgd_t_max
  rds[, closed_at_t_max := beyond]
  rds[beyond, status := "closed"]
  rds[, is_incomplete := status == "open"]
  funnel$tmax <- data.table::data.table(rule = "OPEN_BEYOND_T_MAX_CLOSED", n = sum(beyond),
                                        action = sprintf("open for %d months or more: closed with no further recovery", cfg$lgd_t_max))
  profile <- .lgd_recovery_profile(cf, rds, cfg$lgd_t_max)
  rds[, recovery_extrapolated := 0]
  ext <- NULL
  if (any(rds$is_incomplete)) {
    inc <- rds[is_incomplete == TRUE, list(default_id, product, ead, months_in_default)]
    # vectorised look-up of the profile at (product, month); a product without
    # closed defaults reads the "all" profile
    pkey <- paste(profile$product, profile$month, sep = "\r")
    rho <- function(src, m) {
      v <- profile$cum_recovery[match(paste(src, pmin(m, cfg$lgd_t_max), sep = "\r"), pkey)]
      data.table::fifelse(is.na(v), 0, v)
    }
    src <- data.table::fifelse(inc$product %in% profile$product, inc$product, "all")
    a <- rho(src, inc$months_in_default); b <- rho(src, rep(cfg$lgd_t_max, nrow(inc)))
    ext <- data.table::data.table(default_id = inc$default_id, product = inc$product, months_in_default = inc$months_in_default,
                                  rho_tau = a, rho_t_max = b, expected_further = inc$ead * pmax(0, b - a),
                                  profile_source = src, lambda = 1)
    rds[match(ext$default_id, default_id), recovery_extrapolated := ext$expected_further]
  }
  funnel$inc <- data.table::data.table(rule = "INCOMPLETE_EXTRAPOLATED", n = sum(rds$is_incomplete),
                                       action = "open workouts younger than t_max: expected further recovery added from the product profile")

  # -- costs, loss, bounds ------------------------------------------------- #
  rds[, cost_indirect := .lgd_allocate_costs(rds, indirect_costs, cfg$lgd_cost_allocation)]
  rds[, lgd_raw := (ead - pv_recovery - recovery_extrapolated + pv_cost + pv_drawing + cost_indirect) / ead]
  rds[, lgd_real := lgd_raw]
  neg <- rds$lgd_raw < 0; above <- rds$lgd_raw > 1
  if (isTRUE(cfg$lgd_floor_at_zero)) rds[neg, lgd_real := 0]
  if (isTRUE(cfg$lgd_cap_at_one)) rds[above, lgd_real := 1]
  funnel$neg <- data.table::data.table(rule = "NEGATIVE_LGD_FLOORED", n = sum(neg),
                                       action = if (isTRUE(cfg$lgd_floor_at_zero)) "floored at 0 in lgd_real (raw value kept)" else "kept negative (lgd_floor_at_zero = FALSE)")
  funnel$above <- data.table::data.table(rule = if (isTRUE(cfg$lgd_cap_at_one)) "LGD_ABOVE_ONE_CAPPED" else "LGD_ABOVE_ONE", n = sum(above),
                                         action = if (isTRUE(cfg$lgd_cap_at_one)) "capped at 1 in lgd_real (raw value kept)" else "kept above 1 (lgd_cap_at_one = FALSE)")
  funnel <- data.table::rbindlist(funnel)
  funnel[, kept := nrow(rds)]

  # -- tidy the RDS -------------------------------------------------------- #
  rds[, c("root_id", "rate_base", "pv_recovery_cash", "pv_artificial") := NULL]
  rds[, year := as.integer(format(default_date, "%Y"))]
  front <- c("default_id", "facility_id", "default_date", "year", "ead", "product", "status", "close_date",
             "months_in_default", "discount_rate", "recovery_nominal", "recovery_artificial", "pv_recovery",
             "cost_nominal", "pv_cost", "drawing_nominal", "pv_drawing", "cost_indirect", "recovery_extrapolated",
             "lgd_raw", "lgd_real", "is_cure", "is_incomplete", "closed_at_t_max", "merged_n", "absorbed",
             "n_cashflows", "last_month")
  data.table::setcolorder(rds, c(front, setdiff(names(rds), front)))
  data.table::setorder(rds, default_date, default_id)

  by_product <- rds[, list(n = .N, cure_rate = mean(is_cure), lra = mean(lgd_real),
                           lra_ew = sum(lgd_real * ead) / sum(ead), share_incomplete = mean(is_incomplete)), by = product][order(product)]
  by_year <- rds[, list(n = .N, cure_rate = mean(is_cure), lra = mean(lgd_real),
                        lra_ew = sum(lgd_real * ead) / sum(ead)), by = year][order(year)]
  summary <- list(n = nrow(rds), n_cure = sum(rds$is_cure), cure_rate = mean(rds$is_cure),
                  lra_default_weighted = mean(rds$lgd_real), lra_exposure_weighted = sum(rds$lgd_real * rds$ead) / sum(rds$ead),
                  lra_raw = mean(rds$lgd_raw), share_incomplete = mean(rds$is_incomplete),
                  discount_rate_mean = mean(rds$discount_rate), ead_total = sum(rds$ead),
                  years = length(unique(rds$year)), by_product = by_product, by_year = by_year)
  ledger <- data.table::rbindlist(list(
    .lgd_ledger_row("discounting", sprintf("%s at default + add-on %s, monthly compounding over whole months",
                                           if (is.null(rates)) sprintf("flat rate %s", fmt_pct(cfg$lgd_discount_rate, 2)) else "reference rate",
                                           fmt_pct(cfg$lgd_discount_add_on, 2))),
    .lgd_ledger_row("cure_treatment", "outstanding at the cure date (ead plus drawings, net of cash recovered) as an artificial recovery on the cure date"),
    .lgd_ledger_row("multiple_defaults", sprintf("%d event(s) merged: gap below %d months or overlapping spells", n_merged, cfg$lgd_cure_window)),
    .lgd_ledger_row("incomplete_workouts", sprintf("%d open event(s) extrapolated from the product recovery profile (lambda = 1); %d closed at t_max = %d",
                                                   sum(rds$is_incomplete), sum(beyond), cfg$lgd_t_max)),
    .lgd_ledger_row("cost_allocation", sprintf("indirect costs %s allocated by %s",
                                               if (is.numeric(indirect_costs) && length(indirect_costs) == 1L) format(indirect_costs) else "by product",
                                               cfg$lgd_cost_allocation)),
    .lgd_ledger_row("bounds", sprintf("floor at zero: %s | cap at one: %s", cfg$lgd_floor_at_zero, cfg$lgd_cap_at_one))))
  msg("  RDS: %s events | cure rate %s | LRA %s (default-weighted) %s (exposure-weighted) | incomplete %s (%.2fs)",
      n_fmt(summary$n), fmt_pct(summary$cure_rate), fmt_pct(summary$lra_default_weighted), fmt_pct(summary$lra_exposure_weighted),
      fmt_pct(summary$share_incomplete), as.numeric(difftime(Sys.time(), t0, units = "secs")))
  structure(list(rds = rds[], recovery_profile = profile, extrapolation = ext, funnel = funnel[], summary = summary,
                 ledger = ledger, config = cfg, obs_date = obs_date,
                 cashflows = if (isTRUE(keep_rows)) cf[, list(default_id = root_id, original_id = default_id, date, type, amount, month, discount_rate, pv)] else NULL),
            class = c("scr_workout", "list"))
}

#' @export
print.scr_workout <- function(x, ...) {
  s <- x$summary
  cat(sprintf("<scr_workout> %s default events | %d calendar years | observation date %s\n", n_fmt(s$n), s$years, format(x$obs_date)))
  cat(sprintf("  cure rate %s | incomplete %s | mean discount rate %s\n", fmt_pct(s$cure_rate), fmt_pct(s$share_incomplete),
              fmt_pct(s$discount_rate_mean, 2)))
  cat(sprintf("  long-run average LGD: %s default-weighted | %s exposure-weighted | raw mean %s\n",
              fmt_pct(s$lra_default_weighted), fmt_pct(s$lra_exposure_weighted), fmt_pct(s$lra_raw)))
  bp <- s$by_product
  for (i in seq_len(nrow(bp))) cat(sprintf("  %-12s n %-5d cure %-7s LRA %-7s (ew %s)\n", bp$product[i], bp$n[i],
                                           fmt_pct(bp$cure_rate[i]), fmt_pct(bp$lra[i]), fmt_pct(bp$lra_ew[i])))
  f <- x$funnel[n > 0]
  for (i in seq_len(nrow(f))) cat(sprintf("  funnel %-28s %d\n", f$rule[i], f$n[i]))
  invisible(x)
}

# ============================================================================ #
# model pieces
# ============================================================================ #

#' Drivers ready for both binners: categorical NA becomes "NA", numeric NA is refused
#' @keywords internal
#' @noRd
.lgd_prepare <- function(dt, drivers, fn = "scr_lgd") {
  dt <- data.table::copy(dt)
  for (f in drivers) {
    v <- dt[[f]]
    if (is.numeric(v)) {
      if (anyNA(v)) stop(fn, "(): driver '", f, "' has missing values; impute them before the LGD stages.", call. = FALSE)
      data.table::set(dt, j = f, value = as.double(v))
    } else {
      v <- as.character(v); v[is.na(v)] <- "NA"
      data.table::set(dt, j = f, value = v)
    }
  }
  dt
}

#' glm (or beta regression) on bin statistics with the removal of wrong-signed drivers
#' @keywords internal
#' @noRd
.lgd_fit_signcheck <- function(w, y, feats, engine = c("binomial", "fractional_logit", "beta"), max_abs_coef = 15) {
  engine <- match.arg(engine)
  if (identical(engine, "beta") && !requireNamespace("betareg", quietly = TRUE)) {
    stop("lgd_severity = \"beta\" needs the 'betareg' package; install it or use \"fractional_logit\".", call. = FALSE)
  }
  log <- list(); fs <- feats; phi <- NA_real_
  clamp <- function(p) pmin(pmax(p, 1e-4), 1 - 1e-4)
  repeat {
    if (!length(fs)) { cf <- c("(Intercept)" = stats::qlogis(clamp(mean(y)))); break }
    X <- as.data.frame(w[, paste0(fs, "_woe"), with = FALSE]); names(X) <- fs
    if (identical(engine, "beta")) {
      n <- length(y); X$.y <- (y * (n - 1) + 0.5) / n
      m <- betareg::betareg(.y ~ ., data = X)
      cf <- stats::coef(m, model = "mean"); phi <- unname(stats::coef(m, model = "precision"))
    } else {
      X$.y <- y
      fam <- if (identical(engine, "binomial")) stats::binomial() else stats::quasibinomial()
      m <- suppressWarnings(stats::glm(.y ~ ., family = fam, data = X, model = FALSE, x = FALSE, y = FALSE))
      cf <- stats::coef(m)
    }
    cf[is.na(cf)] <- 0
    b <- cf[fs]
    bad <- fs[b <= 0 | abs(b) > max_abs_coef]
    if (!length(bad)) {
      for (f in fs) log[[f]] <- data.table::data.table(variable = f, coef = unname(b[f]), action = "kept", reason = "OK")
      break
    }
    worst <- fs[which.min(ifelse(fs %in% bad, b, Inf))]
    if (!worst %in% bad) worst <- bad[1]
    log[[worst]] <- data.table::data.table(variable = worst, coef = unname(b[worst]), action = "removed",
                                           reason = if (b[worst] <= 0) "SIGN_REVERSED" else "COEF_TOO_LARGE")
    fs <- setdiff(fs, worst)
  }
  names(cf)[1] <- "(Intercept)"
  sc <- if (length(log)) data.table::rbindlist(log)[order(action, variable)] else
    data.table::data.table(variable = character(), coef = numeric(), action = character(), reason = character())
  list(coef = cf[c("(Intercept)", fs)], features = fs, sign_check = sc, engine = engine, phi = phi)
}

#' Linear predictor from coefficients and `_woe` columns
#' @keywords internal
#' @noRd
.lgd_link <- function(coef, w, feats, n) {
  eta <- rep(unname(coef["(Intercept)"]), n)
  for (f in feats) eta <- eta + unname(coef[f]) * w[[paste0(f, "_woe")]]
  eta
}

#' Somers' D of the prediction with respect to the realised value
#'
#' `(C - D) / (pairs untied on the realised value)`, the pairs counted
#' exactly in `O(n log n)` by the compiled kernel (`.scr_somers()`); `NA`
#' when the prediction is constant. Generalised AUC is `(D + 1) / 2`.
#' @keywords internal
#' @noRd
.lgd_somers <- function(p, r) .scr_somers(p, r, const_p = NA_real_)

#' Loss capture ratio: area above the diagonal of the model curve over the ideal one
#' @keywords internal
#' @noRd
.lgd_lcr <- function(p, r, e) {
  loss <- r * e; tot <- sum(loss); te <- sum(e)
  if (!is.finite(tot) || tot <= 0 || te <= 0) return(NA_real_)
  area <- function(o) {
    x <- c(0, cumsum(e[o]) / te); yv <- c(0, cumsum(loss[o]) / tot)
    sum(diff(x) * (yv[-1] + yv[-length(yv)]) / 2) - 0.5
  }
  a_i <- area(order(-r))
  if (a_i <= 0) return(NA_real_)
  area(order(-p)) / a_i
}

#' Accuracy and discrimination of a continuous prediction, with a bootstrap CI
#' @keywords internal
#' @noRd
.lgd_metrics <- function(pred, real, ead, n_boot = 200L, level = 0.95, seed = NULL, nthread = 1L) {
  ok <- is.finite(pred) & is.finite(real) & is.finite(ead)
  pred <- as.double(pred[ok]); real <- as.double(real[ok]); ead <- as.double(ead[ok]); n <- length(pred)
  na <- list(somers_lo = NA_real_, somers_hi = NA_real_, gauc_lo = NA_real_, gauc_hi = NA_real_,
             lcr_lo = NA_real_, lcr_hi = NA_real_)
  if (n < 3L) return(c(list(n = n, rmse = NA_real_, mae = NA_real_, r2 = NA_real_, spearman = NA_real_, somers_d = NA_real_,
                            gauc = NA_real_, lcr = NA_real_), na, list(n_boot = 0L, level = level)))
  sst <- sum((real - mean(real))^2)
  out <- list(n = n, rmse = sqrt(mean((pred - real)^2)), mae = mean(abs(pred - real)),
              r2 = if (sst > 0) 1 - sum((real - pred)^2) / sst else NA_real_,
              spearman = suppressWarnings(stats::cor(pred, real, method = "spearman")),
              somers_d = .lgd_somers(pred, real))
  out$gauc <- (out$somers_d + 1) / 2
  out$lcr <- .lgd_lcr(pred, real, ead)
  out <- c(out, na, list(n_boot = 0L, level = level))
  if (n_boot >= 2L) {
    .scr_local_seed(seed)
    seeds <- sample.int(.Machine$integer.max, n_boot)
    .scr_rng_guard()   # the replicates re-seed; the user's stream moves by the n_boot draws above only
    # the closure carries the data and nothing else from this frame; its parent
    # is the namespace, which a PSOCK worker loads to reach the compiled kernel
    env <- list2env(list(pred = pred, real = real, ead = ead, n = n, somers = .lgd_somers, lcr = .lgd_lcr),
                    parent = environment(.lgd_metrics))
    fun <- function(sd) {
      set.seed(sd)
      j <- sample.int(n, n, replace = TRUE)
      c(somers(pred[j], real[j]), lcr(pred[j], real[j], ead[j]))
    }
    environment(fun) <- env
    b <- do.call(rbind, .scr_lapply(seeds, fun, nthread = nthread))
    a <- (1 - level) / 2
    qs <- stats::quantile(b[, 1], c(a, 1 - a), na.rm = TRUE, names = FALSE)
    ql <- stats::quantile(b[, 2], c(a, 1 - a), na.rm = TRUE, names = FALSE)
    out$somers_lo <- qs[1]; out$somers_hi <- qs[2]
    out$gauc_lo <- (qs[1] + 1) / 2; out$gauc_hi <- (qs[2] + 1) / 2
    out$lcr_lo <- ql[1]; out$lcr_hi <- ql[2]
    out$n_boot <- as.integer(n_boot)
  }
  out
}

#' Cohort split on the default date: the last share of the dates is the hold-out
#' @keywords internal
#' @noRd
.lgd_split <- function(dates, holdout, seed) {
  n <- length(dates)
  cutoff <- as.Date(stats::quantile(as.numeric(dates), probs = 1 - holdout, type = 1, names = FALSE), origin = "1970-01-01")
  ho <- dates > cutoff
  method <- "cohort"
  if (!any(ho) || all(ho)) {
    .scr_local_seed(seed); ho <- seq_len(n) %in% sample.int(n, max(1L, round(holdout * n))); method <- "random"; cutoff <- as.Date(NA)
  }
  list(holdout = ho, cutoff = cutoff, method = method, n_train = sum(!ho), n_holdout = sum(ho))
}

#' Cure stage: the binary engine on is_cure, WOE, hold-out, glm with sign check
#' @keywords internal
#' @noRd
.lgd_cure_stage <- function(dt_tr, dt_ho, y_tr, y_ho, drivers, cfg) {
  none <- function(note) {
    p <- mean(y_tr)
    list(active = FALSE, fit = NULL, features = character(), binned = character(),
         coef = c("(Intercept)" = stats::qlogis(pmin(pmax(p, 1e-4), 1 - 1e-4))),
         sign_check = data.table::data.table(variable = character(), coef = numeric(), action = character(), reason = character()),
         bins = NULL, holdout = NULL, note = note)
  }
  if (sum(y_tr) < 20L || sum(1 - y_tr) < 20L) return(none("fewer than 20 cures or 20 non-cures on train: constant cure rate"))
  dtr <- data.table::copy(dt_tr)[, is_cure := as.integer(y_tr)]
  fit <- strip_failed_features(fit_binning(dtr, "is_cure", drivers, cfg))
  binned <- names(fit$results)
  if (!length(binned)) return(none("no driver could be binned: constant cure rate"))
  app_tr <- apply_woe(fit, dt_tr, binned, "both"); app_ho <- apply_woe(fit, dt_ho, binned, "both")
  ho <- holdout_check(app_tr, app_ho, y_tr, y_ho, binned, cfg)
  s <- data.table::as.data.table(fit$summary)
  iv <- stats::setNames(s$total_iv, s$feature)
  cand <- binned[iv[binned] >= cfg$iv_min & ho$holdout_ok[match(binned, ho$feature)] %in% TRUE]
  fr <- .lgd_fit_signcheck(app_tr, y_tr, cand, engine = "binomial", max_abs_coef = cfg$max_abs_coef)
  bins <- data.table::rbindlist(lapply(binned, function(f) {
    r <- fit$results[[f]]
    data.table::data.table(variable = f, bin_id = as.integer(r$id), bin = r$bin, woe = r$woe, iv = r$iv,
                           count = r$count, count_pos = r$count_pos, cure_rate = r$count_pos / pmax(1L, r$count),
                           selected = f %in% fr$features)
  }))
  bins <- merge(bins, ho[, list(variable = feature, iv_holdout, psi, psi_flag, holdout_ok, holdout_reason)], by = "variable", sort = FALSE)
  idx <- function(app) data.table::as.data.table(lapply(stats::setNames(fr$features, paste0(fr$features, "_cure")), function(f)
    match(app[[paste0(f, "_bin")]], fit$results[[f]]$bin)))
  list(active = TRUE, fit = fit, features = fr$features, binned = binned, coef = fr$coef, sign_check = fr$sign_check,
       bins = bins[], holdout = ho, note = if (length(fr$features)) "OK" else "no driver survived screening, hold-out and the sign check: constant cure rate",
       woe_tr = app_tr, woe_ho = app_ho, idx_tr = idx(app_tr), idx_ho = idx(app_ho))
}

#' Severity stage: the continuous binner on lgd of non-cures, fractional logit or beta
#' @keywords internal
#' @noRd
.lgd_severity_stage <- function(dt_tr, dt_ho, y_tr, y_ho, drivers, cfg) {
  dtr <- data.table::copy(dt_tr)[, lgd_sev := y_tr]; dho <- data.table::copy(dt_ho)[, lgd_sev := y_ho]
  fit <- .scr_bin_continuous(dtr, "lgd_sev", drivers, min_bins = 2L, max_bins = cfg$max_bins, min_share = cfg$min_bin_pct,
                             min_n = cfg$lgd_min_defaults_bin, monotone = "auto", scale = "mean", nthread = cfg$nthread)
  ho <- .cbins_holdout(fit, dtr, dho, "lgd_sev", alpha = cfg$psi_alpha)
  s <- merge(data.table::as.data.table(fit$summary), ho$summary, by = "feature", sort = FALSE)
  cand <- s[n_bins >= 2L & holdout_ok %in% TRUE, feature]
  w_tr <- .cbins_apply_value(fit, dt_tr, drivers, unbinned = 0); w_ho <- .cbins_apply_value(fit, dt_ho, drivers, unbinned = 0)
  fr <- .lgd_fit_signcheck(w_tr, y_tr, cand, engine = cfg$lgd_severity, max_abs_coef = cfg$max_abs_coef)
  bins <- data.table::rbindlist(lapply(drivers, function(f) {
    r <- fit$results[[f]]
    data.table::data.table(variable = f, bin_id = as.integer(r$id), bin = r$bin, mean = r$mean, eta2_share = r$iv,
                           count = r$count, selected = f %in% fr$features)
  }))
  bins <- merge(bins, ho$bins[, list(variable = feature, bin, n_holdout, mean_holdout)], by = c("variable", "bin"), sort = FALSE)
  bins <- merge(bins, s[, list(variable = feature, eta2, eta2_holdout, psi, psi_flag, holdout_ok, holdout_reason)], by = "variable", sort = FALSE)
  idx <- function(dt) { i <- .cbins_apply_idx(fit, dt, fr$features); data.table::setnames(i, paste0(names(i), "_sev")); i }
  list(fit = fit, features = fr$features, coef = fr$coef, sign_check = fr$sign_check, engine = fr$engine, phi = fr$phi,
       bins = bins[order(variable, bin_id)], summary = s, woe_tr = w_tr, woe_ho = w_ho, idx_tr = idx(dt_tr), idx_ho = idx(dt_ho),
       note = if (length(fr$features)) "OK" else "no driver survived binning, hold-out and the sign check: constant severity")
}

#' Predictions of the two stages for a prepared driver table
#' @keywords internal
#' @noRd
.lgd_predict <- function(x, dt) {
  n <- nrow(dt)
  cu <- x$cure
  p_cure <- if (isTRUE(cu$active) && length(cu$features)) {
    stats::plogis(.lgd_link(cu$coef, apply_woe(cu$fit, dt, cu$features, "woe"), cu$features, n))
  } else rep(stats::plogis(unname(cu$coef["(Intercept)"])), n)
  if (!x$has_cures) p_cure <- rep(0, n)
  sv <- x$severity
  sev <- if (length(sv$features)) {
    stats::plogis(.lgd_link(sv$coef, .cbins_apply_value(sv$fit, dt, sv$features, unbinned = 0), sv$features, n))
  } else rep(stats::plogis(unname(sv$coef["(Intercept)"])), n)
  data.table::data.table(p_cure = p_cure, severity = sev, lgd_pred = p_cure * x$lgd_cure + (1 - p_cure) * sev)
}

#' Pool index of predictions from the pool table
#' @keywords internal
#' @noRd
.lgd_pool_of <- function(pred, pools) {
  k <- nrow(pools)
  if (k <= 1L) return(rep(1L, length(pred)))
  findInterval(pred, pools$pred_hi[-k], left.open = TRUE) + 1L
}

#' Reference value: mean realised LGD of the two worst calendar years, per pool
#' @keywords internal
#' @noRd
.lgd_reference_value <- function(scored, pools) {
  yr <- scored[, list(lgd = mean(lgd_real)), by = list(pool, year = as.integer(format(default_date, "%Y")))]
  vapply(pools$pool, function(p) {
    v <- sort(yr$lgd[yr$pool == p], decreasing = TRUE)
    if (!length(v)) NA_real_ else mean(utils::head(v, 2L))
  }, numeric(1))
}

#' Downturn table for a set of pools (provisional when no periods are given)
#' @keywords internal
#' @noRd
.lgd_downturn_table <- function(pools, scored, method, add_on, periods = NULL, min_n = 10L) {
  # training rows only, like the long-run averages: the hold-out stays independent
  if ("sample" %in% names(scored)) scored <- scored[sample == "train"]
  tb <- pools[, list(pool, n, lra, moc_c)]
  tb[, reference_value := .lgd_reference_value(scored, pools)]
  tb[, dt_type3 := lra + add_on]
  tb[, `:=`(dt_observed = NA_real_, n_downturn = 0L)]
  if (identical(method, "type1") && !is.null(periods)) {
    in_dt <- .lgd_in_periods(scored$default_date, periods)
    obs <- scored[in_dt, list(n_dt = .N, dt_obs = mean(lgd_real)), by = pool]
    m <- match(tb$pool, obs$pool)
    tb[!is.na(m), `:=`(dt_observed = obs$dt_obs[m[!is.na(m)]], n_downturn = obs$n_dt[m[!is.na(m)]])]
  }
  tb[, method_used := method]
  if (identical(method, "type1")) {
    tb[is.na(dt_observed) | n_downturn < min_n, method_used := if (is.null(periods)) "type3_provisional" else "type3_fallback"]
  }
  fallback <- if (identical(method, "none")) tb$lra else tb$dt_type3
  tb[, dt := data.table::fifelse(method_used == "type1", dt_observed, fallback)]
  tb[, lgd_dt := pmin(1, pmax(lra + moc_c, dt + moc_c))]
  tb[, impact := lgd_dt - pmin(1, lra + moc_c)]
  tb[, below_reference := !is.na(reference_value) & lgd_dt < reference_value]
  tb[]
}

#' Which dates fall inside any of the periods
#' @keywords internal
#' @noRd
.lgd_in_periods <- function(dates, periods) {
  p <- data.table::as.data.table(periods)
  if (!all(c("start", "end") %in% names(p))) stop("`periods` needs the columns `start` and `end`.", call. = FALSE)
  p[, `:=`(start = as.Date(start), end = as.Date(end))]
  if (any(is.na(p$start) | is.na(p$end) | p$end < p$start)) stop("`periods`: every row needs start <= end.", call. = FALSE)
  out <- rep(FALSE, length(dates))
  for (i in seq_len(nrow(p))) out <- out | (dates >= p$start[i] & dates <= p$end[i])
  out
}

# ============================================================================ #
# scr_lgd
# ============================================================================ #

#' Two-stage LGD model and pools on the reference data set
#'
#' Fits the standard two-stage structure on the RDS of [scr_workout()]:
#' \deqn{\mathrm{LGD} = P(\mathrm{cure}\mid x)\,\mathrm{LGD}^{\mathrm{cure}} + \big(1 - P(\mathrm{cure}\mid x)\big)\,\mathrm{E}[\mathrm{LGD}\mid \mathrm{no\ cure}, x]}
#' The **cure stage** is a binary model on `is_cure` with the scorecard
#' machinery: optimal binning of the drivers on the training cohorts, WOE,
#' hold-out revalidation with frozen bins and a logistic regression on the
#' WOE columns with the sign check (every coefficient positive). The
#' **severity stage** bins the same drivers against the realised LGD of
#' the non-cures with [scr_bin_continuous()] (bin means, monotone, at least
#' `lgd_min_defaults_bin` defaults per bin, hold-out revalidated) and fits a
#' fractional logit (`glm` with a quasi-binomial family on the bin means)
#' or, with `lgd_severity = "beta"`, a beta regression through the
#' `betareg` package. `LGD^cure` is the mean realised LGD of the cures on
#' train (costs and the discount effect, never zero by decree).
#'
#' The split is by cohort of default: the last `holdout` share of the
#' default dates is the hold-out. Metrics on both samples: RMSE, MAE,
#' R-squared, Spearman rho, Somers' D of the prediction with respect to
#' the realised LGD (generalised AUC `(D + 1) / 2`) with a bootstrap
#' confidence interval, and the loss capture ratio. Pools come from
#' [scr_lgd_pools()]. The object carries a provisional downturn (type 3
#' add-on, or none, by configuration) and no floor until
#' [scr_lgd_downturn()] and [scr_lgd_floor()] run.
#'
#' @param x An [scr_workout()] object.
#' @param drivers Column names of the RDS to use as drivers.
#' @param config A [scr_config()]; keys `lgd_*`, the binning and hold-out
#'   keys of stage 2, `max_abs_coef`, `n_boot`, `ci_level`, `seed`, `nthread`.
#' @param holdout Share of the cohorts held out (by default date).
#' @param date_col Column of the RDS with the default date.
#'
#' @return An object of class `scr_lgd`: `split`, `drivers`, `cure` (fit,
#'   features, coef, sign_check, bins, holdout), `severity` (fit, features,
#'   coef, engine, sign_check, bins), `lgd_cure`, `has_cures`, `scored` (one
#'   row per default: `sample`, `p_cure`, `severity`, `lgd_pred`, `pool`,
#'   `lgd_real`), `bins_idx`, `samples` (predicted vs realised by decile of
#'   the prediction), `metrics`, `pools`, `downturn`, `floors`, `workout`
#'   (the profile and summary of the RDS), `model_card`, `ledger`, `config`.
#'
#' @family irb-lgd
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, n_boot = 20)
#' wo <- scr_workout(scr_demo_lgd, scr_demo_lgd_cashflows, rates = scr_demo_rates, config = cfg)
#' m <- scr_lgd(wo, drivers = c("product", "ltv", "prior_dpd_max", "months_on_book", "region"),
#'              config = cfg)
#' m
#' m$pools[, c("pool", "n", "lra", "lra_ew", "moc_c", "lgd_dt")]
#' m$metrics
#' @export
scr_lgd <- function(x, drivers, config = scr_config(), holdout = 0.3, date_col = "default_date") {
  if (!inherits(x, "scr_workout")) stop("scr_lgd(): `x` must come from scr_workout().", call. = FALSE)
  check_config(config, "scr_lgd")
  cfg <- config
  old <- scr_verbose(isTRUE(cfg$verbose)); on.exit(scr_verbose(old), add = TRUE)
  .scr_num1(holdout, "holdout", lower = 0, upper = 1, open_lower = TRUE)
  rds <- x$rds
  miss <- setdiff(c(drivers, date_col), names(rds))
  if (length(miss)) stop("scr_lgd(): column(s) not found in the RDS: ", lst(miss), call. = FALSE)
  if (!length(drivers)) stop("scr_lgd(): give at least one driver.", call. = FALSE)
  t0 <- Sys.time()
  msg_stage(10, "LGD model: cure x severity, pools")
  dt <- .lgd_prepare(rds, drivers)
  dates <- as.Date(rds[[date_col]])
  sp <- .lgd_split(dates, holdout, cfg$seed)
  tr <- which(!sp$holdout); ho <- which(sp$holdout)
  msg("  split by %s: train %s | hold-out %s%s", sp$method, n_fmt(sp$n_train), n_fmt(sp$n_holdout),
      if (!is.na(sp$cutoff)) sprintf(" (defaults after %s)", format(sp$cutoff)) else "")
  dt_tr <- dt[tr, drivers, with = FALSE]; dt_ho <- dt[ho, drivers, with = FALSE]
  y_cure <- as.integer(rds$is_cure)
  has_cures <- any(y_cure[tr] == 1L)

  # -- cure stage ---------------------------------------------------------- #
  cure <- .lgd_cure_stage(dt_tr, dt_ho, y_cure[tr], y_cure[ho], drivers, cfg)
  msg("  cure stage: %s", if (length(cure$features)) sprintf("%d driver(s): %s", length(cure$features), lst(cure$features)) else cure$note)
  lgd_cure <- if (has_cures) mean(rds$lgd_real[tr][y_cure[tr] == 1L]) else 0

  # -- severity stage ------------------------------------------------------ #
  y_sev <- pmin(1, pmax(0, rds$lgd_real))
  nc_tr <- tr[y_cure[tr] == 0L]; nc_ho <- ho[y_cure[ho] == 0L]
  if (length(nc_tr) < 10L) stop("scr_lgd(): fewer than 10 non-cured defaults on train.", call. = FALSE)
  sev <- .lgd_severity_stage(dt[nc_tr, drivers, with = FALSE], dt[nc_ho, drivers, with = FALSE], y_sev[nc_tr], y_sev[nc_ho], drivers, cfg)
  msg("  severity stage (%s): %s", sev$engine, if (length(sev$features)) sprintf("%d driver(s): %s", length(sev$features), lst(sev$features)) else sev$note)

  obj <- list(cure = cure, severity = sev, lgd_cure = lgd_cure, has_cures = has_cures)
  pr <- .lgd_predict(obj, dt[, drivers, with = FALSE])
  scored <- data.table::data.table(default_id = rds$default_id, sample = data.table::fifelse(sp$holdout, "holdout", "train"),
                                   default_date = dates, product = rds$product, ead = rds$ead,
                                   months_in_default = rds$months_in_default, is_cure = rds$is_cure,
                                   is_incomplete = rds$is_incomplete, lgd_real = rds$lgd_real, pr)
  bins_idx <- data.table::data.table(default_id = rds$default_id, sample = scored$sample)
  if (length(cure$features)) {
    ci <- data.table::rbindlist(list(cure$idx_tr, cure$idx_ho)); ci[, .i := c(tr, ho)]; data.table::setorder(ci, .i); ci[, .i := NULL]
    bins_idx <- cbind(bins_idx, ci)
  }
  if (length(sev$features)) {
    si <- .cbins_apply_idx(sev$fit, dt, sev$features); data.table::setnames(si, paste0(names(si), "_sev"))
    bins_idx <- cbind(bins_idx, si)
  }
  for (nm in c("woe_tr", "woe_ho", "idx_tr", "idx_ho")) { cure[[nm]] <- NULL; sev[[nm]] <- NULL }

  # -- metrics ------------------------------------------------------------- #
  metrics <- data.table::rbindlist(lapply(c("train", "holdout"), function(nm) {
    s <- scored[sample == nm]
    m <- .lgd_metrics(s$lgd_pred, s$lgd_real, s$ead, n_boot = cfg$n_boot, level = cfg$ci_level, seed = cfg$seed, nthread = cfg$nthread)
    data.table::data.table(sample = nm, as.data.frame(m))
  }))
  for (i in seq_len(nrow(metrics))) msg("  %-8s RMSE %.4f  R2 %.3f  Spearman %.3f  gAUC %.3f [%.3f, %.3f]  LCR %.3f", metrics$sample[i],
                                        metrics$rmse[i], metrics$r2[i], metrics$spearman[i], metrics$gauc[i], metrics$gauc_lo[i], metrics$gauc_hi[i], metrics$lcr[i])
  samples <- .lgd_deciles(scored)

  # -- pools, provisional downturn ----------------------------------------- #
  ledger <- data.table::rbindlist(list(
    x$ledger,
    .lgd_ledger_row("split", sprintf("%s split, hold-out %s%s", sp$method, fmt_pct(holdout), if (!is.na(sp$cutoff)) sprintf(", cut-off %s", format(sp$cutoff)) else "")),
    .lgd_ledger_row("cure_stage", cure$note),
    .lgd_ledger_row("severity_stage", sprintf("%s: %s", sev$engine, sev$note)),
    if (nrow(cure$sign_check[action == "removed"])) .lgd_ledger_row("cure_sign_check", lst(cure$sign_check[action == "removed", variable]), "SIGN_REVERSED_OR_TOO_LARGE"),
    if (nrow(sev$sign_check[action == "removed"])) .lgd_ledger_row("severity_sign_check", lst(sev$sign_check[action == "removed", variable]), "SIGN_REVERSED_OR_TOO_LARGE")
  ), use.names = TRUE)
  out <- structure(c(obj, list(
    split = sp[c("method", "cutoff", "n_train", "n_holdout")], drivers = drivers, date_col = date_col,
    scored = scored, bins_idx = bins_idx, samples = samples, metrics = metrics, pools = NULL, downturn = NULL, floors = NULL,
    workout = list(recovery_profile = x$recovery_profile, summary = x$summary, funnel = x$funnel, extrapolation = x$extrapolation, obs_date = x$obs_date),
    ledger = ledger, config = cfg, model_card = NULL)), class = c("scr_lgd", "list"))
  pools <- scr_lgd_pools(out)
  out$ledger <- data.table::rbindlist(list(out$ledger, attr(pools, "ledger")), use.names = TRUE)
  data.table::setattr(pools, "ledger", NULL)
  out$pools <- data.table::copy(pools)
  out$scored[, pool := .lgd_pool_of(lgd_pred, out$pools)]
  dtb <- .lgd_downturn_table(out$pools, out$scored, cfg$lgd_downturn, cfg$lgd_downturn_add_on, periods = NULL)
  out$downturn <- list(table = dtb, periods = NULL, method = cfg$lgd_downturn, add_on = cfg$lgd_downturn_add_on, status = "provisional",
                       reason = NA_character_)
  out$pools[, `:=`(lgd_dt = dtb$lgd_dt, floor = 0, lgd_final = dtb$lgd_dt)]
  out$ledger <- data.table::rbindlist(list(out$ledger, .lgd_ledger_row("downturn",
    sprintf("provisional: %s%s", if (identical(cfg$lgd_downturn, "none")) "no downturn (LRA + MoC)" else sprintf("type 3 add-on %s", fmt_pct(cfg$lgd_downturn_add_on)),
            if (identical(cfg$lgd_downturn, "type1")) "; run scr_lgd_downturn() with the downturn periods" else ""), "DOWNTURN_PENDING")), use.names = TRUE)
  out$model_card <- .lgd_model_card(out)
  msg("  pools: %d | LRA from %s to %s | provisional downturn %s (%.2fs)", nrow(out$pools), fmt_pct(min(out$pools$lra)), fmt_pct(max(out$pools$lra)),
      cfg$lgd_downturn, as.numeric(difftime(Sys.time(), t0, units = "secs")))
  out
}

#' Predicted vs realised by decile of the prediction, per sample
#' @keywords internal
#' @noRd
.lgd_deciles <- function(scored) {
  br <- unique(stats::quantile(scored[sample == "train", lgd_pred], probs = seq(0.1, 0.9, 0.1), names = FALSE))
  scored[, list(n = .N, ead = sum(ead), pred_mean = mean(lgd_pred), real_mean = mean(lgd_real),
                real_ew = sum(lgd_real * ead) / sum(ead), cure_rate = mean(is_cure)),
         by = list(sample, decile = findInterval(lgd_pred, br, left.open = TRUE) + 1L)][order(sample, decile)]
}

#' @keywords internal
#' @noRd
.lgd_model_card <- function(m) {
  ws <- m$workout$summary; cfg <- m$config
  list(
    package = sprintf("scorecraft %s", as.character(utils::packageVersion("scorecraft"))),
    fitted_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    n_defaults = ws$n, n_train = m$split$n_train, n_holdout = m$split$n_holdout, split_method = m$split$method,
    split_cutoff = if (is.na(m$split$cutoff)) "" else format(m$split$cutoff), observation_date = format(m$workout$obs_date),
    cure_rate = ws$cure_rate, share_incomplete = ws$share_incomplete, lra_default_weighted = ws$lra_default_weighted,
    lra_exposure_weighted = ws$lra_exposure_weighted, lgd_cure = m$lgd_cure,
    drivers = paste(m$drivers, collapse = ", "),
    cure_features = paste(m$cure$features, collapse = ", "), cure_engine = "logistic on WOE",
    severity_features = paste(m$severity$features, collapse = ", "), severity_engine = m$severity$engine,
    n_pools = nrow(m$pools), downturn_method = m$downturn$method, downturn_status = m$downturn$status,
    downturn_add_on = m$downturn$add_on, floor_asset_class = m$floors$asset_class %||% "none",
    framework = cfg$framework, discount_add_on = cfg$lgd_discount_add_on, t_max = cfg$lgd_t_max, cure_window = cfg$lgd_cure_window,
    cost_allocation = cfg$lgd_cost_allocation, floor_at_zero = cfg$lgd_floor_at_zero, cap_at_one = cfg$lgd_cap_at_one,
    min_defaults_bin = cfg$lgd_min_defaults_bin,
    gauc_train = m$metrics[sample == "train", gauc], gauc_holdout = m$metrics[sample == "holdout", gauc],
    rmse_holdout = m$metrics[sample == "holdout", rmse], lcr_holdout = m$metrics[sample == "holdout", lcr])
}

#' @export
print.scr_lgd <- function(x, ...) {
  cat(sprintf("<scr_lgd> %s defaults | train %s / hold-out %s (%s split%s) | cure rate %s\n",
              n_fmt(nrow(x$scored)), n_fmt(x$split$n_train), n_fmt(x$split$n_holdout), x$split$method,
              if (!is.na(x$split$cutoff)) sprintf(" after %s", format(x$split$cutoff)) else "", fmt_pct(x$workout$summary$cure_rate)))
  cat(sprintf("  cure stage: %s | severity stage (%s): %s | LGD of a cure %s\n",
              if (length(x$cure$features)) lst(x$cure$features) else "constant", x$severity$engine,
              if (length(x$severity$features)) lst(x$severity$features) else "constant", fmt_pct(x$lgd_cure)))
  for (i in seq_len(nrow(x$metrics))) {
    m <- x$metrics[i]
    cat(sprintf("  %-8s n %-5d RMSE %.4f  R2 %.3f  Spearman %.3f  gAUC %.3f [%.3f, %.3f]  LCR %.3f\n", m$sample, m$n, m$rmse, m$r2,
                m$spearman, m$gauc, m$gauc_lo, m$gauc_hi, m$lcr))
  }
  cat(sprintf("  pools %d | downturn %s (%s) | floor %s\n", nrow(x$pools), x$downturn$method, x$downturn$status,
              if (is.null(x$floors)) "not applied" else sprintf("%s, binding in %s of the defaults", x$floors$asset_class, fmt_pct(x$floors$binding_share))))
  p <- x$pools
  cat("  pool   n     pred        LRA     LRA ew   MoC C    LGD DT   floor    final\n")
  for (i in seq_len(nrow(p))) cat(sprintf("  %-4d %-5d %6.3f  %8s  %8s  %6.3f  %8s  %6.3f  %8s\n", p$pool[i], p$n[i], p$pred_mean[i],
                                          fmt_pct(p$lra[i]), fmt_pct(p$lra_ew[i]), p$moc_c[i], fmt_pct(p$lgd_dt[i]), p$floor[i], fmt_pct(p$lgd_final[i])))
  invisible(x)
}

# ============================================================================ #
# pools, downturn, floors, ELBE
# ============================================================================ #

#' LGD pools from the predicted LGD
#'
#' Cuts the training predictions into `n_pools` quantile bands, merges the
#' bands with fewer than `min_defaults` defaults into the neighbour with
#' the closer long-run average, then merges adjacent bands whose long-run
#' averages break the increasing order (pool-adjacent violators), so that
#' the pools are ordered both in predicted and in realised LGD. Per pool:
#' the default-weighted long-run average (the regulatory estimate), the
#' exposure-weighted one, the standard error, the category-C margin of
#' conservatism (one-sided 95% t interval on the mean) and their sum.
#'
#' @param x An [scr_lgd()] object.
#' @param n_pools Target number of pools; `NULL` uses `lgd_n_pools`.
#' @param min_defaults Minimum defaults per pool; `NULL` uses
#'   `lgd_min_defaults_bin`.
#'
#' @return A `data.table` with one row per pool: `pool`, `pred_lo`,
#'   `pred_hi`, `pred_mean`, `n`, `share`, `ead`, `lra`, `lra_ew`, `sd`,
#'   `se`, `moc_c`, `lra_moc`, `merged_from`.
#'
#' @family irb-lgd
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, n_boot = 20)
#' wo <- scr_workout(scr_demo_lgd, scr_demo_lgd_cashflows, rates = scr_demo_rates, config = cfg)
#' m <- scr_lgd(wo, drivers = c("product", "ltv", "prior_dpd_max"), config = cfg)
#' scr_lgd_pools(m, n_pools = 4)
#' @export
scr_lgd_pools <- function(x, n_pools = NULL, min_defaults = NULL) {
  if (!inherits(x, "scr_lgd")) stop("scr_lgd_pools(): `x` must come from scr_lgd().", call. = FALSE)
  k <- as.integer(n_pools %||% x$config$lgd_n_pools); min_n <- as.integer(min_defaults %||% x$config$lgd_min_defaults_bin)
  if (length(k) != 1L || is.na(k) || k < 1L) stop("scr_lgd_pools(): `n_pools` must be a positive integer.", call. = FALSE)
  if (length(min_n) != 1L || is.na(min_n) || min_n < 1L) stop("scr_lgd_pools(): `min_defaults` must be a positive integer.", call. = FALSE)
  s <- x$scored[sample == "train"]
  pred <- s$lgd_pred; real <- s$lgd_real; ead <- s$ead
  probs <- seq(0, 1, length.out = k + 1L)[-c(1L, k + 1L)]
  br <- unique(stats::quantile(pred, probs = probs, type = 7, names = FALSE))
  # snap every break to the midpoint between the adjacent distinct predictions:
  # a boundary sitting on an observed value would let floating-point noise
  # (R versus SQL) flip the pool of that observation
  u <- sort(unique(pred))
  br <- vapply(br, function(q) { hi <- u[u > q]; if (!length(hi)) q else (max(u[u <= q]) + hi[1]) / 2 }, numeric(1))
  br <- unique(br); br <- br[br > min(pred) & br < max(pred)]
  pre <- findInterval(pred, br, left.open = TRUE) + 1L
  K <- length(br) + 1L
  n_b <- tabulate(pre, nbins = K); s_b <- vapply(seq_len(K), function(b) sum(real[pre == b]), numeric(1))
  ledger <- list()
  g <- .cbin_merge(n_b, s_b, max_bins = K, min_n = min_n, min_share = 0)$group
  n_small <- K - max(g)
  if (n_small > 0) ledger$small <- .lgd_ledger_row("pool_merge", sprintf("%d band(s) below %d defaults merged into the neighbour with the closer LRA", n_small, min_n), "MIN_DEFAULTS")
  n_g <- vapply(seq_len(max(g)), function(i) sum(n_b[g == i]), numeric(1)); s_g <- vapply(seq_len(max(g)), function(i) sum(s_b[g == i]), numeric(1))
  g2 <- .cbin_pava(n_g, s_g, increasing = TRUE)
  if (max(g2) < max(g)) ledger$pava <- .lgd_ledger_row("pool_merge", sprintf("%d adjacent pool(s) merged: long-run average not increasing in the prediction", max(g) - max(g2)), "NOT_MONOTONE")
  g <- g2[g]
  pool <- g[pre]
  edges <- c(-Inf, br, Inf)
  tb <- data.table::rbindlist(lapply(seq_len(max(g)), function(i) {
    b <- which(g == i); w <- pool == i
    n <- sum(w); sdv <- if (n > 1L) stats::sd(real[w]) else NA_real_
    se <- if (n > 1L) sdv / sqrt(n) else NA_real_
    moc <- if (n > 1L && is.finite(se)) stats::qt(0.95, n - 1L) * se else 0
    data.table::data.table(pool = i, pred_lo = edges[min(b)], pred_hi = edges[max(b) + 1L], pred_mean = mean(pred[w]),
                           n = n, share = n / length(pred), ead = sum(ead[w]), lra = mean(real[w]),
                           lra_ew = sum(real[w] * ead[w]) / sum(ead[w]), sd = sdv, se = se, moc_c = moc, lra_moc = mean(real[w]) + moc,
                           merged_from = paste(b, collapse = ";"))
  }))
  data.table::setattr(tb, "ledger", if (length(ledger)) data.table::rbindlist(ledger) else NULL)
  tb
}

#' Downturn LGD per pool
#'
#' Quantifies the downturn per pool from user-supplied downturn periods.
#' `method = "type1"` (observed impact): the default-weighted realised LGD
#' of the training defaults whose default date falls inside the periods; a
#' pool with fewer than ten such defaults falls back to type 3. `method = "type3"`:
#' the long-run average plus `add_on`. `method = "none"`: the long-run
#' average. The reference value (a challenger, not a bound) is the mean of
#' the two worst calendar years of the pool. Both use the training rows
#' only, so the hold-out stays independent evidence. The downturn LGD used for
#' capital is
#' \deqn{\mathrm{LGD}^{DT} = \min\!\big(1,\ \max(\mathrm{LRA} + \mathrm{MoC},\ \mathrm{DT} + \mathrm{MoC})\big)}
#' and the impact `LGD^DT - min(1, LRA + MoC)` is reported per pool.
#'
#' @param x An [scr_lgd()] object.
#' @param periods A table with `start` and `end` dates of the downturn
#'   periods. Required for `"type1"`.
#' @param method `"type1"`, `"type3"` or `"none"`; `NULL` uses `lgd_downturn`.
#' @param add_on Type-3 add-on; `NULL` uses `lgd_downturn_add_on`.
#' @param reason Free text recorded in the ledger, mandatory: the choice of
#'   periods and method is an analyst decision.
#'
#' @return The `scr_lgd` object with `downturn` (`table` per pool:
#'   `lra`, `moc_c`, `dt_observed`, `n_downturn`, `dt_type3`,
#'   `reference_value`, `method_used`, `dt`, `lgd_dt`, `impact`,
#'   `below_reference`; `periods`, `method`, `add_on`, `status`, `reason`)
#'   and the pool columns `lgd_dt` and `lgd_final` updated.
#'
#' @family irb-lgd
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, n_boot = 20)
#' wo <- scr_workout(scr_demo_lgd, scr_demo_lgd_cashflows, rates = scr_demo_rates, config = cfg)
#' m <- scr_lgd(wo, drivers = c("product", "ltv", "prior_dpd_max"), config = cfg)
#' m <- scr_lgd_downturn(m, periods = data.frame(start = as.Date("2022-01-01"),
#'                                                end = as.Date("2023-12-31")),
#'                       reason = "reference rate above 13% in 2022-2023")
#' m$downturn$table
#' @export
scr_lgd_downturn <- function(x, periods = NULL, method = NULL, add_on = NULL, reason = NULL) {
  if (!inherits(x, "scr_lgd")) stop("scr_lgd_downturn(): `x` must come from scr_lgd().", call. = FALSE)
  cfg <- x$config
    overridden <- (!is.null(method) && !identical(method, cfg$lgd_downturn)) ||
    (!is.null(add_on) && !isTRUE(all.equal(add_on, cfg$lgd_downturn_add_on)))
  method <- match.arg(method %||% cfg$lgd_downturn, c("type1", "type3", "none"))
  add_on <- add_on %||% cfg$lgd_downturn_add_on
  .scr_num1(add_on, "add_on", lower = 0, upper = 1)
  if (is.null(reason) || !is.character(reason) || length(reason) != 1L || !nzchar(trimws(reason))) {
    stop("scr_lgd_downturn(): a `reason` is mandatory (the downturn periods and method go in the ledger).", call. = FALSE)
  }
  if (identical(method, "type1") && is.null(periods)) stop("scr_lgd_downturn(): `periods` (start, end) are required for method \"type1\".", call. = FALSE)
  p <- if (is.null(periods)) NULL else data.table::as.data.table(periods)[, list(start = as.Date(start), end = as.Date(end))]
  base <- x$pools[, list(pool, pred_lo, pred_hi, pred_mean, n, share, ead, lra, lra_ew, sd, se, moc_c, lra_moc, merged_from)]
  tb <- .lgd_downturn_table(base, x$scored, method, add_on, periods = p)
  x$downturn <- list(table = tb, periods = p, method = method, add_on = add_on, status = "final", reason = reason %||% NA_character_)
  # copy before `:=`: the pool table is shared by reference with the caller's object
  x$pools <- data.table::copy(x$pools)
  x$pools[, lgd_dt := tb$lgd_dt]
  x$pools[, lgd_final := pmax(lgd_dt, floor)]
  x$ledger <- data.table::rbindlist(list(x$ledger, .lgd_ledger_row("downturn",
    sprintf("%s%s; add-on %s; periods %s; mean impact %s", method, if (overridden) " (overrides the configuration)" else "", fmt_pct(add_on),
            if (is.null(p)) "none" else paste(sprintf("%s to %s", format(p$start), format(p$end)), collapse = ", "),
            fmt_pct(mean(tb$impact))), reason %||% "")), use.names = TRUE)
  x$model_card <- .lgd_model_card(x)
  x
}

#' Input floor on the downturn LGD per pool
#'
#' Applies the LGD input floor of the framework's parameter table, blended
#' between the unsecured and the collateralised floor with the secured
#' share of the exposure:
#' \deqn{\mathrm{floor} = \mathrm{floor}_U\,(1 - s) + \mathrm{floor}_S\,s,\qquad
#'       \mathrm{LGD}^{\mathrm{final}} = \max(\mathrm{LGD}^{DT}, \mathrm{floor})}
#' A missing unsecured floor (residential mortgages, whose floor applies to
#' the whole exposure) uses the collateral floor throughout; an asset class
#' with no floor at all yields a floor of zero.
#'
#' @param x An [scr_lgd()] object.
#' @param params An [scr_irb_params()] object; `NULL` uses the configured
#'   framework. Edits are detected and recorded in the ledger.
#' @param asset_class Row of `params$lgd_floor`; `NULL` uses
#'   `asset_class` of the configuration.
#' @param secured_share Secured share of the exposure in `[0, 1]`: one value
#'   or one per pool. `NULL` means unsecured.
#' @param collateral Column of `params$lgd_floor` for the secured part:
#'   `"real_estate"`, `"financial"`, `"receivables"` or `"other_physical"`.
#'
#' @return The `scr_lgd` object with `floors` (`table` per pool with
#'   `lgd_dt`, `floor_unsecured`, `floor_secured`, `secured_share`,
#'   `floor`, `lgd_final`, `binding`; `asset_class`, `collateral`,
#'   `framework`, `params_modified`, `binding_share`) and the pool columns
#'   `floor` and `lgd_final` updated.
#'
#' @family irb-lgd
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, n_boot = 20)
#' wo <- scr_workout(scr_demo_lgd, scr_demo_lgd_cashflows, rates = scr_demo_rates, config = cfg)
#' m <- scr_lgd(wo, drivers = c("product", "ltv", "prior_dpd_max"), config = cfg)
#' m <- scr_lgd_floor(m, asset_class = "retail_other", secured_share = 0.4)
#' m$floors$table
#' @export
scr_lgd_floor <- function(x, params = NULL, asset_class = NULL, secured_share = NULL,
                          collateral = c("real_estate", "financial", "receivables", "other_physical")) {
  if (!inherits(x, "scr_lgd")) stop("scr_lgd_floor(): `x` must come from scr_lgd().", call. = FALSE)
  collateral <- match.arg(collateral)
  params <- .check_params(params %||% scr_irb_params(x$config$framework), "scr_lgd_floor")
  asset_class <- asset_class %||% x$config$asset_class
  lf <- params$lgd_floor
  i <- match(asset_class, lf$asset_class)
  if (is.na(i)) stop("scr_lgd_floor(): unknown `asset_class` '", asset_class, "'. See scr_irb_params()$lgd_floor.", call. = FALSE)
  fu <- lf$unsecured[i]; fs <- lf[[collateral]][i]
  k <- nrow(x$pools)
  s <- secured_share %||% 0
  if (!length(s) %in% c(1L, k) || any(!is.finite(s)) || any(s < 0 | s > 1)) stop("scr_lgd_floor(): `secured_share` must be in [0, 1], one value or one per pool.", call. = FALSE)
  s <- rep_len(as.double(s), k)
  fl <- if (is.na(fu) && is.na(fs)) rep(0, k) else if (is.na(fu)) rep(fs, k) else if (is.na(fs)) rep(fu, k) else fu * (1 - s) + fs * s
  tb <- data.table::data.table(pool = x$pools$pool, n = x$pools$n, lgd_dt = x$pools$lgd_dt, floor_unsecured = fu, floor_secured = fs,
                               secured_share = s, floor = fl)
  tb[, lgd_final := pmax(lgd_dt, floor)]
  tb[, binding := lgd_final > lgd_dt]
  # copy before `:=`: the pool table is shared by reference with the caller's object
  x$pools <- data.table::copy(x$pools)
  x$pools[, `:=`(floor = fl, lgd_final = pmax(lgd_dt, fl))]
  x$floors <- list(table = tb[], asset_class = asset_class, collateral = collateral, framework = params$framework,
                   params_modified = isTRUE(params$modified), binding_share = sum(tb$n[tb$binding]) / sum(tb$n))
  x$ledger <- data.table::rbindlist(list(x$ledger, .lgd_ledger_row("floor",
    sprintf("%s %s: unsecured %s, %s %s, secured share %s; binding in %s of the defaults%s", params$framework, asset_class,
            if (is.na(fu)) "n/a" else fmt_pct(fu), collateral, if (is.na(fs)) "n/a" else fmt_pct(fs),
            paste(fmt_pct(unique(s)), collapse = "/"), fmt_pct(x$floors$binding_share),
            if (isTRUE(params$modified)) "; params_modified = TRUE" else ""))), use.names = TRUE)
  x$model_card <- .lgd_model_card(x)
  x
}

#' ELBE and in-default LGD on a grid of months since default
#'
#' For every pool and every reference age `tau` of the grid, the expected
#' loss best estimate is the mean realised LGD of the training defaults of
#' the pool that were still in workout at `tau` (so that at `tau = 0` it
#' equals the pool's long-run average), and the in-default LGD adds the
#' unexpected-loss increment
#' \deqn{\Delta^{UL}(\tau) = \max(0,\ \mathrm{LGD}^{DT} - \mathrm{LRA})\;\frac{\rho(T_{\max}) - \rho(\tau - 1)}{\rho(T_{\max})}}
#' read from the recovery profile of the pool's product mix, where
#' \eqn{\rho(\tau - 1)} is the cumulative discounted recovery rate of the
#' months before age `tau` (zero at `tau = 0`): the downturn uplift
#' shrinks as the recoveries come in. The consistency table checks
#' that `lgd_in_default` at `tau = 0` reproduces the pool's `lgd_dt`.
#'
#' @param x An [scr_lgd()] object.
#' @param grid Months since default; `NULL` uses `lgd_elbe_grid`.
#'
#' @return An object of class `scr_elbe`: `table` (`months_since_default`,
#'   `pool`, `n_open`, `share_open`, `recovered_share`, `elbe`, `delta_ul`,
#'   `lgd_in_default`), `consistency` (per pool at `tau = 0`), `grid`, `t_max`.
#'
#' @family irb-lgd
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, n_boot = 20)
#' wo <- scr_workout(scr_demo_lgd, scr_demo_lgd_cashflows, rates = scr_demo_rates, config = cfg)
#' m <- scr_lgd(wo, drivers = c("product", "ltv", "prior_dpd_max"), config = cfg)
#' e <- scr_elbe(m)
#' e
#' @export
scr_elbe <- function(x, grid = NULL) {
  if (!inherits(x, "scr_lgd")) stop("scr_elbe(): `x` must come from scr_lgd().", call. = FALSE)
  grid <- sort(unique(as.integer(grid %||% x$config$lgd_elbe_grid)))
  if (any(grid < 0)) stop("scr_elbe(): the grid must be non-negative months.", call. = FALSE)
  t_max <- x$config$lgd_t_max
  prof <- x$workout$recovery_profile
  s <- x$scored[sample == "train"]
  rho_mix <- function(p, tau) {
    w <- s[pool == p, list(ead = sum(ead)), by = product]
    if (!nrow(w) || tau < 0) return(0)
    v <- vapply(w$product, function(pr) {
      src <- if (pr %in% prof$product) pr else "all"
      r <- prof[product == src & month == min(tau, t_max), cum_recovery]; if (length(r)) r else 0
    }, numeric(1))
    sum(v * w$ead) / sum(w$ead)
  }
  rows <- lapply(x$pools$pool, function(p) {
    lra <- x$pools$lra[x$pools$pool == p]; dt <- x$pools$lgd_dt[x$pools$pool == p]
    r_max <- rho_mix(p, t_max)
    data.table::rbindlist(lapply(grid, function(tau) {
      open <- s[pool == p & (months_in_default > tau | (is_incomplete & months_in_default >= tau) | tau == 0L)]
      n_p <- sum(s$pool == p)
      # recoveries collected before age tau: profile months 0 .. tau - 1 (none at tau = 0)
      rec <- if (r_max > 0) rho_mix(p, tau - 1L) / r_max else 0
      elbe <- if (nrow(open)) mean(open$lgd_real) else NA_real_
      dul <- max(0, dt - lra) * (1 - rec)
      data.table::data.table(months_since_default = tau, pool = p, n_open = nrow(open), share_open = nrow(open) / max(1L, n_p),
                             recovered_share = rec, elbe = elbe, delta_ul = dul, lgd_in_default = elbe + dul)
    }))
  })
  tb <- data.table::rbindlist(rows)[order(months_since_default, pool)]
  z <- tb[months_since_default == min(months_since_default)]
  cons <- data.table::data.table(pool = x$pools$pool, lra = x$pools$lra, lgd_dt = x$pools$lgd_dt,
                                 elbe_0 = z$elbe[match(x$pools$pool, z$pool)], lgd_in_default_0 = z$lgd_in_default[match(x$pools$pool, z$pool)])
  cons[, ok := abs(elbe_0 - lra) < 1e-9 & abs(lgd_in_default_0 - lgd_dt) < 1e-9]
  structure(list(table = tb, consistency = cons[], grid = grid, t_max = t_max), class = c("scr_elbe", "list"))
}

#' @export
print.scr_elbe <- function(x, ...) {
  cat(sprintf("<scr_elbe> %d pools x %d reference ages (months since default: %s) | t_max %d\n",
              length(unique(x$table$pool)), length(x$grid), paste(x$grid, collapse = ", "), x$t_max))
  cat(sprintf("  consistency at tau = 0: %s\n", if (all(x$consistency$ok)) "ELBE equals the LRA and the in-default LGD equals the downturn LGD" else "BROKEN"))
  w <- data.table::dcast(x$table, pool ~ months_since_default, value.var = "lgd_in_default")
  cat("  in-default LGD by pool and age\n")
  cat(sprintf("  %-5s %s\n", "pool", paste(sprintf("%7s", paste0("m", x$grid)), collapse = " ")))
  for (i in seq_len(nrow(w))) cat(sprintf("  %-5d %s\n", w$pool[i], paste(sprintf("%7s", fmt_pct(as.numeric(w[i, -1]))), collapse = " ")))
  invisible(x)
}

# ============================================================================ #
# production: apply, SQL
# ============================================================================ #

#' @rdname scr_apply
#' @export
scr_apply.scr_lgd <- function(x, newdata, what = c("pool", "lgd", "all"), ...) {
  what <- match.arg(what)
  dt <- data.table::as.data.table(newdata)
  miss <- setdiff(x$drivers, names(dt))
  if (length(miss)) stop("newdata lacks the driver column(s): ", lst(miss), call. = FALSE)
  dt <- .lgd_prepare(dt[, x$drivers, with = FALSE], x$drivers, "scr_apply")
  pr <- .lgd_predict(x, dt)
  pool <- .lgd_pool_of(pr$lgd_pred, x$pools)
  out <- data.table::data.table(pool = pool, lgd_lra = x$pools$lra[pool], lgd_dt = x$pools$lgd_dt[pool], lgd_final = x$pools$lgd_final[pool])
  if (identical(what, "pool")) return(out[])
  out <- cbind(pr, out)
  if (identical(what, "lgd")) out <- out[, list(lgd_pred, pool, lgd_lra, lgd_dt, lgd_final)]
  out[]
}

#' @rdname scr_sql
#' @export
scr_sql.scr_lgd <- function(x, table = NULL, dialect = NULL, file = NULL, ...) {
  .scr_local_verbose(x)
  cfg <- x$config
  if (!is.null(table)) cfg$sql_table <- table
  if (!is.null(dialect)) cfg$sql_dialect <- dialect
  .sql_out(.lgd_build_sql(x, cfg), file)
}

#' LGD SQL: drivers CTE, bin statistics of both stages, logits, pool, floor
#' @keywords internal
#' @noRd
.lgd_build_sql <- function(x, cfg) {
  dialect <- match.arg(cfg$sql_dialect, c("ansi", "postgres", "mysql", "mariadb", "sqlserver", "oracle", "spark", "hive",
                                          "databricks", "bigquery", "snowflake", "redshift", "duckdb", "sqlite"))
  keep <- cfg$sql_keep_columns
  cf <- x$cure$features; sf <- x$severity$features
  drivers <- union(cf, sf)
  ind <- function(v) paste0("    ", v)
  sel_lines <- function(cols) { l <- ind(cols); paste0(l, c(rep(",", length(l) - 1L), "")) }
  exprs <- character()
  if (length(cf)) exprs <- c(exprs, .sql_select_exprs(OptimalBinningWoE::obwoe_sql(
    obj = x$cure$fit, table = "base_scr", features = cf, output = "woe", style = "select", dialect = dialect,
    digits = NULL, comment = FALSE, suffix_woe = "_cwoe", bin_separator = cfg$bin_separator)))
  if (length(sf)) {
    if (length(exprs)) exprs[length(exprs)] <- paste0(exprs[length(exprs)], ",")
    exprs <- c(exprs, .sql_select_exprs(OptimalBinningWoE::obwoe_sql(
      obj = x$severity$fit, table = "base_scr", features = sf, output = "woe", style = "select", dialect = dialect,
      digits = NULL, comment = FALSE, suffix_woe = "_swoe", bin_separator = cfg$bin_separator)))
  }
  logit <- function(coef, feats, suffix) {
    terms <- c(.sql_num(unname(coef["(Intercept)"])), vapply(feats, function(f) sprintf("%s * %s%s", .sql_num(unname(coef[f])), f, suffix), character(1)))
    sprintf("1 / (1 + EXP(-(%s)))", paste(terms, collapse = " + "))
  }
  p_cure <- if (x$has_cures) logit(x$cure$coef, cf, "_cwoe") else "0"
  sev <- logit(x$severity$coef, sf, "_swoe")
  k <- nrow(x$pools)
  pool_case <- if (k <= 1L) "1" else sprintf("CASE %s ELSE %d END", paste(sprintf("WHEN lgd_pred <= %s THEN %d", .sql_num(x$pools$pred_hi[-k]), seq_len(k - 1L)), collapse = " "), k)
  by_pool <- function(v) sprintf("CASE pool %s ELSE %s END", paste(sprintf("WHEN %d THEN %s", x$pools$pool, .sql_num(v)), collapse = " "), .sql_num(v[k]))
  greatest <- if (dialect %in% c("sqlite", "sqlserver")) "CASE WHEN lgd_dt >= lgd_floor THEN lgd_dt ELSE lgd_floor END" else "GREATEST(lgd_dt, lgd_floor)"
  base_cols <- c(keep, drivers)
  c("-- =============================================================",
    sprintf("-- scorecraft | LGD model | %d pool(s) | dialect: %s", k, dialect),
    sprintf("-- Generated on %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf("-- cure stage: %s | severity stage (%s): %s | LGD of a cure %s",
            if (length(cf)) paste(cf, collapse = ", ") else "constant", x$severity$engine,
            if (length(sf)) paste(sf, collapse = ", ") else "constant", .sql_num(x$lgd_cure)),
    sprintf("-- downturn %s (%s) | floor %s", x$downturn$method, x$downturn$status, x$floors$asset_class %||% "none"),
    "-- Block 1 (CTE base_scr): the driver columns.",
    "-- Block 2 (CTE woe_scr): bin statistics of both stages, emitted by OptimalBinningWoE::obwoe_sql().",
    "-- Block 3: cure probability, severity, predicted LGD, pool, pool parameters and the floored LGD.",
    "-- =============================================================", "",
    "WITH base_scr AS (", "  SELECT", if (length(base_cols)) sel_lines(base_cols) else "    *", sprintf("  FROM %s", cfg$sql_table), "),",
    "woe_scr AS (", "  SELECT", if (length(keep)) sprintf("    %s,", keep), if (length(exprs)) exprs else "    1 AS one_scr", "  FROM base_scr", "),",
    "link_scr AS (", "  SELECT", "    w.*,", sprintf("    %s AS p_cure,", p_cure), sprintf("    %s AS severity", sev), "  FROM woe_scr w", "),",
    "pred_scr AS (", "  SELECT", "    l.*,", sprintf("    p_cure * %s + (1 - p_cure) * severity AS lgd_pred", .sql_num(x$lgd_cure)), "  FROM link_scr l", "),",
    "pool_scr AS (", "  SELECT", "    p.*,", sprintf("    %s AS pool", pool_case), "  FROM pred_scr p", "),",
    "param_scr AS (", "  SELECT", "    q.*,", sprintf("    %s AS lgd_lra,", by_pool(x$pools$lra)), sprintf("    %s AS lgd_dt,", by_pool(x$pools$lgd_dt)),
    sprintf("    %s AS lgd_floor", by_pool(x$pools$floor)), "  FROM pool_scr q", ")", "",
    "SELECT", if (length(keep)) sprintf("    %s,", keep),
    "    p_cure,", "    severity,", "    lgd_pred,", "    pool,", "    lgd_lra,", "    lgd_dt,", "    lgd_floor,",
    sprintf("    %s AS lgd_final", greatest), "FROM param_scr;") |> .sql_lines()
}

# ============================================================================ #
# validation
# ============================================================================ #

#' Validation battery of an LGD model
#'
#' Runs the three blocks of the usual LGD validation on the hold-out sample
#' (or on `newdata`) against the training reference:
#'
#' * **Calibration.** Per pool and for the portfolio, the one-sided t-test
#'   of realised against estimated LGD (the pool long-run average), where
#'   under-estimation is the failure: `p = 1 - Phi(t)`; the loss shortfall
#'   `1 - sum(LGD_real E) / sum(LGD_pred E)`; the coverage of the realised
#'   mean by the downturn LGD; the regression of realised on predicted.
#' * **Discrimination.** Somers' D / generalised AUC of the prediction with
#'   its bootstrap interval, compared with the training value through
#'   `S = (gAUC_init - gAUC_curr) / sigma_curr`; Spearman rho; the loss
#'   capture ratio; R-squared.
#' * **Stability.** PSI of the pool distribution and of the bins of every
#'   driver of both stages, with the fixed and the n-adjusted threshold.
#' * **Pools.** Homogeneity within a pool (Welch test between the halves of
#'   the pool split at its median prediction; a small p-value means a pool
#'   that still discriminates) and heterogeneity between adjacent pools
#'   (Welch test; a large p-value means pools that do not differ).
#'
#' Traffic lights use the p-value thresholds of `config$pd_lights` (shared
#' with the PD validation) (red below the
#' first, amber below the second) and the fixed PSI thresholds.
#'
#' @param x An [scr_lgd()] object.
#' @param newdata `NULL` (the hold-out), an [scr_workout()] object or a
#'   table with the drivers, `lgd_real`, `ead` and the default date.
#'
#' @return An object of class `scr_lgd_validation`: `calibration` (per
#'   pool), `portfolio`, `discrimination`, `stability` (`pools`, `drivers`),
#'   `homogeneity`, `heterogeneity`, `summary` (test, statistic, p, light),
#'   `sample`, `n`.
#'
#' @family irb-lgd
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, n_boot = 20)
#' wo <- scr_workout(scr_demo_lgd, scr_demo_lgd_cashflows, rates = scr_demo_rates, config = cfg)
#' m <- scr_lgd(wo, drivers = c("product", "ltv", "prior_dpd_max"), config = cfg)
#' v <- scr_lgd_validate(m)
#' v
#' v$calibration
#' @export
scr_lgd_validate <- function(x, newdata = NULL) {
  if (!inherits(x, "scr_lgd")) stop("scr_lgd_validate(): `x` must come from scr_lgd().", call. = FALSE)
  cfg <- x$config
  lights <- cfg$pd_lights
  light_p <- function(p) data.table::fifelse(is.na(p), "grey", data.table::fifelse(p <= lights[1], "red", data.table::fifelse(p <= lights[2], "amber", "green")))
  ref <- x$scored[sample == "train"]; ref_idx <- x$bins_idx[sample == "train"]
  if (is.null(newdata)) {
    cur <- x$scored[sample == "holdout"]; cur_idx <- x$bins_idx[sample == "holdout"]; label <- "holdout"
  } else {
    nd <- if (inherits(newdata, "scr_workout")) newdata$rds else data.table::as.data.table(newdata)
    need <- c(x$drivers, "lgd_real", "ead")
    miss <- setdiff(need, names(nd))
    if (length(miss)) stop("scr_lgd_validate(): newdata lacks column(s): ", lst(miss), call. = FALSE)
    if (any(!is.finite(as.double(nd$lgd_real))) || any(!is.finite(as.double(nd$ead)) | as.double(nd$ead) <= 0)) {
      stop("scr_lgd_validate(): newdata needs a finite `lgd_real` and a positive `ead` on every row.", call. = FALSE)
    }
    dt <- .lgd_prepare(nd[, x$drivers, with = FALSE], x$drivers, "scr_lgd_validate")
    pr <- .lgd_predict(x, dt)
    cur <- data.table::data.table(default_id = if ("default_id" %in% names(nd)) as.character(nd$default_id) else as.character(seq_len(nrow(nd))),
                                  sample = "new", default_date = if (x$date_col %in% names(nd)) as.Date(nd[[x$date_col]]) else as.Date(NA),
                                  ead = as.double(nd$ead), lgd_real = as.double(nd$lgd_real), pr)
    cur[, pool := .lgd_pool_of(lgd_pred, x$pools)]
    cur_idx <- data.table::data.table(default_id = cur$default_id, sample = "new")
    if (length(x$cure$features)) {
      w <- apply_woe(x$cure$fit, dt, x$cure$features, "bin")
      for (f in x$cure$features) data.table::set(cur_idx, j = paste0(f, "_cure"), value = match(w[[paste0(f, "_bin")]], x$cure$fit$results[[f]]$bin))
    }
    if (length(x$severity$features)) {
      si <- .cbins_apply_idx(x$severity$fit, dt, x$severity$features)
      for (f in x$severity$features) data.table::set(cur_idx, j = paste0(f, "_sev"), value = si[[f]])
    }
    label <- "newdata"
  }
  if (!nrow(cur)) stop("scr_lgd_validate(): no observation to validate.", call. = FALSE)
  cur[, lgd_est := x$pools$lra[pool]]
  cur[, lgd_dt_est := x$pools$lgd_dt[pool]]

  # -- calibration per pool and portfolio ---------------------------------- #
  ttest <- function(real, est) {
    n <- length(real); if (n < 2L) return(c(NA_real_, NA_real_))
    s <- stats::sd(real); if (!is.finite(s) || s == 0) return(c(if (mean(real) > est) Inf else -Inf, if (mean(real) > est) 0 else 1))
    t <- (mean(real) - est) / (s / sqrt(n)); c(t, 1 - stats::pnorm(t))
  }
  calib <- cur[, {
    tt <- ttest(lgd_real, lgd_est[1])
    list(n = .N, ead = sum(ead), lgd_est = lgd_est[1], lgd_dt = lgd_dt_est[1], real_mean = mean(lgd_real), real_ew = sum(lgd_real * ead) / sum(ead),
         real_lo = mean(lgd_real) - 1.96 * stats::sd(lgd_real) / sqrt(.N), real_hi = mean(lgd_real) + 1.96 * stats::sd(lgd_real) / sqrt(.N),
         t = tt[1], p = tt[2], dt_covers = mean(lgd_real) <= lgd_dt_est[1])
  }, by = pool][order(pool)]
  calib[, light := light_p(p)]
  tt <- ttest(cur$lgd_real, mean(cur$lgd_est))
  reg <- if (isTRUE(stats::sd(cur$lgd_pred) > 0)) stats::coef(stats::lm(cur$lgd_real ~ cur$lgd_pred)) else c(NA_real_, NA_real_)
  portfolio <- data.table::data.table(
    sample = label, n = nrow(cur), real_mean = mean(cur$lgd_real), est_mean = mean(cur$lgd_est), pred_mean = mean(cur$lgd_pred),
    t = tt[1], p = tt[2], light = light_p(tt[2]),
    loss_shortfall = 1 - sum(cur$lgd_real * cur$ead) / sum(cur$lgd_est * cur$ead),
    dt_coverage = sum(cur$lgd_real * cur$ead) <= sum(cur$lgd_dt_est * cur$ead),
    reg_intercept = unname(reg[1]), reg_slope = unname(reg[2]))

  # -- discrimination ------------------------------------------------------ #
  m <- .lgd_metrics(cur$lgd_pred, cur$lgd_real, cur$ead, n_boot = cfg$n_boot, level = cfg$ci_level, seed = cfg$seed, nthread = cfg$nthread)
  g_init <- x$metrics[sample == "train", gauc]
  sigma <- if (is.finite(m$gauc_hi) && is.finite(m$gauc_lo)) (m$gauc_hi - m$gauc_lo) / (2 * stats::qnorm(1 - (1 - cfg$ci_level) / 2)) else NA_real_
  S <- if (is.finite(sigma) && sigma > 0) (g_init - m$gauc) / sigma else NA_real_
  p_S <- if (is.finite(S)) 1 - stats::pnorm(S) else NA_real_
  disc <- data.table::data.table(sample = label, n = m$n, gauc = m$gauc, gauc_lo = m$gauc_lo, gauc_hi = m$gauc_hi, gauc_init = g_init,
                                 S = S, p = p_S, light = light_p(p_S), somers_d = m$somers_d, spearman = m$spearman, r2 = m$r2,
                                 rmse = m$rmse, mae = m$mae, lcr = m$lcr, lcr_lo = m$lcr_lo, lcr_hi = m$lcr_hi)

  # -- stability ----------------------------------------------------------- #
  psi_row <- function(name, base, comp, levels) {
    r <- scr_psi(base, comp, levels = levels, alpha = cfg$psi_alpha)
    data.table::data.table(item = name, psi = r$psi, flag_fixed = r$flag_fixed, critical = r$critical, flag_adjusted = r$flag_adjusted,
                           light = data.table::fifelse(is.na(r$flag_fixed), "grey", data.table::fifelse(r$flag_fixed == "stable", "green",
                                                       data.table::fifelse(r$flag_fixed == "moderate", "amber", "red"))))
  }
  st_pools <- psi_row("pool", ref$pool, cur$pool, levels = x$pools$pool)
  dcols <- setdiff(names(cur_idx), c("default_id", "sample"))
  st_drivers <- if (length(dcols)) data.table::rbindlist(lapply(dcols, function(cn) {
    lv <- sort(unique(c(ref_idx[[cn]], cur_idx[[cn]]))); lv <- lv[!is.na(lv)]
    psi_row(cn, ref_idx[[cn]], cur_idx[[cn]], levels = lv)
  })) else data.table::data.table(item = character(), psi = numeric(), flag_fixed = character(), critical = numeric(), flag_adjusted = character(), light = character())

  # -- homogeneity within, heterogeneity between ---------------------------- #
  welch <- function(a, b) {
    if (length(a) < 2L || length(b) < 2L) return(NA_real_)
    if (stats::sd(a) == 0 && stats::sd(b) == 0) return(if (mean(a) == mean(b)) 1 else 0)
    tryCatch(stats::t.test(a, b)$p.value, error = function(e) NA_real_)
  }
  homog <- cur[, {
    md <- stats::median(lgd_pred); lo <- lgd_real[lgd_pred <= md]; hi <- lgd_real[lgd_pred > md]
    list(n = .N, n_low = length(lo), n_high = length(hi), mean_low = if (length(lo)) mean(lo) else NA_real_,
         mean_high = if (length(hi)) mean(hi) else NA_real_, p = welch(lo, hi))
  }, by = pool][order(pool)]
  homog[, light := light_p(p)]
  pl <- sort(unique(cur$pool))
  heter <- if (length(pl) >= 2L) data.table::rbindlist(lapply(seq_len(length(pl) - 1L), function(i) {
    a <- cur$lgd_real[cur$pool == pl[i]]; b <- cur$lgd_real[cur$pool == pl[i + 1L]]
    p <- welch(a, b)
    data.table::data.table(pool_a = pl[i], pool_b = pl[i + 1L], n_a = length(a), n_b = length(b), mean_a = mean(a), mean_b = mean(b), p = p,
                           light = data.table::fifelse(is.na(p), "grey", data.table::fifelse(p <= 0.05, "green", data.table::fifelse(p <= 0.10, "amber", "red"))))
  })) else data.table::data.table(pool_a = integer(), pool_b = integer(), n_a = integer(), n_b = integer(), mean_a = numeric(), mean_b = numeric(), p = numeric(), light = character())

  worst <- function(l) if (!length(l) || all(is.na(l))) "grey" else if ("red" %in% l) "red" else if ("amber" %in% l) "amber" else "green"
  summary <- data.table::rbindlist(list(
    data.table::data.table(test = "calibration_portfolio_t", statistic = portfolio$t, p = portfolio$p, light = portfolio$light),
    data.table::data.table(test = "calibration_pools_t", statistic = suppressWarnings(max(calib$t, na.rm = TRUE)), p = suppressWarnings(min(calib$p, na.rm = TRUE)), light = worst(calib$light)),
    data.table::data.table(test = "loss_shortfall", statistic = portfolio$loss_shortfall, p = NA_real_, light = if (portfolio$loss_shortfall < -0.10) "red" else if (portfolio$loss_shortfall < 0) "amber" else "green"),
    data.table::data.table(test = "downturn_coverage", statistic = as.numeric(portfolio$dt_coverage), p = NA_real_, light = if (portfolio$dt_coverage) "green" else "red"),
    data.table::data.table(test = "gauc_vs_initial", statistic = disc$S, p = disc$p, light = disc$light),
    data.table::data.table(test = "psi_pools", statistic = st_pools$psi, p = NA_real_, light = st_pools$light),
    data.table::data.table(test = "psi_drivers", statistic = if (nrow(st_drivers)) suppressWarnings(max(st_drivers$psi, na.rm = TRUE)) else NA_real_, p = NA_real_, light = worst(st_drivers$light)),
    data.table::data.table(test = "homogeneity_within_pools", statistic = NA_real_, p = suppressWarnings(min(homog$p, na.rm = TRUE)), light = worst(homog$light)),
    data.table::data.table(test = "heterogeneity_between_pools", statistic = NA_real_, p = suppressWarnings(max(heter$p, na.rm = TRUE)), light = worst(heter$light))
  ))
  summary[!is.finite(statistic), statistic := NA_real_]; summary[!is.finite(p), p := NA_real_]
  structure(list(calibration = calib[], portfolio = portfolio, discrimination = disc,
                 stability = list(pools = st_pools, drivers = st_drivers), homogeneity = homog[], heterogeneity = heter,
                 summary = summary[], sample = label, n = nrow(cur)), class = c("scr_lgd_validation", "list"))
}

#' @export
print.scr_lgd_validation <- function(x, ...) {
  cat(sprintf("<scr_lgd_validation> sample %s | n %s\n", x$sample, n_fmt(x$n)))
  p <- x$portfolio; d <- x$discrimination
  cat(sprintf("  calibration: realised %s vs estimate %s | t %.2f p %.3f [%s] | loss shortfall %s | downturn covers: %s\n",
              fmt_pct(p$real_mean), fmt_pct(p$est_mean), p$t, p$p, p$light, fmt_pct(p$loss_shortfall), p$dt_coverage))
  cat(sprintf("  discrimination: gAUC %.3f [%.3f, %.3f] vs initial %.3f (S %.2f, p %.3f) [%s] | Spearman %.3f | LCR %.3f\n",
              d$gauc, d$gauc_lo, d$gauc_hi, d$gauc_init, d$S, d$p, d$light, d$spearman, d$lcr))
  s <- x$stability$pools
  cat(sprintf("  stability: pool PSI %.4f (%s; adjusted %s) | drivers: %s\n", s$psi, s$flag_fixed, s$flag_adjusted,
              if (nrow(x$stability$drivers)) paste(sprintf("%s %.3f", x$stability$drivers$item, x$stability$drivers$psi), collapse = ", ") else "none"))
  for (i in seq_len(nrow(x$summary))) cat(sprintf("  %-30s %-7s\n", x$summary$test[i], x$summary$light[i]))
  invisible(x)
}

# ============================================================================ #
# export
# ============================================================================ #

#' @rdname scr_export
#' @param validation For the IRB models (`scr_pd`, `scr_lgd`, `scr_ead`): the
#'   matching validation object ([scr_pd_validate()], [scr_lgd_validate()],
#'   [scr_ead_validate()]); `NULL` runs it on the hold-out where possible.
#' @param elbe For `scr_lgd`: an [scr_elbe()] object; `NULL` computes it.
#' @param tag For `scr_lgd` and `scr_ead`: the file tag (`lgd_<tag>.xlsx`,
#'   `ead_<tag>.xlsx`). `scr_pd` names its files after the target
#'   (`pd_<target>.xlsx`) and `scr_capital` after the framework
#'   (`capital_<framework>.xlsx`).
#' @export
scr_export.scr_lgd <- function(x, dir, stamp = TRUE, validation = NULL, elbe = NULL, tag = "model", ...) {
  .scr_local_verbose(x)
  .need_openxlsx()
  out_dir <- .export_dir(dir, stamp)
  tag <- .file_tag(tag)
  v <- validation %||% scr_lgd_validate(x)
  e <- elbe %||% scr_elbe(x)
  ws <- x$workout$summary
  sheets <- list(
    "RDS_Funnel"            = x$workout$funnel,
    "RDS_Summary"           = .kv_table(ws[setdiff(names(ws), c("by_product", "by_year"))]),
    "RDS_By_Product"        = ws$by_product,
    "RDS_By_Year"           = ws$by_year,
    "Recovery_Profile"      = x$workout$recovery_profile,
    "Cure_Bins"             = x$cure$bins,
    "Cure_Coefficients"     = data.frame(term = names(x$cure$coef), estimate = unname(x$cure$coef), stringsAsFactors = FALSE),
    "Severity_Bins"         = x$severity$bins,
    "Severity_Coefficients" = data.frame(term = names(x$severity$coef), estimate = unname(x$severity$coef), engine = x$severity$engine, stringsAsFactors = FALSE),
    "Sign_Check"            = data.table::rbindlist(list(cbind(stage = "cure", x$cure$sign_check), cbind(stage = "severity", x$severity$sign_check)), use.names = TRUE),
    "Pools"                 = x$pools,
    "Downturn"              = x$downturn$table,
    "Downturn_Periods"      = x$downturn$periods,
    "Floors"                = if (is.null(x$floors)) NULL else x$floors$table,
    "ELBE_Grid"             = e$table,
    "Validation"            = v$summary,
    "Validation_Calibration" = v$calibration,
    "Validation_Discrimination" = v$discrimination,
    "Validation_Stability"  = data.table::rbindlist(list(v$stability$pools, v$stability$drivers), use.names = TRUE),
    "Model_Card"            = .kv_table(x$model_card),
    "Decision_Ledger"       = x$ledger)
  files <- list(xlsx = .scr_write_xlsx(sheets, file.path(out_dir, sprintf("lgd_%s.xlsx", tag))),
                sql = file.path(out_dir, sprintf("sql_lgd_%s.sql", tag)))
  writeLines(scr_sql(x), files$sql)
  for (f in files) msg("  %s", f)
  x$files <- files
  invisible(x)
}

## silence the "no visible binding" NOTE for the columns this file uses inside `[.data.table`
utils::globalVariables(c(
  ".i", "amount", "below_reference", "binding", "close_date", "close_last", "closed_at_t_max", "cost_indirect",
  "cum_recovery", "decile", "default_date", "default_id", "discount_rate", "dt", "dt_observed", "dt_type3", "ead",
  "elbe_0", "eta2", "eta2_holdout", "facility_id", "gauc", "impact", "is_cure", "is_incomplete", "kept", "lcr",
  "lgd_dt", "lgd_dt_est", "lgd_est", "lgd_in_default_0", "lgd_lra", "lgd_pred", "lgd_raw", "lgd_real", "lgd_sev",
  "light", "lra", "lra_ew", "lra_moc", "mean_holdout", "merged_from", "method_used", "moc_c", "month",
  "months_in_default", "months_since_default", "n_downturn", "n_holdout", "ok", "p", "pool", "pred_hi", "pred_lo",
  "pred_mean", "product", "psi_flag", "pv", "pv_artificial", "pv_cost", "pv_drawing", "pv_recovery",
  "pv_recovery_cash", "rate", "rate_base", "recovery_artificial", "recovery_extrapolated", "recovery_nominal",
  "reference_value", "rmse", "root_id", "se", "share", "statistic", "status", "status_last", "year"))

# NSE column names used in data.table expressions of this file
utils::globalVariables(c(
  "drawing_nominal",
  "end",
  "lgd_final",
  "somers",
  "start"
))
