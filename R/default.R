# ============================================================================ #
# default.R - default definition engine and default rates by cohort (IRB)
# ============================================================================ #
# The default flag is the target of every IRB model. It is built from a
# monthly panel with a small state machine: a trigger (days past due above a
# threshold with material arrears, or an unlikeliness-to-pay flag) opens an
# event; the event closes after a probation of consecutive months without a
# trigger. Obligor-level flags propagate to every facility of the obligor.
# scr_default_rate() then turns any flagged panel into one-year default
# rates by cohort, the long-run average and its benchmark.
# ============================================================================ #

#' Build the default flag from a monthly panel
#'
#' Applies the standard definition of default to a panel with one row per
#' unit (`id`) and month (`date`): a unit enters default when `dpd` reaches
#' `default_days` and the arrears are material (both `default_abs` in
#' currency units and `default_rel` as a share of `exposure`, when arrears
#' and exposure are supplied), or when `utp` (unlikeliness to pay) is `TRUE`.
#' It leaves default after `default_probation` consecutive months without a
#' trigger (`default_probation_restructured` when `restructured` was `TRUE`
#' at any point of the event). With `obligor` supplied and
#' `default_level = "obligor"`, a unit whose obligor has more than
#' `default_pulling` of its exposure in default is pulled into default too,
#' and a defaulted obligor defaults all its units.
#'
#' Rows must be monthly; gaps are tolerated (the probation counts observed
#' months). The result keeps the row-level flags: they are the product.
#'
#' @param data A `data.frame` or `data.table`, one row per `id` and `date`.
#' @param id,date Column names of the unit identifier and the month.
#' @param dpd Column name of days past due (integer). Optional when `utp`
#'   is given.
#' @param arrears,exposure Column names of the overdue amount and the total
#'   exposure, both optional; when given, the materiality test applies.
#' @param utp Column name of a logical unlikeliness-to-pay flag, optional.
#' @param restructured Column name of a logical distressed-restructuring
#'   flag, optional.
#' @param obligor Column name of the obligor when `id` is a facility,
#'   optional; enables the pulling effect.
#' @param config A [scr_config()]; keys `default_*`.
#'
#' @return An object of class `scr_default`: `flags` (`id`, `date`,
#'   `default` 0/1, `event_id`, `trigger`, `months_in_default`, `cured`),
#'   `events` (one row per event: `event_id`, `id`, `start`, `end`,
#'   `trigger`, `cured`, `months`), `summary`, `ledger` and `config`.
#'
#' @family irb-parameters
#' @examples
#' d <- scr_default(scr_demo_panel, id = "id", date = "ref_date", dpd = "dpd",
#'                  arrears = "arrears", exposure = "exposure",
#'                  restructured = "restructured",
#'                  config = scr_config(verbose = FALSE))
#' d
#' head(d$events)
#' @export
scr_default <- function(data, id, date, dpd = NULL, arrears = NULL, exposure = NULL, utp = NULL,
                        restructured = NULL, obligor = NULL, config = scr_config()) {
  check_config(config, "scr_default")
  cfg <- config
  old <- scr_verbose(isTRUE(cfg$verbose)); on.exit(scr_verbose(old), add = TRUE)
  dt <- data.table::as.data.table(data)
  need <- c(id, date, dpd, arrears, exposure, utp, restructured, obligor)
  miss <- setdiff(need, names(dt))
  if (length(miss)) stop("scr_default(): column(s) not found: ", lst(miss), call. = FALSE)
  if (is.null(dpd) && is.null(utp)) stop("scr_default(): give `dpd`, `utp` or both.", call. = FALSE)
  t0 <- Sys.time()

  p <- data.table::data.table(id = as.character(dt[[id]]), date = as.Date(dt[[date]]))
  p[, dpd := if (is.null(dpd)) 0L else as.integer(dt[[dpd]])]
  p[, utp := if (is.null(utp)) FALSE else isTRUE_vec(dt[[utp]])]
  p[, restr := if (is.null(restructured)) FALSE else isTRUE_vec(dt[[restructured]])]
  p[, expo := if (is.null(exposure)) NA_real_ else as.double(dt[[exposure]])]
  p[, arr := if (is.null(arrears)) NA_real_ else as.double(dt[[arrears]])]
  p[, obligor := if (is.null(obligor)) id else as.character(dt[[obligor]])]
  if (anyNA(p$id) || anyNA(p$date)) stop("scr_default(): `id` and `date` cannot be missing.", call. = FALSE)
  if (anyDuplicated(p, by = c("id", "date"))) stop("scr_default(): duplicated (id, date) rows.", call. = FALSE)
  data.table::setorder(p, id, date)

  # -- row triggers ------------------------------------------------------- #
  material <- if (!is.null(arrears) && !is.null(exposure)) {
    (p$arr >= cfg$default_abs) & (p$arr / pmax(p$expo, .Machine$double.eps) > cfg$default_rel)
  } else if (!is.null(arrears)) {
    p$arr >= cfg$default_abs
  } else TRUE
  p[, trig_dpd := !is.na(dpd) & dpd >= cfg$default_days & material %in% TRUE]
  p[, trig := trig_dpd | utp %in% TRUE]

  # -- state machine per unit ------------------------------------------- #
  # vectorised over the sorted panel; units with a restructuring (when its
  # probation differs) run the row loop. Results are placed by row position:
  # split() orders the units by the locale collation, setorder() by the C
  # locale, so a positional rbindlist() of split() groups can misalign rows.
  n_ids <- data.table::uniqueN(p$id)
  p[, c("default", "ev", "trigger", "months", "cured") :=
      .default_flags(p$id, p$trig, p$utp, p$restr, cfg$default_probation, cfg$default_probation_restructured, cfg$nthread)]

  # -- pulling effect at obligor level ------------------------------------ #
  pulled <- 0L
  if (!is.null(obligor) && identical(cfg$default_level, "obligor")) {
    p[, share := {
      e <- if (all(is.na(expo))) rep(1, .N) else data.table::fifelse(is.na(expo), 0, expo)
      s <- sum(e); if (s > 0) sum(e[default == 1L]) / s else mean(default == 1L)
    }, by = c("obligor", "date")]
    pull <- p$share > cfg$default_pulling & p$default == 0L
    pulled <- sum(pull)
    if (pulled) {
      p[pull, `:=`(default = 1L, trigger = "pulling")]
      # re-run the probation with the pulled months as triggers
      p[, trig := trig | trigger %in% "pulling"]
      p[, c("default", "ev", "trigger2", "months", "cured") :=
          .default_flags(p$id, p$trig, p$utp, p$restr, cfg$default_probation, cfg$default_probation_restructured, cfg$nthread)]
      p[trigger2 != "" & !(trigger %in% "pulling"), trigger := trigger2]
      p[, trigger2 := NULL]
    }
    p[, share := NULL]
  }

  p[, event_id := data.table::fifelse(ev > 0L, paste0(id, "#", ev), NA_character_)]
  # grouped aggregates GForce can run (no per-event R call), then the trigger
  # a panel with no default month gives an empty table (j is not evaluated on
  # an empty subset, so min()/max() never see a zero-length Date)
  events <- if (any(p$default == 1L)) {
    p[default == 1L, list(id = id[1L], start = min(date), end = max(date), months = .N,
                          cured = max(cured)), by = "event_id"]
  } else {
    data.table::data.table(event_id = character(), id = p$id[0L], start = p$date[0L], end = p$date[0L],
                           months = integer(), cured = p$cured[0L])
  }
  ev_trig <- p[default == 1L & trigger != "", list(trigger = trigger[1L]), by = "event_id"]
  events[, trigger := ev_trig$trigger[match(event_id, ev_trig$event_id)]]
  data.table::setcolorder(events, c("event_id", "id", "start", "end", "trigger", "months", "cured"))
  data.table::setorder(events, id, start)

  flags <- p[, list(id, date, default, event_id, trigger, months_in_default = months, cured)]
  n_ev <- nrow(events)
  by_trig <- if (n_ev) prop.table(table(events$trigger)) else numeric()
  summary <- list(
    n_ids = n_ids, n_rows = nrow(p), n_events = n_ev,
    share_by_trigger = as.list(by_trig),
    median_months_in_default = if (n_ev) stats::median(events$months) else NA_real_,
    share_cured = if (n_ev) mean(events$cured == 1L) else NA_real_,
    n_pulled_rows = pulled,
    rules = list(days = cfg$default_days, abs = cfg$default_abs, rel = cfg$default_rel, level = cfg$default_level,
                 probation = cfg$default_probation, probation_restructured = cfg$default_probation_restructured,
                 pulling = cfg$default_pulling, materiality = !is.null(arrears)))
  ledger <- data.table::data.table(
    action = "default_definition",
    detail = sprintf("dpd >= %d%s%s; probation %d (%d after restructuring); level %s%s",
                     cfg$default_days,
                     if (!is.null(arrears)) sprintf(" and arrears >= %s", format(cfg$default_abs)) else "",
                     if (!is.null(arrears) && !is.null(exposure)) sprintf(" and arrears/exposure > %s", format(cfg$default_rel)) else "",
                     cfg$default_probation, cfg$default_probation_restructured, cfg$default_level,
                     if (!is.null(obligor)) sprintf("; pulling > %s", format(cfg$default_pulling)) else ""),
    date = format(Sys.Date()))
  msg("  default flag: %s units, %s events, %s cured (%.2fs)", n_fmt(n_ids), n_fmt(n_ev),
      if (n_ev) fmt_pct(summary$share_cured) else "-", as.numeric(difftime(Sys.time(), t0, units = "secs")))
  structure(list(flags = flags, events = events, summary = summary, ledger = ledger, config = cfg),
            class = c("scr_default", "list"))
}

#' @keywords internal
#' @noRd
isTRUE_vec <- function(x) {
  if (is.logical(x)) x %in% TRUE else as.integer(x) %in% 1L
}

#' State machine of one unit: default rows, event numbers, trigger, months, cure
#' @keywords internal
#' @noRd
.default_run <- function(trig, utp, restr, probation, probation_restr) {
  n <- length(trig)
  default <- integer(n); ev <- integer(n); trigger <- character(n); months <- integer(n); cured <- integer(n)
  in_def <- FALSE; quiet <- 0L; k <- 0L; m <- 0L; ev_restr <- FALSE
  for (i in seq_len(n)) {
    if (!in_def) {
      if (isTRUE(trig[i])) {
        in_def <- TRUE; k <- k + 1L; m <- 0L; quiet <- 0L; ev_restr <- isTRUE(restr[i])
        trigger[i] <- if (isTRUE(utp[i])) "utp" else "dpd"
      }
    }
    if (in_def) {
      m <- m + 1L
      ev_restr <- ev_restr || isTRUE(restr[i])
      default[i] <- 1L; ev[i] <- k; months[i] <- m
      if (isTRUE(trig[i])) quiet <- 0L else quiet <- quiet + 1L
      need <- if (ev_restr) probation_restr else probation
      if (quiet >= need) { cured[i] <- 1L; in_def <- FALSE; quiet <- 0L }
    }
  }
  list(default = default, ev = ev, trigger = trigger, months = months, cured = cured)
}

#' The state machine of .default_run() over a whole panel sorted by (id, date)
#'
#' With a constant probation `P` the machine has a closed form: a row is in
#' default when the last trigger of its unit is at most `P` observed months
#' before it (`dist <= P`), it cures when `dist == P`, and an event starts on
#' a default row whose predecessor was not in default or had just cured.
#' Units with a restructuring flag (when `probation_restr != probation`)
#' change `P` inside an event and go through the row loop.
#' @keywords internal
#' @noRd
.default_flags <- function(id, trig, utp, restr, probation, probation_restr, nthread = 1L) {
  n <- length(trig)
  if (!n) return(list(integer(), integer(), character(), integer(), integer()))
  i <- seq_len(n)
  first <- c(TRUE, id[-1L] != id[-n])
  gstart <- cummax(data.table::fifelse(first, i, 0L))           # first row of the unit
  lt <- cummax(data.table::fifelse(trig %in% TRUE, i, 0L))      # last trigger so far
  dist <- i - lt
  default <- as.integer(lt >= gstart & dist <= probation)
  cured <- as.integer(default == 1L & dist == probation)
  start <- default == 1L & (first | c(0L, default[-n]) == 0L | c(0L, cured[-n]) == 1L)
  cs <- cumsum(start)
  ev <- data.table::fifelse(default == 1L, cs - c(0L, cs)[gstart], 0L)
  estart <- cummax(data.table::fifelse(start, i, 0L))
  months <- data.table::fifelse(default == 1L, i - estart + 1L, 0L)
  trigger <- data.table::fifelse(start, data.table::fifelse(utp %in% TRUE, "utp", "dpd"), "")
  if (!identical(as.integer(probation_restr), as.integer(probation)) && any(restr %in% TRUE)) {
    rix <- which(id %in% unique(id[restr %in% TRUE]))
    grp <- split(rix, id[rix])
    runs <- .scr_lapply(grp, function(ix) .default_run(trig[ix], utp[ix], restr[ix], probation, probation_restr),
                        nthread = nthread, fork_only = TRUE)
    pos <- unlist(grp, use.names = FALSE)
    r <- data.table::rbindlist(runs)
    default[pos] <- r$default; ev[pos] <- r$ev; trigger[pos] <- r$trigger; months[pos] <- r$months; cured[pos] <- r$cured
  }
  list(default, as.integer(ev), trigger, as.integer(months), cured)
}

#' @export
print.scr_default <- function(x, ...) {
  s <- x$summary
  cat(sprintf("<scr_default> %s units x %s rows | %s events | level %s\n",
              n_fmt(s$n_ids), n_fmt(s$n_rows), n_fmt(s$n_events), s$rules$level))
  cat(sprintf("  rule: dpd >= %d%s | probation %d months (%d after restructuring)%s\n",
              s$rules$days, if (isTRUE(s$rules$materiality)) " with material arrears" else "",
              s$rules$probation, s$rules$probation_restructured,
              if (s$n_pulled_rows > 0) sprintf(" | %s rows pulled by the obligor", n_fmt(s$n_pulled_rows)) else ""))
  if (s$n_events) {
    cat(sprintf("  triggers: %s | median %s months in default | cured %s\n",
                paste(sprintf("%s %s", names(s$share_by_trigger), fmt_pct(unlist(s$share_by_trigger))), collapse = ", "),
                format(s$median_months_in_default), fmt_pct(s$share_cured)))
  }
  invisible(x)
}

# -- default rates by cohort -------------------------------------------------- #

#' One-year default rates by cohort and the long-run average
#'
#' From a flagged panel (an [scr_default()] object or any table with a 0/1
#' default column by unit and month), computes for every cohort start the
#' population of non-defaulted units, the share that defaults within
#' `horizon` months, optionally by `grade` or `segment` (as observed at the
#' cohort start) and exposure-weighted when `exposure` is given. The long-run
#' average is the arithmetic mean of the cohort rates. When the analyst
#' proposes an adjusted value (`lra_adjusted`, for instance after judging
#' that the period lacks bad years), it is benchmarked against the larger of
#' the last five years' mean and the whole period's mean, and a flag records
#' when it sits below that benchmark.
#'
#' @param x An `scr_default` or a `data.frame`/`data.table`.
#' @param id,date,default Column names (ignored for an `scr_default`).
#' @param horizon Months of the default window after the cohort start.
#' @param by Cohort frequency: `"month"`, `"quarter"` or `"year"`.
#' @param grade,segment Optional column names observed at the cohort start.
#' @param exposure Optional column name; adds exposure-weighted rates.
#' @param lra_adjusted Optional adjusted long-run average proposed by the
#'   analyst, in `[0, 1]`; benchmarked and flagged, never applied.
#' @param config A [scr_config()]; only `pd_dr_by` is read (the default of
#'   `by`).
#'
#' @return An object of class `scr_dr`: `table` (cohort rates, by grade or
#'   segment when given), `portfolio` (one row per cohort: `n`, `defaults`,
#'   `dr`), `lra` (`mean`, `weighted_mean`, `recent5_mean`, `benchmark`,
#'   `adjusted`, `flag_below_benchmark`, `min`, `max`, `sd`, `n_cohorts`,
#'   `years`), `horizon` and `by`.
#'
#' @family irb-parameters
#' @examples
#' d <- scr_default(scr_demo_panel, id = "id", date = "ref_date", dpd = "dpd",
#'                  config = scr_config(verbose = FALSE))
#' dr <- scr_default_rate(d, by = "quarter")
#' dr
#' dr$table
#' @export
scr_default_rate <- function(x, id = "id", date = "date", default = "default", horizon = 12L,
                             by = NULL, grade = NULL, segment = NULL, exposure = NULL, lra_adjusted = NULL,
                             config = scr_config()) {
  check_config(config, "scr_default_rate")
  if (!is.null(lra_adjusted)) .scr_num1(lra_adjusted, "lra_adjusted", lower = 0, upper = 1)
  if (!is.numeric(horizon) || length(horizon) != 1L || is.na(horizon) || horizon < 1L) {
    stop("scr_default_rate(): `horizon` must be a positive number of months.", call. = FALSE)
  }
  by <- by %||% config$pd_dr_by
  by <- match.arg(by, c("month", "quarter", "year"))
  horizon <- as.integer(horizon)
  if (inherits(x, "scr_default")) {
    dt <- x$flags[, list(id, date, default)]
    if (!is.null(grade) || !is.null(segment) || !is.null(exposure)) {
      stop("scr_default_rate(): pass a table (not an scr_default) to use `grade`, `segment` or `exposure`.", call. = FALSE)
    }
  } else {
    dt <- data.table::as.data.table(x)
    need <- c(id, date, default, grade, segment, exposure)
    miss <- setdiff(need, names(dt))
    if (length(miss)) stop("scr_default_rate(): column(s) not found: ", lst(miss), call. = FALSE)
    dt <- data.table::data.table(id = as.character(dt[[id]]), date = as.Date(dt[[date]]),
                                 default = as.integer(dt[[default]]),
                                 grade = if (is.null(grade)) NA_character_ else as.character(dt[[grade]]),
                                 segment = if (is.null(segment)) NA_character_ else as.character(dt[[segment]]),
                                 exposure = if (is.null(exposure)) NA_real_ else as.double(dt[[exposure]]))
  }
  if (!"grade" %in% names(dt)) dt[, grade := NA_character_]
  if (!"segment" %in% names(dt)) dt[, segment := NA_character_]
  if (!"exposure" %in% names(dt)) dt[, exposure := NA_real_]
  data.table::setorder(dt, id, date)
  dates <- sort(unique(dt$date))
  mth <- as.integer(format(dates, "%m"))
  starts <- switch(by, month = dates, quarter = dates[mth %in% c(1L, 4L, 7L, 10L)], year = dates[mth == 1L])
  last_ok <- max(dates)
  starts <- starts[.add_months(starts, horizon) <= last_ok]
  if (!length(starts)) stop("scr_default_rate(): no cohort has a complete ", horizon, "-month window.", call. = FALSE)

  keys <- c("grade", "segment")
  # every cohort at once: the population at each start and one rolling join
  # for the first default after it
  # plain logical subset (a data.table auto-index on a big panel costs a sort)
  sel <- dt$default %in% 0L & as.double(dt$date) %in% as.double(starts)
  pop <- dt[sel, list(id, cohort = date, grade, segment, exposure)]
  pop[, d := .dr_outcome(dt[default == 1L, list(id, date)], id, cohort, .add_months_u(cohort, horizon))]
  tab <- pop[, list(n = .N, defaults = sum(d), dr = mean(d),
                    ead = if (all(is.na(exposure))) NA_real_ else sum(exposure, na.rm = TRUE),
                    dr_weighted = if (all(is.na(exposure))) NA_real_ else sum(exposure * d, na.rm = TRUE) / sum(exposure, na.rm = TRUE)),
             by = c("cohort", keys)]
  data.table::setcolorder(tab, c("cohort", keys))
  if (all(is.na(tab$grade))) tab[, grade := NULL]
  if (all(is.na(tab$segment))) tab[, segment := NULL]
  if (all(is.na(tab$ead))) tab[, c("ead", "dr_weighted") := NULL]
  data.table::setorderv(tab, c("cohort", intersect(keys, names(tab))))

  # long-run average on the whole population (grade/segment rows collapsed)
  port <- tab[, list(n = sum(n), defaults = sum(defaults)), by = "cohort"][, dr := defaults / n]
  yrs <- as.numeric(difftime(max(port$cohort), min(port$cohort), units = "days")) / 365.25
  recent <- port[cohort > .add_months(max(cohort), -60L)]
  lra <- list(mean = mean(port$dr), weighted_mean = sum(port$defaults) / sum(port$n),
              recent5_mean = mean(recent$dr), n_cohorts = nrow(port), years = yrs)
  # the benchmark of an adjusted long-run average: the larger of the last
  # five years' mean and the whole period's mean
  lra$benchmark <- max(lra$recent5_mean, lra$mean)
  lra$adjusted <- lra_adjusted
  lra$flag_below_benchmark <- !is.null(lra_adjusted) && lra_adjusted < lra$benchmark - 1e-12
  lra$min <- min(port$dr); lra$max <- max(port$dr); lra$sd <- stats::sd(port$dr)
  structure(list(table = tab, portfolio = port, lra = lra, horizon = horizon, by = by),
            class = c("scr_dr", "list"))
}

#' Shift a Date by whole months (same day when it exists, else month end)
#' @keywords internal
#' @noRd
.add_months <- function(d, k) {
  lt <- as.POSIXlt(d)
  lt$mon <- lt$mon + as.integer(k)
  out <- as.Date(lt)
  # roll back overflowed days (e.g. 31 Jan + 1 month) to the month end
  bad <- as.integer(format(out, "%d")) < as.integer(format(d, "%d")) & !is.na(out)
  if (any(bad)) out[bad] <- out[bad] - as.integer(format(out[bad], "%d"))
  out
}

#' .add_months() evaluated once per distinct date (cohort starts repeat on
#' every row of the population)
#' @keywords internal
#' @noRd
.add_months_u <- function(d, k) {
  u <- unique(d)
  .add_months(u, k)[match(as.double(d), as.double(u))]
}

#' Default within (t0, t1]: 1 when the unit has a default row dated after t0
#' and at most t1
#'
#' One rolling join (next observation carried backward) finds the first
#' default date strictly after `t0` for every (id, t0).
#' @keywords internal
#' @noRd
.dr_outcome <- function(def_rows, id, t0, t1) {
  if (!length(id)) return(integer())
  if (!nrow(def_rows)) return(integer(length(id)))
  x <- unique(data.table::data.table(id = as.character(def_rows$id), t = as.double(def_rows$date)))
  q <- data.table::data.table(id = as.character(id), t = as.double(t0) + 0.5)
  nxt <- x[q, on = c("id", "t"), roll = -Inf, mult = "first", x.t]
  as.integer(!is.na(nxt) & nxt <= as.double(t1))
}

#' @export
print.scr_dr <- function(x, ...) {
  l <- x$lra
  cat(sprintf("<scr_dr> %d %sly cohorts over %.1f years | horizon %d months\n", l$n_cohorts, x$by, l$years, x$horizon))
  cat(sprintf("  default rate: mean %s | weighted %s | min %s | max %s | sd %s\n",
              fmt_pct(l$mean, 2), fmt_pct(l$weighted_mean, 2), fmt_pct(l$min, 2), fmt_pct(l$max, 2), fmt_pct(l$sd, 2)))
  cat(sprintf("  long-run average %s | benchmark %s (max of last-5-years %s and all-years %s)\n",
              fmt_pct(l$mean, 2), fmt_pct(l$benchmark, 2), fmt_pct(l$recent5_mean, 2), fmt_pct(l$mean, 2)))
  if (!is.null(l$adjusted)) {
    cat(sprintf("  adjusted long-run average %s%s\n", fmt_pct(l$adjusted, 2),
                if (isTRUE(l$flag_below_benchmark)) " - BELOW the benchmark, justify and cover with a margin" else " - at or above the benchmark"))
  }
  if (l$years < 5) cat("  note: fewer than five years of cohorts; the average is not a long-run one yet\n")
  invisible(x)
}

# NSE column names used in data.table expressions of this file
utils::globalVariables(c(
  "arr",
  "cured",
  "default",
  "defaults",
  "ev",
  "event_id",
  "expo",
  "restr",
  "start",
  "trig",
  "trig_dpd",
  "trigger",
  "trigger2",
  "x.t"
))
