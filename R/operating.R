# ============================================================================ #
# operating.R - the operating point of a score under constraints
# ============================================================================ #
# The cells of the count table are accumulated from the chosen end of the
# score: every cell boundary is a candidate cut, and every cumulative
# quantity (volume, events, cost, value, capacity) is a cumulative sum over
# cells, so the whole curve costs O(K) after the one pass over the rows.
# ============================================================================ #

#' Operating point of a score under constraints
#'
#' Accumulates the score from one end, one score value at a time, and finds
#' the cut that maximizes the value of the decision under volume, budget,
#' daily capacity and event-rate constraints: how many customers to target,
#' how many alerts to raise, or how many applicants to approve.
#'
#' @section Side:
#'
#' `side = "event"` selects from the event-rich end of the score (targeting
#' under propensity, alerting under fraud, collections under credit);
#' `side = "safe"` accepts from the safe end (approval under credit). The
#' default follows the objective and the direction: propensity selects from
#' the event-rich end (`"event"`); risk with `higher_is_riskier` (fraud)
#' alerts from the event-rich end (`"event"`); risk with `higher_is_safer`
#' (credit) approves from the safe end (`"safe"`).
#'
#' The selected rows are `score >= cut` when the selection starts at the high
#' scores and `score < cut` when it starts at the low ones, the convention of
#' [scr_cutoff()]. Every cut sits between two adjacent distinct scores (or on
#' a bucket edge when the scores were pooled into `max_cells` cells), as in
#' [scr_bands()]; the last row of the curve selects every row (`cut` is
#' `-Inf` or `Inf`).
#'
#' @section Curve:
#'
#' One row per candidate cut, in increasing depth: `cut`, `depth` (share of
#' the volume selected), `n_sel`, `events_sel`, `rate_sel` with its Jeffreys
#' interval `rate_lo`, `rate_hi` (on the Kish effective size under weights),
#' `capture` (share of all events selected), `lift` (`rate_sel` over the
#' overall rate), `marginal_rate` (the event rate of the score value just
#' added, smoothed by pool adjacent violators toward the event-rich end of
#' the score), `cost` (`cost_select * n_sel`), `value` and `feasible`. With
#' a date, `day_q` (the `day_quantile` quantile of the selected volume per
#' day) and `pct_days_over` (share of days above `max_per_day`).
#'
#' The economics:
#' \itemize{
#'   \item side `"event"`: `value = gain_event * events_sel - cost_select *
#'     n_sel`; with a `value` column and no `gain_event`, the sum of the
#'     value of the selected events replaces `gain_event * events_sel`.
#'   \item side `"safe"`: `value = revenue_good * nonevents_sel - loss_bad *
#'     events_sel - cost_select * n_sel`, the cumulative profit of
#'     [scr_strategy()] when `cost_select = 0` (a missing one of
#'     `revenue_good` and `loss_bad` counts as 0).
#' }
#' Without economics (`gain_event`, a `value` column, `revenue_good` or
#' `loss_bad`), `value` is `NA`. Non-events are rows with a known outcome
#' that are not events; rows with a missing outcome count in the volume and
#' the cost only.
#'
#' The curve is thinned to about `n_points` rows evenly spread in depth; the
#' rows of the optimum, the last row meeting each constraint, the deepest
#' row and the rows nearest 1%, 5%, 10%, 20% and 50% are always kept. The
#' optimum and the constraints are evaluated on every cell boundary.
#'
#' @section Constraints and optimum:
#'
#' `max_n` (`n_sel <= max_n`), `max_share` (`depth <= max_share`), `budget`
#' (`cost_select * n_sel <= budget`), `max_per_day` (`day_q <=
#' max_per_day`), `min_rate` (`rate_sel >= min_rate`, side `"event"`) and
#' `max_rate` (`rate_sel <= max_rate`, the event rate among the accepted,
#' side `"safe"`). A row is feasible when it meets every constraint given.
#'
#' The optimum is the feasible row with the highest value (the smallest
#' depth on a tie); without economics, the deepest feasible row. A
#' constraint is `binding` when dropping it alone, the others kept,
#' improves the optimum: a higher value, or a greater depth when there are
#' no economics. A constraint that is slack at the optimum is therefore
#' never named, and a constraint that stops the curve at the row that is
#' the best anyway is not binding either. Constraints that stop the optimum
#' at the same row bind jointly (none improves it alone) and are named
#' together. When no constraint binds, `binding` is `"value"` with
#' economics (no row is worth more than the optimum) and `"end of the
#' curve"` without (every row is selected).
#' `shadow_price` is the marginal value of the next score value beyond the
#' optimum, per additional selected case: `gain_event * marginal_rate -
#' cost_select` on side `"event"` (with a `value` column, the smoothed event
#' value per case of that score value), `revenue_good * (1 -
#' marginal_rate) - loss_bad * marginal_rate - cost_select` on side
#' `"safe"`. It is what one more selected case is worth when a constraint
#' binds; divide it by `cost_select` for the value of one more unit of
#' budget. When no row is feasible, the optimum is `NA` and a warning names
#' the constraints that the first row already breaks.
#'
#' @section Daily capacity:
#'
#' A day is a distinct value of the date column: with a monthly date, read
#' the capacity per month. For every candidate cut, the selected volume of
#' each day is a cumulative sum over a table of counts per day and score
#' value; `day_q` is its quantile across the days (type 7 of
#' [stats::quantile()]). The quantile grows with the depth, so the deepest
#' cut within `max_per_day` is found by bisection. Rows with a missing date
#' count in the curve but not in the daily volumes. A scorecard uses the
#' dates of its scored sample when they exist. A date-time column counts
#' every distinct time as a day: convert it with [as.Date()] first.
#'
#' @inheritParams scr_bands
#' @param x An object from [scr_scorecard()], or a `data.frame` with one row
#'   per scored case (or one row per score value with `counts = TRUE`).
#' @param side `"event"` or `"safe"`; `NULL` follows the objective and the
#'   direction (see the section Side).
#' @param gain_event Gain per selected event (side `"event"`).
#' @param cost_select Cost per selected case (default 0).
#' @param revenue_good,loss_bad Revenue per accepted non-event and loss per
#'   accepted event (side `"safe"`).
#' @param max_n Largest volume selected.
#' @param max_share Largest share of the volume selected, in (0, 1].
#' @param budget Largest cost, `cost_select * n_sel`; needs a positive
#'   `cost_select`.
#' @param max_per_day Largest selected volume per day, read at the
#'   `day_quantile` quantile of the days; needs a date.
#' @param day_quantile Quantile of the daily volume compared with
#'   `max_per_day` (0.9: nine days in ten within capacity).
#' @param date For a data.frame: name of a date column, for the daily
#'   capacity. For a scorecard: a column of the scored sample; `NULL` uses
#'   its `date` column when present.
#' @param min_rate Smallest event rate among the selected (side `"event"`).
#' @param max_rate Largest event rate among the accepted (side `"safe"`).
#' @param sample For a scorecard: the sample the curve is read on
#'   (`"holdout"`). For a data.frame: the name of a column with sample
#'   labels, as in [scr_bands()]; the curve is read on `study`.
#' @param study For a data.frame: the label of the sample the curve is read
#'   on; `NULL` takes the first label other than the first level (the
#'   reference of [scr_bands()]), or the only one.
#' @param n_points About how many rows of the curve to keep.
#' @param level Confidence level of the Jeffreys intervals. For a scorecard,
#'   `NULL` uses `config$study_level` (0.95).
#' @param value Optional column of the value of every case: with no
#'   `gain_event`, the value of the selected events is the gain.
#'
#' @return An object of class `c("scr_operating", "list")`:
#'   \describe{
#'     \item{`curve`}{The thinned curve (see the section Curve).}
#'     \item{`optimum`}{One row: the columns of the curve at the optimum,
#'       `binding`, `shadow_price` and `next_rate` (the marginal rate of the
#'       next score value).}
#'     \item{`constraints`}{One row per constraint given: `constraint`,
#'       `limit`, `at_optimum` (the constrained quantity at the optimum) and
#'       `binding`.}
#'     \item{`side`, `objective`, `direction`, `target`, `sample`, `level`,
#'       `gain_event`, `cost_select`, `revenue_good`, `loss_bad`,
#'       `day_quantile`, `call`}{The settings.}
#'     \item{`economics`, `value_column`}{Whether the curve has a value, and
#'       whether it comes from the `value` column instead of `gain_event`.}
#'     \item{`select_high`}{`TRUE` when the selection starts at the high
#'       scores (`score >= cut`), `FALSE` at the low ones (`score < cut`).}
#'     \item{`n_cells`, `n_days`, `quantized`, `weighted`}{The number of
#'       candidate cuts (score values or pooled cells), the number of
#'       distinct dates (`NA` without a date), whether the scores were pooled
#'       into `max_cells` cells, and whether weights were used.}
#'     \item{`message`}{`NA`, or the explanation of an infeasible or a
#'       loss-making optimum.}
#'   }
#'
#' @references
#' Brown, L. D., Cai, T. T. and DasGupta, A. (2001). Interval estimation for
#' a binomial proportion. *Statistical Science*, 16(2), 101-133.
#' \doi{10.1214/ss/1009213286}
#'
#' Thomas, L. C., Crook, J. and Edelman, D. (2017). *Credit Scoring and Its
#' Applications*, 2nd edition. SIAM. \doi{10.1137/1.9781611974560}
#'
#' @seealso [scr_cutoff()] and [scr_strategy()] for the cut-off sweep and the
#'   strategy table of a scorecard, [scr_claims()] to test statements about
#'   the selected rates.
#' @family score-studies
#' @examples
#' set.seed(1)
#' x <- rnorm(5000)
#' d <- data.frame(score = round(500 + 50 * x),
#'                 y = rbinom(5000, 1, plogis(-1.5 + 1.2 * x)),
#'                 day = as.Date("2026-01-01") + sample(0:29, 5000, TRUE))
#' # targeting under propensity: a gain per responder, a cost per contact, a budget
#' op <- scr_operating(d, objective = "propensity", gain_event = 40, cost_select = 6,
#'                     budget = 6000)
#' op
#' op$optimum[, c("cut", "depth", "n_sel", "rate_sel", "value", "binding", "shadow_price")]
#'
#' # alerting under a daily capacity: at most 40 alerts on nine days in ten
#' scr_operating(d, objective = "propensity", max_per_day = 40, date = "day")$optimum
#' @export
scr_operating <- function(x, ...) UseMethod("scr_operating")

#' @rdname scr_operating
#' @export
scr_operating.scr_scorecard <- function(x, side = NULL, gain_event = NULL, cost_select = 0, revenue_good = NULL,
                                        loss_bad = NULL, max_n = NULL, max_share = NULL, budget = NULL,
                                        max_per_day = NULL, day_quantile = 0.9, date = NULL, min_rate = NULL,
                                        max_rate = NULL, sample = "holdout", n_points = 200L, level = NULL,
                                        max_cells = 1e5, ...) {
  fn <- "scr_operating"
  .study_dots(list(...), fn)
  .study_chr1(sample, "sample", fn)
  s <- x$samples[[sample]]
  if (is.null(s)) stop(fn, "(): sample '", sample, "' not in the scorecard (", lst(names(x$samples)), ").", call. = FALSE)
  if (!is.null(date)) .study_chr1(date, "date", fn)
  by <- date %||% if ("date" %in% names(s)) "date" else NULL
  inp <- .study_input(x, sample = sample, reference = sample, by = by, max_cells = max_cells, fn = fn)
  .operating_fit(inp, sample, side, gain_event, cost_select, revenue_good, loss_bad, max_n, max_share, budget,
                 max_per_day, day_quantile, min_rate, max_rate, level %||% x$config$study_level %||% 0.95,
                 n_points, !is.null(by), sys.call())
}

#' @rdname scr_operating
#' @export
scr_operating.data.frame <- function(x, side = NULL, gain_event = NULL, cost_select = 0, revenue_good = NULL,
                                     loss_bad = NULL, max_n = NULL, max_share = NULL, budget = NULL,
                                     max_per_day = NULL, day_quantile = 0.9, date = NULL, min_rate = NULL,
                                     max_rate = NULL, sample = NULL, n_points = 200L, score = "score", y = "y",
                                     objective = "risk", direction = NULL, weight = NULL, value = NULL,
                                     study = NULL, counts = FALSE, n = "n", events = "events",
                                     value_events = NULL, level = 0.95, max_cells = 1e5, ...) {
  fn <- "scr_operating"
  .study_dots(list(...), fn)
  if (!is.null(date)) .study_chr1(date, "date", fn)
  smp <- "all"
  if (!is.null(sample)) {
    .study_chr1(sample, "sample", fn)
    if (!sample %in% names(x)) stop(fn, "(): column(s) ", sample, " not in `x`.", call. = FALSE)
    sv <- x[[sample]]
    lv <- if (is.factor(sv)) levels(droplevels(sv)) else sort(unique(as.character(sv[!is.na(sv)])))
    if (!is.null(study)) .study_chr1(study, "study", fn)
    smp <- study %||% if (length(lv) > 1L) lv[2] else lv[1]
  } else if (!is.null(study)) {
    stop(fn, "(): `study` needs a `sample` column.", call. = FALSE)
  }
  # only the rows of the sample read are aggregated
  inp <- .study_input(x, score = score, y = y, objective = objective, direction = direction, weight = weight,
                      value = value, value_events = value_events, sample = sample,
                      reference = if (!is.null(sample)) smp, study = if (!is.null(sample)) smp,
                      counts = counts, n = n, events = events, by = date, max_cells = max_cells, fn = fn)
  .operating_fit(inp, smp, side, gain_event, cost_select, revenue_good, loss_bad, max_n, max_share, budget,
                 max_per_day, day_quantile, min_rate, max_rate, level, n_points, !is.null(date), sys.call())
}

#' Optional number in a range; NULL passes
#' @keywords internal
#' @noRd
.op_num <- function(x, name, lower = -Inf, upper = Inf, open_lower = FALSE) {
  if (is.null(x)) return(NULL)
  .scr_num1(x, name, lower = lower, upper = upper, open_lower = open_lower)
}

#' Curve, constraints and optimum on one sample of a count table
#' @keywords internal
#' @noRd
.operating_fit <- function(inp, smp, side, gain_event, cost_select, revenue_good, loss_bad, max_n, max_share,
                           budget, max_per_day, day_quantile, min_rate, max_rate, level, n_points, daily, call) {
  fn <- "scr_operating"
  level <- .study_level(level, fn)
  n_points <- .study_whole(n_points, "n_points", fn, lower = 2)
  riskier <- identical(inp$direction, "higher_is_riskier")
  side <- side %||% if (identical(inp$objective, "propensity") || riskier) "event" else "safe"
  if (!is.character(side) || length(side) != 1L || !side %in% c("event", "safe")) {
    stop(fn, "(): `side` must be NULL, \"event\" or \"safe\".", call. = FALSE)
  }
  ev <- identical(side, "event")
  .op_num(gain_event, "gain_event", lower = 0); .scr_num1(cost_select, "cost_select", lower = 0)
  .op_num(revenue_good, "revenue_good", lower = 0); .op_num(loss_bad, "loss_bad", lower = 0)
  .op_num(max_n, "max_n", lower = 0, open_lower = TRUE)
  .op_num(max_share, "max_share", lower = 0, upper = 1, open_lower = TRUE)
  .op_num(budget, "budget", lower = 0)
  .op_num(max_per_day, "max_per_day", lower = 0)
  .scr_num1(day_quantile, "day_quantile", lower = 0, upper = 1)
  .op_num(min_rate, "min_rate", lower = 0, upper = 1); .op_num(max_rate, "max_rate", lower = 0, upper = 1)
  # every argument belongs to one side
  wrong <- if (ev) c(revenue_good = !is.null(revenue_good), loss_bad = !is.null(loss_bad), max_rate = !is.null(max_rate))
           else c(gain_event = !is.null(gain_event), min_rate = !is.null(min_rate))
  if (any(wrong)) {
    stop(fn, "(): ", lst(names(wrong)[wrong]), " do(es) not apply to side \"", side, "\".", call. = FALSE)
  }
  if (!is.null(budget) && !(cost_select > 0)) stop(fn, "(): `budget` needs a positive `cost_select`.", call. = FALSE)
  if (!is.null(max_per_day) && !daily) stop(fn, "(): `max_per_day` needs a date column (`date`).", call. = FALSE)
  has_value <- isTRUE(inp$meta$has_value)
  econ <- if (ev) !is.null(gain_event) || has_value else !is.null(revenue_good) || !is.null(loss_bad)

  h <- inp$hist
  cells <- .study_cells(h, smp)
  K <- nrow(cells)
  if (!K) stop(fn, "(): sample '", smp, "' has no scored row.", call. = FALSE)
  # selection from the high scores: the event-rich end under higher_is_riskier, the safe end otherwise
  high <- ev == riskier
  o <- if (high) rev(seq_len(K)) else seq_len(K)
  bc <- .study_bcut(cells)
  cut <- if (high) c(bc[rev(seq_len(K - 1L))], -Inf) else c(bc, Inf)
  n <- cells$n[o]; ny <- cells$n_y[o]; e <- cells$e[o]
  N <- sum(n); NY <- sum(ny); E <- sum(e)
  cn <- cumsum(n); cny <- cumsum(ny); ce <- cumsum(e)
  rate <- ifelse(cny > 0, ce / cny, NA_real_)
  neff <- .study_kish(cny, cumsum(cells$w2_y[o]))
  ci <- .study_jeffreys(rate * neff, neff, level)
  R <- if (NY > 0) E / NY else NA_real_
  # marginal rate: pool adjacent violators from the event-poor to the event-rich end of the score
  poor <- if (riskier) seq_len(K) else rev(seq_len(K))
  gp <- .study_pav(cells$e[poor], cells$n_y[poor])
  bn <- as.numeric(rowsum(cells$n_y[poor], gp, reorder = TRUE))
  blk_rate <- ifelse(bn > 0, as.numeric(rowsum(cells$e[poor], gp, reorder = TRUE)) / bn, NA_real_)
  mr <- numeric(K); mr[poor] <- blk_rate[gp]; mr <- mr[o]
  rg <- revenue_good %||% 0; lb <- loss_bad %||% 0
  value <- rep(NA_real_, K); mv <- rep(NA_real_, K)
  if (econ) {
    if (ev && !is.null(gain_event)) {
      value <- gain_event * ce - cost_select * cn
      mv <- gain_event * mr - cost_select
    } else if (ev) {
      # the value of the selected events, and its smoothed amount per case for the margin
      value <- cumsum(cells$ve[o]) - cost_select * cn
      bv <- ifelse(bn > 0, as.numeric(rowsum(cells$ve[poor], gp, reorder = TRUE)) / bn, NA_real_)
      vv <- numeric(K); vv[poor] <- bv[gp]
      mv <- vv[o] - cost_select
    } else {
      value <- rg * (cny - ce) - lb * ce - cost_select * cn
      mv <- rg * (1 - mr) - lb * mr - cost_select
    }
  }
  depth <- cn / N

  # constraints on every candidate cut; a relative slack of 1e-12 keeps a limit met exactly feasible
  tol <- 1e-12
  chk <- list(); lim <- list()
  if (!is.null(max_n)) { chk$max_n <- cn <= max_n * (1 + tol); lim$max_n <- max_n }
  if (!is.null(max_share)) { chk$max_share <- depth <= max_share * (1 + tol); lim$max_share <- max_share }
  if (!is.null(budget)) { chk$budget <- cost_select * cn <= budget * (1 + tol); lim$budget <- budget }
  if (!is.null(min_rate)) { chk$min_rate <- !is.na(rate) & rate >= min_rate * (1 - tol); lim$min_rate <- min_rate }
  if (!is.null(max_rate)) { chk$max_rate <- !is.na(rate) & rate <= max_rate * (1 + tol); lim$max_rate <- max_rate }
  dv <- if (daily) .operating_days(h, smp, cells, high, day_quantile) else NULL
  if (!is.null(max_per_day)) {
    # the daily quantile never falls with the depth: bisection for the deepest cut within capacity
    jd <- 0L; a <- 1L; b <- K
    while (a <= b) {
      m <- (a + b) %/% 2L
      if (dv$q(m) <= max_per_day * (1 + tol)) { jd <- m; a <- m + 1L } else b <- m - 1L
    }
    chk$max_per_day <- seq_len(K) <= jd; lim$max_per_day <- max_per_day
  }
  feasible <- Reduce(`&`, chk, rep(TRUE, K))

  # optimum of a set of feasible rows: the best value (the smallest depth on
  # a tie), else the deepest row; NA when no row is feasible
  best <- function(ok) {
    f <- which(ok)
    if (!length(f)) return(NA_integer_)
    if (!econ) return(max(f))
    v <- value[f]
    f[which(v == max(v))[1]]
  }
  # `b` improves on `a`: worth more, or deeper when there are no economics
  better <- function(b, a) !is.na(b) && if (econ) value[b] > value[a] else b > a
  note <- NA_character_
  jo <- best(feasible)
  if (is.na(jo)) {
    first <- names(chk)[!vapply(chk, `[`, logical(1), 1L)]
    note <- sprintf("no cut meets the constraints: the first score value already breaks %s", lst(first))
    warning(fn, "(): ", note, "; the optimum is NA.", call. = FALSE)
  }
  binding <- NA_character_; shadow <- NA_real_; next_rate <- NA_real_
  if (!is.na(jo)) {
    # a constraint binds when dropping it alone improves the optimum
    bind <- names(chk)[vapply(seq_along(chk), function(i) better(best(Reduce(`&`, chk[-i], rep(TRUE, K))), jo),
                              logical(1))]
    if (!length(bind) && length(chk)) {
      # constraints that stop the optimum at the same row bind jointly: none
      # improves it alone, so name those the nearest better row breaks
      imp <- if (econ) which(value > value[jo]) else which(seq_len(K) > jo)
      if (length(imp)) bind <- names(chk)[!vapply(chk, `[`, logical(1), imp[1])]
    }
    binding <- if (length(bind)) paste(bind, collapse = ", ") else if (econ) "value" else "end of the curve"
    if (jo < K) { shadow <- mv[jo + 1L]; next_rate <- mr[jo + 1L] }
    if (econ && value[jo] < 0) note <- "every feasible cut loses value; selecting nothing is worth 0"
  }

  # thinned curve: even in depth, plus the rows the optimum and the constraints pick
  keep <- if (K <= n_points) seq_len(K) else {
    c(.op_nearest(depth, seq_len(n_points) / n_points), .op_nearest(depth, c(0.01, 0.05, 0.10, 0.20, 0.50)),
      1L, K, jo, vapply(chk, function(v) if (any(v)) max(which(v)) else 1L, integer(1)))
  }
  keep <- sort(unique(keep[!is.na(keep)]))
  curve <- data.table::data.table(
    cut = cut[keep], depth = depth[keep], n_sel = cn[keep], events_sel = ce[keep], rate_sel = rate[keep],
    rate_lo = ci$lo[keep], rate_hi = ci$hi[keep], capture = if (E > 0) ce[keep] / E else NA_real_,
    lift = if (isTRUE(R > 0)) rate[keep] / R else NA_real_, marginal_rate = mr[keep], cost = cost_select * cn[keep],
    value = value[keep],
    feasible = feasible[keep])
  if (daily) {
    M <- dv$vol(keep)
    curve[, `:=`(day_q = dv$q_of(M), pct_days_over = if (is.null(max_per_day) || !dv$D) NA_real_ else
      colMeans(M > max_per_day * (1 + tol)))]
  }
  opt <- if (is.na(jo)) curve[NA_integer_] else curve[match(jo, keep)]
  opt[, `:=`(binding = binding, shadow_price = shadow, next_rate = next_rate)]
  at <- function(cn_) if (is.na(jo)) NA_real_ else switch(cn_, max_n = cn[jo], max_share = depth[jo],
    budget = cost_select * cn[jo], min_rate = rate[jo], max_rate = rate[jo], max_per_day = dv$q(jo))
  bnd <- if (is.na(binding)) character() else strsplit(binding, ", ", fixed = TRUE)[[1]]
  nl <- as.character(names(lim))
  cons <- data.table::data.table(constraint = nl, limit = as.double(unlist(lim)),
                                 at_optimum = vapply(nl, at, numeric(1), USE.NAMES = FALSE), binding = nl %in% bnd)
  structure(list(
    curve = curve, optimum = opt, constraints = cons, side = side, objective = inp$objective,
    direction = inp$direction, target = inp$target, sample = smp, level = level, economics = econ,
    gain_event = gain_event, cost_select = cost_select, revenue_good = revenue_good, loss_bad = loss_bad,
    value_column = ev && is.null(gain_event) && has_value, select_high = high, day_quantile = day_quantile,
    n_days = if (daily) dv$D else NA_integer_, n_cells = K, quantized = isTRUE(inp$meta$quantized),
    weighted = isTRUE(inp$meta$weighted), message = note, call = call),
    class = c("scr_operating", "list"))
}

#' Rows whose depth is nearest each target (the lower one on a tie)
#' @keywords internal
#' @noRd
.op_nearest <- function(depth, targets) {
  K <- length(depth)
  i <- findInterval(targets, depth)
  lo <- pmax(i, 1L); hi <- pmin(i + 1L, K)
  ifelse(abs(depth[hi] - targets) < abs(depth[lo] - targets), hi, lo)
}

#' Selected volume per day at any depth, from the counts per day and score value
#'
#' `vol(js)` returns a D x length(js) matrix of the volume selected on each
#' day by the first `js` cells in selection order; `q(j)` the
#' `day_quantile` of one column. A rolling join on the cumulative volume of
#' every day gives each column in O(D log K).
#' @keywords internal
#' @noRd
.operating_days <- function(h, smp, cells, high, day_quantile) {
  hd <- h[h[["sample"]] == smp & !is.na(h[["group"]])]
  K <- nrow(cells)
  key <- if ("g" %in% names(cells)) "g" else "s"
  asc <- match(hd[[key]], cells[[key]])
  dd <- data.table::data.table(day = hd[["group"]], pos = if (high) K - asc + 1L else asc, n = hd$n)
  data.table::setkeyv(dd, c("day", "pos"))
  dd[, cum := cumsum(n), by = "day"]
  days <- sort(unique(dd$day))
  D <- length(days)
  vol <- function(js) {
    if (!D) return(matrix(numeric(), 0L, length(js)))
    q <- data.table::data.table(day = rep(days, times = length(js)), pos = rep(as.integer(js), each = D))
    v <- dd[q, on = c("day", "pos"), roll = TRUE, x.cum]
    v[is.na(v)] <- 0
    matrix(v, nrow = D)
  }
  q_of <- function(M) if (!D) rep(NA_real_, ncol(M)) else
    apply(M, 2L, stats::quantile, probs = day_quantile, names = FALSE, type = 7)
  list(vol = vol, q_of = q_of, q = function(j) q_of(vol(j)), D = D)
}

#' @export
print.scr_operating <- function(x, ...) {
  ends <- if (x$select_high) "the high scores, score >= cut" else "the low scores, score < cut"
  cat(sprintf("<scr_operating> target \"%s\" | objective %s | %s | side %s (from %s)\n", x$target, x$objective,
              x$direction, x$side, ends))
  eco <- if (!x$economics) "none" else if (identical(x$side, "event")) {
    if (isTRUE(x$value_column)) sprintf("value column, cost %s", format(x$cost_select))
    else sprintf("gain %s per event, cost %s", format(x$gain_event), format(x$cost_select))
  } else sprintf("revenue %s, loss %s, cost %s", format(x$revenue_good %||% 0), format(x$loss_bad %||% 0),
                 format(x$cost_select))
  cat(sprintf("  sample '%s' | %d score values%s | economics: %s\n", x$sample, x$n_cells,
              if (is.na(x$n_days)) "" else sprintf(" | %d dates", x$n_days), eco))
  if (nrow(x$constraints)) {
    cat(sprintf("  constraints: %s\n", paste(sprintf("%s %s", x$constraints$constraint,
                                                     vapply(signif(x$constraints$limit, 6), format, character(1),
                                                            big.mark = ",", scientific = FALSE, trim = TRUE)),
                             collapse = ", ")))
  }
  o <- x$optimum
  if (is.na(o$depth)) {
    cat(sprintf("\nOptimum: NA (%s)\n", x$message))
  } else {
    cat(sprintf("\nOptimum: cut %s | depth %.1f%% (n %s) | rate %s [%s, %s] | capture %s | value %s\n",
                format(o$cut, digits = 7), 100 * o$depth, .study_n(o$n_sel), .claims_pct(o$rate_sel),
                .claims_pct(o$rate_lo), .claims_pct(o$rate_hi), .claims_pct(o$capture),
                if (is.na(o$value)) "-" else format(round(o$value, 2), big.mark = ",")))
    cat(sprintf("  binding: %s | shadow price %s per additional case | next marginal rate %s\n", o$binding,
                if (is.na(o$shadow_price)) "-" else format(signif(o$shadow_price, 4)), .claims_pct(o$next_rate)))
    if (!is.na(x$message)) cat(sprintf("  note: %s\n", x$message))
  }
  cv <- x$curve
  rows <- unique(c(.op_nearest(cv$depth, c(0.01, 0.05, 0.10, 0.20, 0.50, 1)),
                   if (!is.na(o$depth)) match(o$depth, cv$depth)))
  rows <- sort(rows[!is.na(rows)])
  daily <- "day_q" %in% names(cv)
  cat("\nCurve (selected rows)\n")
  cat(sprintf("  %7s %12s %10s %-24s %8s %6s %9s %14s%s %8s\n", "depth", "cut", "n_sel", "rate [lo, hi]", "capture",
              "lift", "marginal", "value", if (daily) "     day_q" else "", "feasible"))
  for (i in rows) {
    r <- cv[i]
    cat(sprintf("  %6.1f%% %12s %10s %-24s %8s %6s %9s %14s%s %8s%s\n", 100 * r$depth, format(r$cut, digits = 7),
                .study_n(r$n_sel), sprintf("%s [%s, %s]", .claims_pct(r$rate_sel), .claims_pct(r$rate_lo),
                                           .claims_pct(r$rate_hi)),
                .claims_pct(r$capture), .study_f(r$lift, "%.2f"), .claims_pct(r$marginal_rate),
                if (is.na(r$value)) "-" else format(round(r$value, 2), big.mark = ","),
                if (daily) sprintf(" %9s", .study_f(r$day_q, "%.1f")) else "", if (r$feasible) "yes" else "no",
                if (!is.na(o$depth) && identical(r$depth, o$depth)) "  <- optimum" else ""))
  }
  invisible(x)
}

#' @rdname scr_export
#' @export
scr_export.scr_operating <- function(x, dir, stamp = TRUE, ...) {
  .study_dots(list(...), "scr_export")
  .need_openxlsx()
  out_dir <- .export_dir(dir, stamp)
  tag <- .file_tag(x$target)
  settings <- .kv_table(list(target = x$target, objective = x$objective, direction = x$direction, side = x$side,
                             sample = x$sample, level = x$level, economics = x$economics, gain_event = x$gain_event,
                             cost_select = x$cost_select, revenue_good = x$revenue_good, loss_bad = x$loss_bad,
                             day_quantile = x$day_quantile, n_days = x$n_days, n_cells = x$n_cells,
                             message = x$message))
  sheets <- lapply(list(Curve = x$curve, Optimum = x$optimum, Constraints = x$constraints, Settings = settings),
                   .study_sheet)
  files <- list(xlsx = .scr_write_xlsx(sheets, file.path(out_dir, sprintf("operating_%s.xlsx", tag))))
  for (f in files) msg("  %s", f)
  x$files <- files
  invisible(x)
}

# data.table column names used without quotes in this file
utils::globalVariables(c("cum", "x.cum", "day_q", "pct_days_over", "binding", "shadow_price", "next_rate"))
