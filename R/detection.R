# ============================================================================ #
# detection.R - time to detection of fraud episodes
# ============================================================================ #
# Only the entities with an event are kept and sorted once by entity and
# time. Within an episode the running maximum of the score rises at a few
# "record" rows, and a threshold is first reached at the record that crosses
# it: every threshold is answered from the record rows, with no pass over
# the rows per threshold. The alert share of all rows comes from one count
# table of the score.
# ============================================================================ #

#' Time to detection of fraud episodes
#'
#' Reads a score on transactions grouped by entity (a card, an account, a
#' device) and measures, for each alert threshold, how many fraud episodes
#' the score detects, how many fraudulent transactions pass before the
#' first alert, how long that takes and how much of the loss is prevented.
#'
#' @section Episodes:
#'
#' An episode is an entity with at least one event row. It starts at its
#' first event row in time order (rows at the same time keep their order in
#' `x`), and the rows of the entity before that start are ignored. At a
#' threshold \eqn{t}, the episode is detected at the first row from its
#' start whose score is on the alert side: `score >= t` under
#' `higher_is_riskier`, `score < t` under `higher_is_safer`. The detecting
#' row may be any row of the entity, an event or not.
#'
#' @section Table:
#'
#' One row per threshold, the strictest (fewest alerts) first:
#' \itemize{
#'   \item `alert_share`: the share of **all** rows of `x` on the alert side
#'     of the threshold, the workload it costs;
#'   \item `episodes`, `detected` and `pct_detected`;
#'   \item `median_events_before` and `mean_events_before`: the event rows
#'     of the episode before the detecting row, over the detected episodes
#'     (0 when the first event row is alerted);
#'   \item `median_time`: the median time from the start of the episode to
#'     the detecting row, in the unit of `time` (days for a `Date`, seconds
#'     for a `POSIXct`);
#'   \item with `amount`: `loss_before`, the amount of the event rows before
#'     the detecting row, plus the whole amount of the episodes never
#'     detected; `loss_total`, the amount of every event row; and
#'     `pct_loss_prevented` = `1 - loss_before / loss_total`. The detecting
#'     row and the rows after it count as prevented.
#' }
#' `alert_shares` are turned into thresholds on all rows, tie-safe as in
#' [scr_bands()]: the boundary between two distinct scores nearest to the
#' share (`target_share` keeps the share asked for, `alert_share` the one
#' realized).
#'
#' Rows with a missing entity or time, or a missing or infinite score, are
#' left out (`n_dropped`); a missing outcome is not an event, and a missing
#' amount counts as 0.
#'
#' @param x A `data.frame` with one row per transaction.
#' @param entity Name of the column that identifies the entity.
#' @param time Name of the time column: a number, a `Date` or a `POSIXct`.
#' @param score,y Column names of the score and of the 0/1 outcome (1 = a
#'   fraudulent row).
#' @param amount Optional column of the amount of each row, for the loss.
#' @param thresholds Score thresholds of the alerts.
#' @param alert_shares Shares of all rows to alert, in (0, 1\], turned into
#'   thresholds; may be given with or instead of `thresholds`.
#' @param direction `"higher_is_riskier"` (default) or `"higher_is_safer"`.
#' @param max_cells Largest number of distinct score values kept exactly in
#'   the count table of the alert shares.
#' @param ... Not used; an unknown argument is an error.
#'
#' @return An object of class `c("scr_detection", "list")`:
#'   \describe{
#'     \item{`table`}{One row per threshold: `threshold`, `target_share`,
#'       `alert_share`, `episodes`, `detected`, `pct_detected`,
#'       `median_events_before`, `mean_events_before`, `median_time`,
#'       `loss_before`, `loss_total` and `pct_loss_prevented`.}
#'     \item{`episodes`}{One row per episode: `entity`, `start`, `rows` (from
#'       the start), `events`, `amount` (of its event rows), `peak_score`
#'       (the highest score from the start; the lowest under
#'       `higher_is_safer`) and, for the strictest threshold, `detected`,
#'       `detect_time`, `time_to_detect`, `events_before` and
#'       `loss_before` (`NA` when the episode is not detected).}
#'     \item{`direction`, `entity`, `time`, `score`, `target`, `amount`,
#'       `n_rows`, `n_dropped`, `n_episode_rows`, `call`}{The settings and
#'       the row counts; `n_episode_rows` is the number of rows of the
#'       episodes, from their starts.}
#'   }
#'
#' @seealso [scr_overlap()] for the rules against the score,
#'   [scr_operating()] for the threshold under a review capacity.
#' @family score-studies
#' @examples
#' local({
#'   set.seed(1)
#'   n <- 30000
#'   d <- data.frame(card = sample(3000, n, TRUE), ts = sort(runif(n, 0, 30)))
#'   # 60 cards are compromised from some day on; their later rows are fraud
#'   hit <- sample(3000, 60)
#'   since <- runif(3000, 5, 25)
#'   d$y <- as.integer(d$card %in% hit & d$ts >= since[d$card])
#'   d$score <- round(100 * plogis(rnorm(n, -2 + 2 * d$y)))
#'   d$amount <- round(rexp(n, 1 / 60), 2)
#'   dt <- scr_detection(d, entity = "card", time = "ts", amount = "amount",
#'                       alert_shares = c(0.01, 0.03, 0.10))
#'   print(dt)
#'   head(dt$episodes)
#' })
#' @export
scr_detection <- function(x, ...) UseMethod("scr_detection")

#' @rdname scr_detection
#' @export
scr_detection.data.frame <- function(x, entity, time, score = "score", y = "y", amount = NULL, thresholds = NULL,
                                     alert_shares = NULL, direction = "higher_is_riskier", max_cells = 1e5, ...) {
  fn <- "scr_detection"
  .study_dots(list(...), fn)
  if (missing(entity) || missing(time)) stop(fn, "(): `entity` and `time` are needed: two column names.", call. = FALSE)
  for (nm in c("entity", "time", "score", "y", "amount")) {
    v <- get(nm)
    if (!is.null(v)) .study_chr1(v, nm, fn)
  }
  if (!is.character(direction) || length(direction) != 1L || !direction %in% c("higher_is_safer", "higher_is_riskier")) {
    stop(fn, "(): `direction` must be \"higher_is_riskier\" or \"higher_is_safer\".", call. = FALSE)
  }
  miss <- setdiff(c(entity, time, score, y, amount), names(x))
  if (length(miss)) stop(fn, "(): column(s) ", lst(miss), " not in `x`.", call. = FALSE)
  for (nm in c(score, amount)) {
    if (!is.numeric(x[[nm]])) stop(fn, "(): column '", nm, "' must be numeric.", call. = FALSE)
  }
  tm <- x[[time]]
  if (!is.numeric(tm) && !inherits(tm, c("Date", "POSIXct"))) {
    stop(fn, "(): column '", time, "' must be a number, a Date or a POSIXct.", call. = FALSE)
  }
  if (is.null(thresholds) && is.null(alert_shares)) stop(fn, "(): give `thresholds`, `alert_shares` or both.", call. = FALSE)
  if (!is.null(thresholds) && (!is.numeric(thresholds) || !length(thresholds) || any(!is.finite(thresholds)))) {
    stop(fn, "(): `thresholds` must be finite numbers.", call. = FALSE)
  }
  if (!is.null(alert_shares) && (!is.numeric(alert_shares) || !length(alert_shares) || anyNA(alert_shares) ||
                                 any(alert_shares <= 0 | alert_shares > 1))) {
    stop(fn, "(): `alert_shares` must be shares in (0, 1].", call. = FALSE)
  }
  ent <- x[[entity]]
  sc <- as.double(x[[score]])
  yy <- .scr_y01(x[[y]], fn)
  ev <- !is.na(yy) & yy == 1L
  ok <- !is.na(ent) & !is.na(tm) & is.finite(sc)
  amt <- NULL
  if (!is.null(amount)) {
    amt <- as.double(x[[amount]]); amt[is.na(amt)] <- 0
    if (any(!is.finite(amt)) || any(amt < 0)) stop(fn, "(): the amounts must be finite and non-negative.", call. = FALSE)
  }
  if (!any(ok)) stop(fn, "(): no row has an entity, a time and a score.", call. = FALSE)
  side <- .study_side(direction)
  high <- side == "high"

  # one count table of the score of all rows: the cuts of the shares and the alert share of every threshold
  all_in <- all(ok)
  h <- .study_hist(if (all_in) sc else sc[ok], integer(sum(ok)), max_cells = max_cells, breaks = thresholds, fn = fn)
  cells <- .study_collapse(h)
  thr <- data.table::data.table(threshold = c(as.double(thresholds), .study_share_cut(cells, sort(unique(alert_shares)), side)),
                                target_share = c(rep(NA_real_, length(thresholds)), sort(unique(alert_shares))))
  thr <- thr[!duplicated(thr$threshold)]
  # the strictest threshold first: the fewest alerts
  thr <- thr[order(if (high) -thr$threshold else thr$threshold)]
  thr[, alert_share := .study_alert_share(cells, threshold, side)]

  ep <- .detection_episodes(ent, tm, sc, ev, amt, ok, thr$threshold, high)
  tab <- data.table::data.table(thr, episodes = ep$E, detected = ep$by_thr$detected,
                                pct_detected = if (ep$E > 0) ep$by_thr$detected / ep$E else NA_real_,
                                median_events_before = ep$by_thr$median_events, mean_events_before = ep$by_thr$mean_events,
                                median_time = ep$by_thr$median_time)
  if (!is.null(amount)) {
    lb <- ep$by_thr$loss_det + (ep$loss_total - ep$by_thr$amt_det)
    tab[, `:=`(loss_before = lb, loss_total = ep$loss_total,
               pct_loss_prevented = if (ep$loss_total > 0) 1 - lb / ep$loss_total else NA_real_)]
  } else {
    tab[, `:=`(loss_before = NA_real_, loss_total = NA_real_, pct_loss_prevented = NA_real_)]
    if (!is.null(ep$episodes)) ep$episodes[, `:=`(amount = NA_real_, loss_before = NA_real_)]
  }
  structure(list(table = tab[], episodes = ep$episodes, direction = direction, entity = entity, time = time,
                 score = score, target = y, amount = amount, n_rows = sum(ok), n_dropped = sum(!ok),
                 n_episode_rows = ep$n_rows, call = sys.call()),
            class = c("scr_detection", "list"))
}

#' Episodes and their first detection at every threshold
#'
#' `thr` is ordered from the strictest threshold. Returns `E` episodes,
#' `by_thr` (per threshold: `detected`, the medians and the mean over the
#' detected episodes, `loss_det` the loss before detection and `amt_det` the
#' whole amount of the detected episodes), `loss_total`, the per-episode
#' table and the number of episode rows.
#' @keywords internal
#' @noRd
.detection_episodes <- function(ent, tm, sc, ev, amt, ok, thr, high) {
  TT <- length(thr)
  zero <- data.table::data.table(detected = rep(0, TT), median_events = NA_real_, mean_events = NA_real_,
                                 median_time = NA_real_, loss_det = 0, amt_det = 0)
  # entities with an event: a filter, then one sort of that subset
  hit <- unique(ent[ok & ev])
  if (!length(hit)) return(list(E = 0L, by_thr = zero, loss_total = 0, episodes = NULL, n_rows = 0L))
  keep <- which(ok & ent %in% hit)
  d <- data.table::data.table(ent = ent[keep], tm = tm[keep], s = if (high) sc[keep] else -sc[keep],
                              ev = as.double(ev[keep]), la = if (is.null(amt)) 0 else amt[keep] * ev[keep])
  data.table::setorderv(d, c("ent", "tm"))
  gid <- data.table::rleid(d$ent)
  n <- nrow(d)
  first <- c(TRUE, gid[-1L] != gid[-n])
  # events up to each row within its entity: a running sum less its value at the start of the entity
  cs <- cumsum(d$ev)
  cev <- cs - (cs - d$ev)[first][gid]
  # the episode starts at the first event row
  d <- d[cev > 0]; gid <- gid[cev > 0]; cev <- cev[cev > 0]
  n <- nrow(d)
  first <- c(TRUE, gid[-1L] != gid[-n])
  E <- gid[n]
  tn <- as.numeric(d$tm)
  since <- tn - tn[first][gid]
  cl <- cumsum(d$la)
  loss_before <- cl - d$la - (cl - d$la)[first][gid]
  ev_before <- cev - d$ev
  ep_amt <- as.double(rowsum(d$la, gid, reorder = TRUE))

  # running maximum of the score within the episode, on dense ranks: an offset
  # per episode makes the global running maximum restart at every episode
  u <- sort(unique(d$s))
  rk <- match(d$s, u) + gid * (length(u) + 1)
  cm <- cummax(rk)
  rec <- which(c(TRUE, cm[-1L] > cm[-n]))
  m1 <- d$s[rec]
  m0 <- c(-Inf, m1[-length(m1)])
  new_ep <- first[rec]
  # thresholds a record is the first to reach: those above the previous record and at or below this one
  # (from the low scores the alert is strict, `score < t`)
  ts <- rev(if (high) thr else -thr)
  cnt_le <- function(v) if (high) findInterval(v, ts) else findInterval(v, ts, left.open = TRUE)
  hi <- cnt_le(m1)
  lo <- ifelse(new_ep, 0L, cnt_le(m0))
  k <- pmax(hi - lo, 0L)
  pr <- rep(seq_along(rec), k)
  pj <- sequence(k, from = lo + 1L)
  # index of the strictest-first threshold: the sorted thresholds end with the strictest
  pj <- TT + 1L - pj
  r <- rec[pr]
  pairs <- data.table::data.table(j = pj, evb = ev_before[r], dt = since[r], lb = loss_before[r], am = ep_amt[gid[r]])
  agg <- pairs[, list(detected = as.double(.N), median_events = as.double(stats::median(evb)), mean_events = mean(evb),
                      median_time = as.double(stats::median(dt)), loss_det = sum(lb), amt_det = sum(am)), keyby = "j"]
  by_thr <- data.table::copy(zero)
  for (cn in names(zero)) data.table::set(by_thr, i = agg$j, j = cn, value = agg[[cn]])

  # one row per episode, with its first detection at the strictest threshold
  last <- c(first[-1L], TRUE)
  peak <- u[cm[last] - gid[last] * (length(u) + 1)]
  strict <- pj == 1L
  det <- rep(NA_integer_, E); det[gid[r[strict]]] <- r[strict]
  episodes <- data.table::data.table(
    entity = d$ent[first], start = d$tm[first], rows = tabulate(gid, E), events = as.double(rowsum(d$ev, gid, reorder = TRUE)),
    amount = ep_amt, peak_score = if (high) peak else -peak, detected = !is.na(det), detect_time = d$tm[det],
    time_to_detect = since[det], events_before = ev_before[det], loss_before = loss_before[det])
  list(E = E, by_thr = by_thr, loss_total = sum(ep_amt), episodes = episodes, n_rows = n)
}

#' @export
print.scr_detection <- function(x, ...) {
  t <- x$table
  cat(sprintf("<scr_detection> entity \"%s\" | time \"%s\" | score \"%s\" (%s) | outcome \"%s\"\n", x$entity, x$time,
              x$score, x$direction, x$target))
  E <- if (nrow(t)) t$episodes[1] else 0
  cat(sprintf("  %s rows%s | %s episodes over %s rows from their starts%s\n", .study_n(x$n_rows),
              if (x$n_dropped > 0) sprintf(" (%s left out)", .study_n(x$n_dropped)) else "", .study_n(E),
              .study_n(x$n_episode_rows),
              if (is.null(x$amount)) "" else sprintf(" | loss %s", format(round(t$loss_total[1], 2), big.mark = ","))))
  cat(sprintf("\n  %12s %11s %9s %8s %15s %13s %12s %15s\n", "threshold", "alert share", "detected", "pct",
              "events before", "(mean)", "median time", "loss prevented"))
  for (i in seq_len(nrow(t))) {
    cat(sprintf("  %12s %11s %9s %8s %15s %13s %12s %15s\n", format(t$threshold[i], digits = 7),
                .study_f(100 * t$alert_share[i], "%.2f%%"), .study_n(t$detected[i]),
                .study_f(100 * t$pct_detected[i], "%.1f%%"), .study_f(t$median_events_before[i], "%.1f"),
                .study_f(t$mean_events_before[i], "%.2f"),
                if (is.na(t$median_time[i])) "-" else format(signif(t$median_time[i], 4)),
                .study_f(100 * t$pct_loss_prevented[i], "%.1f%%")))
  }
  invisible(x)
}

#' @rdname scr_export
#' @export
scr_export.scr_detection <- function(x, dir, stamp = TRUE, ...) {
  .study_dots(list(...), "scr_export")
  .need_openxlsx()
  out_dir <- .export_dir(dir, stamp)
  tag <- .file_tag(x$target)
  settings <- .kv_table(list(entity = x$entity, time = x$time, score = x$score, target = x$target,
                             amount = x$amount, direction = x$direction, n_rows = x$n_rows,
                             n_dropped = x$n_dropped, n_episode_rows = x$n_episode_rows))
  ep <- x$episodes
  # a date-time is written as text: the workbook keeps the value as read
  if (!is.null(ep)) {
    ep <- data.table::copy(ep)
    for (cn in c("start", "detect_time")) if (inherits(ep[[cn]], "POSIXct")) data.table::set(ep, j = cn, value = format(ep[[cn]]))
  }
  sheets <- list(Thresholds = x$table, Episodes = ep, Settings = settings)
  sheets <- lapply(sheets[!vapply(sheets, is.null, logical(1))], .study_sheet)
  files <- list(xlsx = .scr_write_xlsx(sheets, file.path(out_dir, sprintf("detection_%s.xlsx", tag))))
  for (f in files) msg("  %s", f)
  x$files <- files
  invisible(x)
}

# data.table column names used without quotes in this file
utils::globalVariables(c("alert_share", "threshold", "evb", "lb", "am", "loss_before", "loss_total",
                         "pct_loss_prevented", "amount"))
