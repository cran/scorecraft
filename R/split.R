# ============================================================================ #
# split.R - Stage 0: typing, event orientation and the train/hold-out split
# ============================================================================ #

#' Stage 0: type the data and split train and hold-out
#'
#' First stage of the pipeline, callable on its own. Converts the target to
#' 0/1 (resolving `event_level`), types the candidates (numerics become
#' `double`, everything else becomes `character`) and splits train and
#' hold-out **before** any supervised fit.
#'
#' The split prefers out-of-time by `date_col`: it is the only one that
#' tests generalisation to a future period. The cut is made on the
#' **distinct** date values, not by row quantile: it picks the smallest set
#' of most recent periods that already reaches `ratio` of the population.
#' Without a date column, or with a single period, it falls back to random
#' stratified by the target. The date column is never a candidate: it is the
#' key of the split and leaves the contest. A text date column is read as an
#' ISO date (`YYYY-MM-DD`, `YYYY/MM/DD`, `YYYY-MM`) or as all-digit periods
#' (`YYYYMM`); rows with a missing date belong to no period and are left out
#' of both train and hold-out, with a warning in the log.
#'
#' @section Event orientation:
#'
#' `event_level` changes **what is modelled**. Passing `0` makes class 0 the
#' event: the sign of every WOE flips, the emitted SQL changes, the points
#' change. For a text target, the second level in alphabetical order is the
#' event by default, and the choice is always reported. Not to be confused
#' with `config$objective`, which only orients the reading and the points
#' scale.
#'
#' @param data A `data.frame` or `data.table` with the target and the candidates.
#' @param target Name of the target column. Binary: 0/1, logical, or a
#'   two-level factor/character.
#' @param date_col Date column of the out-of-time cut. `NULL` uses a random
#'   stratified split.
#' @param ratio Target hold-out fraction.
#' @param seed Seed of the random split. `NULL` draws from the session's
#'   random stream (reproducible only through a `set.seed()` of your own;
#'   [scr_select()] passes the `seed` of [scr_config()]). A seed is applied
#'   locally: the session's random stream is restored on exit.
#' @param event_level Which target value counts as the event. `NULL` uses
#'   the convention (`1`, or the second alphabetical level).
#' @param drop Columns that are never candidates (identifiers, sibling
#'   targets, free text). They stay in the funnel as `00.config`.
#' @param copy If `TRUE` (default), works on a copy of `data`. `FALSE`
#'   modifies a `data.table` by reference (target and typing), saving
#'   memory; a `data.frame` is always converted, hence copied.
#'
#' @return An `scr_split` object with `data` (typed), `target`, `train_idx`,
#'   `holdout_idx`, `method`, `cutoff`, `date_col` and `cols` (`features`,
#'   `var_num`, `var_cat`, `dropped`, `event`).
#'
#' @family stages
#' @examples
#' sp <- scr_split(scr_demo, "default", date_col = "ref_date", drop = "id")
#' sp
#' length(sp$train_idx); length(sp$holdout_idx)
#' @export
scr_split <- function(data, target, date_col = NULL, ratio = 0.30, seed = NULL,
                      event_level = NULL, drop = character(), copy = TRUE) {
  if (!is.data.frame(data)) stop("`data` must be a data.frame or data.table.", call. = FALSE)
  if (!is.character(target) || length(target) != 1L) {
    stop("`target` must be the name of a single column.", call. = FALSE)
  }
  .scr_num1(ratio, "ratio", lower = 0, upper = 1, open_lower = TRUE)

  dt <- if (data.table::is.data.table(data)) {
    if (isTRUE(copy)) data.table::copy(data) else data
  } else {
    data.table::as.data.table(data)
  }
  if (!is.null(date_col) && !date_col %in% names(dt)) {
    msg("  '%s' does not exist in the table - falling back to a random stratified split.", date_col)
    date_col <- NULL
  }

  # The split comes BEFORE typing, on purpose: after it the date column would
  # be text, and as.numeric() of text returns NA - the out-of-time cut would
  # silently degenerate into a random split.
  y_raw <- .target_as_int(dt[[target]], target, event_level)
  data.table::set(dt, j = target, value = y_raw$y)
  split <- split_train_holdout(dt, target, date_col, ratio, seed)

  drop_cols <- drop
  if (!is.null(date_col) && !date_col %in% drop_cols) drop_cols <- union(drop_cols, date_col)
  cols <- prepare_columns(dt, target, drop_cols)
  cols$event <- y_raw$event

  structure(list(data = dt, target = target, train_idx = split$train_idx,
                 holdout_idx = split$holdout_idx, method = split$method,
                 cutoff = split$cutoff, date_col = date_col, cols = cols,
                 ratio = ratio, seed = seed),
            class = c("scr_split", "list"))
}

#' @export
print.scr_split <- function(x, ...) {
  cat(sprintf("<scr_split> target \"%s\" | %s rows: train %s, hold-out %s\n",
              x$target, n_fmt(nrow(x$data)), n_fmt(length(x$train_idx)), n_fmt(length(x$holdout_idx))))
  cat(sprintf("  method: %s%s\n", x$method,
              if (!is.na(x$cutoff)) sprintf(" (hold-out from %s)", x$cutoff) else ""))
  cat(sprintf("  candidates: %d (%d numeric, %d categorical) | dropped: %d\n",
              length(x$cols$features), length(x$cols$var_num), length(x$cols$var_cat), length(x$cols$dropped)))
  cat(sprintf("  event: class '%s'%s\n", x$cols$event$label,
              if (isTRUE(x$cols$event$inverted)) " (TARGET INVERTED by event_level)" else ""))
  invisible(x)
}

#' Convert the target to a 0/1 integer, resolving event_level
#' @keywords internal
#' @noRd
.target_as_int <- function(y, target, event_level = NULL) {
  if (is.null(y)) stop("Target '", target, "' is missing from the table.", call. = FALSE)
  if (anyNA(y)) stop("Target '", target, "' has NA - a missing target is not accepted.", call. = FALSE)
  raw <- y
  if (is.logical(y)) y <- as.integer(y)
  if (is.numeric(y) && is.object(y)) y <- as.double(y)   # e.g. bit64::integer64

  if (is.factor(y) || is.character(y)) {
    lv <- sort(unique(as.character(y)))
    if (length(lv) != 2L) stop("Target '", target, "' is not binary (", length(lv), " levels).", call. = FALSE)
    ev <- if (is.null(event_level)) lv[2] else as.character(event_level)
    if (!ev %in% lv) {
      stop("event_level '", ev, "' does not exist in target '", target, "'. Levels: ",
           paste(lv, collapse = ", "), ".", call. = FALSE)
    }
    msg("  text target: event = '%s', non-event = '%s'%s", ev, setdiff(lv, ev),
        if (is.null(event_level)) " (alphabetical default)" else " (event_level)")
    return(list(y = as.integer(as.character(y) == ev),
                event = list(label = ev, inverted = FALSE, input_class = class(raw)[1])))
  }
  # checked BEFORE as.integer(): truncation would silently turn 0.5 into 0
  if (!all(y %in% c(0, 1))) {
    stop("Target '", target, "' must be 0/1 (found: ", lst(sort(unique(y))), ").", call. = FALSE)
  }
  y <- as.integer(y)
  label <- "1"; inverted <- FALSE
  if (!is.null(event_level)) {
    el <- suppressWarnings(as.integer(event_level))
    if (is.na(el) || !el %in% c(0L, 1L)) {
      stop("event_level for a 0/1 target must be 0 or 1 (got '", event_level, "').", call. = FALSE)
    }
    if (el == 0L) {
      msg("  event_level = 0: target INVERTED - class 0 becomes the event.")
      y <- 1L - y; label <- "0"; inverted <- TRUE
    }
  }
  list(y = y, event = list(label = label, inverted = inverted, input_class = class(raw)[1]))
}

#' Separate candidates from dropped columns and normalise types (by reference)
#' @keywords internal
#' @noRd
prepare_columns <- function(dt, target, drop_cols = character()) {
  dropped  <- intersect(drop_cols, names(dt))
  features <- setdiff(names(dt), c(target, dropped))
  var_num <- var_cat <- character()
  for (f in features) {
    x <- dt[[f]]
    if (is.numeric(x)) {
      var_num <- c(var_num, f)
      # is.object(): a classed double such as bit64::integer64 (what DBI
      # drivers return for BIGINT) is `double` in storage but its bits are
      # not the value; as.double() dispatches to the method that converts
      if (!is.double(x) || is.object(x)) data.table::set(dt, j = f, value = as.double(x))
    } else {
      var_cat <- c(var_cat, f)
      if (!is.character(x)) data.table::set(dt, j = f, value = as.character(x))
    }
  }
  list(target = target, features = features, var_num = var_num, var_cat = var_cat, dropped = dropped)
}

#' Train/hold-out split: out-of-time by whole periods, or stratified
#' @keywords internal
#' @noRd
split_train_holdout <- function(dt, target, date_col = NULL, ratio = 0.30, seed = NULL) {
  n <- nrow(dt)
  if (!is.null(date_col) && date_col %in% names(dt) && data.table::uniqueN(dt[[date_col]]) > 1L) {
    d <- .date_as_num(dt[[date_col]], date_col)
    n_na <- sum(is.na(d))
    u <- sort(unique(d[!is.na(d)]))
    k <- length(u)
    # rows per period in one pass (a scan per period is O(n * periods))
    share_ge <- rev(cumsum(rev(tabulate(match(d, u), nbins = k) / n)))
    cand  <- which(share_ge >= ratio)
    cut_i <- if (length(cand)) max(cand) else k
    if (cut_i == 1L) cut_i <- min(2L, k)
    tr <- which(d < u[cut_i]); ho <- which(d >= u[cut_i])
    if (length(tr) && length(ho)) {
      cutoff <- as.character(dt[[date_col]][which(d == u[cut_i])[1]])
      msg("  OOT: %d period(s) in train, %d in hold-out (hold-out starts at %s, %.1f%% of rows)",
          cut_i - 1L, k - cut_i + 1L, cutoff, 100 * length(ho) / n)
      if (n_na) {
        msg("  WARNING: %s row(s) with a missing '%s' belong to no period and are left out of both train and hold-out.",
            n_fmt(n_na), date_col)
      }
      return(list(train_idx = tr, holdout_idx = ho, method = "out-of-time", cutoff = cutoff))
    }
    msg("  OOT cut degenerate on '%s' - falling back to random stratified.", date_col)
  }
  .scr_local_seed(seed)
  y  <- dt[[target]]
  tr <- integer(0)
  for (lv in unique(y)) {
    idx <- which(y == lv)
    k   <- min(length(idx), max(1L, floor((1 - ratio) * length(idx))))
    tr  <- c(tr, if (length(idx) == 1L) idx else sample(idx, k))
  }
  tr <- sort(unique(tr))
  list(train_idx = tr, holdout_idx = setdiff(seq_len(n), tr),
       method = "stratified random", cutoff = NA_character_)
}

#' Date column as an ordered number (NA where it cannot be read)
#'
#' A `Date`, `POSIXct` or number is used as is. Text or a factor (what a database often
#' returns) is read as an ISO date (`YYYY-MM-DD` or `YYYY/MM/DD`); failing
#' that, `YYYY-MM` as the first day of the month, and all-digit text
#' (`YYYYMM`, `YYYYMMDD`) as a number; any other factor keeps the order of
#' its levels. Without this, `as.numeric()` of text is all `NA` and the out-of-time cut
#' degenerates into a random split.
#' @keywords internal
#' @noRd
.date_as_num <- function(x, date_col) {
  fac <- if (is.factor(x)) x
  if (is.factor(x)) x <- as.character(x)
  if (!is.character(x)) return(suppressWarnings(as.numeric(x)))
  ok <- !is.na(x) & nzchar(x)
  if (!any(ok)) return(rep(NA_real_, length(x)))
  for (fm in c("%Y-%m-%d", "%Y/%m/%d")) {
    d <- as.numeric(as.Date(substr(x, 1L, 10L), format = fm))
    if (!anyNA(d[ok])) return(d)
  }
  if (all(grepl("^[0-9]{4}[-/][0-9]{2}$", x[ok]))) {
    return(as.numeric(as.Date(paste0(gsub("/", "-", x, fixed = TRUE), "-01"), format = "%Y-%m-%d")))
  }
  if (all(grepl("^[0-9]+$", x[ok]))) return(as.numeric(x))
  if (!is.null(fac)) return(as.numeric(fac))   # a factor keeps the order of its levels
  msg("  '%s' is text that is neither an ISO date nor all digits - it cannot order the periods.", date_col)
  rep(NA_real_, length(x))
}
