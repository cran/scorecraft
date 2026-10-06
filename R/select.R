# ============================================================================ #
# select.R - orchestration of stages 0-4 for one target, and for several
# ============================================================================ #

#' Select variables for the scorecard
#'
#' Shortcut that chains [scr_split()], [scr_triage()], [scr_bin()] and
#' [scr_model()] on a table and a binary target, and returns an object with
#' the shortlist, the complete audit funnel, the gains table and the
#' production SQL of the approved variables. Every stage remains callable on
#' its own for whoever wants more control (hybrid interface).
#'
#' @section Reproducibility:
#'
#' With the same `data`, the same `target` and the same `config$seed`, the
#' result is identical with one or several `nthread`: the seed governs the
#' random split, the cross-validation, the classifier subsample, the trees
#' and the bootstrap, and the binning is deterministic per column.
#'
#' @param data A `data.frame` or `data.table` with the target, the candidates
#'   and, if any, the date column of the out-of-time cut.
#' @param target Name of the target column (0/1, logical, or a two-level
#'   factor/character).
#' @param config An object from [scr_config()].
#' @param drop Columns that are never candidates. They stay in the funnel as
#'   `00.config`.
#' @param date_col Date column of the out-of-time cut. Defaults to
#'   `config$oot_date_col`.
#' @param event_level Which target value counts as the event; see [scr_split()].
#' @param export Directory to write the deliverables to. `NULL` (default)
#'   writes nothing; use [scr_export()] later.
#' @param copy If `TRUE` (default), works on a copy of `data`.
#'
#' @return An object of class `scr_result`. Read it with [scr_selected()],
#'   [scr_funnel()], [scr_gains()], [scr_sql()], [scr_leakage()] and
#'   [summary()]; continue with [scr_scorecard()]; write it with [scr_export()].
#'
#' @seealso [scr_run()] for several targets straight from the database,
#'   [scr_scorecard()] for the next step.
#' @family stages
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   xgb_rounds = 60, n_boot = 20)
#' res <- scr_select(scr_demo, "default", config = cfg, drop = "id",
#'                   date_col = "ref_date")
#' res
#' scr_selected(res)
#' head(scr_funnel(res, only_selected = TRUE))
#' @export
scr_select <- function(data, target, config = scr_config(), drop = character(),
                       date_col = config$oot_date_col, event_level = NULL,
                       export = NULL, copy = TRUE) {
  check_config(config, "scr_select")
  cfg <- config
  cfg$oot_date_col <- date_col
  old <- scr_verbose(isTRUE(cfg$verbose)); on.exit(scr_verbose(old), add = TRUE)
  t0 <- Sys.time()

  msg_stage(0, "ingestion, typing and split")
  split <- scr_split(data, target, date_col = date_col, ratio = cfg$holdout_ratio, seed = cfg$seed,
                     event_level = event_level, drop = drop, copy = copy)
  msg("  split %s: train=%s hold-out=%s%s", split$method, n_fmt(length(split$train_idx)),
      n_fmt(length(split$holdout_idx)), if (!is.na(split$cutoff)) sprintf(" (cut at %s)", split$cutoff) else "")
  msg("  columns: %d | candidates: %d (%d numeric, %d categorical) | dropped: %d",
      ncol(split$data), length(split$cols$features), length(split$cols$var_num),
      length(split$cols$var_cat), length(split$cols$dropped))
  vv <- vocab(cfg)
  msg("  convention: objective=\"%s\" - %s; %s", cfg$objective, vv$target1, vv$points)
  if (length(split$holdout_idx) < 100L)
    msg("  WARNING: hold-out with %d rows - revalidation will be weak.", length(split$holdout_idx))

  triage <- scr_triage(split, cfg)
  bins   <- scr_bin(triage, cfg)
  models <- scr_model(bins, cfg)

  msg_stage(7, "funnel, gains, SQL and summary")
  funnel <- build_funnel(split$cols, triage, bins, models, cfg)
  gains  <- build_gains(bins, models, cfg)
  sel    <- models$consensus$selected
  sql    <- build_sql_woe(bins$fit, triage$ledger, sel, cfg, target)

  y_tr <- bins$y_train; y_ho <- bins$y_holdout
  meta <- list(target = target, event = split$cols$event, objective = cfg$objective,
               score_direction = resolve_direction(cfg),
               n_total = nrow(split$data), n_train = length(split$train_idx),
               n_holdout = length(split$holdout_idx), split_method = split$method,
               split_cutoff = split$cutoff, date_col = date_col,
               event_rate_train = mean(y_tr), event_rate_holdout = mean(y_ho),
               n_cols = ncol(split$data), n_candidates = length(split$cols$features),
               n_after_triage = length(triage$keep), n_binned = length(bins$binned),
               n_after_screening = length(bins$pos_screen), n_after_holdout = length(bins$pos_holdout),
               n_after_prune = length(bins$pool), derived = triage$derived,
               seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  summary_md <- build_summary(meta, funnel, models, cfg)
  msg("Target %s finished in %.1fs - %d variables approved out of %d candidates.",
      target, meta$seconds, length(sel), length(split$cols$features))

  date_vec <- if (!is.null(date_col) && date_col %in% names(split$data)) split$data[[date_col]] else NULL

  res <- structure(
    list(target = target, config = cfg, meta = meta,
         split = split[setdiff(names(split), "data")],
         triage = triage[c("profile", "ledger", "keep", "derived")],
         fit = bins$fit, screen = bins$screen, holdout = bins$holdout, prune = bins$prune,
         pool = bins$pool, derived_excluded = bins$derived_excluded,
         models = models[c("votes", "metrics")], consensus = models$consensus,
         funnel = funnel, gains = gains, sql = sql, summary_md = summary_md,
         data_clean = triage$clean, date = date_vec, files = NULL),
    class = c("scr_result", "list"))
  if (!is.null(export)) res$files <- scr_export(res, export)$files
  res
}

#' Run the selection for several targets straight from the database
#'
#' For each target: fetches the table, runs [scr_select()] and writes the
#' deliverables. A failure on one target is recorded and the loop continues.
#'
#' @section Table convention:
#'
#' `table` accepts the `{target}` placeholder, replaced by the lower-case
#' target name. Without the placeholder, the same table is used for every
#' target.
#'
#' @param con A DBI connection, from [scr_connect()].
#' @param table Table name, with an optional `{target}`.
#' @param targets Vector with the names of the target columns.
#' @param config An object from [scr_config()], used for every target.
#' @param drop Columns that are never candidates.
#' @param date_col Date column of the out-of-time split, passed on to
#'   [scr_select()]. Defaults to `config$oot_date_col`; an explicit `NULL`
#'   forces a random stratified split.
#' @param event_level Passed on to [scr_select()].
#' @param sample_frac Sampling fraction. A scalar or a list named by target.
#' @param max_rows Row cap per target. `NULL` switches it off.
#' @param export Root output directory; each target writes to a subdirectory.
#'
#' @return An `scr_runset` object: a named list of `scr_result` (or, for the
#'   targets that failed, a list with `error`).
#'
#' @seealso [scr_compare()] and [scr_core()] to read the run set.
#' @family portfolio
#' @examplesIf requireNamespace("RSQLite", quietly = TRUE)
#' con <- scr_connect(driver = RSQLite::SQLite(), dbname = ":memory:")
#' d <- scr_demo; d$ref_date <- as.character(d$ref_date)   # SQLite has no Date type
#' DBI::dbWriteTable(con, "dtm", d)
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   xgb_rounds = 60, n_boot = 20)
#' rs <- scr_run(con, "dtm", targets = c("default", "churn"), config = cfg,
#'               drop = c("id", "ref_date", "default", "churn"))
#' rs
#' scr_compare(rs)
#' DBI::dbDisconnect(con)
#' @export
scr_run <- function(con, table, targets, config = scr_config(), drop = character(), date_col = config$oot_date_col,
                    event_level = NULL, sample_frac = 1.0, max_rows = NULL, export = NULL) {
  check_config(config, "scr_run")
  old <- scr_verbose(isTRUE(config$verbose)); on.exit(scr_verbose(old), add = TRUE)
  if (!length(targets)) stop("`targets` is empty.", call. = FALSE)
  frac_of <- function(vt) if (is.list(sample_frac)) sample_frac[[vt]] %||% 1.0 else sample_frac
  res <- vector("list", length(targets)); names(res) <- targets
  for (vt in targets) {
    msg("\n############### TARGET: %s ###############", vt)
    res[[vt]] <- tryCatch({
      tab <- gsub("{target}", tolower(vt), table, fixed = TRUE)
      dt  <- scr_fetch(con, tab, frac_of(vt), config$seed, max_rows)
      cfg <- config; cfg$sql_table <- tab
      out <- if (is.null(export)) NULL else file.path(export, tolower(vt))
      # sibling targets are never candidates for one another
      scr_select(dt, vt, config = cfg, drop = setdiff(drop, vt), date_col = date_col,
                 event_level = event_level, export = out, copy = FALSE)
    }, error = function(e) {
      msg("ERROR on target '%s': %s", vt, conditionMessage(e))
      list(target = vt, error = conditionMessage(e))
    })
    gc(verbose = FALSE)
  }
  ok <- vapply(res, inherits, logical(1), "scr_result")
  msg("\nDone: %d of %d target(s) succeeded.", sum(ok), length(res))
  structure(res, class = c("scr_runset", "list"))
}

#' Compare runs across targets
#'
#' One row per target, with the funnel, the hold-out performance of the best
#' model (with CI) and the warning signs.
#'
#' @param x An object from [scr_run()], or a named list of `scr_result`.
#'
#' @return A `data.table` with one row per successful target.
#'
#' @family portfolio
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   xgb_rounds = 60, n_boot = 20)
#' r1 <- scr_select(scr_demo, "default", config = cfg, drop = c("id", "churn"),
#'                  date_col = "ref_date")
#' r2 <- scr_select(scr_demo, "churn", config = cfg, drop = c("id", "default"),
#'                  date_col = "ref_date")
#' scr_compare(list(default = r1, churn = r2))
#' scr_core(list(default = r1, churn = r2), min_targets = 2)
#' @export
scr_compare <- function(x) {
  ok <- vapply(x, inherits, logical(1), "scr_result")
  if (!any(ok)) stop("no valid result to compare.", call. = FALSE)
  data.table::rbindlist(lapply(names(x)[ok], function(nm) {
    r <- x[[nm]]; m <- r$meta; f <- r$funnel; mt <- r$models$metrics
    best <- if (nrow(mt) && any(is.finite(mt$auc))) mt[which.max(mt$auc)] else NULL
    data.table::data.table(
      target = nm, rows = m$n_total, train = m$n_train, holdout = m$n_holdout,
      event_rate = round(m$event_rate_train, 4), candidates = m$n_candidates, triage = m$n_after_triage,
      screening = m$n_after_screening, holdout_ok = m$n_after_holdout, pool = m$n_after_prune,
      approved = length(scr_selected(r)),
      max_iv_approved = round(suppressWarnings(max(f[approved == TRUE, total_iv], na.rm = TRUE)), 3),
      n_iv_suspect = nrow(f[approved == TRUE & iv_suspect %in% TRUE]),
      best_model = if (!is.null(best)) best$model else NA_character_,
      auc = if (!is.null(best)) round(best$auc, 4) else NA_real_,
      auc_lo = if (!is.null(best)) round(best$auc_lo, 4) else NA_real_,
      auc_hi = if (!is.null(best)) round(best$auc_hi, 4) else NA_real_,
      ks = if (!is.null(best)) round(best$ks, 4) else NA_real_,
      seconds = round(m$seconds), relaxation = r$consensus$meta$relaxation)
  }), fill = TRUE)[]
}

#' Variables that cross several targets
#'
#' Which variables were approved on how many targets. A stable core across
#' targets is the best argument in favour of a variable.
#'
#' @param x An object from [scr_run()], or a named list of `scr_result`.
#' @param min_targets Minimum number of targets to enter the result.
#'
#' @return A `data.table` with `feature`, `n_targets`, `targets` and `mean_rank`.
#'
#' @family portfolio
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   xgb_rounds = 60, n_boot = 20)
#' r1 <- scr_select(scr_demo, "default", config = cfg, drop = c("id", "churn"),
#'                   date_col = "ref_date")
#' r2 <- scr_select(scr_demo, "churn", config = cfg, drop = c("id", "default"),
#'                  date_col = "ref_date")
#' scr_core(list(default = r1, churn = r2), min_targets = 2)
#' @export
scr_core <- function(x, min_targets = 2L) {
  ok <- vapply(x, inherits, logical(1), "scr_result")
  if (!any(ok)) stop("no valid result.", call. = FALSE)
  d <- data.table::rbindlist(lapply(names(x)[ok], function(nm) {
    f <- x[[nm]]$funnel[approved == TRUE]
    data.table::data.table(target = nm, feature = f$feature, pos = f$consensus_rank)
  }))
  out <- d[, .(n_targets = data.table::uniqueN(target), mean_rank = round(mean(pos, na.rm = TRUE), 1),
               targets = paste(sort(unique(target)), collapse = ", ")), by = feature]
  out <- out[n_targets >= min_targets]
  data.table::setorder(out, -n_targets, mean_rank)
  out[]
}
