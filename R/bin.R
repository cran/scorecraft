# ============================================================================ #
# bin.R - Stage 2: optimal binning (parallel by column), screening, hold-out, pruning
# ============================================================================ #
# Almost everything here is delegation to the engine: obwoe() bins,
# obwoe_select() applies the eight admission rules, obwoe_apply()
# materialises the WOE space and obwoe_prune() removes redundancy. What is
# ours: parallelising by column and revalidating on the hold-out with
# FROZEN bins (recomputed IV + PSI), because the engine's screening looks at
# the training rows only.
# ============================================================================ #

#' Stage 2: optimal binning, screening, hold-out revalidation and pruning
#'
#' Fits the bins **on the training rows only** (cut points and WOE use the
#' target), in parallel by column, and applies in sequence:
#'
#' \enumerate{
#'   \item **Screening**, native to the engine: eight admission rules
#'     (`IV_BELOW_MIN`, `IV_SUSPICIOUS`, `NOT_MONOTONIC`, `TOO_FEW_BINS`,
#'     `TOO_MANY_BINS`, `SMALL_BIN`, `DEGENERATE_BIN`, `BINNING_ERROR`).
#'   \item **Hold-out revalidation** with frozen bins: IV recomputed on the
#'     same labels, train/hold-out PSI (the fixed threshold decides; the
#'     n-adjusted one is reported) and the fraction of hold-out without a bin.
#'   \item **Redundancy pruning** by rank correlation on the WOE space,
#'     ranked by hold-out IV. Under `allow_derived_final = FALSE` the derived
#'     flags leave before this step (`derived_excluded`), so that a flag
#'     that cannot be delivered never prunes a real column.
#' }
#'
#' @section Parallelism:
#'
#' Columns are split into `config$nthread` chunks, each chunk is binned by a
#' worker and the fits are merged. The result is identical to the serial one
#' (a test pins this): the engine is deterministic per column.
#'
#' @param triage An object from [scr_triage()].
#' @param config An object from [scr_config()].
#'
#' @return An `scr_bins` object with `fit` (an `obwoe` object), `screen`
#'   (`summary` and `full`), `holdout`, `prune`, `pool` (eligible for the
#'   models), `derived_excluded`, the counts `binned`, `pos_screen` and
#'   `pos_holdout` (the survivors of each gate, in order), the WOE matrices
#'   `woe_train`/`woe_holdout`, the originating `triage` and the `config`.
#'
#' @family stages
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1)
#' sp <- scr_split(scr_demo, "default", date_col = "ref_date", drop = "id")
#' bn <- scr_bin(scr_triage(sp, cfg), cfg)
#' bn
#' head(bn$holdout)
#' @export
scr_bin <- function(triage, config = scr_config()) {
  if (!inherits(triage, "scr_triage")) stop("`triage` must come from scr_triage().", call. = FALSE)
  check_config(config, "scr_bin")
  old <- scr_verbose(isTRUE(config$verbose)); on.exit(scr_verbose(old), add = TRUE)
  cfg <- config
  target <- triage$split$target
  clean  <- triage$clean
  feats  <- setdiff(names(clean), target)
  dt_tr  <- clean[triage$split$train_idx]; dt_ho <- clean[triage$split$holdout_idx]
  y_tr   <- dt_tr[[target]];               y_ho  <- dt_ho[[target]]

  msg_stage(2, "optimal binning and selection (OptimalBinningWoE)")
  fit    <- time_it(sprintf("obwoe (%d features, %d thread(s))", length(feats), cfg$nthread),
                    fit_binning(dt_tr, target, feats, cfg))
  screen <- time_it("obwoe_select", screen_features(fit, cfg))

  binned     <- screen$summary[!(error %in% TRUE), feature]
  pos_screen <- screen$summary[selected %in% TRUE, feature]
  msg("  binned: %d | passed screening: %d", length(binned), length(pos_screen))
  mot <- screen$summary[!(selected %in% TRUE), .N, by = reason][order(-N)]
  for (i in seq_len(nrow(mot))) msg("  screening failed for %-28s %d", mot$reason[i], mot$N[i])

  fit    <- strip_failed_features(fit)
  app_tr <- time_it("obwoe_apply (train)",    apply_woe(fit, dt_tr, binned))
  app_ho <- time_it("obwoe_apply (hold-out)", apply_woe(fit, dt_ho, binned))

  holdout <- time_it("hold-out revalidation (frozen IV + PSI)",
                     holdout_check(app_tr, app_ho, y_tr, y_ho, binned, cfg))
  pos_holdout <- intersect(pos_screen, holdout[holdout_ok %in% TRUE, feature])
  failed_ho   <- setdiff(pos_screen, pos_holdout)
  if (length(failed_ho)) {
    msg("  hold-out failed %d of %d that passed screening", length(failed_ho), length(pos_screen))
    mho <- holdout[feature %in% failed_ho, .N, by = holdout_reason][order(-N)]
    for (i in seq_len(nrow(mho))) msg("    %-46s %d", mho$holdout_reason[i], mho$N[i])
  }

  # The `_bin` columns only serve the revalidation: discarding them is what
  # lets a large table fit in memory.
  bins_tr <- grep("_bin$", names(app_tr), value = TRUE)
  if (length(bins_tr)) {
    app_tr[, (bins_tr) := NULL]
    app_ho[, (intersect(bins_tr, names(app_ho))) := NULL]
  }

  # A derived flag that cannot reach the deliverable leaves BEFORE the
  # pruning: otherwise it could remove, as redundant, a real column that
  # correlates with it, and then be excluded itself - losing both.
  derived_out <- character()
  to_prune <- pos_holdout
  if (!isTRUE(cfg$allow_derived_final)) {
    derived_out <- intersect(pos_holdout, triage$derived)
    if (length(derived_out)) {
      msg("  %d derived flag(s) excluded from the final selection (allow_derived_final = FALSE)", length(derived_out))
      to_prune <- setdiff(pos_holdout, derived_out)
    }
  }

  ranking <- holdout[feature %in% to_prune][order(-iv_holdout), feature]
  prune   <- prune_redundancy(app_tr, to_prune, ranking, cfg)
  if (nrow(prune$dropped)) msg("  redundancy removed %d: %s", nrow(prune$dropped), lst(prune$dropped$feature))
  pool <- prune$keep
  msg("  pool eligible for the models: %d", length(pool))

  structure(list(fit = fit, screen = screen, holdout = holdout, prune = prune, pool = pool,
                 binned = binned, pos_screen = pos_screen, pos_holdout = pos_holdout,
                 derived_excluded = derived_out, woe_train = app_tr, woe_holdout = app_ho,
                 y_train = y_tr, y_holdout = y_ho, triage = triage, config = cfg),
            class = c("scr_bins", "list"))
}

#' @export
print.scr_bins <- function(x, ...) {
  cat(sprintf("<scr_bins> %d binned | screening %d | hold-out %d | pruning %d | pool %d\n",
              length(x$binned), length(x$pos_screen), length(x$pos_holdout), length(x$prune$keep), length(x$pool)))
  m <- x$screen$summary[!(selected %in% TRUE), .N, by = reason][order(-N)]
  for (i in seq_len(nrow(m))) cat(sprintf("  screening: %-26s %d\n", m$reason[i], m$N[i]))
  invisible(x)
}

#' Effective algorithm per variable type
#'
#' `jedi` is universal. A numeric-only (e.g. `ir`) or categorical-only (e.g.
#' `ivb`) algorithm applies where it is valid; the other type falls back to
#' `jedi`, and this is reported.
#' @keywords internal
#' @noRd
.algorithm_for <- function(cfg, type) {
  alg <- cfg$algorithm
  if (identical(alg, "auto")) return("jedi")
  tb  <- OptimalBinningWoE::obwoe_algorithms()
  ok  <- tb[[type]][match(alg, tb$algorithm)]
  if (isTRUE(ok)) alg else "jedi"
}

#' Fit the binning by column chunks in parallel and merge the fits
#' @keywords internal
#' @noRd
fit_binning <- function(dt_train, target, features, cfg) {
  num <- features[vapply(features, function(f) is.numeric(dt_train[[f]]), logical(1))]
  cat <- setdiff(features, num)
  alg_num <- .algorithm_for(cfg, "numerical"); alg_cat <- .algorithm_for(cfg, "categorical")
  if (!identical(cfg$algorithm, "auto")) {
    if (length(num) && alg_num != cfg$algorithm) msg("  algorithm \"%s\" is not numeric: numerics use \"jedi\".", cfg$algorithm)
    if (length(cat) && alg_cat != cfg$algorithm) msg("  algorithm \"%s\" is not categorical: categoricals use \"jedi\".", cfg$algorithm)
  }
  ctrl <- obwoe_control(cfg)
  k <- max(1L, cfg$nthread)
  tasks <- c(
    lapply(.scr_chunks(num, k), function(ch) list(feats = ch, alg = alg_num)),
    lapply(.scr_chunks(cat, k), function(ch) list(feats = ch, alg = alg_cat)))
  tasks <- tasks[vapply(tasks, function(t) length(t$feats) > 0L, logical(1))]

  fits <- .scr_lapply(tasks, function(t) {
    OptimalBinningWoE::obwoe(
      data = as.data.frame(dt_train[, c(target, t$feats), with = FALSE]),
      target = target, feature = t$feats, min_bins = cfg$min_bins, max_bins = cfg$max_bins,
      algorithm = t$alg, control = ctrl)
  }, nthread = k)

  fit <- .merge_obwoe(fits, features)
  s <- data.table::as.data.table(fit$summary)
  errs <- s[error == TRUE, feature]
  if (length(errs)) msg("  binning failed on %d feature(s): %s", length(errs), lst(errs))
  fit
}

#' Merge several obwoe objects into one, in the original feature order
#' @keywords internal
#' @noRd
.merge_obwoe <- function(fits, features) {
  if (length(fits) == 1L) return(fits[[1]])
  out <- fits[[1]]
  out$results <- do.call(c, lapply(fits, `[[`, "results"))
  out$summary <- do.call(rbind, lapply(fits, `[[`, "summary"))
  ord <- match(features, names(out$results))
  ord <- ord[!is.na(ord)]
  out$results <- out$results[ord]
  out$summary <- out$summary[match(names(out$results), out$summary$feature), , drop = FALSE]
  rownames(out$summary) <- NULL
  out$n_features <- length(out$results)
  out$algorithm  <- unique(vapply(fits, function(f) as.character(f$algorithm)[1], character(1)))
  out
}

#' Remove from the fit the features whose binning failed
#' @keywords internal
#' @noRd
strip_failed_features <- function(fit) {
  s   <- data.table::as.data.table(fit$summary)
  bad <- s[error %in% TRUE, feature]
  if (length(bad)) {
    fit$results    <- fit$results[setdiff(names(fit$results), bad)]
    fit$summary    <- as.data.frame(s[!feature %in% bad])
    fit$n_features <- length(fit$results)
  }
  fit
}

#' Native screening: one row per feature, with verdict and reason
#' @keywords internal
#' @noRd
screen_features <- function(fit, cfg) {
  arg <- list(obj = fit, iv_min = cfg$iv_min, iv_max = cfg$iv_max,
              require_monotonic = cfg$require_monotonic, monotonicity = cfg$monotonicity,
              min_bins = cfg$screen_min_bins, max_bins = cfg$max_bins, min_bin_pct = cfg$min_bin_pct,
              allow_degenerate = cfg$allow_degenerate, sort_by = "iv",
              bin_separator = cfg$bin_separator)
  if (!is.null(cfg$screen_top_n)) arg$top_n <- as.integer(cfg$screen_top_n)
  sel  <- data.table::as.data.table(.obwoe_quiet(do.call(OptimalBinningWoE::obwoe_select, arg)))
  full <- data.table::as.data.table(.obwoe_quiet(do.call(OptimalBinningWoE::obwoe_select, c(arg, list(detail = "full")))))
  list(summary = sel, full = full)
}

#' Evaluate an OptimalBinningWoE call without the Rcpp subset-proxy warnings
#'
#' With Rcpp >= 1.1, `OptimalBinningWoE::obwoe_gains_score()` (reached by
#' `obwoe_select(detail = "full")`) reorders its input with the aliasing
#' pattern `x = x[idx]`, and Rcpp's bounds check then warns
#' "subscript out of bounds (index i >= vector size s)" once per element,
#' thousands of times per run. The values returned are correct (checked
#' against a hand computation); the fix belongs upstream. Only that exact
#' message is muffled; every other condition passes through.
#' @keywords internal
#' @noRd
.obwoe_quiet <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("^subscript out of bounds \\(index -?[0-9]+ >= vector size -?[0-9]+\\)$", conditionMessage(w))) {
      invokeRestart("muffleWarning")
    }
  })
}

#' Materialise the WOE space with a frozen fit (subset of features)
#' @keywords internal
#' @noRd
apply_woe <- function(fit, dt, features, what = "both") {
  features <- intersect(features, names(fit$results))
  fit_sub <- fit
  fit_sub$results <- fit_sub$results[features]
  s <- data.table::as.data.table(fit_sub$summary)
  fit_sub$summary <- as.data.frame(s[feature %in% features])
  fit_sub$n_features <- length(features)
  app <- OptimalBinningWoE::obwoe_apply(as.data.frame(dt[, features, with = FALSE]), fit_sub,
                                        keep_original = FALSE)
  data.table::setDT(app)
  cols <- switch(what,
    woe  = paste0(features, "_woe"),
    bin  = paste0(features, "_bin"),
    both = unlist(lapply(features, function(f) c(paste0(f, "_bin"), paste0(f, "_woe")))))
  app[, intersect(cols, names(app)), with = FALSE]
}

#' Revalidate each feature on the hold-out with the FROZEN training bins
#' @keywords internal
#' @noRd
holdout_check <- function(app_train, app_holdout, y_train, y_holdout, features, cfg) {
  rows <- .scr_lapply(features, function(f) {
    cb <- paste0(f, "_bin")
    if (!cb %in% names(app_train) || !cb %in% names(app_holdout)) {
      return(data.table::data.table(
        feature = f, iv_train_bins = NA_real_, iv_holdout = NA_real_, iv_ratio = NA_real_,
        psi = NA_real_, psi_flag = NA_character_, psi_critical = NA_real_,
        psi_flag_adjusted = NA_character_, pct_unbinned = NA_real_, holdout_ok = FALSE,
        holdout_reason = "NO_BIN_APPLIED"))
    }
    b_tr <- as.character(app_train[[cb]]); b_ho <- as.character(app_holdout[[cb]])
    ok_tr <- !is.na(b_tr); ok_ho <- !is.na(b_ho)
    pct_unbinned <- if (length(ok_ho)) 1 - mean(ok_ho) else 1
    iv_tr <- scr_iv(b_tr[ok_tr], y_train[ok_tr])
    iv_ho <- scr_iv(b_ho[ok_ho], y_holdout[ok_ho])
    ratio <- if (iv_tr > 0) iv_ho / iv_tr else NA_real_
    lv    <- union(unique(b_tr[ok_tr]), unique(b_ho[ok_ho]))
    ps    <- if (length(lv) >= 2L && any(ok_ho))
      scr_psi(b_tr[ok_tr], b_ho[ok_ho], levels = lv, alpha = cfg$psi_alpha)
    else list(psi = NA_real_, flag_fixed = NA_character_, critical = NA_real_, flag_adjusted = NA_character_)

    reasons <- character()
    if (is.finite(ratio) && ratio < cfg$iv_ratio_min) reasons <- c(reasons, "IV_DROPS_ON_HOLDOUT")
    if (iv_ho < cfg$iv_min)                           reasons <- c(reasons, "IV_LOW_ON_HOLDOUT")
    if (is.finite(ps$psi) && ps$psi > cfg$psi_max)    reasons <- c(reasons, "PSI_UNSTABLE")
    if (!is.finite(ps$psi))                           reasons <- c(reasons, "PSI_UNDEFINED")
    if (pct_unbinned > cfg$max_unbinned)              reasons <- c(reasons,
      sprintf("UNBINNED_POPULATION_ON_HOLDOUT:%.1f%%", 100 * pct_unbinned))

    data.table::data.table(
      feature = f, iv_train_bins = iv_tr, iv_holdout = iv_ho, iv_ratio = ratio,
      psi = ps$psi, psi_flag = ps$flag_fixed, psi_critical = ps$critical,
      psi_flag_adjusted = ps$flag_adjusted, pct_unbinned = pct_unbinned,
      holdout_ok = length(reasons) == 0L,
      holdout_reason = if (length(reasons)) paste(reasons, collapse = ";") else "OK")
  }, nthread = cfg$nthread, fork_only = TRUE)   # per column the work is cheap: shipping both WOE tables to PSOCK workers is not
  data.table::rbindlist(rows)
}

#' Remove redundancy on the WOE space
#'
#' Greedy pruning with the rule of `OptimalBinningWoE::obwoe_prune()`: the
#' worse-ranked member of the strongest pair at or above `corr_cutoff` leaves,
#' one at a time. For `"pearson"` and `"spearman"` the correlation matrix is
#' computed by the compiled kernel (`.scr_cor_matrix()`: each column ranked
#' once, one BLAS cross-product) and the greedy sweep runs on it
#' (`.scr_prune_matrix()`), with the same result as `obwoe_prune()` at a
#' fraction of the cost on wide tables; the other methods go through
#' `obwoe_prune()`. The WOE columns are read in place, not copied.
#' @keywords internal
#' @noRd
prune_redundancy <- function(app_train, features, ranking, cfg) {
  empty <- data.table::data.table(feature = character(), correlated_with = character(), correlation = numeric())
  if (length(features) < 2L) return(list(keep = features, dropped = empty))
  x <- lapply(paste0(features, "_woe"), function(f) app_train[[f]])
  names(x) <- features
  const <- names(x)[vapply(x, function(v) !is.finite(stats::sd(v)) || stats::sd(v) == 0, logical(1))]
  if (length(const)) x <- x[setdiff(names(x), const)]
  if (length(x) < 2L) {
    return(list(keep = names(x), dropped = if (length(const)) data.table::data.table(
      feature = const, correlated_with = NA_character_, correlation = NA_real_) else empty))
  }
  rk <- intersect(ranking, names(x))
  fast <- cfg$corr_method %in% c("pearson", "spearman") &&
    all(vapply(x, function(v) is.numeric(v) && all(is.finite(v)), logical(1)))
  pr <- if (fast) {
    .scr_prune_matrix(.scr_cor_matrix(x, cfg$corr_method, cfg$nthread), rk, cfg$corr_cutoff)
  } else {
    OptimalBinningWoE::obwoe_prune(as.data.frame(x), ranking = rk, cutoff = cfg$corr_cutoff, method = cfg$corr_method)
  }
  dropped <- data.table::as.data.table(pr$dropped)
  if (nrow(dropped)) data.table::setnames(dropped, names(dropped)[1], "feature")
  if (length(const)) {
    dropped <- data.table::rbindlist(list(dropped, data.table::data.table(
      feature = const, correlated_with = NA_character_, correlation = NA_real_)), use.names = TRUE, fill = TRUE)
  }
  list(keep = pr$keep, dropped = if (nrow(dropped)) dropped else empty)
}
