# ============================================================================ #
# model.R - Stages 3 and 4: classifiers on the WOE space and consensus
# ============================================================================ #
# Different families of evidence over the SAME WOE columns:
#   glmnet   - regularisation (elastic net): who survives the shrinkage
#   xgboost  - boosting: accumulated gain (required Import)
#   lightgbm - alternative boosting, optional (Suggests)
#   ranger   - random forest: permutation importance
# Each model is trained on train and evaluated on the HOLD-OUT; its hold-out
# Gini becomes its weight in the consensus.
# ============================================================================ #

#' Stages 3 and 4: multi-strategy selection and consensus
#'
#' Trains the classifiers enabled in the configuration on the WOE columns of
#' the eligible pool, measures each on the hold-out (AUC/KS/Gini with a
#' bootstrap CI) and combines the votes:
#'
#' \preformatted{
#' consensus_score = mean of the importance rank percentiles, weighted by the
#'                   hold-out Gini of each model
#' votes           = how many models elected the feature (top-K, or non-zero
#'                   coefficient in the elastic net)
#' }
#'
#' The final cut respects `[target_min, target_max]`. If the strict
#' consensus does not reach `target_min`, relaxation happens in **named**,
#' recorded steps (`min_votes` reduced; completed by score), never
#' resurrecting a feature failed by an earlier gate.
#'
#' @param bins An object from [scr_bin()].
#' @param config An object from [scr_config()].
#'
#' @return An `scr_models` object with `votes` (one row per model and
#'   feature), `metrics` (one row per model, with CI), `consensus` (`table`,
#'   `selected`, `meta`) and the originating `bins`.
#'
#' @family stages
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   xgb_rounds = 60, n_boot = 20)
#' sp <- scr_split(scr_demo, "default", date_col = "ref_date", drop = "id")
#' md <- scr_model(scr_bin(scr_triage(sp, cfg), cfg), cfg)
#' md
#' md$consensus$selected
#' @export
scr_model <- function(bins, config = scr_config()) {
  if (!inherits(bins, "scr_bins")) stop("`bins` must come from scr_bin().", call. = FALSE)
  check_config(config, "scr_model")
  old <- scr_verbose(isTRUE(config$verbose)); on.exit(scr_verbose(old), add = TRUE)
  cfg <- config
  pool <- bins$pool
  y_tr <- bins$y_train; y_ho <- bins$y_holdout

  msg_stage(3, "classifiers on the WOE space")
  if (length(pool) >= 2L) {
    i_tr <- subsample_stratified(y_tr, cfg$model_max_rows, cfg$seed)
    i_ho <- subsample_stratified(y_ho, cfg$model_max_rows, cfg$seed)
    if (length(i_tr) < length(y_tr) || length(i_ho) < length(y_ho))
      msg("  stratified subsample for the models: train %s, hold-out %s (cap %s)",
          n_fmt(length(i_tr)), n_fmt(length(i_ho)), n_fmt(cfg$model_max_rows))
    # rows AND columns in one subset: the WOE tables hold every binned
    # feature, the classifiers only need the pool
    wc <- intersect(paste0(pool, "_woe"), names(bins$woe_train))
    models <- run_classifiers(bins$woe_train[i_tr, wc, with = FALSE], bins$woe_holdout[i_ho, wc, with = FALSE],
                              y_tr[i_tr], y_ho[i_ho], pool, cfg)
  } else {
    msg("  pool with %d feature(s): classifiers do not run - consensus degenerates to the IV.", length(pool))
    iv <- bins$holdout[feature %in% pool, iv_holdout]
    models <- list(
      votes = data.table::data.table(model = "(no model)", feature = pool, importance = iv,
                                     score_pct = rep(1, length(pool)), vote = rep(TRUE, length(pool))),
      metrics = data.table::data.table(model = "(no model)", auc = NA_real_, auc_lo = NA_real_,
                                       auc_hi = NA_real_, ks = NA_real_, ks_lo = NA_real_, ks_hi = NA_real_,
                                       gini = NA_real_, gini_lo = NA_real_, gini_hi = NA_real_,
                                       n_vote = length(pool), note = "pool < 2 features: no classifier ran"))
  }

  msg_stage(4, "consensus")
  consensus <- build_consensus(pool, models$votes, models$metrics, bins$holdout, cfg)
  msg("  approved: %d (min_votes used: %d of %d configured) | relaxation: %s",
      length(consensus$selected), consensus$meta$min_votes_used, consensus$meta$min_votes_config,
      consensus$meta$relaxation)
  msg("  %s", lst(consensus$selected, 40))

  structure(list(votes = models$votes, metrics = models$metrics, consensus = consensus,
                 bins = bins, config = cfg), class = c("scr_models", "list"))
}

#' @export
print.scr_models <- function(x, ...) {
  cat(sprintf("<scr_models> %d model(s) | pool %d | approved %d | relaxation: %s\n",
              nrow(x$metrics), x$consensus$meta$n_pool, length(x$consensus$selected), x$consensus$meta$relaxation))
  for (i in seq_len(nrow(x$metrics))) {
    m <- x$metrics[i]
    cat(sprintf("  %-9s AUC %.4f [%.4f, %.4f]  KS %.4f  votes %d\n", m$model, m$auc, m$auc_lo, m$auc_hi, m$ks, m$n_vote))
  }
  invisible(x)
}

#' Run every enabled classifier
#' @keywords internal
#' @noRd
run_classifiers <- function(app_train, app_holdout, y_train, y_holdout, features, cfg) {
  cols  <- paste0(features, "_woe")
  missing <- setdiff(cols, names(app_train))
  if (length(missing)) stop("run_classifiers(): WOE column missing: ", lst(missing), call. = FALSE)
  x_tr <- as.matrix(app_train[, cols, with = FALSE])
  x_ho <- as.matrix(app_holdout[, cols, with = FALSE])
  colnames(x_tr) <- colnames(x_ho) <- features

  adapters <- list()
  if (isTRUE(cfg$use_glmnet))   adapters$glmnet   <- .fit_glmnet
  if (isTRUE(cfg$use_xgboost))  adapters$xgboost  <- .fit_xgboost
  if (isTRUE(cfg$use_ranger))   adapters$ranger   <- .fit_ranger
  if (isTRUE(cfg$use_lightgbm)) adapters$lightgbm <- .fit_lightgbm

  votes <- list(); metrics <- list()
  for (nm in names(adapters)) {
    res <- tryCatch(adapters[[nm]](x_tr, y_train, x_ho, y_holdout, cfg),
                    error = function(e) list(err = conditionMessage(e)))
    if (!is.null(res$err)) { msg("  %-9s unavailable: %s", nm, res$err); next }
    imp <- res$importance
    imp[!is.finite(importance), importance := 0]
    imp[, score_pct := rank_pct(importance)]
    if (is.null(res$vote)) {
      k <- min(cfg$model_top_k, nrow(imp))
      cut_k <- sort(imp$importance, decreasing = TRUE)[k]
      imp[, vote := importance > 0 & importance >= cut_k]
    } else {
      imp[, vote := feature %in% res$vote]
    }
    m <- scr_metrics(res$score, y_holdout, ci = TRUE, n_boot = cfg$n_boot, level = cfg$ci_level,
                     seed = cfg$seed, nthread = cfg$nthread)
    msg("  %-9s hold-out AUC=%.4f [%.4f, %.4f] KS=%.4f | votes=%d%s", nm, m$auc, m$auc_lo, m$auc_hi,
        m$ks, sum(imp$vote), if (nzchar(res$note %||% "")) paste0(" | ", res$note) else "")
    votes[[nm]]   <- data.table::data.table(model = nm, imp)
    metrics[[nm]] <- data.table::data.table(model = nm, auc = m$auc, auc_lo = m$auc_lo, auc_hi = m$auc_hi,
                                            ks = m$ks, ks_lo = m$ks_lo, ks_hi = m$ks_hi, gini = m$gini,
                                            gini_lo = m$gini_lo, gini_hi = m$gini_hi,
                                            n_vote = sum(imp$vote), note = res$note %||% "")
  }
  if (!length(votes)) stop("No classifier ran - a consensus cannot be formed.", call. = FALSE)
  list(votes = data.table::rbindlist(votes), metrics = data.table::rbindlist(metrics))
}

# -- adapters --------------------------------------------------------------- #

#' Elastic net: the regularisation itself is the selector (coefficient != 0)
#' @keywords internal
#' @noRd
.fit_glmnet <- function(x_tr, y_tr, x_ho, y_ho, cfg) {
  if (!requireNamespace("glmnet", quietly = TRUE)) stop("package 'glmnet' is not installed")
  .scr_local_seed(cfg$seed)   # the folds are random; the user's stream is restored on exit
  cv <- glmnet::cv.glmnet(x_tr, y_tr, family = "binomial", alpha = cfg$en_alpha,
                          nfolds = cfg$cv_folds, standardize = TRUE)
  s <- "lambda.1se"; b <- as.numeric(stats::coef(cv, s = s))[-1L]; note <- ""
  if (all(b == 0)) {
    s <- "lambda.min"; b <- as.numeric(stats::coef(cv, s = s))[-1L]
    note <- "lambda.1se zeroed everything; lambda.min used"
  }
  sds <- apply(x_tr, 2L, stats::sd)
  list(importance = data.table::data.table(feature = colnames(x_tr), importance = abs(b) * sds),
       vote  = colnames(x_tr)[b != 0],
       score = as.numeric(stats::predict(cv, newx = x_ho, s = s, type = "link")), note = note)
}

#' XGBoost: accumulated gain, with early stopping on the hold-out
#' @keywords internal
#' @noRd
.fit_xgboost <- function(x_tr, y_tr, x_ho, y_ho, cfg) {
  # xgboost draws its own seed from R's generator when `seed` is not a
  # parameter (>= 2.0) and subsamples with it (1.7): seed it, locally
  .scr_local_seed(cfg$seed)
  p <- list(objective = "binary:logistic", eval_metric = "auc", eta = cfg$xgb_eta,
            max_depth = cfg$xgb_max_depth, subsample = cfg$xgb_subsample,
            colsample_bytree = cfg$xgb_colsample, min_child_weight = cfg$xgb_min_child_weight,
            nthread = cfg$nthread)
  d_tr <- xgboost::xgb.DMatrix(data = x_tr, label = y_tr)
  d_ho <- xgboost::xgb.DMatrix(data = x_ho, label = y_ho)
  # `watchlist` was renamed `evals` in the reworked R interface: ask the
  # function itself rather than guess the version that did it
  arg_eval <- if ("evals" %in% names(formals(xgboost::xgb.train))) "evals" else "watchlist"
  args <- list(params = p, data = d_tr, nrounds = cfg$xgb_rounds,
               early_stopping_rounds = cfg$xgb_early_stopping, verbose = 0L)
  args[[arg_eval]] <- list(valid = d_ho)
  m <- do.call(xgboost::xgb.train, args)
  imp <- data.table::as.data.table(xgboost::xgb.importance(model = m))
  imp <- if (nrow(imp)) imp[, .(feature = Feature, importance = Gain)] else
    data.table::data.table(feature = character(), importance = numeric())
  absent <- setdiff(colnames(x_tr), imp$feature)
  if (length(absent)) imp <- data.table::rbindlist(list(imp, data.table::data.table(feature = absent, importance = 0)))
  # the booster attribute is the 0-based index of the best round (1.7 and
  # >= 2.0 alike); predict() already stops there
  best <- tryCatch(xgboost::xgb.attr(m, "best_iteration"), error = function(e) NULL)
  n_trees <- if (length(best) && is.finite(suppressWarnings(as.numeric(best[1]))))
    as.integer(as.numeric(best[1])) + 1L else m$best_iteration %||% m$niter %||% cfg$xgb_rounds
  list(importance = imp, vote = NULL, score = as.numeric(stats::predict(m, d_ho)),
       note = sprintf("%s trees", n_trees), model = m)
}

#' Random forest: permutation importance
#' @keywords internal
#' @noRd
.fit_ranger <- function(x_tr, y_tr, x_ho, y_ho, cfg) {
  if (!requireNamespace("ranger", quietly = TRUE)) stop("package 'ranger' is not installed")
  df <- as.data.frame(x_tr, check.names = FALSE)
  df[["target__"]] <- factor(y_tr, levels = c(0L, 1L))
  .scr_rng_guard()   # ranger and its predict() draw from the R stream even with `seed`
  m <- ranger::ranger(dependent.variable.name = "target__", data = df, num.trees = cfg$rf_trees,
                      importance = cfg$rf_importance, probability = TRUE, num.threads = cfg$nthread,
                      seed = cfg$seed, respect.unordered.factors = "order")
  imp <- ranger::importance(m)
  pr  <- stats::predict(m, data = as.data.frame(x_ho, check.names = FALSE))$predictions
  list(importance = data.table::data.table(feature = names(imp), importance = as.numeric(imp)),
       vote = NULL, score = as.numeric(pr[, "1"]), note = sprintf("importance=%s", cfg$rf_importance))
}

#' LightGBM: only joins if the package exists on this machine
#' @keywords internal
#' @noRd
.fit_lightgbm <- function(x_tr, y_tr, x_ho, y_ho, cfg) {
  if (!requireNamespace("lightgbm", quietly = TRUE)) stop("package 'lightgbm' is not installed")
  .scr_local_seed(cfg$seed)
  d_tr <- lightgbm::lgb.Dataset(data = x_tr, label = y_tr)
  d_ho <- lightgbm::lgb.Dataset.create.valid(d_tr, data = x_ho, label = y_ho)
  m <- lightgbm::lgb.train(
    params = list(objective = "binary", metric = "auc", learning_rate = cfg$xgb_eta,
                  max_depth = cfg$xgb_max_depth, feature_fraction = cfg$xgb_colsample,
                  bagging_fraction = cfg$xgb_subsample, bagging_freq = 1L,
                  # xgboost's min_child_weight is a sum of hessians, not a row
                  # count: its LightGBM counterpart is min_sum_hessian_in_leaf
                  # (documented alias `min_child_weight`), not min_data_in_leaf
                  min_sum_hessian_in_leaf = cfg$xgb_min_child_weight, num_threads = cfg$nthread,
                  verbosity = -1L, seed = cfg$seed),
    data = d_tr, nrounds = cfg$xgb_rounds, valids = list(valid = d_ho),
    early_stopping_rounds = cfg$xgb_early_stopping, verbose = -1L)
  imp <- data.table::as.data.table(lightgbm::lgb.importance(m))
  imp <- if (nrow(imp)) imp[, .(feature = Feature, importance = Gain)] else
    data.table::data.table(feature = character(), importance = numeric())
  absent <- setdiff(colnames(x_tr), imp$feature)
  if (length(absent)) imp <- data.table::rbindlist(list(imp, data.table::data.table(feature = absent, importance = 0)))
  list(importance = imp, vote = NULL, score = as.numeric(stats::predict(m, x_ho)),
       note = sprintf("%d trees", m$best_iter %||% cfg$xgb_rounds), model = m)
}

# -- consensus -------------------------------------------------------------- #

#' Weighted consensus with named relaxation
#' @keywords internal
#' @noRd
build_consensus <- function(pool, votes, metrics, iv_ref, cfg) {
  if (!length(pool)) {
    return(list(table = data.table::data.table(feature = character(), votes = integer(),
                                               consensus_score = numeric(), iv_holdout = numeric(),
                                               selected = logical(), reason = character(),
                                               consensus_rank = integer()),
                selected = character(),
                meta = list(n_pool = 0L, min_votes_used = cfg$min_votes, min_votes_config = cfg$min_votes,
                            relaxation = "empty pool", scarce = TRUE, cut_top_n = 0L, weights = NULL)))
  }
  w <- if (isTRUE(cfg$weight_by_gini)) metrics[, .(model, weight = pmax(0, gini))] else metrics[, .(model, weight = 1)]
  w[!is.finite(weight), weight := 0]
  if (sum(w$weight) <= 0) w[, weight := 1]

  v  <- merge(votes[feature %in% pool], w, by = "model", all.x = TRUE)
  tb <- v[, .(votes = sum(vote), consensus_score = sum(score_pct * weight) / sum(weight)), by = feature]
  tb <- merge(tb, iv_ref[, .(feature, iv_holdout)], by = "feature", all.x = TRUE)
  data.table::setorder(tb, -consensus_score, -iv_holdout)

  mv    <- as.integer(cfg$min_votes)
  pass  <- tb[votes >= mv, feature]
  relax <- character()
  while (length(pass) < cfg$target_min && mv > 1L) {
    mv    <- mv - 1L
    pass  <- tb[votes >= mv, feature]
    relax <- c(relax, sprintf("min_votes reduced to %d", mv))
  }
  if (length(pass) < cfg$target_min) {
    short  <- cfg$target_min - length(pass)
    extras <- utils::head(setdiff(tb$feature, pass), short)
    if (length(extras)) {
      pass  <- c(pass, extras)
      relax <- c(relax, sprintf("%d completed by consensus score (zero votes)", length(extras)))
    }
  }
  pass <- tb[feature %in% pass, feature]
  cut  <- character()
  if (length(pass) > cfg$target_max) {
    cut  <- pass[(cfg$target_max + 1L):length(pass)]
    pass <- pass[seq_len(cfg$target_max)]
  }
  tb[, selected := feature %in% pass]
  tb[, reason := data.table::fifelse(
    selected & votes >= as.integer(cfg$min_votes), "OK",
    data.table::fifelse(selected, "INCLUDED_BY_RELAXATION",
    data.table::fifelse(feature %in% cut, "NOT_IN_TOP_N", "INSUFFICIENT_VOTES")))]
  tb[, consensus_rank := seq_len(.N)]

  scarce <- length(pass) < cfg$target_min
  if (scarce) {
    msg("  WARNING: consensus returned %d variable(s), below the requested minimum (%d) - the eligible pool has %d. No failed variable was used to fill the number.",
        length(pass), cfg$target_min, length(pool))
  }
  list(table = tb[], selected = pass,
       meta = list(n_pool = length(pool), min_votes_used = mv, min_votes_config = as.integer(cfg$min_votes),
                   relaxation = if (length(relax)) paste(relax, collapse = "; ") else "none",
                   scarce = scarce, cut_top_n = length(cut), weights = w))
}
