# ============================================================================ #
# scorecard.R - Stages 4 and 5: points scorecard, aligned, with a challenger
# ============================================================================ #

#' Stages 4 and 5: points scorecard, aligned to the declared scale
#'
#' Fits a logistic regression on the WOE columns of the shortlist, checks
#' the sign of the coefficients, aligns the logit to the declared scale with
#' [scr_align()] (always) and distributes the points per bin. Measures the
#' score on train and hold-out with a bootstrap CI (always), builds the
#' gains with bands **frozen on train**, the score PSI and the CSI per
#' variable (fixed and n-adjusted thresholds), the
#' calibration and the rank-order diagnostics. Optionally fits a tree
#' challenger on the same WOE columns, aligned to the same scale, with an
#' explicit `supports_scorecard = FALSE`: it compares, it never
#' produces points or reason codes.
#'
#' @section Sign check:
#'
#' The engine's WOE is event-oriented, so every glm coefficient must be
#' positive. A variable with a non-positive coefficient (or above
#' `max_abs_coef` in absolute value) is explaining what another already
#' explained, with the sign reversed; it is removed and the model refitted,
#' one at a time, the most negative first, and each removal is recorded in
#' `sign_check`. The last remaining variable is never removed: it is kept
#' and flagged `NON_POSITIVE_COEF_KEPT_LAST`.
#' The final shortlist of the scorecard (`features`) is what [scr_sql()]
#' covers.
#'
#' @section Points per bin:
#'
#' With `score = a + b * logit` and `logit = alpha + sum(beta_j * woe_ij)`:
#' \deqn{\mathrm{points}_{ij} = b\,\beta_j\,\mathrm{woe}_{ij},\qquad
#'       \mathrm{base} = a + b\,\alpha.}
#' `points_style = "distributed"` spreads `base / k` over each
#' characteristic (Siddiqi, 2006, chapter 6), leaving `base_points = 0`. The
#' exact points stay in `points_raw`; `points` is the rounded version when
#' `points_round = TRUE`. The exact score (`score`) and the whole-points
#' score (`score_points`) are both returned by [scr_apply()] and both
#' emitted by [scr_sql()]. A row that falls in no fitted bin (a category
#' never seen on train, a missing value without a missing bin) gets a WOE of
#' 0 from the binning engine, hence the points of a WOE of 0: `0`, or
#' `base / k` under `"distributed"`.
#'
#' @param x An object from [scr_select()].
#' @param features Variables of the scorecard. Defaults to [scr_selected()].
#' @param base_score,base_odds,pdo,direction The scale; `NULL` uses the
#'   configuration of `x`. See [scr_config()] and [scr_align()].
#' @param align_method `"regression"` or `"direct"`; `NULL` uses the configuration.
#' @param challenger `NULL`, `"xgboost"` or `"lightgbm"`; `NULL` uses the configuration.
#' @param points_style `"base_plus_deviation"` or `"distributed"`; `NULL` uses the configuration.
#' @param n_boot CI resamples; `NULL` uses the configuration.
#' @param seed Seed; `NULL` uses the configuration.
#'
#' @return An `scr_scorecard` object. Main components: `features`, `coef`,
#'   `sign_check`, `alignment` (an `scr_align` object), `points`,
#'   `base_points`, `samples` (train and hold-out: `link`, `prob`, `score`,
#'   `score_points`, `y`, `date`), `metrics`, `gains`, `stability` (`score`
#'   and `variables`), `calibration`, `rank_order`, `challenger`,
#'   `model_card` and `sql`. Also `scale` (`base_score`, `base_odds`, `pdo`,
#'   `factor`, `offset`, `direction`, `odds_orientation`), `breaks` (the
#'   score bands frozen on train), `monitoring_plan` (see
#'   [scr_monitoring_plan()]), `holdout_bins`, `fit` and `ledger` (the frozen
#'   binning and pre-processing that [scr_apply()] and [scr_sql()]
#'   reproduce) and, after a lab commit, `decisions` and `provenance`.
#'
#' @family stages
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   xgb_rounds = 60, n_boot = 20)
#' res <- scr_select(scr_demo, "default", config = cfg, drop = "id",
#'                   date_col = "ref_date")
#' sc <- scr_scorecard(res)
#' sc
#' head(sc$points[, c("variable", "bin", "woe", "points")])
#' sc$metrics
#' sc$alignment
#' @export
scr_scorecard <- function(x, features = NULL, base_score = NULL, base_odds = NULL, pdo = NULL,
                          direction = NULL, align_method = NULL, challenger = NULL,
                          points_style = NULL, n_boot = NULL, seed = NULL) {
  check_result(x, "scr_scorecard")
  cfg <- x$config
  for (nm in c("base_score", "base_odds", "pdo", "align_method", "challenger", "points_style", "n_boot", "seed")) {
    v <- get(nm); if (!is.null(v)) cfg[[nm]] <- v
  }
  if (!is.null(direction)) cfg$direction <- direction
  cfg <- .scr_validate_config(cfg)
  old <- scr_verbose(isTRUE(cfg$verbose)); on.exit(scr_verbose(old), add = TRUE)
  dir <- resolve_direction(cfg)

  features <- features %||% scr_selected(x)
  missing <- setdiff(features, names(x$fit$results))
  if (length(missing)) stop("variable(s) without a fitted binning: ", lst(missing), call. = FALSE)
  if (!length(features)) stop("no variable for the scorecard.", call. = FALSE)

  target <- x$target
  tr <- x$split$train_idx; ho <- x$split$holdout_idx
  dt_tr <- x$data_clean[tr]; dt_ho <- x$data_clean[ho]
  y_tr <- dt_tr[[target]]; y_ho <- dt_ho[[target]]

  msg_stage("4-5", "scorecard and alignment")
  w_tr <- apply_woe(x$fit, dt_tr, features, "both")
  w_ho <- apply_woe(x$fit, dt_ho, features, "both")

  # -- glm with sign check ------------------------------------------------- #
  fitres <- .fit_glm_signcheck(w_tr, y_tr, features, cfg)
  features <- fitres$features
  if (nrow(fitres$sign_check[action != "kept"])) {
    msg("  sign check removed: %s", lst(fitres$sign_check[action != "kept", variable]))
  }
  msg("  scorecard with %d variables | scale %s/%s/%s | %s", length(features), cfg$base_score,
      cfg$base_odds, cfg$pdo, dir)

  link_tr <- .glm_link(fitres$coef, w_tr, features)
  link_ho <- .glm_link(fitres$coef, w_ho, features)

  # -- alignment (always) ------------------------------------------------- #
  al <- scr_align(link_tr, y_tr, base_score = cfg$base_score, base_odds = cfg$base_odds, pdo = cfg$pdo,
                  direction = dir, method = cfg$align_method, n_bands = cfg$align_bands)
  msg("  alignment: score = %.4f + %.4f * logit (%s%s)", al$a, al$b, al$calibration$method,
      if (is.finite(al$calibration$r2)) sprintf(", adj. R2 %.3f", al$calibration$r2) else "")

  # -- points per bin ------------------------------------------------------ #
  pts <- .build_points(x$fit, features, fitres$coef, al, cfg, w_tr)

  # -- scored samples ------------------------------------------------------ #
  date_vec <- x$date
  samples <- list(
    train   = .score_sample(link_tr, y_tr, w_tr, pts, al, if (!is.null(date_vec)) date_vec[tr] else NULL),
    holdout = .score_sample(link_ho, y_ho, w_ho, pts, al, if (!is.null(date_vec)) date_vec[ho] else NULL))

  # -- metrics with CI, frozen gains, stability ---------------------------- #
  hie <- identical(dir, "higher_is_riskier")
  metrics <- data.table::rbindlist(lapply(names(samples), function(nm) {
    s <- samples[[nm]]
    m <- scr_metrics(s$score, s$y, higher_is_event = hie, ci = TRUE, n_boot = cfg$n_boot,
                     level = cfg$ci_level, seed = cfg$seed, nthread = cfg$nthread)
    data.table::data.table(sample = nm, direction = dir, as.data.frame(m))
  }))
  for (i in seq_len(nrow(metrics))) msg("  %-8s AUC %.4f [%.4f, %.4f]  KS %.4f  Gini %.4f", metrics$sample[i],
                                        metrics$auc[i], metrics$auc_lo[i], metrics$auc_hi[i], metrics$ks[i], metrics$gini[i])
  breaks <- .score_breaks(samples$train$score, cfg$score_groups)
  gains  <- data.table::rbindlist(lapply(names(samples), function(nm)
    data.table::data.table(sample = nm, .score_gains(samples[[nm]]$score, samples[[nm]]$y, breaks, dir))))
  stability <- .scorecard_stability(samples, w_tr, w_ho, pts, breaks, cfg)
  # integer bin index of every hold-out row, per variable: what a CSI timeline
  # by vintage needs, at a fraction of the cost of keeping the labels
  holdout_bins <- data.table::as.data.table(lapply(stats::setNames(features, features), function(f)
    match(w_ho[[paste0(f, "_bin")]], pts$table[variable == f, bin])))
  calibration <- .calibration(samples$holdout, breaks, al)
  rank_order  <- .rank_order(gains[sample == "holdout"])

  # -- challenger ---------------------------------------------------------- #
  chall <- NULL
  if (!is.null(cfg$challenger)) {
    chall <- tryCatch(.fit_challenger(cfg$challenger, w_tr, w_ho, y_tr, y_ho, features, cfg, al, samples),
                      error = function(e) { msg("  challenger not fitted: %s", conditionMessage(e)); NULL })
  }

  fit_sub <- x$fit
  fit_sub$results <- fit_sub$results[features]
  fit_sub$summary <- fit_sub$summary[fit_sub$summary$feature %in% features, , drop = FALSE]
  fit_sub$n_features <- length(features)

  sc <- structure(list(
    target = target, features = features, coef = fitres$coef, sign_check = fitres$sign_check,
    alignment = al, direction = dir, odds_orientation = al$odds_orientation,
    scale = list(base_score = cfg$base_score, base_odds = cfg$base_odds, pdo = cfg$pdo,
                 factor = al$factor, offset = al$offset, direction = dir, odds_orientation = al$odds_orientation),
    points = pts$table, base_points = pts$base_points, base_points_raw = pts$base_points_raw,
    holdout_bins = holdout_bins, monitoring_plan = NULL,
    points_style = cfg$points_style, points_round = cfg$points_round,
    samples = samples, metrics = metrics, breaks = breaks, gains = gains, stability = stability,
    calibration = calibration, rank_order = rank_order, challenger = chall,
    fit = fit_sub, ledger = x$triage$ledger, config = cfg, event = x$meta$event,
    date_col = x$meta$date_col, lab = x$lab, decisions = x$lab$ledger,
    provenance = if (!is.null(x$lab)) x$funnel[feature %in% features, .(feature, provenance, manual_reason)] else NULL,
    sql = NULL, model_card = NULL, files = NULL
  ), class = c("scr_scorecard", "list"))
  sc$sql <- build_sql_score(sc)
  sc$monitoring_plan <- scr_monitoring_plan(sc)
  sc$model_card <- .model_card(sc, x)
  sc
}

# -- fit -------------------------------------------------------------------- #

#' glm on WOE with iterative removal of wrong-signed coefficients
#' @keywords internal
#' @noRd
.fit_glm_signcheck <- function(w_tr, y_tr, features, cfg) {
  log <- list(); feats <- features
  repeat {
    X <- as.data.frame(w_tr[, paste0(feats, "_woe"), with = FALSE]); names(X) <- feats
    X$.y <- y_tr
    m <- suppressWarnings(stats::glm(.y ~ ., family = stats::binomial(), data = X,
                                     model = FALSE, x = FALSE, y = FALSE))
    cf <- stats::coef(m)
    cf[is.na(cf)] <- 0
    b <- cf[feats]
    bad <- feats[b <= 0 | abs(b) > cfg$max_abs_coef]
    if (!length(bad) || length(feats) <= 1L) {
      for (f in feats) log[[f]] <- data.table::data.table(variable = f, coef = unname(b[f]), action = "kept",
                                                          reason = if (b[f] > 0) "OK" else "NON_POSITIVE_COEF_KEPT_LAST")
      break
    }
    worst <- feats[which.min(ifelse(feats %in% bad, b, Inf))]
    if (!worst %in% bad) worst <- bad[1]
    log[[worst]] <- data.table::data.table(variable = worst, coef = unname(b[worst]), action = "removed",
                                           reason = if (b[worst] <= 0) "SIGN_REVERSED" else "COEF_TOO_LARGE")
    feats <- setdiff(feats, worst)
  }
  names(cf)[1] <- "(Intercept)"
  list(coef = cf[c("(Intercept)", feats)], features = feats,
       sign_check = data.table::rbindlist(log)[order(action, variable)])
}

#' Logit from coefficients and WOE columns
#' @keywords internal
#' @noRd
.glm_link <- function(coef, w, features) {
  eta <- rep(unname(coef["(Intercept)"]), nrow(w))
  for (f in features) eta <- eta + unname(coef[f]) * w[[paste0(f, "_woe")]]
  eta
}

#' Points table per bin
#' @keywords internal
#' @noRd
.build_points <- function(fit, features, coef, al, cfg, w_tr) {
  k <- length(features)
  base_raw <- al$a + al$b * unname(coef["(Intercept)"])
  distrib  <- identical(cfg$points_style, "distributed")
  n_tr <- nrow(w_tr)
  rows <- lapply(features, function(f) {
    r <- fit$results[[f]]
    beta <- unname(coef[f])
    p_raw <- al$b * beta * r$woe + if (distrib) base_raw / k else 0
    bl <- w_tr[[paste0(f, "_bin")]]
    cnt <- as.integer(tabulate(match(bl, r$bin), nbins = length(r$bin)))
    data.table::data.table(variable = f, bin_id = as.integer(r$id), bin = r$bin, woe = r$woe,
                           count_train = cnt, pct_train = cnt / n_tr, count_fit = r$count,
                           pos_rate = r$count_pos / pmax(1L, r$count), iv = r$iv, coef = beta,
                           points_raw = p_raw,
                           points = if (isTRUE(cfg$points_round)) round(p_raw) else p_raw)
  })
  tb <- data.table::rbindlist(rows)
  base_points_raw <- if (distrib) 0 else base_raw
  list(table = tb[], base_points_raw = base_points_raw,
       base_points = if (isTRUE(cfg$points_round)) round(base_points_raw) else base_points_raw,
       unbinned_points = .unbinned_points(distrib, base_raw, k, cfg$points_round))
}

#' Points of a row that falls in no fitted bin (an unseen category, a
#' missing value without a missing bin)
#'
#' The binning engine gives such a row a WOE of 0 (its `na_woe`), so its
#' exact contribution to the score is the distributed share of the base,
#' `base / k`, under `points_style = "distributed"` and 0 otherwise. The
#' whole-points score must add the same, in R and in SQL, or `score` and
#' `score_points` part ways on exactly those rows.
#' @keywords internal
#' @noRd
.unbinned_points <- function(distrib, base_raw, k, round_pts) {
  v <- if (isTRUE(distrib)) base_raw / k else 0
  if (isTRUE(round_pts)) round(v) else v
}

#' The same value, read from a fitted scorecard
#' @keywords internal
#' @noRd
.sc_unbinned_points <- function(sc) {
  .unbinned_points(identical(sc$points_style, "distributed"),
                   sc$alignment$a + sc$alignment$b * unname(sc$coef["(Intercept)"]),
                   length(sc$features), sc$points_round)
}

#' Score a sample: link, prob, exact score and whole-points score
#' @keywords internal
#' @noRd
.score_sample <- function(link, y, w, pts, al, date = NULL) {
  score <- al$a + al$b * link
  sp <- rep(pts$base_points, length(link))
  for (f in unique(pts$table$variable)) {
    p <- pts$table[variable == f]
    pf <- p$points[match(w[[paste0(f, "_bin")]], p$bin)]
    pf[is.na(pf)] <- pts$unbinned_points
    sp <- sp + pf
  }
  out <- data.table::data.table(link = link, prob = stats::plogis(link), score = score, score_points = sp, y = y)
  if (!is.null(date)) out[, date := date]
  out
}

# -- gains, stability, calibration, rank-order ------------------------------ #

#' @keywords internal
#' @noRd
.score_breaks <- function(score, n_groups) {
  probs <- seq(0, 1, length.out = n_groups + 1L)[-c(1L, n_groups + 1L)]
  unique(c(-Inf, stats::quantile(score, probs = probs, na.rm = TRUE, names = FALSE), Inf))
}

#' Score gains per frozen band, from the risky side to the safe side
#' @keywords internal
#' @noRd
.score_gains <- function(score, y, breaks, direction) {
  band <- cut(score, breaks = breaks, include.lowest = TRUE)
  d <- data.table::data.table(band = band, score = score, y = as.integer(y))[
    , .(n = .N, events = sum(y), min_score = min(score), mean_score = mean(score), max_score = max(score)),
    by = band]
  # riskiest band first: low score under higher_is_safer, high under higher_is_riskier
  d <- if (identical(direction, "higher_is_safer")) d[order(band)] else d[order(-as.integer(band))]
  n_tot <- sum(d$n); e_tot <- sum(d$events); ne_tot <- n_tot - e_tot
  d[, `:=`(id = seq_len(.N), pct = n / n_tot, event_rate = events / n, non_events = n - events)]
  d[, `:=`(cum_pct = cumsum(pct), cum_event_pct = cumsum(events) / max(1, e_tot),
           cum_nonevent_pct = cumsum(non_events) / max(1, ne_tot))]
  d[, `:=`(ks = abs(cum_event_pct - cum_nonevent_pct),
           lift = event_rate / (e_tot / n_tot),
           cum_lift = (cumsum(events) / cumsum(n)) / (e_tot / n_tot),
           odds = (n - events + 0.5) / (events + 0.5))]
  d[, log_odds := log(odds)]
  data.table::setcolorder(d, c("id", "band", "n", "pct", "events", "non_events", "event_rate",
                               "min_score", "mean_score", "max_score", "cum_pct", "cum_event_pct",
                               "cum_nonevent_pct", "ks", "lift", "cum_lift", "odds", "log_odds"))
  d[, band := as.character(band)]
  d[]
}

#' Score PSI and per-variable CSI (with the signed points shift)
#' @keywords internal
#' @noRd
.scorecard_stability <- function(samples, w_tr, w_ho, pts, breaks, cfg) {
  ps <- scr_psi(samples$train$score, samples$holdout$score, breaks = breaks, alpha = cfg$psi_alpha)
  score_tb <- data.table::data.table(sample = "holdout", psi = ps$psi, flag_fixed = ps$flag_fixed,
                                     critical = ps$critical, flag_adjusted = ps$flag_adjusted,
                                     n_base = ps$n_base, n_compare = ps$n_compare)
  vars <- unique(pts$table$variable)
  csi <- data.table::rbindlist(.scr_lapply(vars, function(f) {
    p  <- pts$table[variable == f]
    cb <- paste0(f, "_bin")
    r  <- scr_psi(w_tr[[cb]], w_ho[[cb]], levels = p$bin, alpha = cfg$psi_alpha)
    sh <- .points_shift(r$table$pct_base, r$table$pct_compare, p$points[match(r$table$band, p$bin)])
    data.table::data.table(variable = f, sample = "holdout", csi = r$psi, flag_fixed = r$flag_fixed,
                           critical = r$critical, flag_adjusted = r$flag_adjusted, points_shift = sh)
  }, nthread = cfg$nthread))
  list(score = score_tb, variables = csi[order(-csi)], score_psi = ps)
}

#' Hold-out calibration: intercept/slope, Brier, ECE per frozen band
#' @keywords internal
#' @noRd
.calibration <- function(s, breaks, al) {
  p <- .score_to_prob(al, s$score)
  y <- s$y
  band <- cut(s$score, breaks = breaks, include.lowest = TRUE)
  # ordered by the band factor (numeric order), not by its label: "(1e+03,Inf]"
  # would otherwise sort before "(512,530]"
  tb <- data.table::data.table(band = band, p = p, y = y)[
    , .(n = .N, expected = mean(p), observed = mean(y)), by = band][order(band)]
  tb[, band := as.character(band)]
  tb[, gap := observed - expected]
  ece <- sum(tb$n / sum(tb$n) * abs(tb$gap))
  lo <- suppressWarnings(stats::glm(y ~ stats::qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6)), family = stats::binomial()))
  cf <- stats::coef(lo)
  list(summary = data.table::data.table(sample = "holdout", n = length(y), brier = mean((p - y)^2),
                                        ece = ece, mce = max(abs(tb$gap)),
                                        intercept = unname(cf[1]), slope = unname(cf[2]),
                                        expected_rate = mean(p), observed_rate = mean(y)),
       table = tb[])
}

#' Rank-order breaks between adjacent bands, with a binomial test
#' @keywords internal
#' @noRd
.rank_order <- function(g) {
  g <- data.table::copy(g)[order(id)]
  g[, `:=`(prev_rate = data.table::shift(event_rate), monotone = NA, p_value = NA_real_)]
  for (i in seq_len(nrow(g))[-1]) {
    # expected: the event rate does not rise while moving from the risky to the safe side
    g$monotone[i] <- g$event_rate[i] <= g$prev_rate[i]
    # P(observe >= events | rate of the previous band): a small value is a significant break
    g$p_value[i] <- stats::pbinom(g$events[i] - 1L, g$n[i], g$prev_rate[i], lower.tail = FALSE)
  }
  g[, break_flag := !(monotone %in% TRUE) & p_value < 0.05]
  g[, .(id, band, n, events, event_rate, prev_rate, monotone, p_value, break_flag)]
}

# -- challenger ------------------------------------------------------------- #

#' Tree challenger on the same WOE columns, aligned to the same scale
#' @keywords internal
#' @noRd
.fit_challenger <- function(engine, w_tr, w_ho, y_tr, y_ho, features, cfg, al, samples) {
  cols <- paste0(features, "_woe")
  x_tr <- as.matrix(w_tr[, cols, with = FALSE]); x_ho <- as.matrix(w_ho[, cols, with = FALSE])
  colnames(x_tr) <- colnames(x_ho) <- features
  fitf <- switch(engine, xgboost = .fit_xgboost, lightgbm = .fit_lightgbm)
  r <- fitf(x_tr, y_tr, x_ho, y_ho, cfg)
  p_ho <- pmin(pmax(r$score, 1e-6), 1 - 1e-6)
  p_tr <- pmin(pmax(as.numeric(stats::predict(r$model, if (engine == "xgboost") xgboost::xgb.DMatrix(x_tr) else x_tr)), 1e-6), 1 - 1e-6)
  raw_tr <- stats::qlogis(p_tr); raw_ho <- stats::qlogis(p_ho)
  # Aligned to the SAME declared scale by the same mechanism: this is what
  # makes the challenger comparable point for point with the scorecard.
  al_c <- scr_align(raw_tr, y_tr, base_score = al$base_score, base_odds = al$base_odds, pdo = al$pdo,
                    direction = al$direction, method = cfg$align_method, n_bands = cfg$align_bands)
  sc_ho <- predict(al_c, raw_ho)
  hie <- identical(al$direction, "higher_is_riskier")
  m <- scr_metrics(sc_ho, y_ho, higher_is_event = hie, ci = TRUE, n_boot = cfg$n_boot,
                   level = cfg$ci_level, seed = cfg$seed, nthread = cfg$nthread)
  msg("  challenger %s: hold-out AUC %.4f [%.4f, %.4f] (scorecard: %.4f)", engine, m$auc, m$auc_lo, m$auc_hi,
      scr_metrics(samples$holdout$score, y_ho, higher_is_event = hie, ci = FALSE)$auc)
  swap <- .swapset(samples$holdout$score, sc_ho, y_ho, al$direction, c(0.5, 0.7, 0.9))
  list(engine = engine, supports_scorecard = FALSE, points = "NOT_APPLICABLE_ENGINE",
       reason_codes = "NOT_APPLICABLE_ENGINE", alignment = al_c,
       importance = r$importance[order(-importance)], note = r$note,
       metrics = data.table::data.table(sample = "holdout", direction = al$direction, as.data.frame(m)),
       score_holdout = sc_ho, swapset = swap)
}

#' Champion vs challenger swap set at a constant approval rate
#' @keywords internal
#' @noRd
.swapset <- function(champ, chall, y, direction, rates) {
  safe_side <- function(s, rate) {
    if (identical(direction, "higher_is_safer")) s >= stats::quantile(s, 1 - rate, names = FALSE)
    else s <= stats::quantile(s, rate, names = FALSE)
  }
  data.table::rbindlist(lapply(rates, function(rt) {
    a <- safe_side(champ, rt); b <- safe_side(chall, rt)
    swap_in <- !a & b; swap_out <- a & !b
    data.table::data.table(approval_rate = rt, n_swap_in = sum(swap_in), n_swap_out = sum(swap_out),
                           event_rate_swap_in = if (any(swap_in)) mean(y[swap_in]) else NA_real_,
                           event_rate_swap_out = if (any(swap_out)) mean(y[swap_out]) else NA_real_,
                           event_rate_champion = mean(y[a]), event_rate_challenger = mean(y[b]))
  }))
}

#' @keywords internal
#' @noRd
.model_card <- function(sc, x) {
  m <- x$meta
  list(
    package = sprintf("scorecraft %s", as.character(utils::packageVersion("scorecraft"))),
    fitted_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    target = sc$target, event_level = sc$event$label, event_inverted = isTRUE(sc$event$inverted),
    objective = x$config$objective, direction = sc$direction, odds_orientation = sc$odds_orientation,
    base_score = sc$scale$base_score, base_odds = sc$scale$base_odds, pdo = sc$scale$pdo,
    factor = sc$scale$factor, offset = sc$scale$offset,
    align_method = sc$alignment$calibration$method, align_intercept = sc$alignment$calibration$intercept,
    align_slope = sc$alignment$calibration$slope, align_r2 = sc$alignment$calibration$r2,
    n_train = m$n_train, n_holdout = m$n_holdout, split_method = m$split_method, split_cutoff = m$split_cutoff,
    event_rate_train = m$event_rate_train, event_rate_holdout = m$event_rate_holdout,
    n_features = length(sc$features), features = paste(sc$features, collapse = ", "),
    binning_algorithm = paste(unique(sc$fit$summary$algorithm), collapse = ", "),
    points_style = sc$points_style, points_round = sc$points_round,
    shortlist_source = if (is.null(x$lab) || is.null(x$lab$shortlist$manual)) "consensus" else "manual",
    n_manual_bins = if (is.null(x$lab)) 0L else sum(x$lab$source[sc$features] == "manual", na.rm = TRUE),
    manual_bins = if (is.null(x$lab)) "" else paste(intersect(names(x$lab$source)[x$lab$source == "manual"], sc$features), collapse = ", "),
    forced_in = if (is.null(x$lab)) "" else paste(setdiff(x$lab$shortlist$final, x$lab$shortlist$consensus), collapse = ", "),
    manual_dropped = if (is.null(x$lab)) "" else paste(setdiff(x$lab$shortlist$consensus, x$lab$shortlist$final), collapse = ", "),
    n_decisions = if (is.null(x$lab)) 0L else nrow(x$lab$ledger),
    manual_removed_by_sign_check = if (is.null(x$lab)) "" else paste(intersect(sc$sign_check[action == "removed", variable],
                                                                                setdiff(x$lab$shortlist$final, x$lab$shortlist$consensus)), collapse = ", "),
    challenger = if (is.null(sc$challenger)) "none" else sc$challenger$engine,
    challenger_supports_scorecard = if (is.null(sc$challenger)) NA else FALSE,
    auc_holdout = sc$metrics[sample == "holdout", auc], ks_holdout = sc$metrics[sample == "holdout", ks],
    psi_score_holdout = sc$stability$score$psi
  )
}

#' @export
print.scr_scorecard <- function(x, ...) {
  cat(sprintf("<scr_scorecard> target \"%s\" | %d variables | %s\n", x$target, length(x$features), x$direction))
  cat(sprintf("  scale: %s points at odds %s:1 (%s), PDO %s | alignment %s\n", format(x$scale$base_score),
              format(x$scale$base_odds), x$odds_orientation, format(x$scale$pdo), x$alignment$calibration$method))
  cat(sprintf("  score = %.4f + %.4f * logit | base_points = %s\n", x$alignment$a, x$alignment$b, format(x$base_points)))
  for (i in seq_len(nrow(x$metrics))) {
    m <- x$metrics[i]
    cat(sprintf("  %-8s n %-7s AUC %.4f [%.4f, %.4f]  KS %.4f  Gini %.4f\n", m$sample, n_fmt(m$n), m$auc, m$auc_lo, m$auc_hi, m$ks, m$gini))
  }
  s <- x$stability$score
  cat(sprintf("  score PSI (hold-out): %.4f - fixed: %s | adjusted (%.4f): %s\n", s$psi, s$flag_fixed, s$critical, s$flag_adjusted))
  rm <- x$sign_check[action == "removed"]
  if (nrow(rm)) cat(sprintf("  removed by the sign check: %s\n", lst(rm$variable)))
  if (!is.null(x$challenger)) {
    cm <- x$challenger$metrics
    cat(sprintf("  challenger %s (supports_scorecard = FALSE): AUC %.4f [%.4f, %.4f]\n",
                x$challenger$engine, cm$auc, cm$auc_lo, cm$auc_hi))
  }
  cat("\nPoints (first rows)\n")
  p <- utils::head(x$points, 8)
  for (i in seq_len(nrow(p))) cat(sprintf("  %-28s %-26s %8.3f %7s\n", p$variable[i], substr(p$bin[i], 1, 26), p$woe[i], format(p$points[i])))
  if (nrow(x$points) > 8) cat(sprintf("  ... (+%d rows)\n", nrow(x$points) - 8))
  invisible(x)
}

#' Score gains per frozen band
#'
#' How the score behaves in each band: count, event rate, KS, lift,
#' cumulative capture and the score interval of the band, which is what lets
#' a cut-off be read straight from the table. The bands are the deciles of
#' the score on **train**, applied frozen to the other samples.
#'
#' @param x An object from [scr_scorecard()].
#' @param sample `NULL` (all), `"train"` or `"holdout"`.
#'
#' @return A `data.table` with one row per sample and band, from the riskiest
#'   band to the safest.
#'
#' @family accessors
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   xgb_rounds = 60, n_boot = 20)
#' res <- scr_select(scr_demo, "default", config = cfg, drop = "id",
#'                   date_col = "ref_date")
#' sc <- scr_scorecard(res)
#' scr_score_gains(sc, "holdout")[, .(band, n, event_rate, min_score, max_score, ks)]
#' scr_score_metrics(sc)
#' @export
scr_score_gains <- function(x, sample = NULL) {
  if (!inherits(x, "scr_scorecard")) stop("scr_score_gains() expects an object from scr_scorecard().", call. = FALSE)
  g <- x$gains
  if (!is.null(sample)) {
    idx <- which(g[["sample"]] %in% sample)
    g <- g[idx]
  }
  g[]
}

#' Score metrics per sample, with CI
#'
#' `n`, events, AUC, KS and Gini of the scorecard score on train and
#' hold-out, with a bootstrap confidence interval and the direction used:
#' the AUC is always reported above 0.5 when the score ranks correctly in
#' its own direction.
#'
#' @param x An object from [scr_scorecard()].
#'
#' @return A `data.table` with one row per sample.
#'
#' @family accessors
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   xgb_rounds = 60, n_boot = 20)
#' res <- scr_select(scr_demo, "default", config = cfg, drop = "id",
#'                   date_col = "ref_date")
#' sc <- scr_scorecard(res)
#' scr_score_metrics(sc)
#' @export
scr_score_metrics <- function(x) {
  if (!inherits(x, "scr_scorecard")) stop("scr_score_metrics() expects an object from scr_scorecard().", call. = FALSE)
  x$metrics[]
}
