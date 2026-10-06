## ----include = FALSE----------------------------------------------------------
has_glmnet <- requireNamespace("glmnet", quietly = TRUE)
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
old_dt <- data.table::setDTthreads(2)

## ----setup, message = FALSE---------------------------------------------------
library(scorecraft)
library(data.table)

## ----constants----------------------------------------------------------------
factor <- 20 / log(2)
offset <- 600 - factor * log(50)
c(factor = factor, offset = offset)

## ----simulate-----------------------------------------------------------------
set.seed(2026)
n   <- 4000
raw <- stats::rnorm(n, mean = stats::qlogis(0.12), sd = 1.2)   # event logit
y   <- stats::rbinom(n, 1, stats::plogis(raw))                 # calibrated outcome

## ----direct-------------------------------------------------------------------
al_risk <- scr_align(raw, y, base_score = 600, base_odds = 50, pdo = 20,
                     direction = "higher_is_safer", method = "direct")
al_risk

## ----trap---------------------------------------------------------------------
al_prop_wrong <- scr_align(raw, y, base_score = 600, base_odds = 50, pdo = 20,
                           direction = "higher_is_riskier", method = "direct")
al_prop_right <- scr_align(raw, y, base_score = 600, base_odds = 1 / 50, pdo = 20,
                           direction = "higher_is_riskier", method = "direct")

# what does "600 points" mean under each object? Invert score = a + b * raw
prob_at <- function(al, score) predict(al, (score - al$a) / al$b, type = "prob")
data.table(
  object      = c("risk 50 (safe:event)", "propensity 50 (event:safe)", "propensity 1/50 (event:safe)"),
  orientation = c(al_risk$odds_orientation, al_prop_wrong$odds_orientation, al_prop_right$odds_orientation),
  p_event_at_600 = round(c(prob_at(al_risk, 600), prob_at(al_prop_wrong, 600), prob_at(al_prop_right, 600)), 4)
)

## ----mirror-------------------------------------------------------------------
head(data.table(
  raw   = round(raw, 3),
  p_true = round(stats::plogis(raw), 4),
  score_risk = round(predict(al_risk, raw), 1),
  score_prop = round(predict(al_prop_right, raw), 1),
  p_risk = round(predict(al_risk, raw, type = "prob"), 4),
  p_prop = round(predict(al_prop_right, raw, type = "prob"), 4)
), 5)
all.equal(predict(al_risk, raw) + predict(al_prop_right, raw), rep(1200, n))

## ----regression---------------------------------------------------------------
al_reg <- scr_align(raw, y, base_score = 600, base_odds = 50, pdo = 20,
                    direction = "higher_is_safer")
al_reg

## ----bands--------------------------------------------------------------------
al_reg$calibration$bands
c(intercept = al_reg$calibration$intercept, slope = al_reg$calibration$slope,
  adj_r2 = al_reg$calibration$r2)

## ----miscalibrated------------------------------------------------------------
raw_bad <- 0.6 * raw - 1            # over-timid and shifted
al_bad_direct <- scr_align(raw_bad, y, method = "direct")
al_bad_reg    <- scr_align(raw_bad, y)
al_bad_reg$calibration[c("intercept", "slope", "r2")]
data.table(
  method = c("direct", "regression"),
  mean_predicted = c(mean(predict(al_bad_direct, raw_bad, type = "prob")),
                     mean(predict(al_bad_reg,    raw_bad, type = "prob"))),
  observed = mean(y)
)

## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(eval = has_glmnet)

## ----select-two, results = "hide"---------------------------------------------
cfg_risk <- scr_config(objective = "risk", verbose = FALSE, nthread = 1, use_glmnet = TRUE,
                       use_ranger = FALSE, use_lightgbm = FALSE, xgb_rounds = 60, n_boot = 20)
cfg_prop <- scr_config(objective = "propensity", verbose = FALSE, nthread = 1, use_glmnet = TRUE,
                       use_ranger = FALSE, use_lightgbm = FALSE, xgb_rounds = 60, n_boot = 20)

res_risk <- scr_select(scr_demo, "default", config = cfg_risk, drop = c("id", "churn"),
                       date_col = "ref_date")
res_prop <- scr_select(scr_demo, "churn", config = cfg_prop, drop = c("id", "default"),
                       date_col = "ref_date")

sc_risk <- scr_scorecard(res_risk)
sc_prop <- scr_scorecard(res_prop)

## ----two-alignments-----------------------------------------------------------
sc_risk$alignment
sc_prop$alignment
c(risk = sc_risk$odds_orientation, propensity = sc_prop$odds_orientation)

## ----two-metrics--------------------------------------------------------------
scr_score_metrics(sc_risk)[, .(sample, direction, auc, auc_lo, auc_hi, ks)]
scr_score_metrics(sc_prop)[, .(sample, direction, auc, auc_lo, auc_hi, ks)]

## ----mean-by-outcome----------------------------------------------------------
rbind(
  sc_risk$samples$holdout[, .(scorecard = "default (risk)", n = .N, mean_score = round(mean(score), 1)), by = y],
  sc_prop$samples$holdout[, .(scorecard = "churn (propensity)", n = .N, mean_score = round(mean(score), 1)), by = y]
)[order(scorecard, y)]

## ----wrong-way----------------------------------------------------------------
s <- sc_risk$samples$holdout
c(read_as_safer  = scr_metrics(s$score, s$y, higher_is_event = FALSE, ci = FALSE)$auc,
  read_as_riskier = scr_metrics(s$score, s$y, higher_is_event = TRUE,  ci = FALSE)$auc)

## ----points-------------------------------------------------------------------
al <- sc_risk$alignment
c(a = al$a, b = al$b, alpha = unname(sc_risk$coef["(Intercept)"]),
  base_exact = sc_risk$base_points_raw, base_points = sc_risk$base_points)
sc_risk$points[variable == "vl_score_01", .(bin, woe, coef, points_raw = round(points_raw, 2), points)]

## ----score-points-------------------------------------------------------------
sc_new <- scr_apply(sc_risk, head(scr_demo, 5), what = "points")
var_pts <- sc_new[, paste0(sc_risk$features, "_points"), with = FALSE]
sc_new[, .(score = round(score, 2), score_points,
           check = sc_risk$base_points + rowSums(var_pts), diff = round(score - score_points, 2))]

## ----odds---------------------------------------------------------------------
g <- scr_score_gains(sc_risk)
g[, scale_odds := exp((mean_score - al$offset) / al$factor)]
g[sample == "holdout", .(band, n, events, mean_score = round(mean_score, 1),
                         odds = round(odds, 2), scale_odds = round(scale_odds, 2))]
g[, .(pdo_implied = log(2) / coef(stats::lm(log_odds ~ mean_score, weights = n))[[2]]), by = sample]

## ----challenger, results = "hide"---------------------------------------------
sc_x <- scr_scorecard(res_risk, challenger = "xgboost")

## ----challenger-show----------------------------------------------------------
sc_x$challenger$supports_scorecard
sc_x$challenger$points
sc_x$challenger$alignment
rbind(sc_x$challenger$metrics[, .(model = "xgboost challenger", sample, auc, auc_lo, auc_hi, ks)],
      scr_score_metrics(sc_x)[sample == "holdout", .(model = "scorecard", sample, auc, auc_lo, auc_hi, ks)])

## ----swapset------------------------------------------------------------------
sc_x$challenger$swapset

## ----portfolio----------------------------------------------------------------
runs <- list(default = res_risk, churn = res_prop)
scr_compare(runs)[, .(target, rows, event_rate, approved, best_model, auc, ks, relaxation)]
scr_core(runs, min_targets = 2)

## ----rescale, results = "hide"------------------------------------------------
sc_700 <- scr_scorecard(res_risk, base_score = 700, base_odds = 30, pdo = 25)

## ----rescale-show-------------------------------------------------------------
sc_700$alignment
# the calibration is the same; only factor and offset change
c(intercept_600 = sc_risk$alignment$calibration$intercept, intercept_700 = sc_700$alignment$calibration$intercept,
  slope_600 = sc_risk$alignment$calibration$slope, slope_700 = sc_700$alignment$calibration$slope)

## ----rescale-metrics----------------------------------------------------------
rbind(scr_score_metrics(sc_risk)[, .(scale = "600/50/20", sample, auc, ks, gini)],
      scr_score_metrics(sc_700)[, .(scale = "700/30/25", sample, auc, ks, gini)])
merge(sc_risk$points[, .(variable, bin, points_600 = points)],
      sc_700$points[, .(variable, bin, points_700 = points)],
      by = c("variable", "bin"))[1:6]

## ----rescale-identity---------------------------------------------------------
a600 <- sc_risk$alignment; a700 <- sc_700$alignment
ln_odds <- (sc_risk$samples$holdout$score - a600$offset) / a600$factor
all.equal(a700$offset + a700$factor * ln_odds, sc_700$samples$holdout$score)

## ----model-card---------------------------------------------------------------
mc <- sc_risk$model_card
keys <- c("target", "objective", "direction", "odds_orientation",
          "base_score", "base_odds", "pdo", "factor", "offset",
          "align_method", "align_intercept", "align_slope", "align_r2",
          "n_train", "event_rate_train", "split_method", "split_cutoff",
          "challenger", "challenger_supports_scorecard", "auc_holdout", "ks_holdout")
data.table(field = keys, value = vapply(mc[keys], function(v) format(v, digits = 6), character(1)))

## ----include = FALSE, eval = TRUE---------------------------------------------
data.table::setDTthreads(old_dt)

