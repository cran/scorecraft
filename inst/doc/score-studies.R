## ----include = FALSE----------------------------------------------------------
has_glmnet <- requireNamespace("glmnet", quietly = TRUE)
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", eval = has_glmnet)
old_dt <- data.table::setDTthreads(2)

## ----setup, message = FALSE---------------------------------------------------
library(scorecraft)
library(data.table)
cfg <- scr_config(verbose = FALSE, nthread = 1, use_glmnet = TRUE, use_ranger = FALSE,
                  use_lightgbm = FALSE, xgb_rounds = 60, n_boot = 20)
res <- scr_select(scr_demo, "default", config = cfg, drop = c("id", "churn"),
                  date_col = "ref_date")
credit <- scr_scorecard(res)
fraud <- scr_scorecard(res, direction = "higher_is_riskier")

cfg_churn <- scr_config(objective = "propensity", verbose = FALSE, nthread = 1, use_glmnet = TRUE,
                        use_ranger = FALSE, use_lightgbm = FALSE, xgb_rounds = 60, n_boot = 20)
res_churn <- scr_select(scr_demo, "churn", config = cfg_churn, drop = c("id", "default"),
                        date_col = "ref_date")
churn <- scr_scorecard(res_churn)

## ----bands--------------------------------------------------------------------
b <- scr_bands(credit, n_bands = 10, n_boot = 50, seed = 1)
b

## ----tiers--------------------------------------------------------------------
tr <- scr_tiers(credit, n_tiers = 5, round_to = 5, n_boot = 30, seed = 1)
tr
tr$stability$cuts

## ----include = FALSE----------------------------------------------------------
mono <- tr$summary[sample == "holdout", monotone]

## ----rag----------------------------------------------------------------------
rg <- scr_rag(credit, n_boot = 50, seed = 1)
rg$summary

## ----fraud-bands--------------------------------------------------------------
bf <- scr_bands(fraud, spacing = "tail", tail_probs = c(0.02, 0.05, 0.10, 0.20, 0.50), n_boot = 0)
bf$table[sample == "holdout", .(band, label, pct, rate, lift, capture)]

## ----fraud-tiers--------------------------------------------------------------
tf <- scr_tiers(fraud, n_tiers = 3, labels = c("pass", "review", "block"))
tf$table[sample == "holdout", .(tier, label, score_lo, score_hi, pct, rate)]

## ----fraud-operating----------------------------------------------------------
d_fraud <- data.frame(score = scr_apply(fraud, scr_demo)$score, y = scr_demo$default,
                      month = scr_demo$ref_date)
op_f <- scr_operating(d_fraud, objective = "risk", direction = "higher_is_riskier",
                      max_per_day = 80, date = "month")
op_f

## ----churn-tiers--------------------------------------------------------------
ta <- scr_tiers(churn, method = "anchored", anchors = c(0.15, "overall", 0.5))
ta

## ----claims-------------------------------------------------------------------
claims <- data.frame(
  name = c("high tier churns", "top of the score", "low tier is quiet"),
  label = c("high", NA, "low"), score_lo = c(NA, 480, NA),
  op = c(">=", ">=", "<="), rate = c(0.50, 0.60, 0.15))
cl <- scr_claims(ta, claims)
cl

## ----claims-floor-------------------------------------------------------------
fl <- scr_claims(ta, claims, type = "floor")
fl$table[, .(name, group, score_lo, score_hi, n, rate, bound, verdict)]

## ----churn-operating----------------------------------------------------------
op_c <- scr_operating(churn, gain_event = 100, cost_select = 20, budget = 3000)
op_c$optimum[, .(cut, depth, n_sel, rate_sel, value, binding, shadow_price)]

## ----cross--------------------------------------------------------------------
ho <- scr_demo[res$split$holdout_idx, ]
d_x <- data.frame(credit = scr_apply(credit, ho)$score, churn = scr_apply(churn, ho)$score,
                  default = ho$default, churned = ho$churn)
cx <- scr_score_cross(d_x, "credit", "churn", y_a = "default", y_b = "churned",
                      objective_b = "propensity", cuts_a = tr, cuts_b = ta)
cx

## ----mix-shift----------------------------------------------------------------
ms <- scr_mix_shift(credit, n_bands = 5)
ms

## ----mix-shift-by-------------------------------------------------------------
scr_mix_shift(credit, by = "date", n_bands = 5)$summary[
  , .(group, rate_base, rate_cmp, delta, mix_total, rate_total, psi)]

## ----segments-----------------------------------------------------------------
sg <- scr_segments(credit, ho, segment = "ds_channel")
sg

## ----production---------------------------------------------------------------
head(scr_apply(tr, c(495, 520, 560)))
sql <- scr_sql(tr, table = "scored", dialect = "postgres")
substr(sql[9], 1, 80)
substr(sql[10], 1, 80)

## ----include = FALSE----------------------------------------------------------
data.table::setDTthreads(old_dt)

