# Shared fixtures. The full pipeline on scr_demo takes a few seconds; the
# result and the scorecard are fitted once per session and reused.

.fx <- new.env(parent = emptyenv())

cfg_test <- function(...) {
  base <- list(verbose = FALSE, nthread = 1L, use_ranger = FALSE, use_lightgbm = FALSE,
               xgb_rounds = 40L, n_boot = 10L)
  do.call(scr_config, utils::modifyList(base, list(...)))
}

res_demo <- function() {
  if (is.null(.fx$res)) {
    .fx$res <- scr_select(scr_demo, "default", config = cfg_test(), drop = c("id", "churn"),
                          date_col = "ref_date")
  }
  .fx$res
}

sc_demo <- function() {
  if (is.null(.fx$sc)) .fx$sc <- scr_scorecard(res_demo())
  .fx$sc
}

noise_names  <- function() grep("^vl_noise_", names(scr_demo), value = TRUE)
signal_names <- function() grep("^vl_score_", names(scr_demo), value = TRUE)

# fraud-like: the same model on a mirrored scale (higher score, more risk)
sc_fraud_demo <- function() {
  if (is.null(.fx$sc_fraud)) .fx$sc_fraud <- scr_scorecard(res_demo(), direction = "higher_is_riskier")
  .fx$sc_fraud
}

# propensity: churn is the desirable event, more points = more likely
sc_prop_demo <- function() {
  if (is.null(.fx$sc_prop)) {
    res <- scr_select(scr_demo, "churn", config = cfg_test(objective = "propensity"),
                      drop = c("id", "default"), date_col = "ref_date")
    .fx$sc_prop <- scr_scorecard(res)
  }
  .fx$sc_prop
}

# DeLong standard error of an AUC from the mid-ranks of the rows: the
# row-level reference for the count version used by the package
ref_auc_se <- function(score, y, higher_is_event = TRUE) {
  ok <- !is.na(score) & !is.na(y)
  s <- as.double(score[ok]); y <- as.integer(y[ok])
  if (!higher_is_event) s <- -s
  e <- y == 1L; n1 <- sum(e); n0 <- sum(!e)
  if (n1 < 2L || n0 < 2L) return(NA_real_)
  r <- data.table::frank(s, ties.method = "average")
  v10 <- (r[e] - data.table::frank(s[e], ties.method = "average")) / n0
  v01 <- 1 - (r[!e] - data.table::frank(s[!e], ties.method = "average")) / n1
  sqrt(stats::var(v10) / n1 + stats::var(v01) / n0)
}
