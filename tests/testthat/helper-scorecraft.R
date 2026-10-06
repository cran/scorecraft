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
