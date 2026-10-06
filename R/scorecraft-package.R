#' scorecraft: scorecard engine with alignment, cut-off strategy, IRB risk parameters and production SQL
#'
#' A professional scorecard is born of eight chained stages, numbered 0 to 7
#' in every message and in [scr_config_keys()], and this package exposes each
#' of them as a function of its own, next to the shortcut [scr_select()] that
#' chains them for the common case:
#'
#' \enumerate{
#'   \item **Split** ([scr_split()]): train/hold-out by whole periods
#'     (out-of-time) before any supervised fit.
#'   \item **Triage** ([scr_triage()]): structural filters, decomposition of
#'     sentinels and missing values, exact duplicates. The data leaves with no
#'     `NA`.
#'   \item **Binning and screening** ([scr_bin()]): optimal bins parallelised
#'     by column, eight admission rules, hold-out revalidation with frozen
#'     bins, redundancy pruning.
#'   \item **Multi-strategy selection** ([scr_model()]): elastic net, boosting
#'     and random forest on the WOE space; consensus weighted by hold-out Gini.
#'   \item **Scorecard** ([scr_scorecard()]): logistic regression on the
#'     shortlist, sign check, points per bin, a tree challenger explicitly
#'     without points.
#'   \item **Alignment** ([scr_align()]): log-odds regression on the raw score
#'     composed with the PDO map, with `odds_orientation` recorded. Runs
#'     automatically inside [scr_scorecard()].
#'   \item **Cut-off and strategy** ([scr_cutoff()], [scr_strategy()],
#'     [scr_reject()]): sweep with frozen cuts, bands with marginal expected
#'     profit, honest reject inference through a sensitivity band.
#' }
#'
#' The deliverables ([scr_export()]) are the audit funnel, the gains tables,
#' the production SQL ([scr_sql()]) with R-SQL equivalence verified by test,
#' and four `.xlsx` workbooks. [scr_monitor()] recomputes PSI/CSI on new data,
#' with both the fixed and the sample-size-adjusted threshold, and never
#' schedules anything by itself.
#'
#' @section IRB risk parameters:
#'
#' The IRB (internal ratings-based) layer, stages 8 to 12 of
#' [scr_config_keys()], turns the scorecard into regulatory parameters and
#' keeps the same contracts (one configuration, ledgers, hold-out
#' revalidation, workbooks, production SQL). [scr_irb_params()] holds every regime-specific
#' number as a table selected by preset (`"bcb"`, `"basel3_final"`,
#' `"crr3"`); [scr_default()] builds the default flag from a monthly panel and
#' [scr_default_rate()] the default rates by cohort with the long-run average.
#' PD: [scr_calibrate()] anchors the scorecard to a central tendency,
#' [scr_grades()] cuts the score into rating grades, [scr_moc()] and
#' [scr_pd()] add the margin of conservatism and the floor, and
#' [scr_pd_validate()] runs the calibration, discrimination and stability
#' tests with traffic lights; [scr_master_scale()], [scr_migration()] and
#' [scr_pd_pit_ttc()] support the grade structure, the migration analysis
#' and the point-in-time bridge. LGD: [scr_workout()] discounts recovery cash
#' flows into realised LGD, [scr_lgd()] fits the cure and severity stages
#' and the pools, [scr_lgd_downturn()], [scr_lgd_floor()] and [scr_elbe()]
#' complete the estimate, [scr_lgd_pools()] and [scr_lgd_validate()] close
#' the pools and the validation. EAD: [scr_ead_data()] builds the realised
#' conversion factors from facility snapshots and [scr_ead()] the pools;
#' [scr_ead_downturn()] and [scr_ead_validate()] add the downturn and the
#' validation. [scr_el()], [scr_irb_rw()], [scr_sa_rw()], [scr_capital()],
#' [scr_pd_stress()] and [scr_ecl()] compute expected loss, risk weights,
#' capital and expected credit loss. Binning
#' against a continuous target goes through [scr_bin_continuous()], whose
#' result the engine reproduces in R and in SQL. The regulatory texts behind
#' the presets are listed in [scr_irb_params()]; users are responsible for
#' checking the tables against the texts in force before any regulatory use.
#'
#' @section Parallelism:
#'
#' Column-wise work (binning, hold-out revalidation, CSI) and the bootstrap
#' run on `config$nthread` workers. The backend follows
#' `getOption("scorecraft.parallel")`: `"fork"` on unix by default, `"psock"`
#' on Windows (and selectable anywhere, e.g. for tests), `"serial"` to switch
#' parallelism off. Results are identical across backends.
#'
#' Forked workers are clones of the parent and, because the garbage collector
#' writes to the objects it marks, each one ends up owning a copy of most of
#' the parent heap. On Linux the number of fork workers is therefore capped
#' at `getOption("scorecraft.fork_mem_fraction", 0.75)` of the memory
#' available divided by the resident size of the session, with a message
#' when the cap applies. Set the option to `Inf` to disable it.
#'
#' @section Reading conventions:
#'
#' `objective` declares the vocabulary and the direction of the scale
#' (`"risk"`: more points, safer; `"propensity"`: more points, more likely)
#' and does **not** change what is modelled. `event_level` changes what is
#' modelled. Both are documented in [scr_config()] and [scr_split()].
#'
#' @keywords internal
#' @importFrom data.table := .N .SD as.data.table copy data.table dcast fcase
#'   fifelse rbindlist set setDT setcolorder setnames setorder setorderv uniqueN
#' @importFrom stats coef glm binomial median quantile sd predict qchisq qnorm
#'   plogis qlogis pbinom setNames lm weighted.mean
#' @importFrom utils head modifyList packageVersion
"_PACKAGE"

## data.table: this package knows how to use `[.data.table`
.datatable.aware <- TRUE

## silence the "no visible binding" NOTE for columns used inside `[.data.table`
utils::globalVariables(c(
  ".", "feature", "error", "selected", "reason", "N", "triage_status",
  "triage_reason", "holdout_ok", "holdout_reason", "iv_holdout", "approved",
  "exit_stage", "consensus_rank", "total_iv", "iv_suspect", "screen_reason",
  "screen_selected", "bin_error", "prune_status", "prune_corr_with", "prune_corr",
  "consensus_selected", "consensus_reason", "votes", "consensus_score",
  "vote", "importance", "score_pct", "weight", "model", "gini",
  "kind", "source", "output", "impute_value", "derived_from", "corr_with", "corr",
  "i.corr_with", "i.corr", "i.of", "of", "band", "y", "s", "n1", "n0", "psi",
  "psi_flag", "pct_unbinned", "iv_train_bins", "iv_ratio", "variable", "bin",
  "points", "woe", "points_raw", "n", "events", "score", "cum_pct", "period",
  "sample", "target", "pos", "n_targets", "mean_rank", "bin_id", "type", "id",
  "n_bins", "ks", "auc", "count", "event_rate", "cut", "decision", "n_pos",
  "pct", "csi", "shift", "level", "pct_base", "pct_compare", "value",
  "iv_quick", "in_model", "n_events", "n_total", "keep", "action", "ln_odds",
  "ln_odds_fit", "raw_mean", "w", "non_events", "cum_event_pct",
  "cum_nonevent_pct", "lift", "cum_lift", "odds", "log_odds", "gap", "prev_rate",
  "monotone", "p_value", "break_flag", "date", "score_points", "min_score",
  "max_score", "mean_score", "ep_per_account", "band_profit", "cum_event_rate",
  "cum_profit", "n_dev", "events_dev", "rate_dev", "n_pop", "n_unknown",
  "coverage", "coverage_flag", "multiplier", "rate_unknown", "events_implied",
  "rate_implied", "points_shift", "Feature", "Gain", ".tmp", ".i",
  "n_vote", "count_train", "pct_train", "observed", "expected",
  "psi_critical", "psi_flag_adjusted", "provenance", "manual_reason", "status", "pct_shift",
  "pct_holdout", "current", "metric", "optimal", "manual", "delta", "blocking"
))
