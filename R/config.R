# ============================================================================ #
# config.R - the single pipeline configuration
# ============================================================================ #
# One config object crosses every stage. Presets touch only the tightness of
# the funnel (how many variables survive); everything else is identical.
# ============================================================================ #

#' Pipeline configuration
#'
#' Builds the configuration object that crosses every stage, from
#' [scr_split()] to [scr_export()]. A preset sets the tightness of the
#' selection funnel; any individual key can be overridden through `...`.
#'
#' @section Presets:
#'
#' A preset touches four keys and nothing else: `target_max`, `min_votes`,
#' `corr_cutoff` and `iv_min`.
#'
#' | preset | variables at the end | `min_votes` | `corr_cutoff` | `iv_min` |
#' |---|---:|---:|---:|---:|
#' | `"aggressive"` | 10 to 15 | 3 | 0.60 | 0.03 |
#' | `"moderate"` | 10 to 25 | 2 | 0.70 | 0.02 |
#' | `"lazy"` | 10 to 40 | 1 | 0.80 | 0.02 |
#'
#' Use [scr_presets()] to see the resolved table and [scr_config_keys()] for
#' the full key dictionary.
#'
#' @section Risk or propensity (`objective`):
#'
#' The two literatures use the same mathematics with opposite conventions.
#' In credit and fraud, `target = 1` is the **bad** case and the scorecard is
#' built so that more points mean less risk. In propensity, `target = 1` is
#' the **good** case and the campaign list needs more points to mean a higher
#' chance of engaging.
#'
#' | | `"risk"` (default) | `"propensity"` |
#' |---|---|---|
#' | `target = 1` means | undesirable event | desirable event |
#' | Points scale | more points = safer | more points = more likely |
#' | Derived `direction` | `"higher_is_safer"` | `"higher_is_riskier"` |
#' | `odds_orientation` | `safe:event` | `event:safe` |
#'
#' **`objective` does not touch the selection.** It does not change the
#' modelled target, the cut points, the IV or the shortlist; it acts on the
#' direction of the points scale, on the odds orientation of the alignment
#' and on the vocabulary of the reports. To model the other class as the
#' event, the argument is `event_level` in [scr_split()] and [scr_select()],
#' and that one, unlike this, rewrites everything.
#'
#' @section Binning algorithm (`algorithm`):
#'
#' `"jedi"` is the default and stays exposed side by side with the
#' alternatives, never hidden behind an `"auto"`. Choices with
#' distinct properties: `"ivb"`, `"dp"` and `"sblp"` are provably optimal for
#' categoricals; `"cm"` (ChiMerge) and `"fetb"` have a principled stopping
#' rule; `"ir"` (isotonic) guarantees monotone WOE; `"fast_mdlp"` is the
#' faithful Fayyad-Irani. The full list is in
#' [OptimalBinningWoE::obwoe_algorithms()]. A numeric-only or
#' categorical-only algorithm applies where it is valid and the other type
#' falls back to `"jedi"`.
#'
#' @section The Information Value gate:
#'
#' \describe{
#'   \item{`iv_min`}{Admission floor. Fails with `IV_BELOW_MIN`.}
#'   \item{`iv_max`}{Admission ceiling. Fails with `IV_SUSPICIOUS`. `1.00`
#'     (default) tolerates a legitimately strong predictor and cuts the
#'     absurd; `0.50` is the engine default, calibrated for credit default.}
#'   \item{`iv_suspect`}{Only the threshold of the report warning. Fails nothing.}
#' }
#'
#' The real leakage detector is `allow_degenerate = FALSE`: a bin with no
#' events or no non-events is the symptom that has no innocent explanation.
#'
#' @section Scorecard scale:
#'
#' `base_score`, `base_odds` and `pdo` are one statement: *at `base_score`
#' points the odds are `base_odds`, and every `pdo` points they double*.
#' `base_odds` is always expressed in the orientation `direction` implies
#' (non-event:event under `higher_is_safer`; event:non-event under
#' `higher_is_riskier`), and the alignment object records
#' `odds_orientation` so this is never implicit. The classic 600/50/20 is
#' Siddiqi's (2006) textbook example, not a parameter published by any
#' bureau. `align_method = "regression"` (default) takes the raw score to
#' that scale by regressing empirical log-odds on score bands, which absorbs
#' reweighting, miscalibration and prior shift; `"direct"` assumes the model
#' is calibrated and uses the logit as is.
#'
#' @param preset One of `"moderate"` (default), `"aggressive"` or `"lazy"`.
#' @param ... Overrides of any configuration key, by name. `NULL` means
#'   "keep the preset value", not "delete the key". An unknown name is an
#'   error, on purpose: a silent override leaves dead configuration in the
#'   file.
#'
#' @return An object of class `scr_config`: a named list with every key resolved.
#'
#' @references
#' Siddiqi, N. (2006). *Credit Risk Scorecards: Developing and Implementing
#' Intelligent Credit Scoring*. Wiley.
#'
#' @seealso [scr_select()] to use the configuration, [scr_presets()] to
#'   compare presets, [scr_config_keys()] for the key dictionary.
#' @family configuration
#' @examples
#' cfg <- scr_config()
#' cfg
#'
#' # propensity: more points = more likely to have the event
#' scr_config(objective = "propensity")$objective
#'
#' # NULL keeps the preset value (here, iv_min = 0.03 from aggressive)
#' scr_config("aggressive", iv_min = NULL)$iv_min
#'
#' # a wrong name fails loudly instead of becoming dead configuration
#' try(scr_config(iv_maximum = 1))
#' @export
scr_config <- function(preset = c("moderate", "aggressive", "lazy"), ...) {

  preset <- match.arg(preset)

  base <- list(

    # ---- general --------------------------------------------------------- #
    preset       = preset,
    seed         = 2203L,
    nthread      = 2L,
    verbose      = TRUE,
    objective    = "risk",

    # ---- Stage 0: split -------------------------------------------------- #
    oot_date_col  = NULL,
    holdout_ratio = 0.30,

    # ---- Stage 1: triage ------------------------------------------------- #
    special_values      = -999,
    max_missing         = 0.70,
    near_constant       = 0.995,
    max_cat_levels      = 300L,
    min_iv_quick        = 0.005,
    quick_iv_groups     = 20L,
    check_duplicates    = TRUE,
    special_min_share   = 0.01,
    special_min_woe     = 0.05,
    flag_suffix         = "__sp",
    allow_derived_final = FALSE,

    # ---- Stage 2: binning ------------------------------------------------ #
    min_bins       = 3L,
    max_bins       = 7L,
    algorithm      = "jedi",
    bin_cutoff     = 0.03,
    max_n_prebins  = 30L,
    max_iterations = 1000L,
    bin_separator  = "%;%",

    # ---- Stage 2: screening ---------------------------------------------- #
    iv_min            = 0.02,
    iv_max            = 1.00,
    iv_suspect        = 0.50,
    require_monotonic = "numeric",
    monotonicity      = "weak",
    min_bin_pct       = 0.02,
    screen_min_bins   = 2L,
    allow_degenerate  = FALSE,
    screen_top_n      = NULL,

    # ---- Stage 2: hold-out revalidation ---------------------------------- #
    iv_ratio_min = 0.50,
    psi_max      = 0.25,
    psi_alpha    = 0.05,
    max_unbinned = 0.02,

    # ---- Stage 2: redundancy --------------------------------------------- #
    corr_cutoff = 0.70,
    corr_method = "spearman",

    # ---- Stage 3: classifiers -------------------------------------------- #
    use_glmnet   = TRUE,
    use_xgboost  = TRUE,
    use_ranger   = TRUE,
    use_lightgbm = TRUE,
    cv_folds     = 5L,
    en_alpha     = 0.50,
    xgb_rounds   = 400L,
    xgb_eta      = 0.05,
    xgb_max_depth = 4L,
    xgb_subsample = 0.80,
    xgb_colsample = 0.80,
    xgb_min_child_weight = 20,
    xgb_early_stopping   = 40L,
    rf_trees      = 300L,
    rf_importance = "permutation",
    model_top_k   = 40L,
    model_max_rows = 200000L,

    # ---- Stage 4: consensus ---------------------------------------------- #
    min_votes  = 2L,
    target_min = 10L,
    target_max = 25L,
    weight_by_gini = TRUE,

    # ---- Stages 4-5: scorecard and alignment ----------------------------- #
    base_score    = 600,
    base_odds     = 50,
    pdo           = 20,
    direction     = NULL,
    align_method  = "regression",
    align_bands   = 10L,
    points_style  = "base_plus_deviation",
    points_round  = TRUE,
    challenger    = NULL,
    max_abs_coef  = 15,
    n_boot        = 200L,
    ci_level      = 0.95,
    score_groups  = 10L,

    # ---- Stage 6: cut-off and strategy ----------------------------------- #
    cutoff_n           = 20L,
    reject_multipliers = c(2, 4, 8),
    lab_max_iv_loss    = 0.10,
    lab_min_bin_pct_hard = 0.005,

    # ---- Stage 7: outputs ------------------------------------------------ #
    sql_table        = "your_table",
    sql_dialect      = "ansi",
    sql_output       = "both",
    sql_keep_columns = character(),

    # ---- Stage 8: default definition (IRB) ------------------------------- #
    default_days       = 90L,
    default_abs        = 100,
    default_rel        = 0.01,
    default_level      = "obligor",
    default_probation  = 3L,
    default_probation_restructured = 12L,
    default_pulling    = 0.20,
    asset_class        = "retail_other",

    # ---- Stage 9: PD calibration and grades (IRB) ------------------------ #
    pd_calibration  = "intercept",
    pd_n_grades     = 10L,
    pd_grade_method = "geometric",
    pd_min_obligors = 30L,
    pd_min_defaults = 20L,
    pd_source       = "lra",
    pd_moc_method   = "ci_timeseries",
    pd_moc_level    = 0.95,
    pd_lights       = c(0.01, 0.05),
    pd_dr_by        = "quarter",

    # ---- Stage 10: LGD (IRB) --------------------------------------------- #
    lgd_discount_add_on  = 0.05,
    lgd_discount_rate    = NA_real_,
    lgd_cap_at_one       = FALSE,
    lgd_floor_at_zero    = TRUE,
    lgd_t_max            = 60L,
    lgd_cure_window      = 9L,
    lgd_min_defaults_bin = 100L,
    lgd_cost_allocation  = "ead",
    lgd_downturn         = "type1",
    lgd_downturn_add_on  = 0.15,
    lgd_n_pools          = 7L,
    lgd_severity         = "fractional_logit",
    lgd_elbe_grid        = c(0L, 6L, 12L, 24L, 36L),

    # ---- Stage 11: EAD / CCF (IRB) --------------------------------------- #
    ccf_horizon        = "fixed",
    ccf_horizon_months = 12L,
    ccf_measure        = "auto",
    ccf_u_star         = 0.95,
    ccf_floor_realised = 0,
    ccf_cap_realised   = NA_real_,
    ccf_min_defaults   = 30L,
    ccf_n_pools        = 5L,
    ccf_moc_alpha      = 0.05,
    ccf_downturn       = "type1",
    ccf_sa_ccf         = 0.40,
    post_default_drawings_in = "lgd",

    # ---- Stage 12: expected loss, capital, ECL (IRB) --------------------- #
    framework            = "bcb",
    capital_approach     = "airb",
    capital_target_ratio = 0.08,
    capital_output_floor = TRUE,
    capital_sensitivity  = TRUE,
    ecl_discount         = "eir",
    ecl_stage_dpd        = c(30L, 90L),
    ecl_sicr_ratio       = 2.0,
    ecl_horizon_months   = 12L
  )

  tight <- switch(preset,
    aggressive = list(target_max = 15L, min_votes = 3L, corr_cutoff = 0.60, iv_min = 0.03),
    moderate   = list(target_max = 25L, min_votes = 2L, corr_cutoff = 0.70, iv_min = 0.02),
    lazy       = list(target_max = 40L, min_votes = 1L, corr_cutoff = 0.80, iv_min = 0.02)
  )
  cfg <- utils::modifyList(base, tight)

  over <- list(...)
  over <- over[!vapply(over, is.null, logical(1))]
  dup <- unique(names(over)[duplicated(names(over))])
  if (length(dup)) {
    stop("scr_config(): repeated key(s): ", paste(dup, collapse = ", "),
         ". Pass each key once.", call. = FALSE)
  }
  unknown <- setdiff(names(over), names(cfg))
  if (length(unknown)) {
    stop("scr_config(): unknown key(s): ", paste(unknown, collapse = ", "),
         ". Fix the name - a silent override hides dead configuration.", call. = FALSE)
  }
  cfg <- utils::modifyList(cfg, over)
  .scr_validate_config(cfg)
}

#' Validation of the keys that, when wrong, would only blow up stages later
#' @keywords internal
#' @noRd
.scr_validate_config <- function(cfg) {
  cfg$target_min <- max(1L, as.integer(cfg$target_min))
  cfg$target_max <- min(40L, max(cfg$target_min, as.integer(cfg$target_max)))

  if (!cfg$objective %in% c("risk", "propensity")) {
    stop("`objective` must be \"risk\" or \"propensity\" (got \"", cfg$objective, "\").", call. = FALSE)
  }
  if (!is.null(cfg$direction) && !cfg$direction %in% c("higher_is_safer", "higher_is_riskier")) {
    stop("`direction` must be NULL, \"higher_is_safer\" or \"higher_is_riskier\".", call. = FALSE)
  }
  algs <- OptimalBinningWoE::obwoe_algorithms()$algorithm
  if (!identical(cfg$algorithm, "auto") && !cfg$algorithm %in% algs) {
    stop("unknown `algorithm`: \"", cfg$algorithm, "\". See OptimalBinningWoE::obwoe_algorithms().", call. = FALSE)
  }
  if (!cfg$align_method %in% c("regression", "direct")) {
    stop("`align_method` must be \"regression\" or \"direct\".", call. = FALSE)
  }
  if (!cfg$points_style %in% c("base_plus_deviation", "distributed")) {
    stop("`points_style` must be \"base_plus_deviation\" or \"distributed\".", call. = FALSE)
  }
  if (!is.null(cfg$challenger) && !cfg$challenger %in% c("xgboost", "lightgbm")) {
    stop("`challenger` must be NULL, \"xgboost\" or \"lightgbm\".", call. = FALSE)
  }
  .scr_num1(cfg$pdo, "pdo", lower = 0, open_lower = TRUE)
  .scr_num1(cfg$base_odds, "base_odds", lower = 0, open_lower = TRUE)
  .scr_num1(cfg$base_score, "base_score")
  .scr_num1(cfg$holdout_ratio, "holdout_ratio", lower = 0, upper = 1, open_lower = TRUE)
  if (!cfg$sql_output %in% c("woe", "bin", "both")) {
    stop("`sql_output` must be \"woe\", \"bin\" or \"both\".", call. = FALSE)
  }
  if (!is.numeric(cfg$nthread) || length(cfg$nthread) != 1L || is.na(cfg$nthread)) {
    stop("`nthread` must be a single number.", call. = FALSE)
  }
  cfg$nthread <- max(1L, as.integer(cfg$nthread))

  # ---- scorecard pipeline keys (stages 0-7) --------------------------------- #
  # Checked here because, wrong, they fail minutes later inside a stage with
  # an engine's message - or, worse, do not fail at all.
  .flag <- function(key) {
    if (!is.logical(cfg[[key]]) || length(cfg[[key]]) != 1L || is.na(cfg[[key]])) {
      stop("`", key, "` must be TRUE or FALSE.", call. = FALSE)
    }
  }
  .int <- function(key, lower, upper = Inf) {
    v <- cfg[[key]]
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || (is.finite(v) && v != round(v)) ||
        v < lower || v > upper) {
      stop("`", key, "` must be a single whole number in [", lower, ", ", upper, "].", call. = FALSE)
    }
  }
  .str1 <- function(key) {
    if (!is.character(cfg[[key]]) || length(cfg[[key]]) != 1L || is.na(cfg[[key]]) || !nzchar(cfg[[key]])) {
      stop("`", key, "` must be a single non-empty string.", call. = FALSE)
    }
  }
  .int("seed", -.Machine$integer.max, .Machine$integer.max)
  for (k in c("verbose", "check_duplicates", "allow_derived_final", "allow_degenerate",
              "use_glmnet", "use_xgboost", "use_ranger", "use_lightgbm", "weight_by_gini",
              "points_round")) .flag(k)
  if (!is.null(cfg$oot_date_col)) .str1("oot_date_col")
  if (anyNA(cfg$special_values) || !(is.numeric(cfg$special_values) || !length(cfg$special_values))) {
    stop("`special_values` must be a numeric vector without NA (numeric() for none).", call. = FALSE)
  }
  .str1("flag_suffix")
  .str1("bin_separator")
  .scr_num1(cfg$max_missing, "max_missing", lower = 0, upper = 1)
  .scr_num1(cfg$near_constant, "near_constant", lower = 0, upper = 1, open_lower = TRUE)
  .int("max_cat_levels", 2)
  .scr_num1(cfg$min_iv_quick, "min_iv_quick", lower = 0)
  .int("quick_iv_groups", 2)
  .scr_num1(cfg$special_min_share, "special_min_share", lower = 0, upper = 0.5)
  .scr_num1(cfg$special_min_woe, "special_min_woe", lower = 0)
  .int("min_bins", 2)
  .int("max_bins", cfg$min_bins)
  .int("screen_min_bins", 1)
  .int("max_n_prebins", 2)
  .int("max_iterations", 1)
  .scr_num1(cfg$bin_cutoff, "bin_cutoff", lower = 0, upper = 1, open_lower = TRUE)
  .scr_num1(cfg$iv_min, "iv_min", lower = 0)
  if (!is.numeric(cfg$iv_max) || length(cfg$iv_max) != 1L || is.na(cfg$iv_max) || cfg$iv_max <= cfg$iv_min) {
    stop("`iv_max` must be a single number (Inf allowed) above `iv_min`.", call. = FALSE)
  }
  .scr_num1(cfg$iv_suspect, "iv_suspect", lower = 0)
  if (!is.character(cfg$require_monotonic) || !cfg$require_monotonic %in% c("numeric", "all", "none")) {
    stop("`require_monotonic` must be \"numeric\", \"all\" or \"none\".", call. = FALSE)
  }
  if (!is.character(cfg$monotonicity) || !cfg$monotonicity %in% c("weak", "strict")) {
    stop("`monotonicity` must be \"weak\" or \"strict\".", call. = FALSE)
  }
  .scr_num1(cfg$min_bin_pct, "min_bin_pct", lower = 0, upper = 1)
  if (cfg$min_bin_pct >= 1) stop("`min_bin_pct` must be below 1.", call. = FALSE)
  if (!is.null(cfg$screen_top_n)) .int("screen_top_n", 1)
  .scr_num1(cfg$iv_ratio_min, "iv_ratio_min", lower = 0)
  .scr_num1(cfg$psi_max, "psi_max", lower = 0, open_lower = TRUE)
  .scr_num1(cfg$psi_alpha, "psi_alpha", lower = 0, upper = 1, open_lower = TRUE)
  .scr_num1(cfg$max_unbinned, "max_unbinned", lower = 0, upper = 1)
  .scr_num1(cfg$corr_cutoff, "corr_cutoff", lower = 0, upper = 1, open_lower = TRUE)
  if (!is.character(cfg$corr_method) ||
      !cfg$corr_method %in% c("pearson", "spearman", "kendall", "hoeffding", "distance", "biweight", "pbend")) {
    stop("`corr_method` must be one of \"pearson\", \"spearman\", \"kendall\", \"hoeffding\", ",
         "\"distance\", \"biweight\", \"pbend\".", call. = FALSE)
  }
  .int("cv_folds", 3)   # cv.glmnet refuses fewer than three folds
  .scr_num1(cfg$en_alpha, "en_alpha", lower = 0, upper = 1)
  .int("xgb_rounds", 1)
  .scr_num1(cfg$xgb_eta, "xgb_eta", lower = 0, upper = 1, open_lower = TRUE)
  .int("xgb_max_depth", 1)
  .scr_num1(cfg$xgb_subsample, "xgb_subsample", lower = 0, upper = 1, open_lower = TRUE)
  .scr_num1(cfg$xgb_colsample, "xgb_colsample", lower = 0, upper = 1, open_lower = TRUE)
  .scr_num1(cfg$xgb_min_child_weight, "xgb_min_child_weight", lower = 0)
  .int("xgb_early_stopping", 1)
  .int("rf_trees", 1)
  if (!is.character(cfg$rf_importance) ||
      !cfg$rf_importance %in% c("permutation", "impurity", "impurity_corrected")) {
    stop("`rf_importance` must be \"permutation\", \"impurity\" or \"impurity_corrected\".", call. = FALSE)
  }
  .int("model_top_k", 1)
  if (!is.numeric(cfg$model_max_rows) || length(cfg$model_max_rows) != 1L || is.na(cfg$model_max_rows) ||
      cfg$model_max_rows < 10) {
    stop("`model_max_rows` must be a single number >= 10 (Inf for no cap).", call. = FALSE)
  }
  .int("min_votes", 1)
  .scr_num1(cfg$max_abs_coef, "max_abs_coef", lower = 0, open_lower = TRUE)
  .int("align_bands", 2)
  .int("n_boot", 1)
  .scr_num1(cfg$ci_level, "ci_level", lower = 0, upper = 1, open_lower = TRUE)
  if (cfg$ci_level >= 1) stop("`ci_level` must be below 1.", call. = FALSE)
  .int("score_groups", 2)
  .int("cutoff_n", 2)
  if (!is.numeric(cfg$reject_multipliers) || !length(cfg$reject_multipliers) ||
      any(!is.finite(cfg$reject_multipliers)) || any(cfg$reject_multipliers <= 0)) {
    stop("`reject_multipliers` must be a vector of positive numbers.", call. = FALSE)
  }
  .scr_num1(cfg$lab_max_iv_loss, "lab_max_iv_loss", lower = 0, upper = 1)
  .scr_num1(cfg$lab_min_bin_pct_hard, "lab_min_bin_pct_hard", lower = 0, upper = 1)
  .str1("sql_table")
  .str1("sql_dialect")
  if (!is.character(cfg$sql_keep_columns) || anyNA(cfg$sql_keep_columns)) {
    stop("`sql_keep_columns` must be a character vector.", call. = FALSE)
  }

  # ---- IRB keys ------------------------------------------------------------ #
  .enum <- function(key, allowed) {
    if (!cfg[[key]] %in% allowed) {
      stop("`", key, "` must be one of ", paste0("\"", allowed, "\"", collapse = ", "),
           " (got \"", cfg[[key]], "\").", call. = FALSE)
    }
  }
  .enum("default_level", c("obligor", "facility"))
  .enum("asset_class", .irb_classes)
  .enum("pd_calibration", c("intercept", "logodds_ab", "scaling", "qmm"))
  .enum("pd_grade_method", c("geometric", "quantile", "supplied"))
  .enum("pd_source", c("lra", "mean_pd"))
  .enum("pd_moc_method", c("ci_timeseries", "ci_binomial", "bootstrap"))
  .enum("pd_dr_by", c("month", "quarter", "year"))
  .enum("lgd_cost_allocation", c("ead", "count", "duration"))
  .enum("lgd_downturn", c("type1", "type3", "none"))
  .enum("lgd_severity", c("fractional_logit", "beta"))
  .enum("ccf_horizon", c("fixed", "cohort", "variable"))
  .enum("ccf_measure", c("auto", "ulf", "lf", "eadf"))
  .enum("ccf_downturn", c("type1", "type3", "none"))
  .enum("post_default_drawings_in", c("lgd", "ccf"))
  .enum("framework", c("bcb", "basel3_final", "crr3"))
  .enum("capital_approach", c("airb", "firb"))
  .enum("ecl_discount", c("eir", "none"))
  .scr_num1(cfg$default_rel, "default_rel", lower = 0, upper = 1)
  .scr_num1(cfg$default_pulling, "default_pulling", lower = 0, upper = 1)
  .scr_num1(cfg$pd_moc_level, "pd_moc_level", lower = 0, upper = 1, open_lower = TRUE)
  .scr_num1(cfg$lgd_discount_add_on, "lgd_discount_add_on", lower = 0)
  .scr_num1(cfg$lgd_downturn_add_on, "lgd_downturn_add_on", lower = 0, upper = 1)
  .scr_num1(cfg$ccf_u_star, "ccf_u_star", lower = 0, upper = 1)
  .scr_num1(cfg$ccf_sa_ccf, "ccf_sa_ccf", lower = 0, upper = 1)
  .scr_num1(cfg$capital_target_ratio, "capital_target_ratio", lower = 0, upper = 1)
  if (length(cfg$pd_lights) != 2L || any(!is.finite(cfg$pd_lights)) || cfg$pd_lights[1] >= cfg$pd_lights[2]) {
    stop("`pd_lights` must be two increasing p-value thresholds (red below the first, amber below the second).", call. = FALSE)
  }
  structure(cfg, class = c("scr_config", "list"))
}

#' Direction of the points scale, resolved
#'
#' An explicit `direction` wins; `NULL` derives it from `objective`.
#' @keywords internal
#' @noRd
resolve_direction <- function(cfg) {
  cfg$direction %||%
    if (identical(cfg$objective, "propensity")) "higher_is_riskier" else "higher_is_safer"
}

#' Report vocabulary under the project convention
#' @keywords internal
#' @noRd
vocab <- function(cfg) {
  if (identical(cfg$objective, "risk")) {
    list(event = "event (undesirable)", target1 = "target = 1 is the BAD case",
         points = "more points = lower probability of the event (safer)", score = "risk")
  } else {
    list(event = "event (desirable)", target1 = "target = 1 is the GOOD case",
         points = "more points = higher probability of the event (more likely)", score = "propensity")
  }
}

#' @export
print.scr_config <- function(x, ...) {
  cat(sprintf("<scr_config> preset \"%s\" | objective \"%s\" | seed %s | threads %d\n\n",
              x$preset, x$objective, x$seed, x$nthread))
  line <- function(lab, val) cat(sprintf("  %-22s %s\n", lab, val))
  v <- vocab(x)
  cat("Convention\n")
  line("target = 1", v$target1)
  line("points scale", sprintf("%s [%s]", v$points, resolve_direction(x)))
  cat("\nFunnel\n")
  line("variables at the end", sprintf("%d to %d", x$target_min, x$target_max))
  line("minimum votes", x$min_votes)
  line("admissible IV", sprintf("[%s, %s)  warning at %s", x$iv_min, x$iv_max, x$iv_suspect))
  line("correlation", sprintf("%s (%s)", x$corr_cutoff, x$corr_method))
  cat("\nBinning\n")
  line("bins", sprintf("%d to %d, algorithm \"%s\"", x$min_bins, x$max_bins, x$algorithm))
  line("monotonicity", sprintf("%s (%s)", x$require_monotonic, x$monotonicity))
  line("smallest bin", fmt_pct(x$min_bin_pct))
  cat("\nScorecard\n")
  line("scale", sprintf("%s points at odds %s:1, PDO %s", x$base_score, x$base_odds, x$pdo))
  line("alignment", sprintf("%s (%d bands)", x$align_method, x$align_bands))
  line("challenger", x$challenger %||% "none")
  line("bootstrap CI", sprintf("%d resamples, %.0f%%", x$n_boot, 100 * x$ci_level))
  cat("\nData\n")
  line("sentinels", if (length(x$special_values)) paste(x$special_values, collapse = ", ") else "(none)")
  line("derived at the end", if (isTRUE(x$allow_derived_final)) "yes" else "no (diagnostic only)")
  line("hold-out", sprintf("%s%s", fmt_pct(x$holdout_ratio),
                           if (!is.null(x$oot_date_col)) sprintf(", OOT by \"%s\"", x$oot_date_col) else ""))
  cat("\nModels\n")
  mods <- c("glmnet", "xgboost", "ranger", "lightgbm")[c(x$use_glmnet, x$use_xgboost, x$use_ranger, x$use_lightgbm)]
  line("enabled", paste(mods, collapse = ", "))
  line("row cap", n_fmt(x$model_max_rows))
  invisible(x)
}

#' Selection presets, side by side
#'
#' Returns the funnel keys resolved per preset, to compare before choosing.
#' `target_min` and `iv_max` are shown for context; the presets leave them
#' unchanged, as they do every other configuration key.
#'
#' @return A `data.frame` with one row per preset.
#'
#' @family configuration
#' @examples
#' scr_presets()
#' @export
scr_presets <- function() {
  p <- c("aggressive", "moderate", "lazy")
  do.call(rbind, lapply(p, function(nm) {
    c0 <- scr_config(nm)
    data.frame(preset = nm, target_min = c0$target_min, target_max = c0$target_max,
               min_votes = c0$min_votes, corr_cutoff = c0$corr_cutoff,
               iv_min = c0$iv_min, iv_max = c0$iv_max, stringsAsFactors = FALSE)
  }))
}

#' Dictionary of configuration keys
#'
#' One row per [scr_config()] key, with the stage it acts on, the default
#' value and what it controls.
#'
#' @param stage Optional filter by stage: `0` to `7` for the scorecard
#'   pipeline, `8` to `12` for the IRB models. `NULL` returns everything.
#'
#' @return A `data.frame` with `key`, `stage`, `default` and `description`.
#'
#' @family configuration
#' @examples
#' head(scr_config_keys(), 8)
#' scr_config_keys(stage = 5)
#' @export
scr_config_keys <- function(stage = NULL) {
  d <- rbind(
    .ck("preset", 0, "moderate", "Tightness of the selection funnel"),
    .ck("objective", 0, "risk", "Convention: target=1 is an undesirable (risk) or desirable (propensity) event"),
    .ck("seed", 0, "2203", "Seed of everything random"),
    .ck("nthread", 0, "2", "Parallel workers (binning by column, bootstrap, CSI)"),
    .ck("verbose", 0, "TRUE", "Progress messages"),
    .ck("oot_date_col", 0, "NULL", "Date column of the out-of-time cut"),
    .ck("holdout_ratio", 0, "0.30", "Target hold-out fraction"),
    .ck("special_values", 1, "-999", "Business sentinels"),
    .ck("max_missing", 1, "0.70", "Ceiling of missing+sentinel share for a numeric"),
    .ck("near_constant", 1, "0.995", "Mode dominance that fails"),
    .ck("max_cat_levels", 1, "300", "Maximum categorical cardinality"),
    .ck("min_iv_quick", 1, "0.005", "Noise floor of the coarse IV"),
    .ck("quick_iv_groups", 1, "20", "Quantiles of the coarse IV"),
    .ck("check_duplicates", 1, "TRUE", "Detect identical columns"),
    .ck("special_min_share", 1, "0.01", "Minimum mass to create the flag"),
    .ck("special_min_woe", 1, "0.05", "Minimum |WOE| to create the flag"),
    .ck("flag_suffix", 1, "__sp", "Suffix of the special-population flag"),
    .ck("allow_derived_final", 1, "FALSE", "A derived flag may enter the deliverable"),
    .ck("min_bins", 2, "3", "Minimum bins requested from the algorithm"),
    .ck("max_bins", 2, "7", "Maximum bins requested from the algorithm"),
    .ck("algorithm", 2, "jedi", "Binning algorithm (see obwoe_algorithms())"),
    .ck("bin_cutoff", 2, "0.03", "Minimum fraction per bin inside the algorithm"),
    .ck("max_n_prebins", 2, "30", "Pre-bins; strongly affects numerics"),
    .ck("max_iterations", 2, "1000", "Maximum optimiser iterations"),
    .ck("bin_separator", 2, "%;%", "Separator of merged categories"),
    .ck("iv_min", 2, "0.02", "IV admission floor"),
    .ck("iv_max", 2, "1.00", "IV admission ceiling (leakage)"),
    .ck("iv_suspect", 2, "0.50", "Threshold of the suspicious-IV warning"),
    .ck("require_monotonic", 2, "numeric", "Where monotonicity is required"),
    .ck("monotonicity", 2, "weak", "Monotonicity strictness"),
    .ck("min_bin_pct", 2, "0.02", "Minimum population of the smallest bin"),
    .ck("screen_min_bins", 2, "2", "Structural floor of bins at admission"),
    .ck("allow_degenerate", 2, "FALSE", "Accept a single-class bin"),
    .ck("screen_top_n", 2, "NULL", "Ranking cut at screening"),
    .ck("iv_ratio_min", 2, "0.50", "Minimum hold-out IV / train IV"),
    .ck("psi_max", 2, "0.25", "Maximum fixed PSI between train and hold-out"),
    .ck("psi_alpha", 2, "0.05", "Alpha of the n-adjusted PSI threshold"),
    .ck("max_unbinned", 2, "0.02", "Tolerated fraction of hold-out without a bin"),
    .ck("corr_cutoff", 2, "0.70", "Correlation that defines redundancy"),
    .ck("corr_method", 2, "spearman", "Correlation method (rank)"),
    .ck("use_glmnet", 3, "TRUE", "Elastic net as a voter"),
    .ck("use_xgboost", 3, "TRUE", "XGBoost as a voter (required by the package)"),
    .ck("use_ranger", 3, "TRUE", "Random forest as a voter"),
    .ck("use_lightgbm", 3, "TRUE", "LightGBM if installed"),
    .ck("cv_folds", 3, "5", "Cross-validation folds of glmnet"),
    .ck("en_alpha", 3, "0.50", "Elastic net alpha"),
    .ck("xgb_rounds", 3, "400", "Maximum boosting trees"),
    .ck("xgb_eta", 3, "0.05", "Boosting learning rate"),
    .ck("xgb_max_depth", 3, "4", "Maximum tree depth"),
    .ck("xgb_subsample", 3, "0.80", "Row fraction per tree"),
    .ck("xgb_colsample", 3, "0.80", "Column fraction per tree"),
    .ck("xgb_min_child_weight", 3, "20", "Minimum weight per leaf"),
    .ck("xgb_early_stopping", 3, "40", "Rounds without hold-out improvement"),
    .ck("rf_trees", 3, "300", "Random forest trees"),
    .ck("rf_importance", 3, "permutation", "Random forest importance type"),
    .ck("model_top_k", 3, "40", "Top-K of each model that counts as a vote"),
    .ck("model_max_rows", 3, "200000", "Row cap for the classifiers"),
    .ck("min_votes", 4, "2", "Minimum votes to enter the shortlist"),
    .ck("target_min", 4, "10", "Floor of variables at the end"),
    .ck("target_max", 4, "25", "Ceiling of variables at the end"),
    .ck("weight_by_gini", 4, "TRUE", "Weight votes by hold-out Gini"),
    .ck("base_score", 5, "600", "Reference score of the scale"),
    .ck("base_odds", 5, "50", "Odds at the reference score (in the orientation of direction)"),
    .ck("pdo", 5, "20", "Points to double the odds"),
    .ck("direction", 5, "NULL", "Scale direction; NULL derives it from objective"),
    .ck("align_method", 5, "regression", "regression (bands + ln(odds) ~ raw) or direct"),
    .ck("align_bands", 5, "10", "Bands of the alignment regression"),
    .ck("points_style", 5, "base_plus_deviation", "Points: base + deviation, or distributed intercept"),
    .ck("points_round", 5, "TRUE", "Round points per bin"),
    .ck("challenger", 5, "NULL", "Tree challenger (xgboost/lightgbm), without points"),
    .ck("max_abs_coef", 5, "15", "Maximum absolute glm coefficient tolerated"),
    .ck("n_boot", 5, "200", "Bootstrap CI resamples (always computed)"),
    .ck("ci_level", 5, "0.95", "CI level"),
    .ck("score_groups", 5, "10", "Bands of the score gains"),
    .ck("cutoff_n", 6, "20", "Candidate cuts of the sweep"),
    .ck("reject_multipliers", 6, "2, 4, 8", "Sensitivity band of reject inference"),
    .ck("lab_max_iv_loss", 6, "0.10", "Coarse classing: advisory hold-out IV loss vs optimal"),
    .ck("lab_min_bin_pct_hard", 6, "0.005", "Coarse classing: blocking floor of a manual bin share"),
    .ck("sql_table", 7, "your_table", "Source table in the generated SQL"),
    .ck("sql_dialect", 7, "ansi", "Dialect of the generated SQL"),
    .ck("sql_output", 7, "both", "Emit WOE, BIN or both"),
    .ck("sql_keep_columns", 7, "character(0)", "Key columns carried untransformed"),
    .ck("default_days", 8, "90", "Days past due that trigger a default"),
    .ck("default_abs", 8, "100", "Absolute materiality threshold on arrears (currency units)"),
    .ck("default_rel", 8, "0.01", "Relative materiality threshold: arrears over exposure"),
    .ck("default_level", 8, "obligor", "Level of the default flag: obligor or facility"),
    .ck("default_probation", 8, "3", "Months without a trigger before leaving default"),
    .ck("default_probation_restructured", 8, "12", "Probation after a distressed restructuring"),
    .ck("default_pulling", 8, "0.20", "Share of defaulted exposure that pulls the whole obligor into default"),
    .ck("asset_class", 8, "retail_other", "Asset class of the floors and the risk-weight function (one key for PD, LGD and capital)"),
    .ck("pd_calibration", 9, "intercept", "Calibration to the central tendency: intercept, logodds_ab, scaling, qmm"),
    .ck("pd_n_grades", 9, "10", "Number of rating grades when derived"),
    .ck("pd_grade_method", 9, "geometric", "Grade construction: geometric, quantile, supplied"),
    .ck("pd_min_obligors", 9, "30", "Minimum obligors per grade before merging"),
    .ck("pd_min_defaults", 9, "20", "Minimum defaults per grade before merging"),
    .ck("pd_source", 9, "lra", "Grade PD: long-run default rate or mean of individual PDs"),
    .ck("pd_moc_method", 9, "ci_timeseries", "Estimation-error margin of conservatism: ci_timeseries, ci_binomial, bootstrap"),
    .ck("pd_moc_level", 9, "0.95", "One-sided confidence level of the estimation-error margin"),
    .ck("pd_lights", 9, "0.01, 0.05", "Traffic-light p-value thresholds: red at or below the first, amber at or below the second (PD and LGD validation)"),
    .ck("pd_dr_by", 9, "quarter", "Cohort frequency of the default-rate series"),
    .ck("lgd_discount_add_on", 10, "0.05", "Add-on over the reference rate in the workout discount rate"),
    .ck("lgd_discount_rate", 10, "NA", "Flat annual discount rate when no rates table is given"),
    .ck("lgd_cap_at_one", 10, "FALSE", "Cap realised LGD at 1"),
    .ck("lgd_floor_at_zero", 10, "TRUE", "Floor realised LGD at 0 for the averages (raw value kept)"),
    .ck("lgd_t_max", 10, "60", "Maximum workout period in months; older open defaults close with no further recovery"),
    .ck("lgd_cure_window", 10, "9", "Months within which two defaults of one facility are one event"),
    .ck("lgd_min_defaults_bin", 10, "100", "Minimum defaults per LGD bin or pool"),
    .ck("lgd_cost_allocation", 10, "ead", "Allocation key of indirect costs: ead, count, duration"),
    .ck("lgd_downturn", 10, "type1", "Downturn quantification: type1 (observed), type3 (add-on), none"),
    .ck("lgd_downturn_add_on", 10, "0.15", "Type-3 downturn add-on"),
    .ck("lgd_n_pools", 10, "7", "Target number of LGD pools"),
    .ck("lgd_severity", 10, "fractional_logit", "Severity engine: fractional_logit or beta"),
    .ck("lgd_elbe_grid", 10, "0, 6, 12, 24, 36", "Months since default of the in-default grid"),
    .ck("ccf_horizon", 11, "fixed", "Reference-date approach: fixed, cohort, variable"),
    .ck("ccf_horizon_months", 11, "12", "Months between reference date and default under the fixed horizon"),
    .ck("ccf_measure", 11, "auto", "Realised measure: auto (ULF below u_star, LF above), ulf, lf, eadf"),
    .ck("ccf_u_star", 11, "0.95", "Utilisation above which the limit factor replaces the CCF"),
    .ck("ccf_floor_realised", 11, "0", "Floor on the realised CCF (NA keeps the raw value)"),
    .ck("ccf_cap_realised", 11, "NA", "Cap on the realised CCF (NA: none)"),
    .ck("ccf_min_defaults", 11, "30", "Minimum defaults per CCF bin or pool"),
    .ck("ccf_n_pools", 11, "5", "Target number of CCF pools"),
    .ck("ccf_moc_alpha", 11, "0.05", "One-sided alpha of the estimation-error margin on the CCF"),
    .ck("ccf_downturn", 11, "type1", "Downturn quantification of the CCF: type1, type3, none"),
    .ck("ccf_sa_ccf", 11, "0.40", "Standardised CCF used for the own-estimate floor"),
    .ck("post_default_drawings_in", 11, "lgd", "Where drawings after default are booked: lgd or ccf"),
    .ck("framework", 12, "bcb", "Parameter preset: bcb, basel3_final, crr3"),
    .ck("capital_approach", 12, "airb", "Capital approach: airb (own LGD, floored) or firb (supervisory LGD, fixed maturity)"),
    .ck("capital_target_ratio", 12, "0.08", "Capital ratio applied to the reported RWA"),
    .ck("capital_output_floor", 12, "TRUE", "Compute the standardised comparison and the output floor"),
    .ck("capital_sensitivity", 12, "TRUE", "Run the fixed sensitivity grid"),
    .ck("ecl_discount", 12, "eir", "Discount the expected credit loss at the effective interest rate, or not"),
    .ck("ecl_stage_dpd", 12, "30, 90", "Days past due that move an exposure to stage 2 and stage 3"),
    .ck("ecl_sicr_ratio", 12, "2", "PD deterioration ratio (now over origination) that signals stage 2"),
    .ck("ecl_horizon_months", 12, "12", "Horizon of the 12-month expected credit loss")
  )
  if (!is.null(stage)) d <- d[d$stage %in% stage, , drop = FALSE]
  rownames(d) <- NULL
  d
}

#' @keywords internal
#' @noRd
.ck <- function(key, stage, default, description) {
  data.frame(key = key, stage = stage, default = default, description = description,
             stringsAsFactors = FALSE)
}

#' Binning control derived from the configuration
#'
#' A single point of truth: the Stage 2 fit and everything that reapplies
#' the bins read from here, which guarantees the cut points are the same.
#' @keywords internal
#' @noRd
obwoe_control <- function(cfg) {
  OptimalBinningWoE::control.obwoe(
    bin_cutoff     = cfg$bin_cutoff,
    max_n_prebins  = cfg$max_n_prebins,
    max_iterations = cfg$max_iterations,
    bin_separator  = cfg$bin_separator,
    verbose        = FALSE
  )
}

#' @keywords internal
#' @noRd
check_config <- function(config, fn) {
  if (!inherits(config, "scr_config")) {
    stop(sprintf("%s(): `config` must come from scr_config().", fn), call. = FALSE)
  }
  invisible(TRUE)
}
