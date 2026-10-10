# ============================================================================ #
# rag.R - red / amber / green lights of a score against its reference
# ============================================================================ #

#' Thresholds of the red / amber / green lights
#'
#' The editable table read by [scr_rag()]: one row per metric with its
#' thresholds and the rule that turns a value into a light. Edit a value, or
#' drop a row to leave a metric out, and pass the table as `plan`.
#'
#' @section Rules:
#'
#' Every rule follows one convention: a light is amber or red only when the
#' confidence interval shows the metric beyond the threshold (for a p-value,
#' when the test rejects); with too few events it is `"grey"`.
#'
#' \describe{
#'   \item{`ci`}{With `higher_better`: green when the upper bound of the
#'     interval reaches `green`, red when it stays below `red`, amber in
#'     between; the mirror image on the lower bound when lower is better.}
#'   \item{`threshold`}{The same on the value alone.}
#'   \item{`p_value`}{Red at or below `red`, amber at or below `green`, green
#'     above.}
#'   \item{`interval`}{Two-sided: green when the interval meets
#'     `[green, green_hi]`, red when it lies entirely outside
#'     `[red, red_hi]`, amber otherwise.}
#'   \item{`psi`}{Effect and significance: red when the index reaches `red`
#'     **and** exceeds the n-adjusted critical value (Yurdakul and Naranjo,
#'     2020), amber when it reaches `green` and exceeds it, green otherwise.
#'     Significance alone never colors a light, which on a large sample
#'     would flag every negligible shift.}
#'   \item{`count`}{Green at or below `green`, red at or above `red` (never
#'     when `red` is `NA`), amber in between.}
#'   \item{`none`}{Reported without a light.}
#' }
#'
#' Every threshold is a convention of this package, documented in `note`,
#' except where a source is cited there; adjust them to the validation
#' policy in force.
#'
#' @section Calibration under risk and propensity:
#'
#' Under `objective = "risk"` every calibration check is one-sided:
#' over-prediction (more expected than observed events) is prudent and only
#' under-prediction is penalized. `oe_ratio` then follows the rule `ci` with
#' lower better, lit on the lower bound of its interval: green when it is at
#' or below 1.10, red above 1.25, amber in between; a conservative model
#' (O/E well below 1) stays green. `band_calibration` tests every band
#' against under-prediction only. Under `"propensity"` both directions
#' count: `oe_ratio` follows the two-sided rule `interval` (green when the
#' interval meets `[0.90, 1.10]`, red when it lies outside `[0.80, 1.25]`)
#' and the band tests are two-sided.
#'
#' @param objective `"risk"` or `"propensity"`: the Gini ratio thresholds
#'   are looser under propensity (0.90 / 0.80 against 0.95 / 0.90), and the
#'   calibration checks are one-sided under risk, two-sided under propensity
#'   (see the section Calibration under risk and propensity).
#'
#' @return A `data.frame` with `family`, `metric`, `green`, `red`,
#'   `green_hi`, `red_hi` (upper edges of the two-sided rule), `higher_better`,
#'   `rule` and `note`.
#'
#' @references
#' European Central Bank (2019). *Instructions for reporting the validation
#' results of internal models: IRB Pillar I models for credit risk*. ECB
#' Banking Supervision.
#'
#' Yurdakul, B. and Naranjo, J. (2020). Statistical properties of the
#' population stability index. *Journal of Risk Model Validation*, 14(4),
#' 89-100.
#'
#' @family score-studies
#' @examples
#' plan <- scr_rag_plan()
#' plan[, c("family", "metric", "green", "red", "rule")]
#' # a stricter policy on the score PSI
#' plan$green[plan$metric == "score_psi"] <- 0.05
#' @export
scr_rag_plan <- function(objective = c("risk", "propensity")) {
  objective <- match.arg(objective)
  prop <- identical(objective, "propensity")
  data.frame(
    family = c("discrimination", "discrimination", "discrimination", "calibration", "calibration",
               "stability", "stability", "variables", "variables", "variables"),
    metric = c("gini_ratio", "auc_change_p", "ks", "oe_ratio", "band_calibration", "score_psi", "rank_order",
               "csi", "iv_ratio", "woe_sign_flip"),
    # O/E: one-sided under risk (only under-prediction is penalized), two-sided under propensity
    green = c(if (prop) 0.90 else 0.95, 0.05, NA, if (prop) 0.90 else 1.10, 0, 0.10, 0, 0.10, 0.80, 0),
    red = c(if (prop) 0.80 else 0.90, 0.01, NA, if (prop) 0.80 else 1.25, 2, 0.25, 2, 0.25, 0.50, NA),
    green_hi = c(NA, NA, NA, if (prop) 1.10 else NA, NA, NA, NA, NA, NA, NA),
    red_hi = c(NA, NA, NA, if (prop) 1.25 else NA, NA, NA, NA, NA, NA, NA),
    higher_better = c(TRUE, TRUE, NA, if (prop) NA else FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE),
    rule = c("ci", "p_value", "none", if (prop) "interval" else "ci", "count", "psi", "count", "psi", "threshold",
             "count"),
    note = c(
      "Gini(study) / Gini(reference) with a bootstrap interval, lit on its upper bound; package convention",
      "one-sided S-test of AUC(reference) - AUC(study), two-sample DeLong standard error (the ECB (2019) instructions take the initial AUC as fixed)",
      "KS of the study sample, reported without a light",
      if (prop) "observed / expected events with a Jeffreys interval, two-sided; package convention"
      else "observed / expected events with a Jeffreys interval, lit on its lower bound: only under-prediction (O/E above 1) is penalized; package convention",
      "bands whose Jeffreys p-value is at or below 0.01 (one-sided against under-prediction under risk, two-sided under propensity); package convention",
      "PSI over bands frozen on the reference: 0.10 / 0.25 market convention combined with the n-adjusted critical value (Yurdakul and Naranjo, 2020)",
      "adjacent bands with a Holm-adjusted rank-order reversal (one-sided Fisher exact test, p < 0.05); package convention",
      "CSI per variable over its frozen bins, same rule as the score PSI",
      "IV(study) / IV(reference) per variable; package convention",
      "WOE sign flips in bins holding at least 5% of the study volume; never red; package convention"),
    stringsAsFactors = FALSE)
}

#' Red / amber / green lights of a score against its reference
#'
#' Reads a study sample against the reference the score was developed on
#' and lights every check red, amber or green, in four families:
#' discrimination, calibration, stability and, for a scorecard, the
#' variables. With `by`, one set of lights per period or segment, from the
#' same count table.
#'
#' @section Checks:
#'
#' \describe{
#'   \item{Discrimination}{`gini_ratio`, Gini(study) / Gini(reference), with
#'     a bootstrap interval from independent resamples of both samples (the
#'     count bootstrap of [scr_bands()], exact up to `boot_cells` distinct
#'     scores and an approximation above); `auc_change_p`, the one-sided
#'     p-value of the S-test
#'     \deqn{S = (AUC_{ref} - AUC_{study}) / \sqrt{se_{study}^2 + se_{ref}^2},}
#'     with both DeLong standard errors computed from the counts per score
#'     value. This is the two-sample form: both AUCs are estimated here, on
#'     independent samples. The ECB (2019) instructions take the initial AUC
#'     as fixed (only the current AUC's standard error), which suits a
#'     development AUC taken from documentation; applied to an estimated
#'     reference AUC, that form rejects too often. `ks` is reported without
#'     a light.}
#'   \item{Calibration}{Needs an expected probability: the alignment of a
#'     scorecard, or the `prob` column of a data.frame; otherwise both
#'     lights are `"grey"` ("no expected probability"). `oe_ratio`, observed over
#'     expected events, with the Jeffreys interval of the observed rate;
#'     `band_calibration`, the number of bands whose Jeffreys test rejects
#'     the band mean expected probability. Under risk both checks are
#'     one-sided, against under-prediction only (over-prediction is
#'     prudent); under propensity they are two-sided (see
#'     [scr_rag_plan()]). Every band is also listed, without a light.}
#'   \item{Stability}{`score_psi` over bands frozen on the reference, lit by
#'     effect and significance; `rank_order`, the number of Holm-significant
#'     reversals between adjacent bands.}
#'   \item{Variables}{For a scorecard read on its hold-out against train:
#'     the CSI of every variable over its frozen bins (same rule as the
#'     PSI), `iv_ratio` = IV(study) / IV(reference), and `woe_sign_flip`,
#'     the number of bins holding at least 5% of the study volume whose WOE
#'     changes sign (amber at most).}
#' }
#'
#' One convention holds for every light: a light is amber or red only when
#' the confidence interval shows the metric beyond the threshold (for a
#' p-value, when the test rejects); with too few events it is `"grey"`. The
#' thresholds and rules are in [scr_rag_plan()].
#'
#' A group with fewer than `min_events` events, or non-events, in the study
#' sample gets `"grey"` lights, with the values still reported. Within a
#' family the worst light wins (red, then amber, then green; `"grey"` only
#' when nothing is lit). The overall light is the worst of discrimination
#' and calibration; stability and variables can raise it to amber, never to
#' red, and an overall without any lit discrimination or calibration check
#' is `"grey"`. The `reason` of the summary says how the overall light was
#' formed, for example "calibration not tested" when no expected
#' probability is available.
#'
#' With `by`, a group of the study sample is compared with the same group
#' of the reference when the reference has it (a segment), and with the
#' whole reference otherwise (a new period). Without a sample column, every
#' group is compared with the whole data. The groups are listed in the order
#' of their labels, those of a numeric column in numeric order.
#'
#' @inheritParams scr_bands
#' @param x An object from [scr_scorecard()], [scr_bands()] or [scr_tiers()],
#'   or a `data.frame` with one row per scored case.
#' @param plan Thresholds table, as returned by [scr_rag_plan()]; `NULL`
#'   uses the defaults for the objective.
#' @param sample For a scorecard: the study sample (`"holdout"`). For a
#'   data.frame: the name of a column with sample labels, or `NULL`. For a
#'   score study: the study samples, `NULL` for every sample but the
#'   reference.
#' @param by Name of a column holding periods or segments (for a scorecard,
#'   a column of its scored samples such as `"date"`).
#' @param level Confidence level of the intervals. For a scorecard or a
#'   score study, `NULL` uses `config$study_level` or the level of the study.
#' @param min_events Fewest events, and non-events, of the study group for a
#'   lit result.
#' @param n_bands Bands frozen on the reference for the PSI, the rank order
#'   and the band calibration. For a scorecard, `NULL` uses
#'   `config$score_groups`. A score study uses its own cuts.
#' @param n_boot Bootstrap resamples of the Gini ratio interval. For a
#'   scorecard, `NULL` uses `config$n_boot`.
#' @param prob For a data.frame: optional column with the expected event
#'   probability of every case.
#'
#' @return An object of class `scr_rag`:
#'   \describe{
#'     \item{`table`}{One row per check: `sample`, `group`, `family`,
#'       `metric`, `level` (`"score"`, a variable name or a band label),
#'       `value`, `lo`, `hi`, `benchmark` (the reference value or critical
#'       value it is read against), `light` (`"green"`, `"amber"`, `"red"`,
#'       `"grey"` or `"none"`) and `reason`.}
#'     \item{`summary`}{One row per sample and group: the light of
#'       `discrimination`, `calibration`, `stability`, `variables` and
#'       `overall`, and the `reason` of the overall light.}
#'     \item{`plan`}{The thresholds used.}
#'     \item{`objective`, `direction`, `level`, `min_events`, `reference`,
#'       `study`, `by`, `cuts`, `target`, `call`}{The settings.}
#'   }
#'
#' @references
#' Brown, L. D., Cai, T. T. and DasGupta, A. (2001). Interval estimation for
#' a binomial proportion. *Statistical Science*, 16(2), 101-133.
#' \doi{10.1214/ss/1009213286}
#'
#' DeLong, E. R., DeLong, D. M. and Clarke-Pearson, D. L. (1988). Comparing
#' the areas under two or more correlated receiver operating characteristic
#' curves: a nonparametric approach. *Biometrics*, 44(3), 837-845.
#'
#' European Central Bank (2019). *Instructions for reporting the validation
#' results of internal models: IRB Pillar I models for credit risk*. ECB
#' Banking Supervision.
#'
#' Siddiqi, N. (2006). *Credit Risk Scorecards: Developing and Implementing
#' Intelligent Credit Scoring*. Wiley.
#'
#' Yurdakul, B. and Naranjo, J. (2020). Statistical properties of the
#' population stability index. *Journal of Risk Model Validation*, 14(4),
#' 89-100.
#'
#' @family score-studies
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   use_lightgbm = FALSE, xgb_rounds = 40, n_boot = 20)
#' res <- scr_select(scr_demo, "default", config = cfg, drop = c("id", "churn"),
#'                   date_col = "ref_date")
#' sc <- scr_scorecard(res)
#' rg <- scr_rag(sc)
#' rg
#' rg$summary
#'
#' # a data.frame with a sample column and an expected probability
#' d <- rbind(data.frame(sample = "dev", sc$samples$train[, c("score", "y", "prob")]),
#'            data.frame(sample = "new", sc$samples$holdout[, c("score", "y", "prob")]))
#' scr_rag(d, prob = "prob", sample = "sample", n_boot = 50)$summary
#' @export
scr_rag <- function(x, ...) UseMethod("scr_rag")

#' @rdname scr_rag
#' @export
scr_rag.scr_scorecard <- function(x, plan = NULL, sample = "holdout", reference = "train", by = NULL,
                                  level = NULL, min_events = 20L, n_bands = NULL, n_boot = NULL, seed = NULL,
                                  max_cells = 1e5, boot_cells = 1e4, ...) {
  .study_dots(list(...), "scr_rag")
  cfg <- x$config
  .study_chr1(sample, "sample", "scr_rag")
  inp <- .study_input(x, sample = sample, reference = reference, by = by, max_cells = max_cells, fn = "scr_rag")
  n_bands <- .study_whole(n_bands %||% cfg$score_groups %||% 10L, "n_bands", "scr_rag", lower = 1)
  cuts <- .study_cuts(.study_cells(inp$hist, reference), n_bands, "uniform", NULL, .study_side(x$direction))$cuts
  # variable checks need the frozen bins of the hold-out against the training counts
  vars <- if (identical(reference, "train") && identical(sample, "holdout") && !is.null(x$holdout_bins)) {
    grp <- if (is.null(by)) NULL else as.character(x$samples$holdout[[by]])
    function(g) .rag_vars_sc(x, if (is.null(grp)) seq_len(nrow(x$samples$holdout)) else which(grp %in% g),
                             cfg$psi_alpha %||% 0.05)
  } else function(g) NULL
  .rag_run(inp, cuts, plan, level %||% cfg$study_level %||% 0.95, min_events, n_boot %||% cfg$n_boot,
           seed %||% cfg$seed, cfg$psi_alpha %||% 0.05, vars, by, sys.call(), var_scope = TRUE,
           boot_cells = boot_cells)
}

#' @rdname scr_rag
#' @export
scr_rag.data.frame <- function(x, score = "score", y = "y", prob = NULL, objective = "risk", direction = NULL,
                               weight = NULL, sample = NULL, reference = NULL, study = NULL, by = NULL,
                               plan = NULL, level = 0.95, min_events = 20L, n_bands = 10L, n_boot = 200L,
                               seed = NULL, max_cells = 1e5, boot_cells = 1e4, ...) {
  .study_dots(list(...), "scr_rag")
  inp <- .study_input(x, score = score, y = y, objective = objective, direction = direction, weight = weight,
                      prob = prob, sample = sample, reference = reference, study = study, by = by,
                      max_cells = max_cells, fn = "scr_rag")
  n_bands <- .study_whole(n_bands, "n_bands", "scr_rag", lower = 1)
  ref_cells <- .study_cells(inp$hist, inp$reference)
  if (!nrow(ref_cells)) stop("scr_rag(): the reference sample '", inp$reference, "' has no scored row.", call. = FALSE)
  cuts <- .study_cuts(ref_cells, n_bands, "uniform", NULL, .study_side(inp$direction))$cuts
  .rag_run(inp, cuts, plan, level, min_events, n_boot, seed, 0.05, function(g) NULL, by, sys.call(),
           boot_cells = boot_cells)
}

#' @rdname scr_rag
#' @export
scr_rag.scr_study <- function(x, plan = NULL, sample = NULL, level = NULL, min_events = 20L, n_boot = 200L,
                              seed = NULL, boot_cells = 1e4, ...) {
  .study_dots(list(...), "scr_rag")
  study <- sample %||% setdiff(x$samples, x$reference)
  if (!length(study)) study <- x$reference
  .study_chr(study, "sample", "scr_rag")
  bad <- setdiff(study, x$samples)
  if (length(bad)) stop("scr_rag(): sample(s) ", lst(bad), " not in the study.", call. = FALSE)
  inp <- list(hist = x$hist, reference = x$reference, samples = unique(c(x$reference, study)),
              objective = x$objective, direction = x$direction, target = x$target,
              meta = list(has_prob = "ep" %in% names(x$hist)))
  .rag_run(inp, x$cuts, plan, level %||% x$level, min_events, n_boot, seed, 0.05, function(g) NULL, NULL,
           sys.call(), boot_cells = boot_cells)
}

#' Read and check a thresholds table
#' @keywords internal
#' @noRd
.rag_read_plan <- function(plan, objective) {
  if (is.null(plan)) return(data.table::as.data.table(scr_rag_plan(objective)))
  if (!is.data.frame(plan)) stop("scr_rag(): `plan` must be a data.frame from scr_rag_plan().", call. = FALSE)
  pl <- data.table::as.data.table(plan)
  need <- c("family", "metric", "green", "red", "higher_better", "rule")
  miss <- setdiff(need, names(pl))
  if (length(miss)) stop("scr_rag(): `plan` lacks the column(s) ", lst(miss), ".", call. = FALSE)
  if (!"green_hi" %in% names(pl)) pl[, green_hi := NA_real_]
  if (!"red_hi" %in% names(pl)) pl[, red_hi := NA_real_]
  if (!"note" %in% names(pl)) pl[, note := ""]
  ref <- scr_rag_plan(objective)
  key <- paste(pl$family, pl$metric)
  bad <- setdiff(key, paste(ref$family, ref$metric))
  if (length(bad)) stop("scr_rag(): unknown check(s) in `plan`: ", lst(bad), ".", call. = FALSE)
  if (anyDuplicated(key)) stop("scr_rag(): a check appears twice in `plan`.", call. = FALSE)
  for (cn in c("green", "red", "green_hi", "red_hi")) {
    if (!is.numeric(pl[[cn]]) && !all(is.na(pl[[cn]]))) stop("scr_rag(): `plan$", cn, "` must be numeric.", call. = FALSE)
    data.table::set(pl, j = cn, value = as.double(pl[[cn]]))
  }
  ok_rule <- c("ci", "threshold", "p_value", "interval", "psi", "count", "none")
  if (!all(pl$rule %in% ok_rule)) stop("scr_rag(): `plan$rule` must be one of ", lst(ok_rule), ".", call. = FALSE)
  pl[, higher_better := as.logical(higher_better)]
  # a lit rule needs its thresholds
  lit <- pl$rule != "none"
  if (any(lit & is.na(pl$green)) || any(pl$rule %in% c("ci", "threshold", "p_value", "interval", "psi") & is.na(pl$red)) ||
      any(pl$rule == "interval" & (is.na(pl$green_hi) | is.na(pl$red_hi)))) {
    stop("scr_rag(): every lit check in `plan` needs its thresholds.", call. = FALSE)
  }
  pl[]
}

#' Light and reason of one value under one plan row
#' @keywords internal
#' @noRd
.rag_rule <- function(row, value, lo = NA_real_, hi = NA_real_, bench = NA_real_) {
  rule <- row$rule
  if (identical(rule, "none")) return(list(light = "none", reason = "reported"))
  if (is.na(value)) return(list(light = "grey", reason = "not computable"))
  g <- row$green; r <- row$red
  f <- function(v) if (is.na(v)) "NA" else .g3(as.double(v))
  switch(rule,
    ci = , threshold = {
      hb <- !isFALSE(row$higher_better)
      if (is.na(lo) || identical(rule, "threshold")) lo <- value
      if (is.na(hi) || identical(rule, "threshold")) hi <- value
      light <- .rag_light(value, lo, hi, g, r, hb)
      # the bound that can show a deviation: the upper one when higher is better
      bnd <- if (rule == "ci") (if (hb) "upper bound" else "lower bound") else "value"
      b <- if (hb) hi else lo
      reason <- if (hb) switch(light, green = sprintf("%s %s >= %s", bnd, f(b), f(g)),
                                red = sprintf("%s %s < %s", bnd, f(b), f(r)),
                                sprintf("%s %s below %s, not below %s", bnd, f(b), f(g), f(r)))
                else switch(light, green = sprintf("%s %s <= %s", bnd, f(b), f(g)),
                            red = sprintf("%s %s > %s", bnd, f(b), f(r)),
                            sprintf("%s %s above %s, not above %s", bnd, f(b), f(g), f(r)))
      list(light = light, reason = reason)
    },
    p_value = {
      light <- if (value <= r) "red" else if (value <= g) "amber" else "green"
      list(light = light, reason = sprintf("p = %s (red at or below %s, amber at or below %s)", f(value), f(r), f(g)))
    },
    interval = {
      if (is.na(lo)) lo <- value
      if (is.na(hi)) hi <- value
      light <- if (hi >= g && lo <= row$green_hi) "green" else if (hi < r || lo > row$red_hi) "red" else "amber"
      list(light = light, reason = switch(light,
        green = sprintf("interval meets [%s, %s]", f(g), f(row$green_hi)),
        red = sprintf("interval outside [%s, %s]", f(r), f(row$red_hi)),
        sprintf("interval misses [%s, %s] but meets [%s, %s]", f(g), f(row$green_hi), f(r), f(row$red_hi))))
    },
    psi = {
      sig <- !is.na(bench) && value > bench
      light <- if (value >= r && sig) "red" else if (value >= g && sig) "amber" else "green"
      size <- if (value >= r) sprintf(">= %s", f(r)) else if (value >= g) sprintf(">= %s", f(g)) else sprintf("< %s", f(g))
      list(light = light, reason = sprintf("index %s %s, %s the n-adjusted critical value %s", f(value), size,
                                           if (sig) "above" else "not above", f(bench)))
    },
    count = {
      light <- if (value <= g) "green" else if (!is.na(r) && value >= r) "red" else "amber"
      list(light = light, reason = sprintf("%d found (green at or below %s%s)", as.integer(value), f(g),
                                           if (is.na(r)) "" else sprintf(", red at %s or more", f(r))))
    })
}

#' Run the checks on every study sample and group of a count table
#' @keywords internal
#' @noRd
.rag_run <- function(inp, cuts, plan, level, min_events, n_boot, seed, psi_alpha, vars, by, call, var_scope = FALSE,
                     boot_cells = 1e4) {
  fn <- "scr_rag"
  level <- .study_level(level, fn)
  min_events <- .study_whole(min_events, "min_events", fn)
  n_boot <- .study_whole(n_boot, "n_boot", fn)
  boot_cells <- .study_boot_cells(boot_cells, fn)
  pl <- .rag_read_plan(plan, inp$objective)
  h <- inp$hist
  reference <- inp$reference
  study <- setdiff(inp$samples, reference)
  self <- !length(study)
  if (self) study <- reference
  grouped <- "group" %in% names(h)
  .scr_local_seed(seed)
  ref_all <- .study_cells(h, reference)
  units <- if (grouped) unique(h[h[["sample"]] %in% study, list(sample, group)])[order(sample, group)] else
    data.table::data.table(sample = study, group = "all")
  # the groups of a numeric `by` in numeric order (their labels sort "10" before "9")
  glv <- inp$group_levels
  if (grouped && !is.null(glv)) {
    pos <- order(match(units$sample, unique(units$sample)), match(units$group, glv))
    units <- units[pos]
  }
  # the reference results are reused across groups compared with the same reference
  ref_memo <- new.env(parent = emptyenv())
  rows <- list(); summ <- list()
  for (i in seq_len(nrow(units))) {
    s <- units$sample[i]; g <- units$group[i]
    hs <- .study_cells(h, s, if (grouped) g)
    same_grp <- grouped && !self && any(h[["sample"]] == reference & h[["group"]] %in% g)
    rkey <- if (same_grp) paste0("g:", g) else "all"
    hr <- if (same_grp) .study_cells(h, reference, g) else ref_all
    vr <- vars(if (grouped) g else NULL)
    out <- .rag_unit(hs, hr, cuts, pl, inp$objective, inp$direction, level, min_events, n_boot, psi_alpha,
                     isTRUE(inp$meta$has_prob), vr, var_scope, ref_memo, rkey, boot_cells)
    out$table[, `:=`(sample = s, group = g)]
    rows[[i]] <- out$table
    summ[[i]] <- data.table::data.table(sample = s, group = g, data.table::as.data.table(out$summary))
  }
  tab <- data.table::rbindlist(rows, use.names = TRUE)
  data.table::setcolorder(tab, c("sample", "group", "family", "metric", "level", "value", "lo", "hi", "benchmark",
                                 "light", "reason"))
  structure(list(table = tab[], summary = data.table::rbindlist(summ), plan = as.data.frame(pl),
                 objective = inp$objective, direction = inp$direction, level = level, min_events = min_events,
                 reference = reference, study = study, by = by, cuts = cuts, target = inp$target, call = call),
            class = c("scr_rag", "list"))
}

#' All checks of one study group against its reference
#' @keywords internal
#' @noRd
.rag_unit <- function(hs, hr, cuts, pl, objective, direction, level, min_events, n_boot, psi_alpha, has_prob,
                      vrows, var_scope, ref_memo, rkey, boot_cells = 1e4) {
  rows <- list()
  add <- function(family, metric, lev, value, lo = NA_real_, hi = NA_real_, bench = NA_real_, light = NULL,
                  reason = NULL, note = "") {
    # a check left out of the plan is not computed
    i <- match(metric, pl[["metric"]])
    if (is.na(i)) return(invisible(NULL))
    row <- pl[i]
    lr <- if (is.null(light)) .rag_rule(row, value, lo, hi, bench) else list(light = light, reason = reason)
    if (!is.na(note) && nzchar(note)) lr$reason <- paste0(lr$reason, "; ", note)
    rows[[length(rows) + 1L]] <<- data.table::data.table(
      family = family, metric = metric, level = lev, value = as.double(value), lo = as.double(lo),
      hi = as.double(hi), benchmark = as.double(bench), light = lr$light, reason = lr$reason)
  }
  E <- sum(hs$e_raw); NE <- sum(hs$n_y_raw) - E

  # -- discrimination ----------------------------------------------------- #
  want_boot <- any(pl$metric == "gini_ratio")
  nb <- if (want_boot) n_boot else 0L
  ds <- .study_discrimination(hs, direction, nb, level, keep = TRUE, boot_cells = boot_cells)
  dr <- ref_memo[[rkey]]
  if (is.null(dr)) {
    dr <- .study_discrimination(hr, direction, nb, level, keep = TRUE, boot_cells = boot_cells)
    ref_memo[[rkey]] <- dr
  }
  ratio <- if (is.finite(dr$gini) && dr$gini > 0) ds$gini / dr$gini else NA_real_
  rlo <- rhi <- NA_real_
  if (!is.na(ratio) && length(ds$boot_auc) && length(dr$boot_auc)) {
    m <- min(length(ds$boot_auc), length(dr$boot_auc))
    den <- 2 * dr$boot_auc[seq_len(m)] - 1
    rb <- (2 * ds$boot_auc[seq_len(m)] - 1) / den
    rb <- rb[den > 0]
    if (length(rb) >= 2L) {
      q <- stats::quantile(rb, c((1 - level) / 2, 1 - (1 - level) / 2), names = FALSE)
      rlo <- q[1]; rhi <- q[2]
    }
  }
  add("discrimination", "gini_ratio", "score", ratio, rlo, rhi, dr$gini)
  st <- .rag_auc_change(hs, hr, direction)
  add("discrimination", "auc_change_p", "score", st$p, bench = st$auc_ref)
  add("discrimination", "ks", "score", ds$ks, ds$ks_lo, ds$ks_hi, dr$ks)

  # -- calibration -------------------------------------------------------- #
  tb <- .study_band_table(hs, cuts, .study_ref_shares(hr, cuts), level, direction, psi_alpha)
  if (!has_prob) {
    add("calibration", "oe_ratio", "score", NA_real_, light = "grey", reason = "no expected probability")
    add("calibration", "band_calibration", "score", NA_real_, light = "grey", reason = "no expected probability")
  } else {
    O <- sum(hs$e); X <- sum(hs$ep); NY <- sum(hs$n_y)
    neff <- .study_kish(NY, sum(hs$w2_y))
    ci <- .study_jeffreys(if (NY > 0) O / NY * neff else NA_real_, neff, level)
    oe <- if (X > 0) O / X else NA_real_
    add("calibration", "oe_ratio", "score", oe, ci$lo * NY / X, ci$hi * NY / X, X)
    raw <- attr(tb, "raw")
    pd <- ifelse(raw$n_y > 0, raw$ep / raw$n_y, NA_real_)
    nb <- .study_kish(raw$n_y, raw$w2_y)
    xb <- tb$rate * nb
    ok <- !is.na(pd) & !is.na(xb) & nb > 0
    p1 <- rep(NA_real_, length(pd))
    # posterior probability under the Jeffreys prior that the true rate is below the expected one
    p1[ok] <- stats::pbeta(pd[ok], xb[ok] + 0.5, nb[ok] - xb[ok] + 0.5)
    pb <- if (identical(objective, "propensity")) 2 * pmin(p1, 1 - p1) else p1
    fail <- !is.na(pb) & pb <= 0.01
    for (b in seq_len(nrow(tb))) {
      add("calibration", "band_calibration", tb$label[b], tb$rate[b], tb$rate_lo[b], tb$rate_hi[b], pd[b],
          light = "none", reason = if (is.na(pb[b])) "no outcome in the band" else
            sprintf("p = %s%s", .g3(pb[b]), if (fail[b]) ", rejected" else ""))
    }
    add("calibration", "band_calibration", "score", if (any(ok)) sum(fail) else NA_real_)
  }

  # -- stability ---------------------------------------------------------- #
  ps <- attr(tb, "psi")
  add("stability", "score_psi", "score", if (is.null(ps)) NA_real_ else ps$psi, bench = if (is.null(ps)) NA_real_ else ps$critical)
  pr <- tb$p_reversal_adj
  add("stability", "rank_order", "score", if (all(is.na(pr))) NA_real_ else sum(pr < 0.05, na.rm = TRUE))

  # -- variables ---------------------------------------------------------- #
  if (!is.null(vrows)) {
    for (k in seq_len(nrow(vrows))) {
      add("variables", vrows$metric[k], vrows$level[k], vrows$value[k], bench = vrows$benchmark[k],
          note = vrows$reason[k])
    }
  } else if (var_scope) {
    for (m in c("csi", "iv_ratio", "woe_sign_flip")) {
      add("variables", m, "score", NA_real_, light = "grey", reason = "bin counts are kept for train against holdout only")
    }
  }

  tab <- if (length(rows)) data.table::rbindlist(rows) else
    data.table::data.table(family = character(), metric = character(), level = character(), value = numeric(),
                           lo = numeric(), hi = numeric(), benchmark = numeric(), light = character(), reason = character())
  # too few events: every light of the group is "grey", the values stay
  if (E < min_events || NE < min_events) {
    lit <- tab$light != "none"
    tab[lit, `:=`(light = "grey", reason = sprintf("too few events (%s events, %s non-events; %d needed)",
                                                   n_fmt(E), n_fmt(NE), min_events))]
  }
  fam <- function(f) .rag_worst(tab$light[tab$family == f])
  sm <- list(discrimination = fam("discrimination"), calibration = fam("calibration"),
             stability = fam("stability"), variables = if (any(tab$family == "variables")) fam("variables") else NA_character_)
  prim <- .rag_worst(c(sm$discrimination, sm$calibration))
  sec <- .rag_worst(c(sm$stability, sm$variables))
  raised <- identical(prim, "green") && sec %in% c("red", "amber")
  sm$overall <- if (raised) "amber" else prim
  # how the overall light was formed
  why <- if (E < min_events || NE < min_events) "too few events" else if (identical(prim, "grey")) {
    "no discrimination or calibration check lit"
  } else c(if (identical(sm$discrimination, "grey")) "discrimination not tested",
           if (identical(sm$calibration, "grey")) "calibration not tested",
           if (raised) sprintf("raised to amber by %s",
                               paste(c("stability", "variables")[c(sm$stability, sm$variables) %in% c("red", "amber")],
                                     collapse = " and ")))
  sm$reason <- paste(why, collapse = "; ")
  list(table = tab, summary = sm)
}

#' One-sided S-test of an AUC drop between two independent samples
#'
#' \eqn{S = (AUC_{ref} - AUC_{study}) / \sqrt{se_{study}^2 + se_{ref}^2}},
#' both standard errors by DeLong from the counts per score value, and
#' `p = 1 - Phi(S)`. Taking the reference AUC as fixed (`se_ref = 0`), as
#' the ECB instructions do for a documented initial AUC, rejects too often
#' when that AUC is itself estimated on a sample.
#' @keywords internal
#' @noRd
.rag_auc_change <- function(hs, hr, direction) {
  one <- function(cells) {
    c1 <- cells$e; c0 <- pmax(cells$n_y - cells$e, 0)
    # cells are in ascending score; a higher score means more events only under higher_is_riskier
    if (!identical(direction, "higher_is_riskier")) { c1 <- rev(c1); c0 <- rev(c0) }
    ok <- sum(c1) > 0 && sum(c0) > 0
    list(auc = if (ok) .auc_ks_counts(c1, c0)$auc else NA_real_,
         se = .study_delong_counts(c1, c0, sum(cells$e_raw), sum(cells$n_y_raw - cells$e_raw)))
  }
  s <- one(hs); r <- one(hr)
  se <- sqrt(s$se^2 + r$se^2)
  z <- if (is.finite(se) && se > 0 && is.finite(s$auc) && is.finite(r$auc)) (r$auc - s$auc) / se else NA_real_
  list(p = if (is.na(z)) NA_real_ else stats::pnorm(z, lower.tail = FALSE), z = z, se = se,
       auc_study = s$auc, auc_ref = r$auc)
}

#' Variable checks of a scorecard on hold-out rows against its training bins
#' @keywords internal
#' @noRd
.rag_vars_sc <- function(x, rows, psi_alpha) {
  y <- x$samples$holdout$y[rows]
  out <- list()
  for (f in x$features) {
    pt <- x$points[x$points$variable == f]
    r <- x$fit$results[[f]]
    B <- nrow(pt)
    idx <- x$holdout_bins[[f]][rows]
    cnt <- tabulate(idx, nbins = B)
    ev <- tabulate(idx[!is.na(y) & y == 1L], nbins = B)
    ps <- .psi_counts(pt$count_train, cnt, pt$bin, psi_alpha, c(0.10, 0.25))
    m <- match(pt$bin, r$bin)
    re <- as.double(r$count_pos[m])
    rn <- if (!is.null(r$count_neg)) as.double(r$count_neg[m]) else as.double(r$count[m]) - re
    # the same smoothing on both samples, so the ratio compares like with like
    bw_r <- .band_woe(re, rn); bw_s <- .band_woe(ev, cnt - ev)
    iv <- function(bw) { t <- (bw$pct_event - bw$pct_nonevent) * bw$log_odds; if (all(is.na(t))) NA_real_ else sum(t, na.rm = TRUE) }
    iv_r <- iv(bw_r); iv_s <- iv(bw_s)
    big <- sum(cnt) > 0 & cnt / max(1, sum(cnt)) >= 0.05
    flip <- big & !is.na(bw_r$log_odds) & !is.na(bw_s$log_odds) & sign(bw_r$log_odds) * sign(bw_s$log_odds) < 0
    out[[f]] <- data.table::data.table(
      metric = c("csi", "iv_ratio", "woe_sign_flip"), level = f,
      value = c(ps$psi, if (is.finite(iv_r) && iv_r > 0) iv_s / iv_r else NA_real_,
                if (sum(cnt) > 0 && !all(is.na(bw_s$log_odds))) sum(flip) else NA_real_),
      benchmark = c(ps$critical, iv_r, NA_real_),
      reason = c("", "", if (any(flip)) paste("bins", paste(pt$bin[flip], collapse = ", ")) else ""))
  }
  data.table::rbindlist(out)
}

#' Workbook sheets of a set of lights
#' @keywords internal
#' @noRd
.rag_sheets <- function(x) list(Lights = x$table, Light_Summary = x$summary, Light_Plan = x$plan)

#' @rdname scr_export
#' @export
scr_export.scr_rag <- function(x, dir, stamp = TRUE, ...) {
  .study_dots(list(...), "scr_export")
  .need_openxlsx()
  out_dir <- .export_dir(dir, stamp)
  tag <- .file_tag(x$target)
  settings <- .kv_table(list(target = x$target, objective = x$objective, direction = x$direction,
                             reference = x$reference, study = x$study, by = x$by %||% "none", level = x$level,
                             min_events = x$min_events))
  sheets <- lapply(c(.rag_sheets(x), list(Settings = settings)), .study_sheet)
  files <- list(xlsx = .scr_write_xlsx(sheets, file.path(out_dir, sprintf("rag_%s.xlsx", tag))))
  for (f in files) msg("  %s", f)
  x$files <- files
  invisible(x)
}

#' @export
print.scr_rag <- function(x, ...) {
  cat(sprintf("<scr_rag> target \"%s\" | objective %s | %s | reference '%s'\n", x$target, x$objective,
              x$direction, x$reference))
  s <- x$summary
  cat(sprintf("  %-10s %-12s %-15s %-12s %-10s %-10s %-8s\n", "sample", "group", "discrimination", "calibration",
              "stability", "variables", "overall"))
  for (i in seq_len(nrow(s))) {
    why <- if (is.null(s$reason) || !nzchar(s$reason[i])) "" else paste0("  (", s$reason[i], ")")
    cat(sprintf("  %-10s %-12s %-15s %-12s %-10s %-10s %-8s%s\n", substr(s$sample[i], 1, 10), substr(s$group[i], 1, 12),
                s$discrimination[i], s$calibration[i], s$stability[i],
                if (is.na(s$variables[i])) "-" else s$variables[i], toupper(s$overall[i]), why))
  }
  t <- x$table[!x$table$light %in% c("green", "none")]
  if (nrow(t)) {
    cat("\nNot green\n")
    for (i in seq_len(nrow(t))) {
      cat(sprintf("  %-6s %-10s %-12s %-15s %-18s %10s  %s\n", t$light[i], substr(t$sample[i], 1, 10),
                  substr(t$group[i], 1, 12), t$metric[i], substr(t$level[i], 1, 18), .study_f(t$value[i], "%.4f"),
                  t$reason[i]))
    }
  }
  invisible(x)
}

# data.table column names used without quotes in this file
utils::globalVariables(c("green_hi", "red_hi", "note", "higher_better", "reason", "group", "sample", "light"))
