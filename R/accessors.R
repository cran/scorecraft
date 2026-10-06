# ============================================================================ #
# accessors.R - the reading interface of the selection result
# ============================================================================ #

#' @keywords internal
#' @noRd
check_result <- function(x, fn) {
  if (!inherits(x, "scr_result")) stop(sprintf("%s() expects an object from scr_select().", fn), call. = FALSE)
  invisible(TRUE)
}

#' Variables approved for the scorecard
#'
#' The final shortlist, in consensus order (the first is the strongest). It
#' is exactly the list [scr_sql()] covers and [scr_scorecard()] fits.
#'
#' After [scr_classing_apply()] the result carries two lists: the automatic
#' consensus and the analyst's final choice; `which` picks one, and the
#' default is the final one so that every downstream function follows the
#' analyst's decision.
#'
#' @param x An object from [scr_select()].
#' @param which `"final"` (default), `"consensus"` or `"manual"` (`NULL`
#'   when no manual choice was made).
#'
#' @return A character vector of column names.
#'
#' @family accessors
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   xgb_rounds = 60, n_boot = 20)
#' res <- scr_select(scr_demo, "default", config = cfg, drop = "id",
#'                   date_col = "ref_date")
#' scr_selected(res)
#' @export
scr_selected <- function(x, which = c("final", "consensus", "manual")) {
  check_result(x, "scr_selected")
  which <- match.arg(which)
  if (is.null(x$lab)) return(if (identical(which, "manual")) NULL else x$consensus$selected)
  x$lab$shortlist[[which]]
}

#' Audit funnel: every input variable and its fate
#'
#' The central deliverable. One row per input column, plus one per derived
#' variable, with the descriptive profile, the verdict of every gate, the
#' votes of every model and `exit_stage`, the exact stage at which the
#' variable failed. No candidate disappears from the report.
#'
#' @section Values of `exit_stage`:
#'
#' \describe{
#'   \item{`00.config`}{Never competed: it was in `drop`.}
#'   \item{`01.triage`}{Constant, near-constant, high cardinality, exact
#'     duplicate, missing share above the ceiling, or no signal in the coarse IV.}
#'   \item{`02.binning`}{The binning algorithm failed on this column.}
#'   \item{`03.screening`}{Failed one of the eight admission rules.}
#'   \item{`04.holdout`}{IV dropped out of sample, unstable PSI, or part of
#'     the hold-out falls in no bin.}
#'   \item{`05.correlation`}{Redundant with a better-ranked variable.}
#'   \item{`05b.derived_excluded`}{Passed everything, but is a column the
#'     pipeline created and `allow_derived_final = FALSE`.}
#'   \item{`06.consensus`}{Not enough votes, or outside the top-N.}
#'   \item{`07.approved`}{Entered the shortlist.}
#'   \item{`08.manual_drop`}{In the automatic consensus, but removed by the
#'     analyst in the coarse classing lab ([scr_classing_choose()]).}
#' }
#'
#' @param x An object from [scr_select()].
#' @param only_selected If `TRUE`, returns only the approved ones.
#' @param cols `"essentials"` (default) gives a lean view; `"all"` gives
#'   everything; or pass a vector of names.
#'
#' @return A `data.table` ordered with the approved first.
#'
#' @family accessors
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   xgb_rounds = 60, n_boot = 20)
#' res <- scr_select(scr_demo, "default", config = cfg, drop = "id",
#'                   date_col = "ref_date")
#' head(scr_funnel(res, only_selected = TRUE))
#' table(scr_funnel(res, cols = "all")$exit_stage)
#' @export
scr_funnel <- function(x, only_selected = FALSE, cols = "essentials") {
  check_result(x, "scr_funnel")
  d <- x$funnel
  if (isTRUE(only_selected)) d <- d[approved == TRUE]
  if (identical(cols, "essentials")) {
    ess <- c("feature", "derived_from", "type", "approved", "exit_stage", "consensus_rank",
             "consensus_score", "votes", "n_bins", "total_iv", "iv_holdout", "ks", "psi",
             "psi_flag_adjusted", "iv_suspect", "triage_reason", "screen_reason",
             "holdout_reason", "prune_corr_with")
    d <- d[, intersect(ess, names(d)), with = FALSE]
  } else if (!identical(cols, "all")) {
    missing <- setdiff(cols, names(d))
    if (length(missing)) stop("column(s) not in the funnel: ", lst(missing), call. = FALSE)
    d <- d[, cols, with = FALSE]
  }
  d[]
}

#' Gains table, at bin level
#'
#' One row per variable and bin, with counts, event rate, WOE, IV, lift,
#' cumulative KS, precision and recall, plus the hold-out IV and the PSI of
#' the variable.
#'
#' @param x An object from [scr_select()].
#' @param only_selected If `TRUE` (default), only the approved variables.
#'
#' @return A `data.table` at bin level.
#'
#' @family accessors
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   xgb_rounds = 60, n_boot = 20)
#' res <- scr_select(scr_demo, "default", config = cfg, drop = "id",
#'                   date_col = "ref_date")
#' g <- scr_gains(res)
#' g[feature == scr_selected(res)[1], .(bin, count, pos_rate, woe, iv)]
#' @export
scr_gains <- function(x, only_selected = TRUE) {
  check_result(x, "scr_gains")
  d <- x$gains
  if (isTRUE(only_selected)) d <- d[approved == TRUE]
  d[]
}

#' Leakage and suspicious-strength audit
#'
#' Separates what the pipeline failed for excessive strength
#' (`IV_SUSPICIOUS`, `DEGENERATE_BIN`) from what it admitted but deserves a
#' second look (IV above `config$iv_suspect`). A bin with no events or no
#' non-events is the symptom with no innocent explanation: the variable
#' determines the outcome on part of the population.
#'
#' @param x An object from [scr_select()].
#' @param threshold Warning threshold. `NULL` uses `config$iv_suspect`.
#'
#' @return An `scr_leakage` object (a list with `barred`, `degenerate`,
#'   `approved_suspect`), with a print method.
#'
#' @family accessors
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   xgb_rounds = 60, n_boot = 20)
#' res <- scr_select(scr_demo, "default", config = cfg, drop = "id",
#'                   date_col = "ref_date")
#' scr_leakage(res)
#' @export
scr_leakage <- function(x, threshold = NULL) {
  check_result(x, "scr_leakage")
  lim <- threshold %||% x$config$iv_suspect
  f <- x$funnel
  barred <- f[!is.na(total_iv) & grepl("IV_SUSPICIOUS", screen_reason)][order(-total_iv)]
  degen  <- f[!is.na(total_iv) & grepl("DEGENERATE_BIN", screen_reason)][order(-total_iv)]
  appr   <- f[approved == TRUE & !is.na(total_iv) & total_iv >= lim][order(-total_iv)]
  sel <- function(d) d[, intersect(c("feature", "type", "n_bins", "total_iv", "iv_holdout", "ks",
                                     "min_bin_pct", "n_degenerate_bins", "screen_reason"), names(d)), with = FALSE]
  structure(list(target = x$target, threshold = lim, iv_max = x$config$iv_max,
                 barred = sel(barred), degenerate = sel(degen), approved_suspect = sel(appr),
                 n_approved = length(scr_selected(x))), class = c("scr_leakage", "list"))
}

#' @export
print.scr_leakage <- function(x, ...) {
  cat(sprintf("<scr_leakage> target \"%s\"\n", x$target))
  cat(sprintf("  admission ceiling (iv_max): %s | warning threshold: %.2f\n\n", x$iv_max, x$threshold))
  cat(sprintf("Barred for IV above the ceiling: %d\n", nrow(x$barred)))
  if (nrow(x$barred)) {
    d <- utils::head(x$barred, 10)
    for (i in seq_len(nrow(d))) cat(sprintf("  %-46s IV %8.3f  KS %.3f\n", d$feature[i], d$total_iv[i], d$ks[i]))
  }
  if (nrow(x$degenerate)) {
    cat(sprintf("\nWith a degenerate bin (no events or no non-events): %d\n", nrow(x$degenerate)))
    d <- utils::head(x$degenerate, 10)
    for (i in seq_len(nrow(d))) cat(sprintf("  %-46s IV %8.3f\n", d$feature[i], d$total_iv[i]))
  }
  cat(sprintf("\nApproved with IV >= %.2f: %d of %d\n", x$threshold, nrow(x$approved_suspect), x$n_approved))
  if (nrow(x$approved_suspect)) {
    d <- x$approved_suspect
    for (i in seq_len(nrow(d))) cat(sprintf("  %-46s IV %8.3f  IV_hold %7.3f\n", d$feature[i], d$total_iv[i], d$iv_holdout[i]))
    cat("\n  Are these variables available at decision time, or measured in the\n",
        " same window as the outcome? If contemporaneous, they belong in `drop`.\n", sep = "")
  }
  invisible(x)
}
