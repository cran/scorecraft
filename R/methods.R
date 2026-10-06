# ============================================================================ #
# methods.R - S3 methods of the selection result and the run set
# ============================================================================ #

#' Result of a selection
#'
#' Object returned by [scr_select()]. The methods below are the supported
#' way of inspecting the result in the console; to extract data, use the
#' accessors ([scr_selected()], [scr_funnel()], [scr_gains()]).
#'
#' @param x,object An `scr_result` object.
#' @param ... Ignored, present for compatibility with the generic.
#'
#' @return `print()` and `plot()` return `x` invisibly; `summary()` returns
#'   an `scr_summary` object; `as.data.frame()` returns the funnel.
#'
#' @name scr_result
#' @family accessors
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   xgb_rounds = 60, n_boot = 20)
#' res <- scr_select(scr_demo, "default", config = cfg, drop = "id",
#'                   date_col = "ref_date")
#' res                      # print: the funnel in one screen
#' summary(res)             # full text report
#' head(as.data.frame(res)) # the funnel as a data.frame
#' plot(res)
NULL

#' @rdname scr_result
#' @export
print.scr_result <- function(x, ...) {
  m <- x$meta
  cat(sprintf("<scr_result> target \"%s\"\n", x$target))
  cat(sprintf("  %s rows (train %s / hold-out %s) | split %s%s\n",
              n_fmt(m$n_total), n_fmt(m$n_train), n_fmt(m$n_holdout), m$split_method,
              if (!is.na(m$split_cutoff)) sprintf(" at %s", m$split_cutoff) else ""))
  cat(sprintf("  event: %s on train, %s on hold-out | %.1fs\n",
              fmt_pct(m$event_rate_train, 2), fmt_pct(m$event_rate_holdout, 2), m$seconds))
  cat(sprintf("  convention: %s%s\n\n",
              if (identical(m$objective, "risk")) "risk (target=1 is the bad case)" else "propensity (target=1 is the good case)",
              if (isTRUE(m$event$inverted)) " | TARGET INVERTED (event_level = 0)" else ""))

  et <- c(m$n_candidates, m$n_after_triage, m$n_binned, m$n_after_screening, m$n_after_holdout,
          m$n_after_prune, length(x$consensus$selected))
  lab <- c("candidates", "1. triage", "2. binning", "3. screening", "4. hold-out", "5. correlation", "6. consensus")
  if (!is.null(x$lab)) { et <- c(et, length(scr_selected(x))); lab <- c(lab, "7. manual") }
  cat("Funnel\n")
  for (i in seq_along(et)) {
    bar <- strrep("#", max(0L, round(28 * et[i] / max(et[1], 1))))
    cat(sprintf("  %-14s %5d %s\n", lab[i], et[i], bar))
  }
  sel <- scr_selected(x)
  cat(sprintf("\nApproved: %d\n", length(sel)))
  if (length(sel)) {
    f <- x$funnel[approved == TRUE][order(consensus_rank)]
    n <- min(5L, nrow(f))
    for (i in seq_len(n)) cat(sprintf("  %2d. %-44s IV %6.3f  KS %.3f\n", f$consensus_rank[i], f$feature[i], f$total_iv[i], f$ks[i]))
    if (nrow(f) > n) cat(sprintf("  ... (+%d) - scr_selected() for the list\n", nrow(f) - n))
  }
  mt <- x$models$metrics
  if (nrow(mt)) {
    cat("\nModels (hold-out)\n")
    for (i in seq_len(nrow(mt))) cat(sprintf("  %-9s AUC %.4f [%.4f, %.4f]  KS %.4f\n", mt$model[i], mt$auc[i], mt$auc_lo[i], mt$auc_hi[i], mt$ks[i]))
  }
  warn <- character()
  if (isTRUE(x$consensus$meta$scarce)) {
    warn <- c(warn, sprintf("consensus returned %d, below the minimum of %d (eligible pool: %d)",
                            length(sel), x$config$target_min, x$consensus$meta$n_pool))
  }
  if (!identical(x$consensus$meta$relaxation, "none")) warn <- c(warn, sprintf("relaxation applied: %s", x$consensus$meta$relaxation))
  n_susp <- nrow(x$funnel[approved == TRUE & iv_suspect %in% TRUE])
  if (n_susp) warn <- c(warn, sprintf("%d approved with IV >= %.2f - run scr_leakage()", n_susp, x$config$iv_suspect))
  if (length(x$derived_excluded)) warn <- c(warn, sprintf("%d derived flag(s) outside the deliverable by policy (allow_derived_final)", length(x$derived_excluded)))
  if (length(warn)) { cat("\nWarnings\n"); for (w in warn) cat(sprintf("  - %s\n", w)) }
  if (!is.null(x$lab)) {
    cat(sprintf("\nCoarse classing: %d manual bin(s), %d forced in, %d dropped, %d decision(s) by %s - scr_decisions()\n",
                sum(x$lab$source == "manual"), length(setdiff(x$lab$shortlist$final, x$lab$shortlist$consensus)),
                length(setdiff(x$lab$shortlist$consensus, x$lab$shortlist$final)), nrow(x$lab$ledger), x$lab$author))
  }
  invisible(x)
}

#' @rdname scr_result
#' @export
summary.scr_result <- function(object, ...) {
  structure(list(target = object$target, text = object$summary_md, meta = object$meta),
            class = c("scr_summary", "list"))
}

#' @export
print.scr_summary <- function(x, ...) {
  cat(x$text, sep = "\n"); cat("\n")
  invisible(x)
}

#' @rdname scr_result
#' @export
as.data.frame.scr_result <- function(x, ...) {
  as.data.frame(scr_funnel(x))
}

#' @rdname scr_result
#' @export
plot.scr_result <- function(x, ...) {
  m <- x$meta
  v <- c(m$n_candidates, m$n_after_triage, m$n_binned, m$n_after_screening, m$n_after_holdout,
         m$n_after_prune, length(scr_selected(x)))
  nm <- c("input", "triage", "binning", "screening", "hold-out", "correlation", "approved")
  op <- graphics::par(mar = c(7, 4, 3, 1)); on.exit(graphics::par(op), add = TRUE)
  bp <- graphics::barplot(v, names.arg = nm, las = 2, main = sprintf("Selection funnel - %s", x$target),
                          ylab = "variables", border = NA,
                          col = c(rep("grey80", length(v) - 1L), "grey35"), ylim = c(0, max(v) * 1.12))
  graphics::text(bp, v, labels = v, pos = 3, cex = 0.85, xpd = TRUE)
  invisible(x)
}

# -- Run set ---------------------------------------------------------------- #

#' Set of runs, one per target
#'
#' Object returned by [scr_run()]: a named list of `scr_result`, plus the
#' errors of the targets that failed. Use [scr_compare()] for the comparison
#' table and [scr_core()] for the variables that cross several targets.
#'
#' @param x An `scr_runset` object.
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#'
#' @examplesIf requireNamespace("RSQLite", quietly = TRUE)
#' con <- scr_connect(driver = RSQLite::SQLite(), dbname = ":memory:")
#' d <- scr_demo[, c("default", "ds_region", "ds_band", "vl_score_01",
#'                   "vl_score_02", "vl_score_05", "vl_hist_01")]
#' DBI::dbWriteTable(con, "dtm", d)
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   use_lightgbm = FALSE, xgb_rounds = 40, n_boot = 10)
#' rs <- scr_run(con, "dtm", targets = "default", config = cfg)
#' rs
#' names(rs)
#' DBI::dbDisconnect(con)
#' @name scr_runset
#' @family portfolio
#' @export
print.scr_runset <- function(x, ...) {
  ok <- vapply(x, inherits, logical(1), "scr_result")
  cat(sprintf("<scr_runset> %d target(s): %d succeeded, %d failed\n\n", length(x), sum(ok), sum(!ok)))
  if (any(ok)) {
    cat(sprintf("  %-12s %8s %8s %8s %8s\n", "target", "rows", "approved", "AUC", "KS"))
    for (nm in names(x)[ok]) {
      r <- x[[nm]]; a <- r$models$metrics
      best <- if (nrow(a) && any(is.finite(a$auc))) a[which.max(a$auc)] else NULL
      cat(sprintf("  %-12s %8s %8d %8s %8s\n", nm, n_fmt(r$meta$n_total), length(scr_selected(r)),
                  if (!is.null(best)) sprintf("%.4f", best$auc) else "-",
                  if (!is.null(best)) sprintf("%.4f", best$ks) else "-"))
    }
  }
  if (any(!ok)) {
    cat("\n  Errors\n")
    for (nm in names(x)[!ok]) cat(sprintf("  %-12s %s\n", nm, x[[nm]]$error))
  }
  invisible(x)
}
