# ============================================================================ #
# funnel.R - audit funnel, gains per bin and executive summary
# ============================================================================ #
# The central deliverable is the FUNNEL: one row per input variable (plus the
# derived ones), with the verdict of every gate and the exact stage each one
# died at. No candidate disappears from the report. One row per run, no
# run_id.
# ============================================================================ #

#' Consolidate the fate of every input variable
#' @keywords internal
#' @noRd
build_funnel <- function(cols, triage, bins, models, cfg, selected = models$consensus$selected, lab = NULL) {
  fun <- data.table::copy(triage$profile)
  if (length(cols$dropped)) {
    fun <- data.table::rbindlist(list(fun, data.table::data.table(
      feature = cols$dropped, derived_from = NA_character_, type = "dropped",
      n_missing = NA_integer_, pct_missing = NA_real_, n_special = NA_integer_, pct_special = NA_real_,
      n_distinct = NA_integer_, pct_mode = NA_real_, iv_quick = NA_real_, woe_special = NA_real_,
      decomposition = "-", triage_status = "drop", triage_reason = "DROP_LIST")),
      use.names = TRUE, fill = TRUE)
  }

  # -- screening ----------------------------------------------------------- #
  sc <- data.table::copy(bins$screen$summary)
  keep_cols <- intersect(c("feature", "n_bins", "total_iv", "iv_class", "ks", "gini", "auc", "max_lift",
                           "min_bin_pct", "n_degenerate_bins", "monotonic", "monotonic_direction",
                           "quality", "selected", "reason", "reason_desc", "error", "error_msg"), names(sc))
  sc <- sc[, keep_cols, with = FALSE]
  map <- c(selected = "screen_selected", reason = "screen_reason", reason_desc = "screen_desc",
           error = "bin_error", error_msg = "bin_error_msg")
  found <- intersect(names(map), names(sc))
  if (length(found)) data.table::setnames(sc, found, unname(map[found]))
  fun <- merge(fun, sc, by = "feature", all.x = TRUE)
  fun[, iv_suspect := !is.na(total_iv) & total_iv >= cfg$iv_suspect]

  # -- hold-out ------------------------------------------------------------ #
  fun <- merge(fun, bins$holdout[, .(feature, iv_holdout, iv_ratio, psi, psi_flag, psi_critical,
                                     psi_flag_adjusted, pct_unbinned, holdout_ok, holdout_reason)],
               by = "feature", all.x = TRUE)

  # -- redundancy ---------------------------------------------------------- #
  fun[, `:=`(prune_status = NA_character_, prune_corr_with = NA_character_, prune_corr = NA_real_)]
  if (nrow(bins$prune$dropped)) {
    pd <- data.table::as.data.table(bins$prune$dropped)
    pd <- data.table::data.table(feature = pd$feature,
                                 corr_with = if ("correlated_with" %in% names(pd)) pd$correlated_with else NA_character_,
                                 corr = if ("correlation" %in% names(pd)) pd$correlation else NA_real_)
    fun[pd, on = "feature", `:=`(prune_status = "removed", prune_corr_with = i.corr_with, prune_corr = i.corr)]
  }
  fun[feature %in% bins$prune$keep, prune_status := "kept"]

  # -- model votes --------------------------------------------------------- #
  votes <- models$votes
  if (nrow(votes)) {
    wide <- data.table::dcast(votes, feature ~ model, value.var = "vote")
    data.table::setnames(wide, setdiff(names(wide), "feature"), paste0("vote_", setdiff(names(wide), "feature")))
    fun <- merge(fun, wide, by = "feature", all.x = TRUE)
    imp <- data.table::dcast(votes, feature ~ model, value.var = "score_pct")
    data.table::setnames(imp, setdiff(names(imp), "feature"), paste0("rankpct_", setdiff(names(imp), "feature")))
    fun <- merge(fun, imp, by = "feature", all.x = TRUE)
  }

  # -- consensus ----------------------------------------------------------- #
  ct <- models$consensus$table
  if (nrow(ct)) {
    fun <- merge(fun, ct[, .(feature, votes, consensus_score, consensus_rank,
                             consensus_selected = selected, consensus_reason = reason)],
                 by = "feature", all.x = TRUE)
  } else {
    fun[, `:=`(votes = NA_integer_, consensus_score = NA_real_, consensus_rank = NA_integer_,
               consensus_selected = NA, consensus_reason = NA_character_)]
  }

  # -- where each one died (sequential assignment: the first rule wins) --- #
  sel <- selected
  fun[, exit_stage := NA_character_]
  fun[is.na(exit_stage) & triage_reason == "DROP_LIST",                    exit_stage := "00.config"]
  fun[is.na(exit_stage) & triage_status == "drop",                         exit_stage := "01.triage"]
  fun[is.na(exit_stage) & (bin_error %in% TRUE | is.na(screen_selected)), exit_stage := "02.binning"]
  fun[is.na(exit_stage) & !(screen_selected %in% TRUE),                   exit_stage := "03.screening"]
  fun[is.na(exit_stage) & !(holdout_ok %in% TRUE),                        exit_stage := "04.holdout"]
  fun[is.na(exit_stage) & prune_status %in% "removed",                    exit_stage := "05.correlation"]
  fun[is.na(exit_stage) & feature %in% bins$derived_excluded,             exit_stage := "05b.derived_excluded"]
  fun[is.na(exit_stage) & !(consensus_selected %in% TRUE),                exit_stage := "06.consensus"]
  fun[feature %in% sel, exit_stage := "07.approved"]
  fun[is.na(exit_stage), exit_stage := "06.consensus"]
  fun[, approved := feature %in% sel]
  # -- coarse classing provenance (who decided), next to exit_stage (where) - #
  fun[, `:=`(provenance = "auto", manual_reason = NA_character_)]
  if (!is.null(lab)) {
    pv <- lab$provenance
    fun[feature %in% names(pv), provenance := unname(pv[feature])]
    rs <- lab$reasons
    if (length(rs)) fun[feature %in% names(rs), manual_reason := vapply(feature, function(f) rs[[f]] %||% NA_character_, character(1))]
    fun[feature %in% setdiff(lab$shortlist$consensus, sel), exit_stage := "08.manual_drop"]
  }

  front <- c("feature", "derived_from", "type", "approved", "exit_stage", "consensus_rank",
             "consensus_score", "votes", "total_iv", "iv_holdout", "ks")
  data.table::setcolorder(fun, c(intersect(front, names(fun)), setdiff(names(fun), front)))
  data.table::setorderv(fun, c("approved", "consensus_rank", "exit_stage", "total_iv"),
                        c(-1L, 1L, 1L, -1L), na.last = TRUE)
  fun[]
}

#' Gains table (bin level) of every binned feature
#' @keywords internal
#' @noRd
build_gains <- function(bins, models, cfg, selected = models$consensus$selected) {
  g <- data.table::copy(bins$screen$full)
  keep_cols <- intersect(c("feature", "type", "bin_id", "bin", "bin_lower", "bin_upper", "n_categories",
                           "categories", "count", "count_perc", "pos", "neg", "pos_rate", "woe", "woe_model",
                           "iv", "total_iv", "lift", "ks_bin", "cum_pos_perc", "cum_neg_perc", "precision",
                           "recall", "f1_score", "kl_divergence", "js_divergence", "iv_class", "ks", "gini",
                           "auc", "monotonic", "quality", "selected", "reason"), names(g))
  g <- g[, keep_cols, with = FALSE]
  if ("selected" %in% names(g)) data.table::setnames(g, "selected", "screen_selected")
  if ("reason" %in% names(g))   data.table::setnames(g, "reason", "screen_reason")
  g <- merge(g, bins$holdout[, .(feature, iv_holdout, psi, psi_flag)], by = "feature", all.x = TRUE)
  ct <- models$consensus$table
  if (nrow(ct)) g <- merge(g, ct[, .(feature, consensus_rank, votes, consensus_score)], by = "feature", all.x = TRUE)
  g[, approved := feature %in% selected]
  data.table::setcolorder(g, c("feature", "approved", "bin_id", "bin"))
  data.table::setorderv(g, c("approved", "total_iv", "feature", "bin_id"), c(-1L, -1L, 1L, 1L), na.last = TRUE)
  g[]
}

# -- Executive summary ------------------------------------------------------ #

#' @keywords internal
#' @noRd
build_summary <- function(meta, funnel, models, cfg, lab = NULL) {
  consensus <- models$consensus
  top <- funnel[approved == TRUE][order(consensus_rank)]
  vv <- vocab(cfg)
  lines <- c(
    sprintf("# scorecraft - %s", meta$target), "",
    sprintf("- Rows: %s (train %s / hold-out %s) - split %s%s", n_fmt(meta$n_total), n_fmt(meta$n_train),
            n_fmt(meta$n_holdout), meta$split_method,
            if (!is.na(meta$split_cutoff)) sprintf(" (cut: %s)", meta$split_cutoff) else ""),
    sprintf("- Event rate: train %s | hold-out %s", fmt_pct(meta$event_rate_train, 2), fmt_pct(meta$event_rate_holdout, 2)),
    sprintf("- Convention (objective = \"%s\"): %s. Score: %s.", cfg$objective, vv$target1, vv$points),
    if (isTRUE(meta$event$inverted))
      "- **Target inverted** by `event_level = 0`: class 0 of the table is the modelled event."
    else sprintf("- Class modelled as the event: `%s`", meta$event$label %||% "1"),
    sprintf("- Preset: %s | seed: %d | variables target: %d to %d | algorithm: %s",
            cfg$preset, cfg$seed, cfg$target_min, cfg$target_max, cfg$algorithm),
    "", "## Funnel", "",
    "| Stage | In | Survived |", "|---|---:|---:|",
    sprintf("| 0. Columns in the table | %d | %d candidates |", meta$n_cols, meta$n_candidates),
    sprintf("| 1. Descriptive triage | %d | %d |", meta$n_candidates, meta$n_after_triage),
    sprintf("| 2. Binning | %d | %d |", meta$n_after_triage, meta$n_binned),
    sprintf("| 3. Screening | %d | %d |", meta$n_binned, meta$n_after_screening),
    sprintf("| 4. Hold-out revalidation | %d | %d |", meta$n_after_screening, meta$n_after_holdout),
    sprintf("| 5. Redundancy | %d | %d |", meta$n_after_holdout, meta$n_after_prune),
    sprintf("| 6. Consensus (%d models) | %d | **%d approved** |", nrow(models$metrics), meta$n_after_prune,
            length(consensus$selected)),
    "", sprintf("Derived flags created at Stage 1: %d", length(meta$derived)),
    sprintf("Consensus relaxation: %s", consensus$meta$relaxation),
    "", "## Models (hold-out)", "",
    "| Model | AUC | 95% CI | KS | Gini | Votes | Weight | Note |", "|---|---:|---|---:|---:|---:|---:|---|")
  weights <- consensus$meta$weights
  for (i in seq_len(nrow(models$metrics))) {
    m <- models$metrics[i]
    w <- if (!is.null(weights)) weights[model == m$model, weight] else NA_real_
    lines <- c(lines, sprintf("| %s | %.4f | [%.4f, %.4f] | %.4f | %.4f | %d | %.3f | %s |",
                              m$model, m$auc, m$auc_lo, m$auc_hi, m$ks, m$gini, m$n_vote,
                              if (length(w)) w[1] else NA_real_, m$note))
  }
  lines <- c(lines, "", sprintf("## Approved variables (%d)", nrow(top)), "",
             "| # | Variable | Type | Bins | IV train | IV hold-out | KS | PSI | Votes | Score |",
             "|---:|---|---|---:|---:|---:|---:|---:|---:|---:|")
  for (i in seq_len(nrow(top))) {
    r <- top[i]
    lines <- c(lines, sprintf("| %d | %s%s | %s | %s | %.4f | %.4f | %.4f | %.3f | %d | %.3f |",
      r$consensus_rank, r$feature, if (!is.na(r$derived_from)) sprintf(" (from %s)", r$derived_from) else "",
      r$type, as.character(r$n_bins), r$total_iv, r$iv_holdout, r$ks, r$psi, as.integer(r$votes), r$consensus_score))
  }
  susp <- funnel[approved == TRUE & iv_suspect %in% TRUE][order(-total_iv)]
  if (nrow(susp)) {
    lines <- c(lines, "", sprintf("## [WARNING] Suspicious IV among the approved (%d of %d)", nrow(susp), nrow(top)), "",
      sprintf("IV >= %.2f is Siddiqi's \"suspicious\" band. Check that these variables are available at decision time and are not contemporaneous with the event.", cfg$iv_suspect),
      "", "| Variable | IV train | IV hold-out | KS | Bins |", "|---|---:|---:|---:|---:|")
    for (i in seq_len(nrow(susp))) {
      r <- susp[i]
      lines <- c(lines, sprintf("| %s | %.4f | %.4f | %.4f | %s |", r$feature, r$total_iv, r$iv_holdout, r$ks, as.character(r$n_bins)))
    }
  }
  der <- funnel[exit_stage == "05b.derived_excluded"][order(-iv_holdout)]
  if (nrow(der)) {
    lines <- c(lines, "", sprintf("## Derived flags outside the deliverable (%d)", nrow(der)), "",
      "They passed every gate but are columns the pipeline creates, not columns of the table (allow_derived_final = FALSE). The MISSING/sentinel state of the source column carries this signal.",
      "", "| Derived | Source | % sentinel | IV train | IV hold-out | KS |", "|---|---|---:|---:|---:|---:|")
    for (i in seq_len(nrow(der))) {
      r <- der[i]
      lines <- c(lines, sprintf("| %s | %s | %s | %.4f | %.4f | %.4f |", r$feature, r$derived_from,
                                fmt_pct(r$pct_special), r$total_iv, r$iv_holdout, r$ks))
    }
  }
  rep_f <- utils::head(funnel[approved == FALSE & !is.na(total_iv)][order(-total_iv)], 10L)
  if (nrow(rep_f)) {
    lines <- c(lines, "", "## Highest IVs failed (for manual review)", "", "| Variable | IV | Stage | Reason |", "|---|---:|---|---|")
    for (i in seq_len(nrow(rep_f))) {
      r <- rep_f[i]
      why <- paste(stats::na.omit(c(
        if (!identical(r$screen_reason, "OK")) r$screen_reason,
        if (!isTRUE(r$holdout_ok)) r$holdout_reason,
        if (identical(r$prune_status, "removed")) sprintf("REDUNDANT_WITH:%s(%.2f)", r$prune_corr_with, r$prune_corr),
        if (!is.na(r$consensus_reason) && r$consensus_reason != "OK") r$consensus_reason)), collapse = "; ")
      lines <- c(lines, sprintf("| %s | %.4f | %s | %s |", r$feature, r$total_iv, r$exit_stage, why))
    }
  }
  if (!is.null(lab)) {
    lg <- lab$ledger
    lines <- c(lines, "", "## Manual interventions (coarse classing)", "",
               sprintf("Final shortlist: %d variables (consensus %d; forced in %d; dropped %d). Author: %s.",
                       length(lab$shortlist$final), length(lab$shortlist$consensus),
                       length(setdiff(lab$shortlist$final, lab$shortlist$consensus)),
                       length(setdiff(lab$shortlist$consensus, lab$shortlist$final)), lab$author), "",
               "| # | Variable | Action | Bins | IV hold-out before -> after | Verdict | Reason |", "|---:|---|---|---|---|---|---|")
    for (i in seq_len(nrow(lg))) {
      lines <- c(lines, sprintf("| %d | %s | %s | %s | %s | %s | %s |", lg$seq[i], lg$variable[i], lg$action[i],
                                if (is.na(lg$n_bins_after[i])) "" else sprintf("%s -> %d", if (is.na(lg$n_bins_before[i])) "?" else lg$n_bins_before[i], lg$n_bins_after[i]),
                                if (is.na(lg$iv_holdout_after[i])) "" else sprintf("%.4f -> %.4f", lg$iv_holdout_before[i], lg$iv_holdout_after[i]),
                                if (is.na(lg$verdict[i])) "" else lg$verdict[i], lg$reason[i]))
    }
  }
  if (isTRUE(consensus$meta$scarce)) {
    lines <- c(lines, "", sprintf(
      "> **Warning**: the consensus returned %d variables, below the requested minimum (%d). The eligible pool had %d. No variable failed by a gate was included to fill the number.",
      length(consensus$selected), cfg$target_min, consensus$meta$n_pool))
  }
  paste(lines, collapse = "\n")
}
