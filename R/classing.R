# ============================================================================ #
# classing.R - coarse classing lab: manual binning and manual variable choice
# ============================================================================ #
# A manual bin is a hand-built entry with the same contract the binning
# engine reads (cutpoints for numerics, "a%;%b" labels for categoricals), so
# scr_apply(), scr_sql() and scr_scorecard() consume it through the very same
# code path as an optimal bin. Every decision carries a reason and lands in an
# append-only ledger; the automatic artefacts are frozen alongside for audit.
# Value semantics throughout: every verb returns the updated lab.
# ============================================================================ #

#' Coarse classing lab: manual binning and manual variable choice
#'
#' Opens a lab on an [scr_select()] result. Inside it the analyst inspects
#' the optimal bins of any binned variable ([scr_classing_view()]), proposes
#' new breaks or groupings ([scr_classing_propose()]), reads the comparison
#' against the optimal bins, accepts or discards each proposal with a
#' mandatory reason ([scr_classing_accept()], [scr_classing_discard()]),
#' chooses the final variable list ([scr_classing_choose()]) and commits
#' everything to a new `scr_result` ([scr_classing_apply()]) that the rest
#' of the pipeline consumes unchanged: [scr_scorecard()], [scr_apply()],
#' [scr_sql()], [scr_export()].
#'
#' @section Contract of a manual bin:
#'
#' A manual bin is recomputed **on the training rows only** (hold-out rows
#' can never define a bin), with the engine's own WOE formula
#' (`ln(%event / %non-event)`, event-oriented, so glm coefficients stay
#' positive), then revalidated on the hold-out with the bins frozen (IV, PSI
#' with both thresholds, unbinned share) and screened with the eight
#' engine rules, so the lab and the pipeline can never disagree. Numeric
#' intervals are right-closed, `(a, b]`, exactly as the engine and its SQL.
#' Re-declaring the optimal cut points of a numeric reproduces the engine's
#' WOE exactly; for a categorical the engine applies a small internal
#' smoothing of its own, so the raw log-ratio of the lab differs from it in
#' the third decimal.
#'
#' @section What is never allowed silently:
#'
#' An empty bin, a degenerate bin (no events or no non-events, unless
#' `laplace > 0`), a bin below `lab_min_bin_pct_hard`, a manual IV crossing
#' `iv_max` (the lab must not manufacture leakage), a category left
#' unassigned, a missing reason. Those block the proposal (`BLOCKED`);
#' accepting one needs `override = TRUE`, and the override is itself a
#' ledger row.
#'
#' @param x An object from [scr_select()].
#' @param features Variables the lab covers. Default: every variable that
#'   reached binning (`names(x$fit$results)`), so a variable failed by
#'   screening can be rebinned and forced in with a reason.
#' @param laplace Smoothing added to the bin counts when recomputing WOE.
#'   `0` (default) is exactly the engine's formula.
#' @param max_iv_loss Advisory threshold: a manual bin whose hold-out IV
#'   falls more than this fraction below the optimal one raises
#'   `IV_LOSS_VS_OPTIMAL`. `NULL` uses `config$lab_max_iv_loss`.
#' @param author Free text recorded in the ledger.
#'
#' @return An `scr_classing` object (the lab), with a print method that
#'   summarises the session: variables touched, before/after IV, verdicts,
#'   reasons, pending proposals and the final choice.
#'
#' @family classing
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   xgb_rounds = 60, n_boot = 20)
#' res <- scr_select(scr_demo, "default", config = cfg, drop = "id",
#'                   date_col = "ref_date")
#' lab <- scr_coarse_classing(res)
#' lab
#' scr_classing_view(lab, "ds_region")
#' p <- scr_classing_propose(lab, "ds_region",
#'                           groups = list(edge = c("NORTH", "SOUTH"),
#'                                         core = c("EAST", "WEST", "CENTRE")))
#' p
#' lab <- scr_classing_accept(lab, p, reason = "edge/core is what pricing uses")
#' lab <- scr_classing_choose(lab, drop = "vl_score_10",
#'                            reason = "not available at decision time")
#' lab
#' res2 <- scr_classing_apply(lab)
#' scr_selected(res2)
#' scr_selected(res2, which = "consensus")
#' sc <- scr_scorecard(res2)
#' sc$model_card$binning_algorithm
#' @export
scr_coarse_classing <- function(x, features = NULL, laplace = 0, max_iv_loss = NULL,
                                author = Sys.info()[["user"]]) {
  check_result(x, "scr_coarse_classing")
  .scr_num1(laplace, "laplace", lower = 0, upper = 5)
  binned <- names(x$fit$results)
  features <- features %||% binned
  unknown <- setdiff(features, binned)
  if (length(unknown)) {
    stop("scr_coarse_classing(): variable(s) never reached binning: ", lst(unknown),
         ". See scr_funnel(x)$exit_stage.", call. = FALSE)
  }
  lab <- structure(list(
    result = x, target = x$target, features = features, laplace = laplace,
    max_iv_loss = max_iv_loss %||% x$config$lab_max_iv_loss %||% 0.10,
    author = as.character(author %||% "unknown"), opened_at = Sys.time(),
    optimal = x$fit$results[features],
    current = x$fit$results[features],
    source = stats::setNames(rep("optimal", length(features)), features),
    accepted = list(), proposals = list(), n_proposals = 0L, ids = .lab_counter(),
    ledger = .empty_ledger(),
    choice = list(keep = NULL, drop = character(), force = character(), reasons = list()),
    checks_cache = .lab_cache()
  ), class = c("scr_classing", "list"))
  lab
}

#' A counter shared by every copy of a lab, so that proposal ids never collide
#' even when several proposals are made on the same lab value.
#' @keywords internal
#' @noRd
.lab_counter <- function() { e <- new.env(parent = emptyenv()); e$n <- 0L; e }

#' Memo of the checks of the optimal entries, shared by every copy of a lab
#' (they depend only on the frozen result, never on the lab's decisions)
#' @keywords internal
#' @noRd
.lab_cache <- function() new.env(parent = emptyenv())

#' @keywords internal
#' @noRd
.empty_ledger <- function() {
  data.table::data.table(
    seq = integer(), at = as.POSIXct(character()), author = character(), variable = character(),
    action = character(), proposal_id = character(), instruction = character(),
    n_bins_before = integer(), n_bins_after = integer(), iv_train_before = numeric(),
    iv_train_after = numeric(), iv_holdout_before = numeric(), iv_holdout_after = numeric(),
    psi_after = numeric(), verdict = character(), warnings = character(), reason = character())
}

#' @keywords internal
#' @noRd
.ledger_add <- function(lab, variable, action, reason = NA_character_, proposal_id = NA_character_,
                        instruction = NA_character_, before = NULL, after = NULL, verdict = NA_character_,
                        warnings = NA_character_) {
  g <- function(o, nm) if (is.null(o) || is.null(o[[nm]]) || !length(o[[nm]])) NA_real_ else as.numeric(o[[nm]][1])
  row <- data.table::data.table(
    seq = nrow(lab$ledger) + 1L, at = Sys.time(), author = lab$author, variable = variable,
    action = action, proposal_id = proposal_id, instruction = instruction,
    n_bins_before = as.integer(g(before, "n_bins")), n_bins_after = as.integer(g(after, "n_bins")),
    iv_train_before = g(before, "total_iv"), iv_train_after = g(after, "total_iv"),
    iv_holdout_before = g(before, "iv_holdout"), iv_holdout_after = g(after, "iv_holdout"),
    psi_after = g(after, "psi"), verdict = verdict, warnings = warnings, reason = reason)
  lab$ledger <- data.table::rbindlist(list(lab$ledger, row), use.names = TRUE, fill = TRUE)
  lab
}

#' @keywords internal
#' @noRd
check_lab <- function(lab, fn) {
  if (!inherits(lab, "scr_classing")) stop(sprintf("%s() expects a lab from scr_coarse_classing().", fn), call. = FALSE)
  invisible(TRUE)
}

#' @keywords internal
#' @noRd
.check_reason <- function(reason, fn) {
  if (is.null(reason) || !is.character(reason) || !length(reason) || any(is.na(reason)) ||
      any(nchar(trimws(reason)) < 5L)) {
    stop(fn, "(): `reason` is mandatory (at least 5 characters per entry) - it is the audit trail.", call. = FALSE)
  }
  trimws(reason)
}

# -- building a manual entry ------------------------------------------------ #

#' Engine-style numeric label
#' @keywords internal
#' @noRd
.num_label <- function(lo, hi) {
  f <- function(v) if (is.infinite(v)) (if (v < 0) "-Inf" else "+Inf") else sprintf("%.6f", v)
  sprintf("(%s;%s]", f(lo), f(hi))
}

#' WOE and IV per bin from counts (engine formula, optional Laplace)
#' @keywords internal
#' @noRd
.woe_from_counts <- function(pos, neg, laplace = 0) {
  k <- length(pos)
  dp <- (pos + laplace) / (sum(pos) + k * laplace)
  dn <- (neg + laplace) / (sum(neg) + k * laplace)
  woe <- log(dp / dn)
  list(woe = woe, iv = (dp - dn) * woe)
}

#' Hand-built numeric entry, on the training rows
#' @keywords internal
#' @noRd
.manual_entry_num <- function(f, xv, y, breaks, laplace = 0) {
  breaks <- as.double(breaks)
  if (!length(breaks)) stop("variable '", f, "': at least one break is needed (two bins).", call. = FALSE)
  if (any(!is.finite(breaks))) stop("variable '", f, "': breaks must be finite.", call. = FALSE)
  if (is.unsorted(breaks, strictly = TRUE)) stop("variable '", f, "': breaks must be strictly increasing.", call. = FALSE)
  edges <- c(-Inf, breaks, Inf)
  labels <- vapply(seq_len(length(edges) - 1L), function(i) .num_label(edges[i], edges[i + 1L]), character(1))
  if (anyDuplicated(labels)) stop("variable '", f, "': two breaks collapse to the same label at six decimals.", call. = FALSE)
  idx <- cut(xv, breaks = edges, right = TRUE, labels = FALSE)
  k <- length(labels)
  pos <- tabulate(idx[y == 1L], nbins = k); neg <- tabulate(idx[y == 0L], nbins = k)
  w <- .woe_from_counts(pos, neg, laplace)
  list(id = as.numeric(seq_len(k)), bin = labels, woe = w$woe, iv = w$iv, count = as.integer(pos + neg),
       count_pos = as.integer(pos), count_neg = as.integer(neg), cutpoints = breaks,
       converged = TRUE, iterations = 0L, feature = f, type = "numerical", algorithm = "manual",
       manual = list(kind = "numeric", breaks = breaks, laplace = laplace))
}

#' Hand-built categorical entry, on the training rows
#' @keywords internal
#' @noRd
.manual_entry_cat <- function(f, xv, y, groups, other_to = NULL, sep = "%;%", laplace = 0) {
  if (!is.list(groups) || length(groups) < 2L) stop("variable '", f, "': `groups` must be a list of at least two character vectors.", call. = FALSE)
  groups <- lapply(groups, function(g) as.character(g))
  all_cats <- unlist(groups, use.names = FALSE)
  dup <- unique(all_cats[duplicated(all_cats)])
  if (length(dup)) stop("variable '", f, "': category in two groups: ", lst(dup), call. = FALSE)
  lv <- unique(xv)
  unknown <- setdiff(all_cats, lv)
  if (length(unknown)) stop("variable '", f, "': category not seen on train: ", lst(unknown), call. = FALSE)
  if (any(grepl(sep, all_cats, fixed = TRUE))) stop("variable '", f, "': a category contains the bin separator '", sep, "'.", call. = FALSE)
  left <- setdiff(lv, all_cats)
  is_other <- rep(FALSE, length(groups))
  if (length(left)) {
    if (is.null(other_to)) {
      stop("variable '", f, "': training categories left unassigned: ", lst(left),
           ". List them in a group or name the catch-all bin with `other_to`.", call. = FALSE)
    }
    ot <- .resolve_bin_ref(other_to, names(groups), length(groups), f)
    groups[[ot]] <- c(groups[[ot]], left)
    is_other[ot] <- TRUE
  } else if (!is.null(other_to)) {
    is_other[.resolve_bin_ref(other_to, names(groups), length(groups), f)] <- TRUE
  }
  labels <- vapply(groups, function(g) paste(g, collapse = sep), character(1))
  if (anyDuplicated(labels)) stop("variable '", f, "': duplicated group.", call. = FALSE)
  map <- stats::setNames(rep(seq_along(groups), lengths(groups)), unlist(groups, use.names = FALSE))
  idx <- unname(map[xv])
  k <- length(groups)
  pos <- tabulate(idx[y == 1L], nbins = k); neg <- tabulate(idx[y == 0L], nbins = k)
  w <- .woe_from_counts(pos, neg, laplace)
  list(id = as.numeric(seq_len(k)), bin = unname(labels), woe = w$woe, iv = w$iv, count = as.integer(pos + neg),
       count_pos = as.integer(pos), count_neg = as.integer(neg), total_iv = sum(w$iv),
       converged = TRUE, iterations = 0L, feature = f, type = "categorical", algorithm = "manual",
       manual = list(kind = "categorical", groups = unname(groups), labels = names(groups) %||% rep(NA_character_, k),
                     is_other = is_other, laplace = laplace))
}

#' @keywords internal
#' @noRd
.resolve_bin_ref <- function(ref, nms, k, f) {
  if (is.character(ref)) {
    i <- match(ref, nms)
    if (is.na(i)) stop("variable '", f, "': no group named '", ref, "'.", call. = FALSE)
    return(i)
  }
  i <- as.integer(ref)
  if (is.na(i) || i < 1L || i > k) stop("variable '", f, "': bin reference ", ref, " is out of range (1..", k, ").", call. = FALSE)
  i
}

#' Groups currently encoded in a categorical entry
#' @keywords internal
#' @noRd
.entry_groups <- function(entry, sep) strsplit(entry$bin, sep, fixed = TRUE)

#' Build a one- or many-feature obwoe object around hand-built entries
#' @keywords internal
#' @noRd
.mini_fit <- function(x, entries) {
  fit <- x$fit
  fit$results <- entries
  fit$summary <- data.frame(
    feature = names(entries), type = vapply(entries, `[[`, character(1), "type"),
    algorithm = vapply(entries, `[[`, character(1), "algorithm"),
    n_bins = vapply(entries, function(e) length(e$bin), integer(1)),
    total_iv = vapply(entries, function(e) sum(e$iv), numeric(1)),
    converged = TRUE, iterations = 0L, error = FALSE, stringsAsFactors = FALSE)
  fit$n_features <- length(entries)
  fit
}

# -- checks and comparison -------------------------------------------------- #

#' Train/hold-out data of the lab, restricted to the columns `cols`
#'
#' Only the variables a check reads are copied: subsetting the rows of the
#' whole `data_clean` (thousands of columns) for every proposal was the
#' dominant cost of a lab session.
#' @keywords internal
#' @noRd
.lab_data <- function(lab, cols) {
  x <- lab$result
  dc <- x$data_clean
  tr <- x$split$train_idx; ho <- x$split$holdout_idx
  cols <- intersect(cols, names(dc))
  list(tr = dc[tr, cols, with = FALSE], ho = dc[ho, cols, with = FALSE],
       y_tr = dc[[x$target]][tr], y_ho = dc[[x$target]][ho])
}

#' Full quality assessment of one entry: engine screening, hold-out, bin tables, lab codes
#' @keywords internal
#' @noRd
.classing_checks <- function(lab, f, entry, optimal_checks = NULL) {
  x <- lab$result; cfg <- x$config; d <- .lab_data(lab, f)
  fit <- .mini_fit(x, stats::setNames(list(entry), f))
  # a degenerate bin has infinite WOE under laplace = 0; the engine's monotonicity
  # test cannot read it, so screening sees a clamped copy (the bin is blocked anyway)
  fit_s <- fit
  if (any(!is.finite(entry$woe))) {
    w <- entry$woe; w[!is.finite(w)] <- sign(w[!is.finite(w)]) * 20
    fit_s$results[[f]]$woe <- w
    fit_s$results[[f]]$iv[!is.finite(fit_s$results[[f]]$iv)] <- 0
  }
  sc <- screen_features(fit_s, cfg)$summary[1]
  w_tr <- apply_woe(fit, d$tr, f, "both"); w_ho <- apply_woe(fit, d$ho, f, "both")
  ho <- holdout_check(w_tr, w_ho, d$y_tr, d$y_ho, f, cfg)[1]
  b_tr <- w_tr[[paste0(f, "_bin")]]; b_ho <- w_ho[[paste0(f, "_bin")]]
  k <- length(entry$bin)
  n_tr <- tabulate(match(b_tr, entry$bin), nbins = k); n_ho <- tabulate(match(b_ho, entry$bin), nbins = k)
  e_tr <- tabulate(match(b_tr[d$y_tr == 1L], entry$bin), nbins = k); e_ho <- tabulate(match(b_ho[d$y_ho == 1L], entry$bin), nbins = k)
  w_h <- .woe_from_counts(e_ho, n_ho - e_ho, lab$laplace)
  bins <- data.table::data.table(
    variable = f, bin_id = seq_len(k), bin = entry$bin, n = n_tr, pct = n_tr / max(1L, sum(n_tr)), events = e_tr,
    event_rate = e_tr / pmax(1L, n_tr), woe = entry$woe, iv = entry$iv,
    n_holdout = n_ho, pct_holdout = n_ho / max(1L, sum(n_ho)), events_holdout = e_ho,
    event_rate_holdout = e_ho / pmax(1L, n_ho), woe_holdout = w_h$woe)
  bins[, pct_shift := pct_holdout - pct]

  codes <- character(); blocking <- character()
  if (any(n_tr == 0L)) blocking <- c(blocking, "EMPTY_BIN")
  degen <- any((e_tr == 0L | e_tr == n_tr) & n_tr > 0L)
  if (degen) { if (lab$laplace > 0 || isTRUE(cfg$allow_degenerate)) codes <- c(codes, "DEGENERATE_BIN") else blocking <- c(blocking, "DEGENERATE_BIN") }
  if (is.finite(sc$total_iv) && sc$total_iv >= cfg$iv_max) blocking <- c(blocking, "IV_SUSPICIOUS")
  hard <- cfg$lab_min_bin_pct_hard %||% 0.005
  if (any(n_tr > 0L & n_tr / sum(n_tr) < hard)) blocking <- c(blocking, "SMALL_BIN_HARD")
  screen_codes <- if (isTRUE(sc$selected) || identical(sc$reason, "OK")) character() else unlist(strsplit(as.character(sc$reason), ";", fixed = TRUE))
  screen_codes <- setdiff(screen_codes, c("IV_SUSPICIOUS", "DEGENERATE_BIN", "BINNING_ERROR"))
  codes <- c(codes, screen_codes)
  if (!identical(ho$holdout_reason, "OK")) codes <- c(codes, unlist(strsplit(ho$holdout_reason, ";", fixed = TRUE)))
  if (!is.null(optimal_checks) && is.finite(optimal_checks$iv_holdout) && is.finite(ho$iv_holdout) &&
      ho$iv_holdout < (1 - lab$max_iv_loss) * optimal_checks$iv_holdout) codes <- c(codes, "IV_LOSS_VS_OPTIMAL")
  if (is.finite(sc$total_iv) && sc$total_iv >= cfg$iv_suspect && sc$total_iv < cfg$iv_max) codes <- c(codes, "IV_SUSPECT")
  # a ratio over a train IV near zero is unstable: reported, never blocking
  if (is.finite(ho$iv_ratio) && is.finite(sc$total_iv) && sc$total_iv < cfg$iv_min) codes <- c(codes, "IV_RATIO_UNSTABLE")
  codes <- unique(codes)
  verdict <- if (length(blocking)) "BLOCKED" else if (length(codes)) "REVIEW" else "ACCEPTABLE"
  summary <- data.table::data.table(
    variable = f, type = entry$type, algorithm = entry$algorithm, n_bins = k,
    total_iv = sc$total_iv, iv_holdout = ho$iv_holdout, iv_ratio = ho$iv_ratio,
    ks = sc$ks %||% NA_real_, gini = sc$gini %||% NA_real_, psi = ho$psi, psi_flag = ho$psi_flag,
    psi_critical = ho$psi_critical, psi_flag_adjusted = ho$psi_flag_adjusted, pct_unbinned = ho$pct_unbinned,
    min_bin_pct = min(n_tr) / max(1L, sum(n_tr)), largest_bin_pct = max(n_tr) / max(1L, sum(n_tr)),
    n_degenerate_bins = sum((e_tr == 0L | e_tr == n_tr) & n_tr > 0L),
    monotonic = isTRUE(sc$monotonic), screen_selected = isTRUE(sc$selected),
    screen_reason = as.character(sc$reason), holdout_ok = isTRUE(ho$holdout_ok), holdout_reason = ho$holdout_reason,
    verdict = verdict, blocking = paste(blocking, collapse = ";"), warnings = paste(codes, collapse = ";"))
  list(summary = summary, bins = bins, screen = sc, holdout = ho, verdict = verdict,
       blocking = blocking, warnings = codes)
}

#' Checks of the optimal entry, cached in the lab
#' @keywords internal
#' @noRd
.optimal_checks <- function(lab, f) {
  memo <- lab$checks_cache
  if (is.environment(memo) && !is.null(memo[[f]])) return(memo[[f]])
  ck <- .classing_checks(lab, f, lab$optimal[[f]])
  if (is.environment(memo)) assign(f, ck, envir = memo)
  ck
}

#' Side-by-side comparison of two check summaries
#' @keywords internal
#' @noRd
.classing_compare <- function(opt, man, cur = NULL) {
  pick <- function(s) c(n_bins = s$n_bins, iv_train = s$total_iv, iv_holdout = s$iv_holdout, iv_ratio = s$iv_ratio,
                        ks = s$ks, psi = s$psi, min_bin_pct = s$min_bin_pct, largest_bin_pct = s$largest_bin_pct,
                        n_degenerate = s$n_degenerate_bins, monotonic = as.numeric(s$monotonic))
  o <- pick(opt$summary); m <- pick(man$summary)
  d <- data.table::data.table(metric = names(o), optimal = unname(o), manual = unname(m), delta = unname(m - o))
  if (!is.null(cur)) d[, current := unname(pick(cur$summary))]
  d[]
}

# -- view ------------------------------------------------------------------- #

#' Inspect the current bins of a variable in the lab
#'
#' Prints the current bins (optimal, or the accepted manual ones) with
#' train and hold-out side by side and a text bar chart of the event rate,
#' or, without `variable`, one line per variable of the lab.
#'
#' @param lab An object from [scr_coarse_classing()].
#' @param variable A variable name, or `NULL` for the overview.
#'
#' @return Invisibly, the bins table (`variable` given) or the overview table.
#'
#' @seealso [scr_coarse_classing()] for a complete session, from lab to
#'   scorecard.
#' @family classing
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   use_lightgbm = FALSE, xgb_rounds = 40, n_boot = 10)
#' d <- scr_demo[, c("default", "ref_date", "ds_region", "ds_band", "vl_score_01",
#'                   "vl_score_02", "vl_score_05", "vl_score_10", "vl_hist_01")]
#' res <- scr_select(d, "default", config = cfg, date_col = "ref_date")
#' lab <- scr_coarse_classing(res)
#' scr_classing_view(lab)
#' scr_classing_view(lab, "ds_region")
#' @export
scr_classing_view <- function(lab, variable = NULL) {
  check_lab(lab, "scr_classing_view")
  if (is.null(variable)) {
    ov <- data.table::rbindlist(lapply(lab$features, function(f) {
      ck <- if (identical(lab$source[[f]], "manual")) lab$accepted[[f]]$checks else .optimal_checks(lab, f)
      data.table::data.table(variable = f, source = lab$source[[f]], n_bins = ck$summary$n_bins,
                             iv_train = ck$summary$total_iv, iv_holdout = ck$summary$iv_holdout,
                             psi = ck$summary$psi, verdict = ck$verdict,
                             in_shortlist = f %in% .lab_final(lab))
    }))
    cat(sprintf("<scr_classing> target \"%s\" | %d variables\n", lab$target, nrow(ov)))
    cat(sprintf("  %-28s %-8s %5s %9s %9s %7s %-11s %s\n", "variable", "source", "bins", "IV train", "IV hold", "PSI", "verdict", "shortlist"))
    for (i in seq_len(nrow(ov))) cat(sprintf("  %-28s %-8s %5d %9.4f %9.4f %7.4f %-11s %s\n", ov$variable[i], ov$source[i], ov$n_bins[i],
                                             ov$iv_train[i], ov$iv_holdout[i], ov$psi[i], ov$verdict[i], if (ov$in_shortlist[i]) "yes" else "-"))
    return(invisible(ov))
  }
  if (!variable %in% lab$features) stop("'", variable, "' is not in the lab.", call. = FALSE)
  entry <- lab$current[[variable]]
  ck <- if (identical(lab$source[[variable]], "manual")) lab$accepted[[variable]]$checks else .optimal_checks(lab, variable)
  .print_bins(variable, entry, ck, lab$source[[variable]], lab$result$config$bin_separator)
  invisible(ck$bins)
}

#' @keywords internal
#' @noRd
.print_bins <- function(f, entry, ck, source, sep) {
  s <- ck$summary; b <- ck$bins
  cat(sprintf("<scr_classing> %s (%s) | current: %s | %d bins | train IV %.4f, hold-out IV %.4f (ratio %.2f)\n",
              f, entry$type, if (source == "manual") "manual" else paste0("optimal (", entry$algorithm, ")"),
              s$n_bins, s$total_iv, s$iv_holdout, s$iv_ratio))
  cat(sprintf("  monotone: %s | min bin %.1f%% | PSI %.4f (%s) | KS %.3f | degenerate bins: %d | verdict: %s%s\n",
              if (isTRUE(s$monotonic)) "yes" else "no", 100 * s$min_bin_pct, s$psi, s$psi_flag_adjusted, s$ks,
              s$n_degenerate_bins, s$verdict, if (nzchar(s$warnings)) paste0(" [", s$warnings, "]") else ""))
  cat(sprintf("  %3s  %-30s %7s %6s %7s %6s %8s %7s | %7s %6s %7s %8s\n", "id", "bin", "n", "%", "events", "rate", "WOE", "IV", "n.hold", "%", "rate", "WOE.hold"))
  for (i in seq_len(nrow(b))) {
    lab_i <- gsub(sep, " | ", b$bin[i], fixed = TRUE)
    cat(sprintf("  %3d  %-30s %7s %5.1f%% %7d %5.1f%% %8.3f %7.3f | %7s %5.1f%% %6.1f%% %8.3f\n", b$bin_id[i], substr(lab_i, 1, 30),
                n_fmt(b$n[i]), 100 * b$pct[i], b$events[i], 100 * b$event_rate[i], b$woe[i], b$iv[i],
                n_fmt(b$n_holdout[i]), 100 * b$pct_holdout[i], 100 * b$event_rate_holdout[i], b$woe_holdout[i]))
  }
  mx <- max(c(b$event_rate, b$event_rate_holdout), na.rm = TRUE)
  cat("  event rate by bin (train | hold-out)\n")
  for (i in seq_len(nrow(b))) {
    bar <- function(v) strrep("#", if (is.finite(mx) && mx > 0) round(18 * v / mx) else 0L)
    cat(sprintf("  %3d  %-18s %5.1f%% | %-18s %5.1f%%\n", b$bin_id[i], bar(b$event_rate[i]), 100 * b$event_rate[i],
                bar(b$event_rate_holdout[i]), 100 * b$event_rate_holdout[i]))
  }
  invisible(NULL)
}

# -- propose ---------------------------------------------------------------- #

#' Propose manual bins for a variable
#'
#' Exactly one of `breaks`, `groups`, `merge`, `split` or `reset` per call;
#' `missing_to` and `other_to` compose with a categorical instruction. The
#' instruction is resolved against the current bins into an absolute
#' specification, the WOE is refitted on the training rows, the bins are
#' applied frozen to the hold-out, and the comparison against the optimal
#' bins is printed. The proposal is a value: nothing changes in the lab
#' until [scr_classing_accept()].
#'
#' @param lab An object from [scr_coarse_classing()].
#' @param variable A variable of the lab.
#' @param breaks Numeric: interior cut points, `(-Inf, b1], (b1, b2], ...`.
#' @param groups Categorical: a list of character vectors, one per bin;
#'   names are display labels. Every training category must be assigned,
#'   or `other_to` must name the catch-all bin.
#' @param merge Bin ids to merge (adjacent for numerics).
#' @param split `c(id, at)`: split numeric bin `id` at `at`.
#' @param missing_to Categorical: fold the `"MISSING"` category into this bin.
#' @param other_to Categorical: the bin that receives every training
#'   category not listed in `groups`.
#' @param reset `TRUE` proposes a return to the optimal bins.
#'
#' @return An `scr_classing_proposal` with `id`, `variable`, `instruction`
#'   (the resolved instruction as text), `entry` (the hand-built bins),
#'   `checks`, `optimal` (the checks of the optimal bins), `compare`,
#'   `verdict` (`ACCEPTABLE`, `REVIEW` or `BLOCKED`), `warnings` and
#'   `blocking`.
#'
#' @seealso [scr_coarse_classing()] for a complete session, from lab to
#'   scorecard.
#' @family classing
#' @export
scr_classing_propose <- function(lab, variable, breaks = NULL, groups = NULL, merge = NULL, split = NULL,
                                 missing_to = NULL, other_to = NULL, reset = FALSE) {
  check_lab(lab, "scr_classing_propose")
  if (!variable %in% lab$features) stop("'", variable, "' is not in the lab.", call. = FALSE)
  given <- c(breaks = !is.null(breaks), groups = !is.null(groups), merge = !is.null(merge),
             split = !is.null(split), reset = isTRUE(reset))
  if (sum(given) != 1L && !(sum(given) == 0L && !is.null(missing_to))) {
    stop("scr_classing_propose(): give exactly one of `breaks`, `groups`, `merge`, `split`, `reset` (or `missing_to` alone).", call. = FALSE)
  }
  x <- lab$result; cfg <- x$config; sep <- cfg$bin_separator
  d <- .lab_data(lab, variable)
  cur <- lab$current[[variable]]
  xv <- d$tr[[variable]]; y <- d$y_tr
  is_num <- identical(cur$type, "numerical")

  if (isTRUE(reset)) {
    entry <- lab$optimal[[variable]]
    instr <- "reset to optimal"
  } else if (is_num) {
    if (!is.null(groups) || !is.null(missing_to) || !is.null(other_to)) stop("'", variable, "' is numeric: use `breaks`, `merge` or `split`.", call. = FALSE)
    cp <- cur$cutpoints
    if (!is.null(breaks)) {
      new <- as.double(breaks)
      if (anyDuplicated(new)) stop("`breaks` for '", variable, "' contains duplicates.", call. = FALSE)
      new <- sort(new); instr <- sprintf("breaks = c(%s)", paste(format(new, trim = TRUE), collapse = ", "))
    }
    if (!is.null(merge)) {
      ids <- sort(as.integer(merge))
      if (length(ids) < 2L || any(diff(ids) != 1L) || min(ids) < 1L || max(ids) > length(cp) + 1L)
        stop("`merge` must name at least two adjacent bin ids of '", variable, "'.", call. = FALSE)
      new <- cp[-(ids[1]:(ids[length(ids)] - 1L))]; instr <- sprintf("merge = c(%s)", paste(ids, collapse = ", "))
    }
    if (!is.null(split)) {
      if (length(split) != 2L) stop("`split` must be c(id, at).", call. = FALSE)
      id <- as.integer(split[1]); at <- as.double(split[2]); edges <- c(-Inf, cp, Inf)
      if (id < 1L || id > length(cp) + 1L || !(at > edges[id] && at < edges[id + 1L]))
        stop("`split`: ", at, " is not strictly inside bin ", id, " of '", variable, "'.", call. = FALSE)
      new <- sort(c(cp, at)); instr <- sprintf("split = c(%d, %s)", id, format(at))
    }
    if (!length(new)) stop("the instruction leaves a single bin for '", variable, "'.", call. = FALSE)
    entry <- .manual_entry_num(variable, xv, y, new, lab$laplace)
  } else {
    if (!is.null(breaks) || !is.null(split)) stop("'", variable, "' is categorical: use `groups`, `merge` or `missing_to`.", call. = FALSE)
    g <- .entry_groups(cur, sep); labels <- rep(NA_character_, length(g))
    ot <- other_to
    if (!is.null(groups)) { g <- lapply(groups, as.character); labels <- names(groups) %||% rep(NA_character_, length(g)); instr <- .fmt_groups(groups) }
    else if (!is.null(merge)) {
      ids <- sort(unique(as.integer(merge)))
      if (length(ids) < 2L || min(ids) < 1L || max(ids) > length(g)) stop("`merge` must name at least two bin ids of '", variable, "'.", call. = FALSE)
      merged <- unlist(g[ids]); g <- c(list(merged), g[-ids]); instr <- sprintf("merge = c(%s)", paste(ids, collapse = ", "))
      labels <- rep(NA_character_, length(g))
    } else instr <- ""
    if (!is.null(missing_to)) {
      g <- lapply(g, function(v) setdiff(v, "MISSING"))
      i <- .resolve_bin_ref(missing_to, labels, length(g), variable)
      g[[i]] <- c(g[[i]], "MISSING"); g <- g[lengths(g) > 0L]
      instr <- trimws(paste(instr, sprintf("missing_to = %s", format(missing_to))))
    }
    if (!is.null(ot)) instr <- trimws(paste(instr, sprintf("other_to = %s", format(ot))))
    names(g) <- if (all(is.na(labels))) NULL else labels
    entry <- .manual_entry_cat(variable, xv, y, g, other_to = ot, sep = sep, laplace = lab$laplace)
  }
  opt <- .optimal_checks(lab, variable)
  ck <- .classing_checks(lab, variable, entry, optimal_checks = opt$summary)
  curck <- if (identical(lab$source[[variable]], "manual")) lab$accepted[[variable]]$checks else NULL
  lab$ids$n <- lab$ids$n + 1L
  id <- sprintf("P%03d", lab$ids$n)
  structure(list(id = id, variable = variable, instruction = instr, entry = entry, checks = ck,
                 optimal = opt, compare = .classing_compare(opt, ck, curck), verdict = ck$verdict,
                 warnings = ck$warnings, blocking = ck$blocking, target = lab$target,
                 lab_opened_at = lab$opened_at, at = Sys.time()),
            class = c("scr_classing_proposal", "list"))
}

#' @keywords internal
#' @noRd
.fmt_groups <- function(groups) {
  nm <- names(groups) %||% rep("", length(groups))
  paste0("groups = list(", paste(vapply(seq_along(groups), function(i) {
    body <- paste0("c(", paste0("\"", groups[[i]], "\"", collapse = ", "), ")")
    if (nzchar(nm[i])) paste0(nm[i], " = ", body) else body
  }, character(1)), collapse = ", "), ")")
}

#' @export
print.scr_classing_proposal <- function(x, ...) {
  cat(sprintf("<scr_classing_proposal> %s %s | %s | %s\n", x$id, x$variable, x$instruction, format(x$at, "%Y-%m-%d %H:%M")))
  cp <- x$compare
  hdr <- if ("current" %in% names(cp)) sprintf("  %-18s %10s %10s %10s %10s\n", "", "optimal", "current", "manual", "delta")
         else sprintf("  %-18s %10s %10s %10s\n", "", "optimal", "manual", "delta")
  cat(hdr)
  for (i in seq_len(nrow(cp))) {
    fmt <- function(v) if (cp$metric[i] %in% c("n_bins", "n_degenerate", "monotonic")) sprintf("%10.0f", v) else sprintf("%10.4f", v)
    if ("current" %in% names(cp)) cat(sprintf("  %-18s %s %s %s %s\n", cp$metric[i], fmt(cp$optimal[i]), fmt(cp$current[i]), fmt(cp$manual[i]), fmt(cp$delta[i])))
    else cat(sprintf("  %-18s %s %s %s\n", cp$metric[i], fmt(cp$optimal[i]), fmt(cp$manual[i]), fmt(cp$delta[i])))
  }
  b <- x$checks$bins
  cat("  manual bins (train | hold-out)\n")
  for (i in seq_len(nrow(b))) cat(sprintf("  %3d  %-30s %7s %5.1f%% %5.1f%% %7.3f | %7s %5.1f%% %5.1f%% %7.3f\n", b$bin_id[i], substr(gsub("%;%", " | ", b$bin[i], fixed = TRUE), 1, 30),
                                          n_fmt(b$n[i]), 100 * b$pct[i], 100 * b$event_rate[i], b$woe[i], n_fmt(b$n_holdout[i]), 100 * b$pct_holdout[i], 100 * b$event_rate_holdout[i], b$woe_holdout[i]))
  if (length(x$warnings)) { cat("  Warnings\n"); for (w in x$warnings) cat(sprintf("    - %s\n", w)) }
  if (length(x$blocking)) { cat("  Blocking\n"); for (w in x$blocking) cat(sprintf("    - %s\n", w)) }
  cat(sprintf("  Verdict: %s%s\n", x$verdict, switch(x$verdict,
    ACCEPTABLE = " - no warning raised.", REVIEW = " - advisory warnings only; accept with a reason or discard.",
    BLOCKED = " - accept only with override = TRUE (recorded in the ledger).")))
  invisible(x)
}

# -- accept / discard ------------------------------------------------------- #

#' Accept or discard a proposal
#'
#' `reason` is mandatory. Accepting replaces the variable's current bins;
#' the previous accepted proposal is marked `superseded`. A `BLOCKED`
#' proposal needs `override = TRUE`, and the override is itself a ledger row.
#'
#' @param lab An object from [scr_coarse_classing()].
#' @param proposal An object from [scr_classing_propose()].
#' @param reason Free text, at least 5 characters.
#' @param override Accept a `BLOCKED` proposal.
#'
#' @return The updated lab, invisibly.
#'
#' @seealso [scr_coarse_classing()] for a complete session, from lab to
#'   scorecard.
#' @family classing
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   use_lightgbm = FALSE, xgb_rounds = 40, n_boot = 10)
#' d <- scr_demo[, c("default", "ref_date", "ds_region", "ds_band", "vl_score_01",
#'                   "vl_score_02", "vl_score_05", "vl_score_10", "vl_hist_01")]
#' res <- scr_select(d, "default", config = cfg, date_col = "ref_date")
#' lab <- scr_coarse_classing(res)
#' p <- scr_classing_propose(lab, "ds_region",
#'                           groups = list(edge = c("NORTH", "SOUTH"),
#'                                         core = c("EAST", "WEST", "CENTRE")))
#' lab <- scr_classing_accept(lab, p, reason = "edge/core is what pricing uses")
#' scr_classing_view(lab, "ds_region")
#' # a second proposal on the same variable, rejected with its reason
#' p2 <- scr_classing_propose(lab, "ds_region",
#'                            groups = list(c("NORTH", "SOUTH", "EAST"), c("WEST", "CENTRE")))
#' lab <- scr_classing_discard(lab, p2, reason = "no business rationale for this grouping")
#' scr_decisions(lab)
#' @export
scr_classing_accept <- function(lab, proposal, reason, override = FALSE) {
  .scr_local_verbose(lab)
  check_lab(lab, "scr_classing_accept")
  .check_proposal(lab, proposal, "scr_classing_accept")
  reason <- .check_reason(reason, "scr_classing_accept")
  if (identical(proposal$verdict, "BLOCKED") && !isTRUE(override)) {
    stop("scr_classing_accept(): proposal ", proposal$id, " is BLOCKED (", paste(proposal$blocking, collapse = ", "),
         "). Pass override = TRUE to accept it anyway; the override is recorded.", call. = FALSE)
  }
  f <- proposal$variable
  before <- if (identical(lab$source[[f]], "manual")) lab$accepted[[f]]$checks$summary else .optimal_checks(lab, f)$summary
  if (!is.null(lab$accepted[[f]])) {
    lab <- .ledger_add(lab, f, "supersede", proposal_id = lab$accepted[[f]]$proposal_id, reason = "superseded by a new accepted proposal")
  }
  is_reset <- identical(proposal$entry$algorithm, lab$optimal[[f]]$algorithm) && identical(proposal$entry$bin, lab$optimal[[f]]$bin)
  lab$current[[f]] <- proposal$entry
  lab$source[[f]] <- if (is_reset) "optimal" else "manual"
  lab$accepted[[f]] <- if (is_reset) NULL else list(proposal_id = proposal$id, instruction = proposal$instruction,
                                                   entry = proposal$entry, checks = proposal$checks, reason = reason,
                                                   at = Sys.time(), override = isTRUE(override), verdict = proposal$verdict)
  lab$n_proposals <- lab$ids$n
  lab$proposals[[proposal$id]] <- "accepted"
  if (identical(proposal$verdict, "BLOCKED")) {
    lab <- .ledger_add(lab, f, "override", proposal_id = proposal$id, reason = reason,
                       warnings = paste(proposal$blocking, collapse = ";"), verdict = proposal$verdict)
  }
  lab <- .ledger_add(lab, f, if (is_reset) "restore" else "accept", proposal_id = proposal$id, instruction = proposal$instruction,
                     before = before, after = proposal$checks$summary, verdict = proposal$verdict,
                     warnings = paste(proposal$warnings, collapse = ";"), reason = reason)
  msg("  %s: %s accepted (%s) - %d bins, hold-out IV %.4f", f, proposal$id, proposal$verdict,
      proposal$checks$summary$n_bins, proposal$checks$summary$iv_holdout)
  invisible(lab)
}

#' @rdname scr_classing_accept
#' @export
scr_classing_discard <- function(lab, proposal, reason) {
  check_lab(lab, "scr_classing_discard")
  .check_proposal(lab, proposal, "scr_classing_discard")
  reason <- .check_reason(reason, "scr_classing_discard")
  lab$n_proposals <- lab$ids$n
  lab$proposals[[proposal$id]] <- "discarded"
  lab <- .ledger_add(lab, proposal$variable, "discard", proposal_id = proposal$id, instruction = proposal$instruction,
                     after = proposal$checks$summary, verdict = proposal$verdict,
                     warnings = paste(proposal$warnings, collapse = ";"), reason = reason)
  invisible(lab)
}

#' @keywords internal
#' @noRd
.check_proposal <- function(lab, p, fn) {
  if (!inherits(p, "scr_classing_proposal")) stop(fn, "(): `proposal` must come from scr_classing_propose().", call. = FALSE)
  if (!identical(p$target, lab$target) || !identical(p$lab_opened_at, lab$opened_at)) {
    stop(fn, "(): the proposal was made on a different lab.", call. = FALSE)
  }
  if (!is.null(lab$proposals[[p$id]])) stop(fn, "(): proposal ", p$id, " was already ", lab$proposals[[p$id]], ".", call. = FALSE)
  invisible(TRUE)
}

# -- variable choice -------------------------------------------------------- #

#' Choose the final variable list manually
#'
#' The final list is `(consensus shortlist + force) - drop`, then
#' intersected with `keep` when given. `force` is allowed only for
#' variables that reached binning; a variable failed for `IV_SUSPICIOUS`
#' (the leakage ceiling) and a derived `__sp` flag under
#' `allow_derived_final = FALSE` are refused unless `override = TRUE`.
#' `reason` is one string for every variable named, or a character vector
#' named by variable.
#'
#' @param lab An object from [scr_coarse_classing()].
#' @param keep Variables to keep (restricts the final list).
#' @param drop Variables to remove from the final list.
#' @param force Variables to add to the final list.
#' @param reason Mandatory when `drop` or `force` is given.
#' @param override Allow a refused `force`.
#'
#' @return The updated lab, invisibly.
#'
#' @seealso [scr_coarse_classing()] for a complete session, from lab to
#'   scorecard.
#' @family classing
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   use_lightgbm = FALSE, xgb_rounds = 40, n_boot = 10)
#' d <- scr_demo[, c("default", "ref_date", "ds_region", "ds_band", "vl_score_01",
#'                   "vl_score_02", "vl_score_05", "vl_score_10", "vl_hist_01")]
#' res <- scr_select(d, "default", config = cfg, date_col = "ref_date")
#' lab <- scr_coarse_classing(res)
#' lab <- scr_classing_choose(lab, drop = "vl_score_10",
#'                            reason = "not available at decision time")
#' lab
#' @export
scr_classing_choose <- function(lab, keep = NULL, drop = NULL, force = NULL, reason = NULL, override = FALSE) {
  check_lab(lab, "scr_classing_choose")
  x <- lab$result; cfg <- x$config
  named <- c(drop, force)
  if (length(named)) {
    reason <- .check_reason(reason, "scr_classing_choose")
    rs <- if (length(reason) == 1L) stats::setNames(rep(reason, length(named)), named) else reason
    miss <- setdiff(named, names(rs))
    if (length(miss)) stop("scr_classing_choose(): no reason for ", lst(miss), ".", call. = FALSE)
  }
  binned <- names(x$fit$results)
  bad <- setdiff(c(keep, drop, force), binned)
  if (length(bad)) stop("scr_classing_choose(): not a binned variable: ", lst(bad), ".", call. = FALSE)
  for (f in force) {
    fr <- x$funnel[feature == f]
    if (nrow(fr) && grepl("IV_SUSPICIOUS", fr$screen_reason %||% "") && !isTRUE(override))
      stop("scr_classing_choose(): '", f, "' failed for IV_SUSPICIOUS (leakage ceiling); raise `iv_max` or pass override = TRUE.", call. = FALSE)
    if (f %in% x$triage$derived && !isTRUE(cfg$allow_derived_final) && !isTRUE(override))
      stop("scr_classing_choose(): '", f, "' is a derived flag and allow_derived_final = FALSE; pass override = TRUE.", call. = FALSE)
    lab <- .ledger_add(lab, f, "force", reason = rs[[f]], verdict = if (isTRUE(override)) "OVERRIDE" else NA_character_)
  }
  for (f in drop) lab <- .ledger_add(lab, f, "drop", reason = rs[[f]])
  if (!is.null(keep)) {
    lab$choice$keep <- keep
    for (f in keep) lab <- .ledger_add(lab, f, "keep", reason = if (!is.null(reason)) reason[1] else "kept explicitly")
  }
  lab$choice$drop  <- union(lab$choice$drop, drop)
  lab$choice$force <- union(setdiff(lab$choice$force, drop), force)
  lab$choice$drop  <- setdiff(lab$choice$drop, force)
  for (f in named) lab$choice$reasons[[f]] <- rs[[f]]
  if (!length(.lab_final(lab))) stop("scr_classing_choose(): the final list would be empty.", call. = FALSE)
  invisible(lab)
}

#' Final shortlist implied by the lab
#' @keywords internal
#' @noRd
.lab_final <- function(lab) {
  fin <- setdiff(union(lab$result$consensus$selected, lab$choice$force), lab$choice$drop)
  if (!is.null(lab$choice$keep)) fin <- intersect(fin, lab$choice$keep)
  fin
}

# -- spec: long format, read and import ------------------------------------- #

#' Classing specification as a long table, with its file round trip
#'
#' One row per bin of every variable in the lab, optimal and manual, with
#' the authoritative columns a reviewer may edit (`lower`/`upper` for
#' numerics, `categories`/`is_other` for categoricals, `reason`) and
#' context columns that are regenerated on read. Open ends are written as
#' `NA`. [scr_classing_read()] validates a file back into a spec and
#' [scr_classing_import()] turns every variable whose bins differ from the
#' lab's current ones into a proposal, so a spreadsheet edit never enters
#' silently.
#'
#' @param lab An object from [scr_coarse_classing()], or an `scr_result`
#'   returned by [scr_classing_apply()].
#' @param file For `scr_classing_spec()`, an optional `.csv` or `.xlsx` path
#'   to write the table to. For `scr_classing_read()`, the path to read.
#'   For `scr_classing_import()`, a path or an `scr_classing_spec` object.
#'
#' @return A `data.frame` of class `scr_classing_spec`.
#'
#' @seealso [scr_coarse_classing()] for a complete session, from lab to
#'   scorecard.
#' @family classing
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   use_lightgbm = FALSE, xgb_rounds = 40, n_boot = 10)
#' d <- scr_demo[, c("default", "ref_date", "ds_region", "ds_band", "vl_score_01",
#'                   "vl_score_02", "vl_score_05", "vl_score_10", "vl_hist_01")]
#' res <- scr_select(d, "default", config = cfg, date_col = "ref_date")
#' lab <- scr_coarse_classing(res)
#' p <- scr_classing_propose(lab, "ds_region",
#'                           groups = list(edge = c("NORTH", "SOUTH"),
#'                                         core = c("EAST", "WEST", "CENTRE")))
#' lab <- scr_classing_accept(lab, p, reason = "edge/core is what pricing uses")
#' sp <- scr_classing_spec(lab)
#' sp
#' # round trip through a file: a fresh lab receives the manual bins as proposals
#' f <- tempfile(fileext = ".csv")
#' scr_classing_spec(lab, file = f)
#' props <- scr_classing_import(scr_coarse_classing(res), scr_classing_read(f))
#' names(props)
#' unlink(f)
#' @export
scr_classing_spec <- function(lab, file = NULL) {
  .scr_local_verbose(lab)
  if (inherits(lab, "scr_result")) {
    if (is.null(lab$lab)) stop("this scr_result carries no classing.", call. = FALSE)
    sp <- lab$lab$spec
  } else {
    check_lab(lab, "scr_classing_spec")
    sep <- lab$result$config$bin_separator
    sp <- data.table::rbindlist(lapply(lab$features, function(f) {
      e <- lab$current[[f]]; man <- lab$accepted[[f]]
      ck <- if (!is.null(man)) man$checks else .optimal_checks(lab, f)
      b <- ck$bins; k <- nrow(b)
      num <- identical(e$type, "numerical")
      edges <- if (num) c(-Inf, e$cutpoints, Inf) else NULL
      data.table::data.table(
        target = lab$target, variable = f, type = if (num) "numeric" else "categorical", bin_id = seq_len(k),
        bin_label = e$bin,
        lower = if (num) ifelse(is.infinite(edges[-length(edges)]), NA_real_, edges[-length(edges)]) else NA_real_,
        upper = if (num) ifelse(is.infinite(edges[-1]), NA_real_, edges[-1]) else NA_real_,
        categories = if (num) NA_character_ else e$bin,
        is_other = if (num || is.null(e$manual)) FALSE else e$manual$is_other,
        source = if (is.null(man)) "optimal" else "manual",
        proposal_id = man$proposal_id %||% NA_character_, author = if (is.null(man)) NA_character_ else lab$author,
        decided_at = if (is.null(man)) NA_character_ else format(man$at, "%Y-%m-%d %H:%M:%S"),
        reason = man$reason %||% NA_character_,
        n_train = b$n, pct_train = b$pct, events_train = b$events, event_rate_train = b$event_rate,
        woe_train = b$woe, iv_train = b$iv, n_holdout = b$n_holdout, pct_holdout = b$pct_holdout,
        event_rate_holdout = b$event_rate_holdout, woe_holdout = b$woe_holdout)
    }))
    sp <- as.data.frame(sp, stringsAsFactors = FALSE)
    class(sp) <- c("scr_classing_spec", "data.frame")
  }
  if (!is.null(file)) {
    if (grepl("\\.xlsx$", file, ignore.case = TRUE)) {
      .need_openxlsx(); .scr_write_xlsx(list(Coarse_Classing = sp), file)
    } else utils::write.csv(sp, file, row.names = FALSE, na = "")
    msg("classing spec written to %s", file)
    return(invisible(sp))
  }
  sp
}

#' @export
print.scr_classing_spec <- function(x, ...) {
  d <- as.data.frame(x)
  cat(sprintf("<scr_classing_spec> %d bins | %d variables (%d manual)\n", nrow(d), length(unique(d$variable)),
              length(unique(d$variable[d$source == "manual"]))))
  print(utils::head(d[, intersect(c("variable", "type", "bin_id", "bin_label", "lower", "upper", "categories", "is_other", "source", "reason"), names(d))], 12), row.names = FALSE)
  if (nrow(d) > 12) cat(sprintf("  ... (+%d rows)\n", nrow(d) - 12))
  invisible(x)
}

#' @rdname scr_classing_spec
#' @param sep Bin separator of the `categories` column (the configuration's
#'   `bin_separator`). It is validated (no empty category, no category in
#'   two bins) and recorded on the spec, so that [scr_classing_import()]
#'   refuses a spec read with a different separator.
#' @export
scr_classing_read <- function(file, sep = "%;%") {
  if (!is.character(sep) || length(sep) != 1L || is.na(sep) || !nzchar(sep)) {
    stop("scr_classing_read(): `sep` must be a single non-empty string.", call. = FALSE)
  }
  if (!is.character(file) || length(file) != 1L || !file.exists(file)) {
    stop("scr_classing_read(): `file` must be the path of an existing .csv or .xlsx file.", call. = FALSE)
  }
  # the text columns are read as text, and only an empty cell is missing: a
  # category "NA" or "01" must come back as written, not as NA or 1
  txt <- c("variable", "type", "bin_label", "categories", "reason")
  d <- if (grepl("\\.xlsx$", file, ignore.case = TRUE)) {
    .need_openxlsx(); openxlsx::read.xlsx(file, sheet = 1, na.strings = "")
  } else {
    hdr <- names(utils::read.csv(file, nrows = 1L, check.names = FALSE))
    cc <- stats::setNames(rep("character", length(intersect(txt, hdr))), intersect(txt, hdr))
    utils::read.csv(file, stringsAsFactors = FALSE, na.strings = "", colClasses = cc, check.names = FALSE)
  }
  for (cn in intersect(c("bin_label", "categories"), names(d))) {
    v <- as.character(d[[cn]])
    # undo the formula-injection guard of the xlsx writer ("'-1" was "-1")
    d[[cn]] <- sub("^'(?=[=+@-])", "", v, perl = TRUE)
  }
  need <- c("variable", "type", "bin_id")
  miss <- setdiff(need, names(d))
  if (length(miss)) stop("scr_classing_read(): missing column(s): ", lst(miss), call. = FALSE)
  if (!"lower" %in% names(d)) d$lower <- NA_real_
  if (!"upper" %in% names(d)) d$upper <- NA_real_
  if (!"categories" %in% names(d)) d$categories <- NA_character_
  if (!"is_other" %in% names(d)) d$is_other <- FALSE
  if (!"reason" %in% names(d)) d$reason <- NA_character_
  to_num <- function(v) { v <- as.character(v); v[v %in% c("-Inf", "Inf", "+Inf")] <- NA; suppressWarnings(as.numeric(v)) }
  d$lower <- to_num(d$lower); d$upper <- to_num(d$upper)
  d$is_other <- as.logical(d$is_other); d$is_other[is.na(d$is_other)] <- FALSE
  d$bin_id <- as.integer(d$bin_id)
  errs <- character()
  for (v in unique(d$variable)) {
    r <- d[d$variable == v, , drop = FALSE]; r <- r[order(r$bin_id), , drop = FALSE]
    if (!identical(r$bin_id, seq_len(nrow(r)))) errs <- c(errs, sprintf("%s: bin_id must be 1..k without gaps", v))
    if (nrow(r) < 2L) errs <- c(errs, sprintf("%s: fewer than two bins", v))
    if (identical(r$type[1], "numeric")) {
      up <- r$upper[-nrow(r)]
      if (anyNA(up) || is.unsorted(up, strictly = TRUE)) errs <- c(errs, sprintf("%s: `upper` must be finite and strictly increasing except on the last bin", v))
      lo <- r$lower[-1]
      if (!isTRUE(all.equal(lo, up))) errs <- c(errs, sprintf("%s: `lower` of bin i must equal `upper` of bin i-1 (contiguity)", v))
      # the outer bins are open: a bound written there would be silently ignored
      if (!is.na(r$lower[1]) || !is.na(r$upper[nrow(r)])) {
        errs <- c(errs, sprintf("%s: the first `lower` and the last `upper` must be empty (open ends, -Inf/+Inf)", v))
      }
    } else {
      if (anyNA(r$categories) || any(!nzchar(r$categories))) errs <- c(errs, sprintf("%s: every categorical bin needs `categories`", v))
      else {
        cs <- as.character(r$categories)
        if (any(startsWith(cs, sep) | endsWith(cs, sep) | grepl(paste0(sep, sep), cs, fixed = TRUE))) {
          errs <- c(errs, sprintf("%s: empty category around the separator '%s'", v, sep))
        }
        cats <- unlist(strsplit(cs, sep, fixed = TRUE))
        dup <- unique(cats[duplicated(cats)])
        if (length(dup)) errs <- c(errs, sprintf("%s: category in two bins (%s)", v, lst(dup, 5)))
      }
      if (sum(r$is_other) > 1L) errs <- c(errs, sprintf("%s: only one bin may be `is_other`", v))
    }
  }
  if (length(errs)) stop("scr_classing_read(): invalid spec\n  - ", paste(errs, collapse = "\n  - "), call. = FALSE)
  attr(d, "sep") <- sep
  class(d) <- c("scr_classing_spec", "data.frame")
  d
}

#' @rdname scr_classing_spec
#' @return `scr_classing_import()` returns a named list of proposals (one
#'   per variable whose bins differ from the lab's current ones), each to be
#'   accepted or discarded.
#' @export
scr_classing_import <- function(lab, file) {
  .scr_local_verbose(lab)
  check_lab(lab, "scr_classing_import")
  sep <- lab$result$config$bin_separator
  spec <- if (is.character(file)) scr_classing_read(file, sep = sep) else file
  if (!inherits(spec, "scr_classing_spec")) stop("scr_classing_import(): `file` must be a path or an scr_classing_spec.", call. = FALSE)
  spec_sep <- attr(spec, "sep")
  if (!is.null(spec_sep) && !identical(spec_sep, sep)) {
    stop("scr_classing_import(): the spec was read with separator '", spec_sep, "' but the lab uses '", sep, "'.", call. = FALSE)
  }
  out <- list()
  for (v in intersect(unique(spec$variable), lab$features)) {
    r <- spec[spec$variable == v, , drop = FALSE]; r <- r[order(r$bin_id), , drop = FALSE]
    cur <- lab$current[[v]]
    if (identical(cur$type, "numerical")) {
      br <- r$upper[-nrow(r)]
      if (isTRUE(all.equal(br, cur$cutpoints))) next
      out[[v]] <- scr_classing_propose(lab, v, breaks = br)
    } else {
      g <- strsplit(as.character(r$categories), sep, fixed = TRUE)
      if (identical(vapply(g, paste, character(1), collapse = sep), cur$bin)) next
      ot <- if (any(r$is_other)) which(r$is_other)[1] else NULL
      out[[v]] <- scr_classing_propose(lab, v, groups = g, other_to = ot)
    }
    out[[v]]$imported_reason <- r$reason[1]
  }
  if (!length(out)) msg("  nothing to import: every variable matches the lab's current bins.")
  out
}

# -- apply ------------------------------------------------------------------ #

#' Commit the lab into a new selection result
#'
#' Returns a new `scr_result` in which the accepted manual entries replace
#' the optimal ones inside `fit` (the automatic fit is frozen as
#' `fit_auto`), the screening and hold-out rows of those variables are
#' recomputed with the very same pipeline functions, the final shortlist is
#' the one implied by [scr_classing_choose()], and the funnel, gains, SQL
#' and summary are rebuilt with a `provenance` column. `scr_selected()` on
#' the result returns the final list (`which = "consensus"` still gives the
#' automatic one). The ledger travels with the result and into
#' [scr_scorecard()] and [scr_export()]. The input result is not modified.
#'
#' @param lab An object from [scr_coarse_classing()].
#'
#' @return An `scr_result` with a `lab` component (`ledger`, `spec`,
#'   `shortlist`, `source`).
#'
#' @seealso [scr_coarse_classing()] for a complete session, from lab to
#'   scorecard.
#' @family classing
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   use_lightgbm = FALSE, xgb_rounds = 40, n_boot = 10)
#' d <- scr_demo[, c("default", "ref_date", "ds_region", "ds_band", "vl_score_01",
#'                   "vl_score_02", "vl_score_05", "vl_score_10", "vl_hist_01")]
#' res <- scr_select(d, "default", config = cfg, date_col = "ref_date")
#' lab <- scr_coarse_classing(res)
#' p <- scr_classing_propose(lab, "ds_region",
#'                           groups = list(edge = c("NORTH", "SOUTH"),
#'                                         core = c("EAST", "WEST", "CENTRE")))
#' lab <- scr_classing_accept(lab, p, reason = "edge/core is what pricing uses")
#' res2 <- scr_classing_apply(lab)
#' scr_selected(res2)
#' scr_decisions(res2)
#' @export
scr_classing_apply <- function(lab) {
  check_lab(lab, "scr_classing_apply")
  x <- lab$result; cfg <- x$config
  old <- scr_verbose(isTRUE(cfg$verbose)); on.exit(scr_verbose(old), add = TRUE)
  res <- x
  res$fit_auto <- x$fit_auto %||% x$fit
  manual <- names(lab$accepted)
  final <- .lab_final(lab)
  msg_stage("C", "coarse classing applied")
  if (length(manual)) {
    for (f in manual) {
      res$fit$results[[f]] <- lab$accepted[[f]]$entry
      i <- match(f, res$fit$summary$feature)
      res$fit$summary$algorithm[i] <- "manual"; res$fit$summary$n_bins[i] <- length(lab$accepted[[f]]$entry$bin)
      res$fit$summary$total_iv[i] <- sum(lab$accepted[[f]]$entry$iv); res$fit$summary$error[i] <- FALSE
    }
    res$fit$algorithm <- union(as.character(x$fit$algorithm), "manual")
    fit_m <- .mini_fit(x, res$fit$results[manual])
    scr <- screen_features(fit_m, cfg)
    res$screen$summary <- data.table::rbindlist(list(x$screen$summary[!feature %in% manual], scr$summary), use.names = TRUE, fill = TRUE)
    res$screen$full    <- data.table::rbindlist(list(x$screen$full[!feature %in% manual], scr$full), use.names = TRUE, fill = TRUE)
    d <- .lab_data(lab, manual)
    w_tr <- apply_woe(fit_m, d$tr, manual, "both"); w_ho <- apply_woe(fit_m, d$ho, manual, "both")
    ho <- holdout_check(w_tr, w_ho, d$y_tr, d$y_ho, manual, cfg)
    res$holdout <- data.table::rbindlist(list(x$holdout[!feature %in% manual], ho), use.names = TRUE, fill = TRUE)
    msg("  %d variable(s) with manual bins: %s", length(manual), lst(manual))
  }
  provenance <- stats::setNames(rep("auto", length(names(res$fit$results))), names(res$fit$results))
  provenance[manual] <- "manual:rebin"
  forced <- setdiff(final, x$consensus$selected); dropped <- setdiff(x$consensus$selected, final)
  provenance[forced]  <- ifelse(forced %in% manual, "manual:rebin+add", "manual:add")
  provenance[dropped] <- ifelse(dropped %in% manual, "manual:rebin+drop", "manual:drop")
  reasons <- lab$choice$reasons
  for (f in manual) if (is.null(reasons[[f]])) reasons[[f]] <- lab$accepted[[f]]$reason
  lab_slot <- list(ledger = lab$ledger, spec = scr_classing_spec(lab), source = lab$source,
                   provenance = provenance, reasons = reasons, author = lab$author,
                   applied_at = Sys.time(), laplace = lab$laplace,
                   shortlist = list(consensus = x$consensus$selected,
                                    manual = if (length(forced) || length(dropped) || !is.null(lab$choice$keep)) final else NULL,
                                    final = final))
  res$lab <- lab_slot
  res <- .lab_refresh(res)
  msg("  final shortlist: %d (consensus %d, forced %d, dropped %d)", length(final), length(x$consensus$selected),
      length(forced), length(dropped))
  res
}

#' Rebuild funnel, gains, SQL, summary and meta after a lab commit
#' @keywords internal
#' @noRd
.lab_refresh <- function(res) {
  cfg <- res$config
  bins <- list(screen = res$screen, holdout = res$holdout, prune = res$prune, derived_excluded = res$derived_excluded)
  models <- list(votes = res$models$votes, consensus = res$consensus, metrics = res$models$metrics)
  final <- res$lab$shortlist$final
  res$funnel <- build_funnel(res$split$cols, res$triage, bins, models, cfg, selected = final, lab = res$lab)
  res$gains  <- build_gains(bins, models, cfg, selected = final)
  res$sql    <- build_sql_woe(res$fit, res$triage$ledger, final, cfg, res$target, provenance = .provenance_line(res$lab))
  res$meta$n_manual_bins <- sum(res$lab$source == "manual")
  res$meta$n_forced_in <- length(setdiff(final, res$consensus$selected))
  res$meta$n_manual_dropped <- length(setdiff(res$consensus$selected, final))
  res$meta$lab <- TRUE
  res$summary_md <- build_summary(res$meta, res$funnel, models, cfg, lab = res$lab)
  res
}

#' One-line provenance statement for SQL headers and model cards
#' @keywords internal
#' @noRd
.provenance_line <- function(lab_slot) {
  if (is.null(lab_slot)) return(NULL)
  man <- names(lab_slot$source)[lab_slot$source == "manual"]
  forced <- setdiff(lab_slot$shortlist$final, lab_slot$shortlist$consensus)
  dropped <- setdiff(lab_slot$shortlist$consensus, lab_slot$shortlist$final)
  sprintf("Provenance: %d manually binned (%s), %d forced in (%s), %d dropped (%s) - see the decision ledger",
          length(man), lst(man, 5), length(forced), lst(forced, 5), length(dropped), lst(dropped, 5))
}

#' Decision ledger of a lab, a result or a scorecard
#'
#' Returns the append-only ledger of manual decisions: every proposal
#' accepted, discarded or superseded, every forced or dropped variable and
#' every override, each with its reason.
#'
#' @param x An `scr_classing` lab, an `scr_result` from
#'   [scr_classing_apply()] or an `scr_scorecard` fitted on one.
#'
#' @return A `data.table`, one row per decision (append-only), or an empty
#'   one when no manual decision exists.
#'
#' @seealso [scr_coarse_classing()] for a complete session, from lab to
#'   scorecard.
#' @family classing
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   use_lightgbm = FALSE, xgb_rounds = 40, n_boot = 10)
#' d <- scr_demo[, c("default", "ref_date", "ds_region", "ds_band", "vl_score_01",
#'                   "vl_score_02", "vl_score_05", "vl_score_10", "vl_hist_01")]
#' res <- scr_select(d, "default", config = cfg, date_col = "ref_date")
#' lab <- scr_coarse_classing(res)
#' p <- scr_classing_propose(lab, "ds_region",
#'                           groups = list(edge = c("NORTH", "SOUTH"),
#'                                         core = c("EAST", "WEST", "CENTRE")))
#' lab <- scr_classing_accept(lab, p, reason = "edge/core is what pricing uses")
#' scr_decisions(lab)
#' scr_decisions(scr_classing_apply(lab))
#' scr_decisions(res)   # no manual decision: an empty ledger
#' @export
scr_decisions <- function(x) {
  if (inherits(x, "scr_classing")) return(x$ledger[])
  if (inherits(x, "scr_result")) return((x$lab$ledger %||% .empty_ledger())[])
  if (inherits(x, "scr_scorecard")) return((x$decisions %||% .empty_ledger())[])
  stop("scr_decisions() expects a lab, an scr_result or an scr_scorecard.", call. = FALSE)
}

#' @export
print.scr_classing <- function(x, ...) {
  n_acc <- sum(unlist(x$proposals) == "accepted"); n_dis <- sum(unlist(x$proposals) == "discarded")
  cat(sprintf("<scr_classing> target \"%s\" | opened %s by %s | %d variables | %d proposals: %d accepted, %d discarded\n",
              x$target, format(x$opened_at, "%Y-%m-%d %H:%M"), x$author, length(x$features), x$ids$n, n_acc, n_dis))
  if (length(x$accepted)) {
    cat(sprintf("  %-26s %-8s %7s %17s %17s %-11s %s\n", "variable", "action", "bins", "IV train", "IV hold-out", "verdict", "reason"))
    for (f in names(x$accepted)) {
      a <- x$accepted[[f]]; o <- .optimal_checks(x, f)$summary; m <- a$checks$summary
      cat(sprintf("  %-26s %-8s %2d->%-3d %8.4f->%-8.4f %8.4f->%-8.4f %-11s %s\n", f, "accepted", o$n_bins, m$n_bins,
                  o$total_iv, m$total_iv, o$iv_holdout, m$iv_holdout, a$verdict, substr(a$reason, 1, 40)))
    }
  }
  dis <- x$ledger[action == "discard"]
  for (i in seq_len(nrow(dis))) cat(sprintf("  %-26s %-8s %s: %s\n", dis$variable[i], "discard", dis$proposal_id[i], substr(dis$reason[i], 1, 40)))
  fin <- .lab_final(x)
  cat(sprintf("  final choice: %d variables | consensus %d | force: %s | drop: %s\n", length(fin), length(x$result$consensus$selected),
              lst(x$choice$force, 5), lst(x$choice$drop, 5)))
  n_over <- nrow(x$ledger[action == "override"])
  if (n_over) cat(sprintf("  %d override(s) on record - see scr_decisions(lab)\n", n_over))
  invisible(x)
}

#' @rdname scr_export
#' @export
scr_export.scr_classing <- function(x, dir, stamp = TRUE, ...) {
  .scr_local_verbose(x)
  .need_openxlsx()
  out_dir <- .export_dir(dir, stamp)
  tag <- .file_tag(x$target)
  checks <- data.table::rbindlist(lapply(x$features, function(f) {
    ck <- if (identical(x$source[[f]], "manual")) x$accepted[[f]]$checks else .optimal_checks(x, f)
    cbind(source = x$source[[f]], ck$summary)
  }), fill = TRUE)
  bins <- data.table::rbindlist(lapply(x$features, function(f) {
    ck <- if (identical(x$source[[f]], "manual")) x$accepted[[f]]$checks else .optimal_checks(x, f)
    cbind(source = x$source[[f]], ck$bins)
  }), fill = TRUE)
  files <- list(xlsx = .scr_write_xlsx(list(
    "01_Spec" = scr_classing_spec(x), "02_Bins" = bins, "03_Checks" = checks, "04_Ledger" = x$ledger),
    file.path(out_dir, sprintf("classing_%s.xlsx", tag))))
  msg("  %s", files$xlsx)
  x$files <- files
  invisible(x)
}
