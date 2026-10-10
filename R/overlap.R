# ============================================================================ #
# overlap.R - rules against the score: what each one catches
# ============================================================================ #
# One grouped pass counts the rows per pattern of rule flags and score
# alert. Every set (a rule, the score, both, either alone, any rule) is a
# sum over that small table, so the cost after the pass grows with the
# number of patterns and rules, not with the rows.
# ============================================================================ #

#' Overlap of rules and a score
#'
#' Compares a set of rules (0/1 flags) with the alerts of a score at one
#' cut: what each rule catches, what the score also catches, what only one
#' of them catches, and which rules the score makes redundant. Built for
#' fraud, where expert rules and a model alert on the same transactions.
#'
#' @section Score alerts:
#'
#' The score alerts the rows on its event-rich side: `score >= cut` under
#' `higher_is_riskier` and `score < cut` under `higher_is_safer`, the
#' convention of [scr_cutoff()]. Give `cut`, or `alert_share`, the share of
#' the rows to alert: the cut is then the boundary between two distinct
#' scores nearest to that share (tie-safe, as in [scr_bands()]), and the
#' share realized is reported in the summary (`share_score`).
#'
#' @section Rules:
#'
#' Per rule, the rows it flags are split by the score alert:
#' \itemize{
#'   \item `n_rule`, `events_rule`, `precision_rule`: every row the rule
#'     flags;
#'   \item `n_both`, `events_both`, `precision_both`: flagged by the rule
#'     and alerted by the score;
#'   \item `n_rule_only`, `events_rule_only`, `precision_rule_only`: flagged
#'     by the rule, not alerted by the score;
#'   \item `n_score_only`, `events_score_only`, `precision_score_only`:
#'     alerted by the score, not flagged by the rule;
#'   \item `caught_by_score` = `events_both / events_rule`, the share of the
#'     events of the rule that the score alerts too, and `retire_candidate`
#'     = `caught_by_score >= retire_at` (`NA` for a rule without events):
#'     the score already catches what the rule catches.
#' }
#' A precision is the event rate of the set, events over rows with a known
#' outcome, with its Jeffreys interval (`_lo`, `_hi`; on the Kish effective
#' size under weights). With `value`, `value_rule`, `value_both`,
#' `value_rule_only` and `value_score_only` are the sums of the value over
#' the events of each set, and `caught_by_score_value` the share in value.
#' A missing rule flag counts as not flagged.
#'
#' @section Summary:
#'
#' `recall_score`, `recall_rules` (any rule) and `recall_any` (the score or
#' any rule) are the shares of all events alerted; `incr_score` =
#' `recall_any - recall_rules` is what the score adds to the rules and
#' `incr_rules` = `recall_any - recall_score` what the rules add to the
#' score. The `value_` columns are the same shares of the event value. The
#' alert volumes are `n_score`, `n_rules` and `n_any`, each with its share
#' of all rows and its precision.
#'
#' Rows with a missing or infinite score, or a zero weight, are left out
#' (`n_dropped`); rows with a missing outcome count in the volumes only.
#'
#' @section Cost:
#'
#' The rows are counted once per distinct pattern of rule flags and score
#' alert, and every set is a sum over that table. Time and memory after the
#' pass grow with the number of distinct patterns times the number of rules
#' (`n_patterns` is reported). A few dozen rules that seldom fire together
#' give a small table; many dense, unrelated rules can give nearly one
#' pattern per row, in which case pass the rules in smaller groups.
#'
#' @inheritParams scr_bands
#' @param x A `data.frame` with one row per case (a transaction).
#' @param rules Names of the rule columns: 0/1 numbers or logicals.
#' @param alert_share Share of the rows the score alerts, in (0, 1\].
#' @param cut Instead of `alert_share`: the score cut of the alerts.
#' @param value Optional column of a value per case (the amount).
#' @param direction `"higher_is_riskier"` (default, a fraud score) or
#'   `"higher_is_safer"`; `NULL` derives it from `objective`.
#' @param retire_at Share of the events of a rule caught by the score from
#'   which the rule is a candidate for retirement.
#' @param level Confidence level of the Jeffreys intervals.
#'
#' @return An object of class `c("scr_overlap", "list")`:
#'   \describe{
#'     \item{`table`}{One row per rule (see the section Rules).}
#'     \item{`summary`}{One row: `n`, `events`, `n_score`, `share_score`,
#'       `precision_score`, `n_rules`, `share_rules`, `precision_rules`,
#'       `n_any`, `share_any`, `precision_any`, `recall_score`,
#'       `recall_rules`, `recall_any`, `incr_score`, `incr_rules` and, with
#'       `value`, `value_events`, `value_recall_score`,
#'       `value_recall_rules`, `value_recall_any`, `value_incr_score` and
#'       `value_incr_rules`.}
#'     \item{`cut`, `alert_share`, `rules`, `retire_at`, `level`,
#'       `objective`, `direction`, `score`, `target`, `value`, `n`,
#'       `n_rows`, `n_dropped`, `n_patterns`, `weighted`, `call`}{The cut
#'       and the settings; `n_patterns` is the number of distinct patterns
#'       of rule flags and score alert.}
#'   }
#'
#' @references
#' Brown, L. D., Cai, T. T. and DasGupta, A. (2001). Interval estimation for
#' a binomial proportion. *Statistical Science*, 16(2), 101-133.
#' \doi{10.1214/ss/1009213286}
#'
#' @seealso [scr_operating()] for the cut of the alerts under a capacity,
#'   [scr_detection()] for the time to detection, [scr_score_cross()] for
#'   two scores on the same rows.
#' @family score-studies
#' @examples
#' local({
#'   set.seed(1)
#'   n <- 20000
#'   x <- rnorm(n)
#'   fraud <- rbinom(n, 1, plogis(-5 + 1.5 * x))
#'   d <- data.frame(score = round(100 * plogis(x + rnorm(n, sd = 0.5))), y = fraud,
#'                   amount = round(rexp(n, 1 / 80), 2),
#'                   # one rule the score covers, one that sees something else
#'                   rule_velocity = as.integer(x > 1.6),
#'                   rule_new_device = rbinom(n, 1, ifelse(fraud == 1, 0.3, 0.01)))
#'   ov <- scr_overlap(d, rules = c("rule_velocity", "rule_new_device"), alert_share = 0.05,
#'                     value = "amount")
#'   print(ov)
#'   ov$table[, c("rule", "n_rule", "precision_rule", "caught_by_score", "retire_candidate")]
#' })
#' @export
scr_overlap <- function(x, ...) UseMethod("scr_overlap")

#' @rdname scr_overlap
#' @export
scr_overlap.data.frame <- function(x, score = "score", y = "y", rules, alert_share = NULL, cut = NULL, value = NULL,
                                   objective = "risk", direction = "higher_is_riskier", retire_at = 0.95,
                                   level = 0.95, weight = NULL, max_cells = 1e5, ...) {
  fn <- "scr_overlap"
  .study_dots(list(...), fn)
  for (nm in c("score", "y", "value", "weight")) {
    v <- get(nm)
    if (!is.null(v)) .study_chr1(v, nm, fn)
  }
  if (missing(rules) || !is.character(rules) || !length(rules) || anyNA(rules) || anyDuplicated(rules)) {
    stop(fn, "(): `rules` must be the distinct names of the rule columns.", call. = FALSE)
  }
  rd <- .study_one_direction(objective, direction, NULL, "risk", fn)
  miss <- setdiff(c(score, y, rules, value, weight), names(x))
  if (length(miss)) stop(fn, "(): column(s) ", lst(miss), " not in `x`.", call. = FALSE)
  for (nm in c(score, value, weight)) {
    if (!is.numeric(x[[nm]])) stop(fn, "(): column '", nm, "' must be numeric.", call. = FALSE)
  }
  if (is.null(alert_share) == is.null(cut)) stop(fn, "(): give `alert_share` or `cut`, one of them.", call. = FALSE)
  if (!is.null(alert_share)) .study_num1(alert_share, "alert_share", fn, lower = 0, upper = 1, open_lower = TRUE)
  if (!is.null(cut) && (!is.numeric(cut) || length(cut) != 1L || is.na(cut))) {
    stop(fn, "(): `cut` must be a single number.", call. = FALSE)
  }
  .study_num1(retire_at, "retire_at", fn, lower = 0, upper = 1, open_lower = TRUE)
  level <- .study_level(level, fn)
  sc <- as.double(x[[score]])
  yy <- .scr_y01(x[[y]], fn)
  ok <- is.finite(sc)
  wtd <- !is.null(weight)
  if (wtd) {
    w <- as.double(x[[weight]])
    if (anyNA(w) || any(!is.finite(w)) || any(w < 0)) stop(fn, "(): the weights must be finite and non-negative.", call. = FALSE)
    # a row with zero weight does not belong to the population
    ok <- ok & w > 0
  }
  all_in <- all(ok)
  rows <- function(v) if (all_in) v else v[ok]
  sc <- rows(sc); yy <- rows(yy)
  if (wtd) w <- rows(w)
  if (!length(sc)) stop(fn, "(): no row has a score (and a positive weight).", call. = FALSE)
  side <- .study_side(rd$direction)
  # the cut of the alerts: tie-safe from the share, on the count table of the score
  if (is.null(cut)) {
    cut <- .study_share_cut(.study_hist(sc, integer(length(sc)), w = if (wtd) w, max_cells = max_cells, fn = fn),
                            alert_share, side)
  }
  cut <- as.double(cut)

  # one grouped pass: rows per pattern of rule flags and score alert
  R <- length(rules)
  rn <- paste0("r__", seq_len(R))
  cols <- list(a = as.integer(if (side == "high") sc >= cut else sc < cut))
  for (k in seq_len(R)) {
    v <- x[[rules[k]]]
    if (!is.logical(v) && !(is.numeric(v) && all(v %in% c(0, 1, NA)))) {
      stop(fn, "(): rule '", rules[k], "' must be a 0/1 or logical column.", call. = FALSE)
    }
    v <- as.integer(rows(v))
    # a missing flag is a rule that did not fire
    v[is.na(v)] <- 0L
    cols[[rn[k]]] <- v
  }
  kn <- !is.na(yy); ev <- kn & yy == 1L
  ww <- if (wtd) w else 1
  cols$k <- ww * kn; cols$e <- ww * ev; cols$q <- ww * ww * kn
  agg <- if (wtd) alist(n = sum(w), n_y = sum(k), e = sum(e), w2 = sum(q)) else
    alist(n = .N, n_y = sum(k), e = sum(e), w2 = sum(q))
  if (wtd) cols$w <- w
  if (!is.null(value)) {
    vv <- as.double(rows(x[[value]])); vv[is.na(vv)] <- 0
    cols$ve <- vv * ww * ev
    agg <- c(agg, alist(ve = sum(ve)))
  }
  dt <- data.table::setDT(cols)
  j <- as.call(c(as.name("list"), agg))
  G <- dt[, eval(j), keyby = c("a", rn)]
  cnt <- setdiff(names(G), c("a", rn))
  X <- as.matrix(G[, cnt, with = FALSE]); storage.mode(X) <- "double"
  M <- as.matrix(G[, rn, with = FALSE]); storage.mode(M) <- "double"
  A <- as.double(G$a)

  # every set is a sum of pattern counts: the rule, the rule and the score, and what is left of each
  rule <- crossprod(M, X)          # R x counts
  both <- crossprod(M * A, X)
  alert <- colSums(X * A)
  tot <- colSums(X)
  sets <- list(rule = rule, both = both, rule_only = pmax(rule - both, 0),
               score_only = pmax(matrix(alert, R, length(cnt), byrow = TRUE, dimnames = dimnames(rule)) - both, 0))
  tab <- data.table::data.table(rule = rules)
  for (s in names(sets)) {
    r <- .overlap_rate(sets[[s]], level)
    for (cn in names(r)) data.table::set(tab, j = sprintf(cn, s), value = r[[cn]])
  }
  # share of the events of the rule that the score alerts too
  caught <- unname(ifelse(rule[, "e"] > 0, both[, "e"] / rule[, "e"], NA_real_))
  data.table::set(tab, j = "caught_by_score", value = caught)
  data.table::set(tab, j = "retire_candidate", value = caught >= retire_at)
  if (!is.null(value)) {
    for (s in names(sets)) data.table::set(tab, j = paste0("value_", s), value = unname(sets[[s]][, "ve"]))
    data.table::set(tab, j = "caught_by_score_value",
                    value = unname(ifelse(rule[, "ve"] != 0, both[, "ve"] / rule[, "ve"], NA_real_)))
  }

  # summary: the score, any rule, and the two together
  any_r <- as.double(rowSums(M) > 0)
  rules_ <- colSums(X * any_r)
  union_ <- colSums(X * pmax(any_r, A))
  rec <- function(v, cn) if (tot[[cn]] > 0) v[[cn]] / tot[[cn]] else NA_real_
  prec <- function(v) if (v[["n_y"]] > 0) v[["e"]] / v[["n_y"]] else NA_real_
  summ <- data.table::data.table(
    n = tot[["n"]], events = tot[["e"]], n_score = alert[["n"]], share_score = alert[["n"]] / tot[["n"]],
    precision_score = prec(alert), n_rules = rules_[["n"]], share_rules = rules_[["n"]] / tot[["n"]],
    precision_rules = prec(rules_), n_any = union_[["n"]], share_any = union_[["n"]] / tot[["n"]],
    precision_any = prec(union_), recall_score = rec(alert, "e"), recall_rules = rec(rules_, "e"),
    recall_any = rec(union_, "e"))
  summ[, `:=`(incr_score = recall_any - recall_rules, incr_rules = recall_any - recall_score)]
  if (!is.null(value)) {
    summ[, `:=`(value_events = tot[["ve"]], value_recall_score = rec(alert, "ve"),
                value_recall_rules = rec(rules_, "ve"), value_recall_any = rec(union_, "ve"))]
    summ[, `:=`(value_incr_score = value_recall_any - value_recall_rules,
                value_incr_rules = value_recall_any - value_recall_score)]
  }
  structure(list(table = tab[], summary = summ[], cut = cut, alert_share = alert_share, rules = rules,
                 retire_at = retire_at, level = level, objective = rd$objective, direction = rd$direction,
                 score = score, target = y, value = value, n = tot[["n"]], n_rows = length(sc),
                 n_dropped = sum(!ok), n_patterns = nrow(G), weighted = wtd, call = sys.call()),
            class = c("scr_overlap", "list"))
}

#' Volume, events and precision (with a Jeffreys interval) of one set per rule
#'
#' `S` has one row per rule and the count columns `n`, `n_y`, `e` and `w2`.
#' The names are templates: `%s` takes the name of the set.
#' @keywords internal
#' @noRd
.overlap_rate <- function(S, level) {
  ny <- unname(S[, "n_y"]); e <- unname(S[, "e"])
  rate <- ifelse(ny > 0, e / ny, NA_real_)
  neff <- .study_kish(ny, unname(S[, "w2"]))
  ci <- .study_jeffreys(rate * neff, neff, level)
  out <- list(as.double(S[, "n"]), as.double(e), rate, ci$lo, ci$hi)
  names(out) <- c("n_%s", "events_%s", "precision_%s", "precision_%s_lo", "precision_%s_hi")
  out
}

#' @export
print.scr_overlap <- function(x, ...) {
  s <- x$summary
  cat(sprintf("<scr_overlap> outcome \"%s\" | score \"%s\" (%s) | %d rules\n", x$target, x$score, x$direction,
              length(x$rules)))
  cat(sprintf("  %s rows%s | %s events | the score alerts %s %s: %s rows (%s)\n", .study_n(x$n_rows),
              if (x$n_dropped > 0) sprintf(" (%s left out)", .study_n(x$n_dropped)) else "", .study_n(s$events),
              if (identical(x$direction, "higher_is_riskier")) "score >=" else "score <", format(x$cut, digits = 7),
              .study_n(s$n_score), .study_f(100 * s$share_score, "%.2f%%")))
  cat(sprintf("\n  %-10s %10s %8s %10s %8s\n", "alerts", "n", "share", "precision", "recall"))
  for (k in c("score", "rules", "any")) {
    cat(sprintf("  %-10s %10s %8s %10s %8s\n", switch(k, score = "score", rules = "any rule", any = "either"),
                .study_n(s[[paste0("n_", k)]]), .study_f(100 * s[[paste0("share_", k)]], "%.2f%%"),
                .claims_pct(s[[paste0("precision_", k)]]), .claims_pct(s[[paste0("recall_", k)]])))
  }
  cat(sprintf("  incremental recall: the score over the rules %s, the rules over the score %s\n",
              .mix_pp(s$incr_score), .mix_pp(s$incr_rules)))
  if (!is.null(x$value)) {
    cat(sprintf("  in value (\"%s\"): recall %s score, %s rules, %s either | incremental %s and %s\n", x$value,
                .claims_pct(s$value_recall_score), .claims_pct(s$value_recall_rules),
                .claims_pct(s$value_recall_any), .mix_pp(s$value_incr_score), .mix_pp(s$value_incr_rules)))
  }
  t <- x$table
  cat(sprintf("\nRules (retire candidate: the score catches at least %s of the events of the rule)\n",
              paste0(.g3(100 * x$retire_at), "%")))
  cat(sprintf("  %-22s %9s %9s %10s %9s %11s %12s %7s  %s\n", "rule", "n", "events", "precision", "n_both",
              "n_rule_only", "ev_rule_only", "caught", "retire"))
  for (i in seq_len(nrow(t))) {
    cat(sprintf("  %-22s %9s %9s %10s %9s %11s %12s %7s  %s\n", substr(t$rule[i], 1, 22), .study_n(t$n_rule[i]),
                .study_n(t$events_rule[i]), .claims_pct(t$precision_rule[i]), .study_n(t$n_both[i]),
                .study_n(t$n_rule_only[i]), .study_n(t$events_rule_only[i]), .claims_pct(t$caught_by_score[i]),
                if (is.na(t$retire_candidate[i])) "-" else if (t$retire_candidate[i]) "yes" else "no"))
  }
  invisible(x)
}

#' @rdname scr_export
#' @export
scr_export.scr_overlap <- function(x, dir, stamp = TRUE, ...) {
  .study_dots(list(...), "scr_export")
  .need_openxlsx()
  out_dir <- .export_dir(dir, stamp)
  tag <- .file_tag(x$target)
  settings <- .kv_table(list(score = x$score, target = x$target, value = x$value, objective = x$objective,
                             direction = x$direction, cut = x$cut, alert_share = x$alert_share,
                             retire_at = x$retire_at, level = x$level, rules = x$rules, n = x$n,
                             n_rows = x$n_rows, n_dropped = x$n_dropped, weighted = x$weighted))
  sheets <- lapply(list(Rules = x$table, Summary = x$summary, Settings = settings), .study_sheet)
  files <- list(xlsx = .scr_write_xlsx(sheets, file.path(out_dir, sprintf("overlap_%s.xlsx", tag))))
  for (f in files) msg("  %s", f)
  x$files <- files
  invisible(x)
}

# data.table column names used without quotes in this file
utils::globalVariables(c("ve", "recall_any",
                         "recall_rules", "recall_score", "incr_score", "incr_rules", "value_recall_any",
                         "value_recall_rules", "value_recall_score", "value_incr_score", "value_incr_rules",
                         "value_events"))
