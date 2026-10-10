# ============================================================================ #
# tiers.R - Likert tiers of a score (3, 5 or 7 levels of risk or propensity)
# ============================================================================ #

#' Tiers of a score: a few labeled levels of risk or propensity
#'
#' Groups the score into a small number of contiguous tiers (typically 3, 5
#' or 7, labeled from "very low" to "very high") **fitted on the reference
#' sample** and read on the reference and the study samples. Tiers are a
#' policy and communication device; they are not a rating scale in the
#' sense of the internal ratings-based approach (see [scr_grades()] for
#' that).
#'
#' @section Method:
#'
#' The reference scores are cut into at most `max_bins` pre-bins of equal
#' share (tie-safe, as in [scr_bands()]); adjacent pre-bins are then pooled
#' (pool adjacent violators) until the event rate rises monotonically
#' toward the event-rich side of the score. Every tier is a run of these
#' blocks, so the tier rates are monotone on the reference.
#'
#' `method = "optimal"` searches every segmentation of the blocks into
#' `n_tiers` contiguous tiers with an exact dynamic program and keeps the
#' one that maximizes the binomial log-likelihood
#' \deqn{\sum_t \left[e_t \ln p_t + (n_t - e_t) \ln(1 - p_t)\right],
#'       \qquad p_t = e_t / n_t}
#' (`criterion = "deviance"`, the same as minimizing the deviance) or the
#' information value (`criterion = "iv"`), subject to: every tier holds at
#' least `min_pct` of the volume; at least `min_events` events **and**
#' `min_events` non-events; and every pair of adjacent tiers is distinct by
#' a one-sided Fisher exact test at `alpha / (L - 1)`, where \eqn{L} is the
#' tier count being tried. The rates and the objective use the weighted
#' counts; the event constraints and the Fisher test use the unweighted
#' counts. When no segmentation into `n_tiers` tiers meets the constraints,
#' `n_tiers - 2`, `n_tiers - 4`, ... down to 2 are tried (an odd count stays
#' odd while possible); when even 2 tiers are infeasible, the last resort is
#' a single tier. Every attempt is recorded in `ledger` and a warning is
#' raised. The tiers are never relabeled silently: the labels follow the
#' number of tiers achieved.
#'
#' The previous tier enters the state of the program, so its cost is
#' \eqn{O(L M^3)} time and \eqn{L M^2} memory for \eqn{L} tiers and \eqn{M}
#' blocks (compiled code); the Fisher tests are cached when \eqn{M \le 200}
#' and recomputed above. `n_tiers` is capped at 9 and `max_bins` (hence
#' \eqn{M}) at 500, and the stability study refits once per resample, so
#' `n_boot` multiplies the cost.
#'
#' `method = "anchored"` places one cut per value of `anchors` (event-rate
#' thresholds; `"overall"` stands for the reference event rate): the cut
#' sits before the first block, from the low-rate side, whose event rate
#' (its lower Jeffreys bound with `conservative = TRUE`) reaches the anchor.
#' `n_tiers` is then `length(anchors) + 1` and the argument is ignored.
#' `method = "quantile"` cuts the reference into `n_tiers` tiers of equal
#' share, the baseline.
#'
#' `round_to` rounds every cut to the nearest multiple (a policy-friendly
#' cut-off, such as 500 or 520 points); the tiers are then re-evaluated on
#' every sample with the rounded cuts, and both sets of cuts are kept. When
#' the scores were pooled into cells (`max_cells`), the count table is
#' rebuilt with the rounded cuts as forced cell edges, so the reported
#' tiers agree exactly with [scr_apply()] and [scr_sql()].
#' Tiers are left-closed: `score >= cut` is the upper side, the convention
#' of [scr_cutoff()].
#'
#' With `n_boot > 0`, the stability of the fit is measured by redrawing the
#' event and non-event counts of the reference pre-bins from a multinomial
#' law and refitting: per cut, the median and interquartile range of the
#' refitted cut in score units, and the agreement rate (the share of the
#' reference volume that keeps its tier).
#'
#' @section Labels:
#'
#' Tiers are numbered by event rate, lowest first: 2 tiers are labeled
#' `"low"`, `"high"`; 3 tiers `"low"`, `"medium"`, `"high"`; 4 tiers
#' `"low"`, `"medium low"`, `"medium high"`, `"high"`; 5 tiers `"very low"`,
#' `"low"`, `"medium"`, `"high"`, `"very high"`; 6 tiers `"very low"`,
#' `"low"`, `"medium low"`, `"medium high"`, `"high"`, `"very high"`; 7
#' tiers add `"extremely low"` and `"extremely high"` to the five; 1, 8 and 9
#' tiers are numbered `"T1"`, `"T2"`, ... The labels
#' describe the event rate, so under `objective = "risk"` they read as
#' risk and under `"propensity"` as propensity (`measure`). The table lists
#' the event-richest tier first.
#'
#' For production, every label also exists with its order in front
#' (`tier_label`): `"01."` for the event-richest tier, the first row of the
#' table, then `"02."`, ... down to the tier with the lowest event rate,
#' under every objective and direction, and for labels given in `labels`
#' too. With five tiers of a credit score, tier 5 is `"01.very high"` and
#' tier 1 is `"05.very low"`: the number `tier` rises with the event rate,
#' the prefix sorts from the highest rate down. The prefix is zero-padded to
#' two digits. [scr_apply()] and [scr_sql()] return these numbered labels
#' (`numbered = FALSE` gives the plain ones), so their output joins to the
#' table by `tier_label`.
#'
#' @inheritParams scr_bands
#' @param n_tiers Number of tiers requested (2 to 9).
#' @param method `"optimal"`, `"anchored"` or `"quantile"`.
#' @param criterion For `"optimal"`: `"deviance"` (binomial log-likelihood)
#'   or `"iv"` (information value).
#' @param anchors For `"anchored"`: event-rate thresholds, as numbers in
#'   (0, 1) or the string `"overall"` (the reference event rate); a
#'   character vector may mix both.
#' @param conservative For `"anchored"`: compare the anchors with the lower
#'   Jeffreys bound of the block rate instead of the rate.
#' @param min_pct Smallest volume share of a tier. For a scorecard, `NULL`
#'   uses `config$tier_min_pct` (0.05).
#' @param min_events Fewest events, and fewest non-events, of a tier. For a
#'   scorecard, `NULL` uses `config$tier_min_events` (20).
#' @param alpha Significance level of the adjacency test (divided by
#'   `n_tiers - 1`) and of `all_distinct`.
#' @param max_bins Pre-bins of the search (2 to 500). For a scorecard, `NULL`
#'   uses `config$tier_max_bins` (100).
#' @param labels Optional labels, one per tier achieved, in ascending order
#'   of the event rate. They get the order prefix of `tier_label` too.
#' @param round_to Optional positive number: cuts are rounded to its
#'   multiples.
#' @param n_boot Bootstrap resamples of the stability study (`0`, the
#'   default, skips it).
#' @param level Confidence level of the Jeffreys intervals. For a scorecard,
#'   `NULL` uses `config$study_level` (0.95).
#' @param seed Seed of the stability bootstrap. A number is local to the
#'   call (the user's random stream is restored on exit); `NULL` draws from
#'   the user's stream and advances it. For a scorecard, `NULL` uses
#'   `config$seed`.
#'
#' @return An object of class `c("scr_study_tiers", "scr_study", "list")`:
#'   \describe{
#'     \item{`table`}{One row per sample and tier, event-richest tier first:
#'       `sample`, `tier`, `label`, `tier_label` (the label with its order
#'       in front, `"01."` for the event-richest tier; see the section
#'       Labels), `score_lo`, `score_hi`, `n`, `pct`,
#'       `events`, `rate`, `rate_lo`, `rate_hi`, `lift`, `pct_event`,
#'       `pct_nonevent`, `woe`, `p_adjacent` (one-sided Fisher exact test
#'       that the tier has a higher event rate than the next lower tier) and
#'       `p_adjacent_adj` (Holm).}
#'     \item{`summary`}{One row per sample: `sample`, `n`, `events`, `rate`,
#'       `n_tiers`, `monotone` (the rates rise with the tier), `all_distinct`
#'       (every `p_adjacent_adj` below `alpha`), `iv` and `psi` (of the tier
#'       mix against the reference).}
#'     \item{`cuts`, `cuts_raw`}{The cuts in use (rounded when `round_to` is
#'       given) and the fitted ones.}
#'     \item{`labels`, `tier_labels`}{The tier labels, plain and numbered,
#'       lowest rate first.}
#'     \item{`codes`, `code_labels`}{Tier number and plain label of every
#'       interval in ascending score order, used by [scr_apply()] and
#'       [scr_sql()].}
#'     \item{`method`, `criterion`, `measure`, `n_tiers_requested`,
#'       `n_tiers`}{The fit.}
#'     \item{`ledger`}{One row per step: `step`, `n_tiers`, `status` and
#'       `detail`.}
#'     \item{`stability`}{`NULL`, or with `n_boot > 0` a list: `cuts` (per
#'       cut: `score`, `median`, `q25`, `q75`, `iqr` and `n_same`, the
#'       resamples that gave the same number of tiers), `agreement`,
#'       `same_count` and `n_boot`.}
#'     \item{`objective`, `direction`, `level`, `alpha`, `min_pct`,
#'       `min_events`, `reference`, `samples`, `target`, `hist`, `call`}{
#'       The settings and the count table.}
#'   }
#'
#' @references
#' Brown, L. D., Cai, T. T. and DasGupta, A. (2001). Interval estimation for
#' a binomial proportion. *Statistical Science*, 16(2), 101-133.
#' \doi{10.1214/ss/1009213286}
#'
#' Siddiqi, N. (2006). *Credit Risk Scorecards: Developing and Implementing
#' Intelligent Credit Scoring*. Wiley.
#'
#' Yurdakul, B. and Naranjo, J. (2020). Statistical properties of the
#' population stability index. *Journal of Risk Model Validation*, 14(4),
#' 89-100.
#'
#' @seealso [scr_bands()], [scr_rag()]; [scr_apply()] and [scr_sql()] assign
#'   the tiers in production.
#' @family score-studies
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   use_lightgbm = FALSE, xgb_rounds = 40, n_boot = 20)
#' res <- scr_select(scr_demo, "default", config = cfg, drop = c("id", "churn"),
#'                   date_col = "ref_date")
#' sc <- scr_scorecard(res)
#' tr <- scr_tiers(sc, n_tiers = 5)
#' tr
#' tr$table[sample == "holdout", .(tier, label, score_lo, score_hi, pct, rate)]
#'
#' # policy cuts on round numbers, and the tiers assigned to new scores
#' tr10 <- scr_tiers(sc, n_tiers = 3, round_to = 10)
#' tr10$cuts
#' head(scr_apply(tr10, c(480, 530, 600)))
#'
#' # anchors on the event rate: below, around and above the overall rate
#' scr_tiers(sc, method = "anchored", anchors = c(0.08, "overall", 0.25))$summary
#' @export
scr_tiers <- function(x, ...) UseMethod("scr_tiers")

#' @rdname scr_tiers
#' @export
scr_tiers.scr_scorecard <- function(x, n_tiers = 5L, method = c("optimal", "anchored", "quantile"),
                                    criterion = c("deviance", "iv"), anchors = NULL, conservative = FALSE,
                                    level = NULL, min_pct = NULL, min_events = NULL, alpha = 0.05,
                                    max_bins = NULL, labels = NULL, round_to = NULL, n_boot = 0L, seed = NULL,
                                    sample = "holdout", reference = "train", max_cells = 1e5, ...) {
  .study_dots(list(...), "scr_tiers")
  cfg <- x$config
  inp <- .study_input(x, sample = sample, reference = reference, max_cells = max_cells, fn = "scr_tiers")
  .tiers_fit(inp, n_tiers, match.arg(method), match.arg(criterion), anchors, conservative,
             level %||% cfg$study_level %||% 0.95, min_pct %||% cfg$tier_min_pct %||% 0.05,
             min_events %||% cfg$tier_min_events %||% 20L, alpha, max_bins %||% cfg$tier_max_bins %||% 100L,
             labels, round_to, n_boot, seed %||% cfg$seed, sys.call(),
             rebuild = function(br) .study_input(x, sample = sample, reference = reference, max_cells = max_cells,
                                                 breaks = br, fn = "scr_tiers"))
}

#' @rdname scr_tiers
#' @export
scr_tiers.data.frame <- function(x, score = "score", y = "y", objective = "risk", direction = NULL,
                                 weight = NULL, sample = NULL, reference = NULL, study = NULL, counts = FALSE,
                                 n = "n", events = "events", n_tiers = 5L,
                                 method = c("optimal", "anchored", "quantile"), criterion = c("deviance", "iv"),
                                 anchors = NULL, conservative = FALSE, level = 0.95, min_pct = 0.05,
                                 min_events = 20L, alpha = 0.05, max_bins = 100L, labels = NULL, round_to = NULL,
                                 n_boot = 0L, seed = NULL, max_cells = 1e5, ...) {
  .study_dots(list(...), "scr_tiers")
  inp <- .study_input(x, score = score, y = y, objective = objective, direction = direction, weight = weight,
                      sample = sample, reference = reference, study = study, counts = counts, n = n,
                      events = events, max_cells = max_cells, fn = "scr_tiers")
  .tiers_fit(inp, n_tiers, match.arg(method), match.arg(criterion), anchors, conservative, level, min_pct,
             min_events, alpha, max_bins, labels, round_to, n_boot, seed, sys.call(),
             rebuild = function(br) .study_input(x, score = score, y = y, objective = objective, direction = direction,
                                                 weight = weight, sample = sample, reference = reference, study = study,
                                                 counts = counts, n = n, events = events, max_cells = max_cells,
                                                 breaks = br, fn = "scr_tiers"))
}

#' Default tier labels by count, lowest event rate first
#'
#' Words for 2 to 7 tiers (an even count has no "medium"); 1, 8 and 9 tiers
#' are numbered "T1", "T2", ...
#' @keywords internal
#' @noRd
.tier_labels <- function(L) {
  switch(as.character(L),
         "2" = c("low", "high"),
         "3" = c("low", "medium", "high"),
         "4" = c("low", "medium low", "medium high", "high"),
         "5" = c("very low", "low", "medium", "high", "very high"),
         "6" = c("very low", "low", "medium low", "medium high", "high", "very high"),
         "7" = c("extremely low", "very low", "low", "medium", "high", "very high", "extremely high"),
         paste0("T", seq_len(L)))
}

#' Tier labels prefixed by their order, "01" for the event-richest tier
#'
#' `labels` are in tier order (lowest event rate first), so the order runs
#' the other way. The prefix is zero-padded to at least two digits, which
#' makes the labels sort from the highest event rate to the lowest.
#' @keywords internal
#' @noRd
.tier_numbered <- function(labels) {
  L <- length(labels)
  if (!L) return(character())
  paste0(formatC(rev(seq_len(L)), width = max(2L, nchar(L)), flag = "0"), ".", labels)
}

#' Stack pool-adjacent-violators: block of every input so that e / n does
#' not decrease along the index
#'
#' Adjacent blocks are pooled while the earlier one has the higher rate; a
#' block without known outcomes (`n = 0`) is pooled with its neighbor.
#' Rates are compared by cross-multiplication, so no division by zero.
#' @keywords internal
#' @noRd
.study_pav <- function(e, n) {
  m <- length(e)
  if (!m) return(integer())
  be <- bn <- numeric(m); first <- integer(m); top <- 0L
  for (i in seq_len(m)) {
    top <- top + 1L; be[top] <- e[i]; bn[top] <- n[i]; first[top] <- i
    while (top > 1L && (bn[top - 1L] <= 0 || bn[top] <= 0 || be[top - 1L] * bn[top] > be[top] * bn[top - 1L])) {
      be[top - 1L] <- be[top - 1L] + be[top]; bn[top - 1L] <- bn[top - 1L] + bn[top]
      top <- top - 1L
    }
  }
  findInterval(seq_len(m), first[seq_len(top)])
}

#' Reference pre-bins, oriented so that index 1 is the event-poorest
#'
#' `bcut[j]` is the score cut between oriented pre-bins `j` and `j + 1`.
#' @keywords internal
#' @noRd
.tier_prebins <- function(cells, max_bins, side) {
  pc <- .study_cuts(cells, max_bins, "uniform", NULL, side)$cuts
  P <- length(pc) + 1L
  S <- .study_sum(cells, findInterval(cells$s, pc) + 1L, P, c("n", "n_y", "e", "w2_y", "n_y_raw", "e_raw"))
  ori <- if (identical(side, "high")) seq_len(P) else rev(seq_len(P))
  list(vol = S$n[ori], n = S$n_y[ori], e = S$e[ori], w2 = S$w2_y[ori], nr = S$n_y_raw[ori], er = S$e_raw[ori],
       bcut = if (identical(side, "high")) pc else rev(pc))
}

#' PAV blocks of the pre-bins, with the last pre-bin of every block
#' @keywords internal
#' @noRd
.tier_blocks <- function(pre) {
  g <- .study_pav(pre$e, pre$n)
  blk <- lapply(pre[c("vol", "n", "e", "w2", "nr", "er")], function(v) as.numeric(rowsum(v, g, reorder = TRUE)))
  blk$last <- cumsum(tabulate(g))
  blk
}

#' @keywords internal
#' @noRd
.tier_led <- function(step, n_tiers, status, detail) {
  data.table::data.table(step = step, n_tiers = as.integer(n_tiers), status = status, detail = detail)
}

#' Fit the tiers on oriented pre-bins (or, for the main quantile fit, on the cells)
#'
#' Returns the ascending cuts, the ledger rows and the requested count.
#' @keywords internal
#' @noRd
.tier_fit <- function(pre, cells, method, n_tiers, criterion, anchors, conservative, level, min_pct,
                      min_events, alpha, side) {
  led <- list()
  if (identical(method, "quantile")) {
    if (!is.null(cells)) {
      cs <- .study_cuts(cells, n_tiers, "uniform", NULL, side)
      cuts <- cs$cuts
    } else {
      # equal shares over the pre-bins, counted from the event-rich side
      P <- length(pre$vol)
      j <- .study_bounds(rev(pre$vol), seq_len(n_tiers - 1L) / n_tiers)
      cuts <- sort(pre$bcut[P - j])
    }
    st <- if (length(cuts) + 1L < n_tiers) "fewer" else "fitted"
    led[[1]] <- .tier_led("quantile", length(cuts) + 1L, st,
                          if (st == "fewer") sprintf("%d tiers requested; tied scores leave %d", n_tiers, length(cuts) + 1L)
                          else "equal-share tiers")
    return(list(cuts = cuts, ledger = led, n_req = n_tiers))
  }
  blk <- .tier_blocks(pre)
  M <- length(blk$e)
  led[[1]] <- .tier_led("prebins", M, "pooled",
                        sprintf("%d pre-bins, %d monotone blocks after pooling adjacent violators", length(pre$e), M))
  if (identical(method, "anchored")) {
    a <- .tier_anchor_values(anchors, sum(pre$e) / sum(pre$n))
    rate <- blk$e / blk$n
    crit <- if (isTRUE(conservative)) .study_jeffreys(rate * .study_kish(blk$n, blk$w2), .study_kish(blk$n, blk$w2), level)$lo else rate
    bnd <- integer()
    for (k in seq_along(a)) {
      j <- which(crit >= a[k])[1]
      if (is.na(j)) {
        led[[length(led) + 1L]] <- .tier_led("anchor", NA, "no cut", sprintf("anchor %s is above every block rate", .g3(a[k])))
      } else if (j == 1L) {
        led[[length(led) + 1L]] <- .tier_led("anchor", NA, "no cut", sprintf("anchor %s is at or below the lowest block rate", .g3(a[k])))
      } else {
        if (blk$last[j - 1L] %in% bnd) {
          led[[length(led) + 1L]] <- .tier_led("anchor", NA, "merged", sprintf("anchor %s gives the same cut as a lower anchor", .g3(a[k])))
        }
        bnd <- c(bnd, blk$last[j - 1L])
      }
    }
    cuts <- sort(unique(pre$bcut[unique(bnd)]))
    n_req <- length(a) + 1L
    led[[length(led) + 1L]] <- .tier_led("anchored", length(cuts) + 1L, if (length(cuts) + 1L < n_req) "fewer" else "fitted",
                                         sprintf("%d anchor(s): %s", length(a), paste(.g3(a), collapse = ", ")))
    return(list(cuts = cuts, ledger = led, n_req = n_req))
  }
  # optimal: n_tiers, then n_tiers - 2, ... down to 2
  crit <- if (identical(criterion, "iv")) 1L else 0L
  cand <- unique(c(seq(n_tiers, 2L, by = -2L), 2L))
  bnd <- integer()
  for (L in cand) {
    r <- cpp_tier_dp(blk$e, blk$n, blk$vol, blk$er, blk$nr, as.integer(L), crit, as.double(min_pct),
                     as.double(min_events), as.double(alpha))
    if (isTRUE(r$feasible)) {
      bnd <- blk$last[r$ends[-L]]
      led[[length(led) + 1L]] <- .tier_led("optimal", L, "fitted",
                                           sprintf("%s %.4f over %d blocks", criterion, r$objective, M))
      break
    }
    led[[length(led) + 1L]] <- .tier_led("optimal", L, "infeasible",
      sprintf("no %d tiers of the %d blocks with share >= %s, events and non-events >= %s, adjacent tiers distinct at %s",
              L, M, .g3(min_pct), format(min_events), .g3(alpha / (L - 1L))))
  }
  if (!length(bnd)) led[[length(led) + 1L]] <- .tier_led("optimal", 1L, "single", "no feasible segmentation: one tier")
  list(cuts = sort(pre$bcut[bnd]), ledger = led, n_req = n_tiers)
}

#' Anchor values: numbers, or "overall" for the reference event rate
#' @keywords internal
#' @noRd
.tier_anchor_values <- function(anchors, overall) {
  if (is.null(anchors) || !length(anchors)) stop("scr_tiers(): method = \"anchored\" needs `anchors`.", call. = FALSE)
  if (is.numeric(anchors)) a <- as.double(anchors) else if (is.character(anchors)) {
    a <- suppressWarnings(as.double(anchors))
    a[!is.na(anchors) & anchors == "overall"] <- overall
  } else stop("scr_tiers(): `anchors` must be numeric or character.", call. = FALSE)
  if (anyNA(a) || any(a <= 0) || any(a >= 1)) {
    stop("scr_tiers(): `anchors` must be event rates in (0, 1) or \"overall\".", call. = FALSE)
  }
  sort(unique(a))
}

#' Tiers fitted on the reference, read on every sample
#' @keywords internal
#' @noRd
.tiers_fit <- function(inp, n_tiers, method, criterion, anchors, conservative, level, min_pct, min_events,
                       alpha, max_bins, labels, round_to, n_boot, seed, call, rebuild = NULL) {
  fn <- "scr_tiers"
  n_tiers <- .study_whole(n_tiers, "n_tiers", fn, lower = 2)
  # the program costs O(L M^3): keep L and M in a range that answers in seconds
  if (n_tiers > 9L) stop(fn, "(): `n_tiers` must be at most 9 (got ", n_tiers, ").", call. = FALSE)
  level <- .study_level(level, fn)
  .scr_num1(min_pct, "min_pct", lower = 0, upper = 0.5)
  min_events <- .study_whole(min_events, "min_events", fn)
  .scr_num1(alpha, "alpha", lower = 0, upper = 1, open_lower = TRUE)
  max_bins <- .study_whole(max_bins, "max_bins", fn, lower = 2)
  if (max_bins > 500L) stop(fn, "(): `max_bins` must be at most 500 (got ", max_bins, ").", call. = FALSE)
  n_boot <- .study_whole(n_boot, "n_boot", fn)
  if (!is.logical(conservative) || length(conservative) != 1L || is.na(conservative)) {
    stop(fn, "(): `conservative` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(round_to)) .scr_num1(round_to, "round_to", lower = 0, open_lower = TRUE)
  if (!is.null(labels) && (!is.character(labels) || anyNA(labels))) {
    stop(fn, "(): `labels` must be a character vector without NA.", call. = FALSE)
  }
  h <- inp$hist
  ref_cells <- .study_cells(h, inp$reference)
  if (!nrow(ref_cells)) stop(fn, "(): the reference sample '", inp$reference, "' has no scored row.", call. = FALSE)
  side <- .study_side(inp$direction)
  pre <- .tier_prebins(ref_cells, max_bins, side)
  fit <- .tier_fit(pre, ref_cells, method, n_tiers, criterion, anchors, conservative, level, min_pct,
                   min_events, alpha, side)
  ledger <- fit$ledger
  cuts_raw <- fit$cuts
  cuts <- cuts_raw
  warn <- character()
  if (length(cuts_raw) + 1L < fit$n_req) {
    warn <- sprintf("%d tiers requested, %d achieved", fit$n_req, length(cuts_raw) + 1L)
  }
  if (!is.null(round_to) && length(cuts)) {
    cuts <- sort(unique(round(cuts_raw / round_to) * round_to))
    st <- if (length(cuts) < length(cuts_raw)) "merged" else "rounded"
    ledger[[length(ledger) + 1L]] <- .tier_led("round", length(cuts) + 1L, st,
      sprintf("cuts rounded to multiples of %s%s", format(round_to),
              if (st == "merged") sprintf("; %d cuts coincide after rounding", length(cuts_raw) - length(cuts)) else ""))
    if (st == "merged") warn <- c(warn, sprintf("rounding to %s merged tiers (%d left)", format(round_to), length(cuts) + 1L))
  }
  L <- length(cuts) + 1L
  labs <- if (is.null(labels)) .tier_labels(L) else {
    if (length(labels) != L) {
      stop(fn, "(): `labels` has ", length(labels), " entries but ", L, " tiers were achieved (see the ledger).", call. = FALSE)
    }
    labels
  }
  if (length(warn)) warning(fn, "(): ", paste(warn, collapse = "; "), " - see `ledger`.", call. = FALSE)
  # tier number of every ascending interval: 1 is the lowest event rate
  codes <- if (side == "high") seq_len(L) else rev(seq_len(L))
  # production labels carry their order: "01" is the event-richest tier
  num <- .tier_numbered(labs)
  # pooled cells can straddle a rounded cut: count again with the rounded
  # cuts as forced edges, so the tables agree with scr_apply() and the SQL
  if (!is.null(round_to) && length(cuts) && isTRUE(inp$meta$quantized) && is.function(rebuild)) {
    h <- rebuild(cuts)$hist
    ledger[[length(ledger) + 1L]] <- .tier_led("round", L, "recounted",
      "scores pooled into cells: tables counted again with the rounded cuts as forced cell edges")
  }
  ref_sh <- .study_ref_shares(.study_cells(h, inp$reference), cuts)
  tabs <- list(); summ <- list()
  for (nm in inp$samples) {
    cells <- .study_cells(h, nm)
    tb <- .study_band_table(cells, cuts, if (!identical(nm, inp$reference)) ref_sh, level, inp$direction)
    o <- attr(tb, "interval"); raw <- attr(tb, "raw")
    tier <- codes[o]
    # each tier against the next lower one (the next row)
    pa <- c(.study_fisher_gt(raw$e[-L], raw$n[-L], raw$e[-1L], raw$n[-1L]), NA_real_)[seq_len(L)]
    tt <- data.table::data.table(sample = nm, tier = tier, label = labs[tier], tier_label = num[tier], tb[, list(score_lo, score_hi, n, pct, events, rate, rate_lo, rate_hi, lift, pct_event, pct_nonevent, woe)],
                                 p_adjacent = pa, p_adjacent_adj = .study_holm(pa))
    tabs[[nm]] <- tt
    ps <- attr(tb, "psi")
    r <- tt$rate[order(tt$tier)]; r <- r[!is.na(r)]
    pv <- tt$p_adjacent_adj[!is.na(tt$p_adjacent_adj)]
    NY <- sum(cells$n_y)
    summ[[nm]] <- data.table::data.table(
      sample = nm, n = sum(cells$n), events = sum(cells$e), rate = if (NY > 0) sum(cells$e) / NY else NA_real_,
      n_tiers = L, monotone = if (length(r) >= 2L) all(diff(r) >= 0) else NA,
      all_distinct = if (length(pv)) all(pv < alpha) else NA,
      iv = if (all(is.na(tb$iv))) NA_real_ else sum(tb$iv, na.rm = TRUE),
      psi = if (is.null(ps)) NA_real_ else ps$psi)
  }
  stab <- if (n_boot > 0L) .tier_stability(pre, ref_cells, cuts_raw, method, fit$n_req, criterion, anchors,
                                           conservative, level, min_pct, min_events, alpha, side, n_boot, seed) else NULL
  structure(list(
    table = data.table::rbindlist(tabs), summary = data.table::rbindlist(summ), cuts = cuts, cuts_raw = cuts_raw,
    labels = labs, tier_labels = num, codes = codes, code_labels = labs[codes], method = method,
    criterion = if (identical(method, "optimal")) criterion else NA_character_,
    measure = inp$objective, n_tiers_requested = fit$n_req, n_tiers = L,
    ledger = data.table::rbindlist(ledger), stability = stab, objective = inp$objective, direction = inp$direction,
    level = level, alpha = alpha, min_pct = min_pct, min_events = min_events, round_to = round_to,
    reference = inp$reference, samples = inp$samples, target = inp$target, quantized = inp$meta$quantized,
    weighted = inp$meta$weighted, sql_table = inp$sql_table, sql_dialect = inp$sql_dialect, hist = h, call = call),
    class = c("scr_study_tiers", "scr_study", "list"))
}

#' Bootstrap stability of the tier cuts
#'
#' Every resample draws the event and non-event counts of the reference
#' pre-bins from one multinomial law (the weighted shares, with the
#' unweighted count of rows with a known outcome as the size) and refits
#' with the same settings. Cuts are matched by rank when the resample gives
#' the same number of tiers.
#' @keywords internal
#' @noRd
.tier_stability <- function(pre, cells, cuts, method, n_req, criterion, anchors, conservative, level, min_pct,
                            min_events, alpha, side, n_boot, seed) {
  .scr_local_seed(seed)
  P <- length(pre$e)
  prob <- c(pre$e, pmax(pre$n - pre$e, 0))
  size <- round(sum(pre$nr))
  K <- length(cuts)
  codes_of <- function(L) if (side == "high") seq_len(L) else rev(seq_len(L))
  t_main <- codes_of(K + 1L)[findInterval(cells$s, cuts) + 1L]
  vol <- sum(cells$n)
  cb <- matrix(NA_real_, n_boot, max(K, 1L))
  same <- agree <- numeric(n_boot)
  if (size >= 1 && sum(prob) > 0) {
    for (b in seq_len(n_boot)) {
      d <- as.double(stats::rmultinom(1L, size, prob))
      eb <- d[seq_len(P)]; nb <- eb + d[P + seq_len(P)]
      pb <- list(vol = nb, n = nb, e = eb, w2 = nb, nr = nb, er = eb, bcut = pre$bcut)
      fb <- .tier_fit(pb, NULL, method, n_req, criterion, anchors, conservative, level, min_pct, min_events, alpha, side)
      kb <- length(fb$cuts)
      same[b] <- kb == K
      if (kb == K && K > 0L) cb[b, ] <- fb$cuts
      t_b <- codes_of(kb + 1L)[findInterval(cells$s, fb$cuts) + 1L]
      agree[b] <- sum(cells$n[t_b == t_main]) / vol
    }
  }
  q <- if (K > 0L) apply(cb, 2L, stats::quantile, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE) else NULL
  tab <- data.table::data.table(cut = seq_len(K), score = cuts,
                                median = if (K) q[2, ] else numeric(), q25 = if (K) q[1, ] else numeric(),
                                q75 = if (K) q[3, ] else numeric(), iqr = if (K) q[3, ] - q[1, ] else numeric(),
                                n_same = rep(sum(same), K))
  list(cuts = tab, agreement = mean(agree), same_count = mean(same), n_boot = n_boot)
}

#' @export
print.scr_study_tiers <- function(x, ...) {
  cat(sprintf("<scr_study_tiers> target \"%s\" | measure %s | %s\n", x$target, x$measure, x$direction))
  study <- setdiff(x$samples, x$reference)
  if (!length(study)) study <- x$reference
  cat(sprintf("  %s%s | %d tiers requested, %d achieved | fitted on '%s'%s\n", x$method,
              if (is.na(x$criterion)) "" else sprintf(" (%s)", x$criterion), x$n_tiers_requested, x$n_tiers,
              x$reference, if (!is.null(x$round_to)) sprintf(" | cuts rounded to %s", format(x$round_to)) else ""))
  cat(sprintf("  cuts: %s\n", if (length(x$cuts)) paste(format(x$cuts), collapse = ", ") else "(none)"))
  .study_print_summary(x$summary, x$level)
  for (nm in study) {
    t <- x$table[x$table$sample == nm]
    cat(sprintf("\nTiers on '%s' (event-richest first)\n", nm))
    # the numbered label, as scr_apply() and scr_sql() assign it (objects fitted before it existed print the plain one)
    lab <- if (is.null(t$tier_label)) t$label else t$tier_label
    cat(sprintf("  %4s %-18s %-24s %7s %8s %-19s %8s\n", "tier", "label", "score", "pct", "rate",
                sprintf("[%.0f%% CI]", 100 * x$level), "p_adj"))
    for (i in seq_len(nrow(t))) {
      cat(sprintf("  %4d %-18s %-24s %7s %8s %-19s %8s\n", t$tier[i], substr(lab[i], 1, 18),
                  substr(sprintf("[%s, %s)", format(t$score_lo[i]), format(t$score_hi[i])), 1, 24),
                  .study_f(100 * t$pct[i], "%.1f%%"), .study_f(100 * t$rate[i], "%.2f%%"),
                  if (is.na(t$rate_lo[i])) "" else sprintf("[%.2f%%, %.2f%%]", 100 * t$rate_lo[i], 100 * t$rate_hi[i]),
                  .study_f(t$p_adjacent_adj[i], "%.3f")))
    }
  }
  notes <- x$ledger[x$ledger$status %in% c("infeasible", "fewer", "merged", "no cut", "single")]
  if (nrow(notes)) {
    cat("\nLedger\n")
    for (i in seq_len(nrow(notes))) cat(sprintf("  %s: %s\n", notes$step[i], notes$detail[i]))
  }
  if (!is.null(x$stability)) {
    cat(sprintf("\nStability (%d resamples): tier agreement %.1f%%, same tier count %.1f%%\n", x$stability$n_boot,
                100 * x$stability$agreement, 100 * x$stability$same_count))
  }
  invisible(x)
}

# data.table column names used without quotes in this file
utils::globalVariables(c("score_lo", "score_hi", "pct", "events", "rate_lo", "rate_hi", "pct_event",
                         "pct_nonevent", "woe"))
