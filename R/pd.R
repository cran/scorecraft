# ============================================================================ #
# pd.R - IRB probability of default: master scale, calibration of the
#        alignment to a central tendency, rating grades on the score, margin
#        of conservatism, floors, migration, validation battery, PIT/TTC
#        bridge, production SQL and the workbook
# ============================================================================ #
# The scorecard owns a calibrated event probability through its alignment
# (ln(odds) = I + S * raw, then the PDO map). Everything here stays on that
# contract: scr_calibrate() produces a SECOND alignment (I*, S*) anchored to
# the long-run default rate, never touching the points; scr_grades() cuts
# the production score into intervals whose PD is monotone; scr_moc() and
# scr_pd() add the margin of conservatism and the floor of the framework;
# scr_pd_validate() runs the standard test battery on a cohort panel. The
# grade is a CASE on the score, so the SQL of the scorecard extends by one
# CTE and the R and SQL paths agree by test.
# ============================================================================ #

# -- master scale ----------------------------------------------------------- #

#' Master scale of PD grades
#'
#' A grade structure with geometric midpoints and geometric-mean boundaries:
#' \deqn{PD_k = PD_1 \, r^{k-1},\quad r = (PD_K / PD_1)^{1/(K-1)},\quad
#'       \mathrm{bound}_k = \sqrt{PD_k \, PD_{k+1}},}
#' so that every grade doubles (or multiplies by `r`) the PD of the one
#' before. With `method = "supplied"` the table comes from the user: a
#' numeric vector of midpoints (boundaries derived as the geometric means)
#' or a `data.frame` with `pd_lo` and `pd_hi` (and optionally `pd_mid`,
#' `label`). Grade 1 is always the safest.
#'
#' @param pd_min,pd_max PD midpoints of the first and the last grade.
#' @param n_grades Number of grades.
#' @param method `"geometric"` (default) or `"supplied"`.
#' @param grades For `"supplied"`: a numeric vector of midpoints or a
#'   `data.frame` with `pd_lo` and `pd_hi`.
#' @param labels Optional grade labels (default `"1"`, `"2"`, ...).
#'
#' @return A `data.table` of class `scr_master_scale` with `grade`, `label`,
#'   `pd_lo`, `pd_mid`, `pd_hi`, and the attributes `ratio` (the geometric
#'   ratio between consecutive midpoints) and `method`.
#'
#' @family irb-pd
#' @examples
#' ms <- scr_master_scale(0.0005, 0.25, n_grades = 8)
#' ms
#' scr_master_scale(method = "supplied", grades = c(0.001, 0.01, 0.05, 0.20))
#' @export
scr_master_scale <- function(pd_min = 0.0003, pd_max = 0.30, n_grades = 10L,
                             method = c("geometric", "supplied"), grades = NULL, labels = NULL) {
  method <- match.arg(method)
  if (identical(method, "geometric")) {
    .scr_num1(pd_min, "pd_min", lower = 0, upper = 1, open_lower = TRUE)
    .scr_num1(pd_max, "pd_max", lower = 0, upper = 1, open_lower = TRUE)
    n_grades <- as.integer(n_grades)
    if (n_grades < 2L) stop("scr_master_scale(): `n_grades` must be at least 2.", call. = FALSE)
    if (pd_max <= pd_min) stop("scr_master_scale(): `pd_max` must exceed `pd_min`.", call. = FALSE)
    r <- (pd_max / pd_min)^(1 / (n_grades - 1L))
    mid <- pd_min * r^(seq_len(n_grades) - 1L)
    bound <- sqrt(mid[-n_grades] * mid[-1L])
    lo <- c(0, bound); hi <- c(bound, 1)
  } else {
    if (is.null(grades)) stop("scr_master_scale(): `grades` is needed with method = \"supplied\".", call. = FALSE)
    if (is.numeric(grades)) {
      mid <- sort(as.double(grades))
      if (length(mid) < 2L || any(mid <= 0 | mid >= 1)) stop("scr_master_scale(): midpoints must be in (0, 1), at least two.", call. = FALSE)
      bound <- sqrt(mid[-length(mid)] * mid[-1L])
      lo <- c(0, bound); hi <- c(bound, 1)
    } else {
      g <- as.data.frame(grades, stringsAsFactors = FALSE)
      if (!all(c("pd_lo", "pd_hi") %in% names(g))) stop("scr_master_scale(): a supplied table needs `pd_lo` and `pd_hi`.", call. = FALSE)
      g <- g[order(g$pd_lo), , drop = FALSE]
      lo <- as.double(g$pd_lo); hi <- as.double(g$pd_hi)
      if (any(hi <= lo) || any(lo < 0) || any(hi > 1)) stop("scr_master_scale(): every grade needs 0 <= pd_lo < pd_hi <= 1.", call. = FALSE)
      if (any(abs(lo[-1L] - hi[-length(hi)]) > 1e-12)) stop("scr_master_scale(): grades must be contiguous (pd_lo of a grade equal to pd_hi of the previous one).", call. = FALSE)
      mid <- if ("pd_mid" %in% names(g)) as.double(g$pd_mid) else
        ifelse(lo > 0 & hi < 1, sqrt(lo * hi), (lo + hi) / 2)
      if (is.null(labels) && "label" %in% names(g)) labels <- as.character(g$label)
    }
    n_grades <- length(mid)
    r <- if (n_grades > 1L) exp(mean(diff(log(mid)))) else NA_real_
  }
  labels <- labels %||% as.character(seq_len(n_grades))
  if (length(labels) != n_grades) stop("scr_master_scale(): `labels` must have one entry per grade.", call. = FALSE)
  out <- data.table::data.table(grade = seq_len(n_grades), label = as.character(labels), pd_lo = lo, pd_mid = mid, pd_hi = hi)
  data.table::setattr(out, "ratio", r); data.table::setattr(out, "method", method)
  data.table::setattr(out, "class", c("scr_master_scale", "data.table", "data.frame"))
  out
}

#' @export
print.scr_master_scale <- function(x, ...) {
  cat(sprintf("<scr_master_scale> %d grades (%s)%s\n", nrow(x), attr(x, "method"),
              if (is.finite(attr(x, "ratio"))) sprintf(" | ratio between midpoints %.3f", attr(x, "ratio")) else ""))
  cat(sprintf("  %-6s %-8s %10s %10s %10s\n", "grade", "label", "pd_lo", "pd_mid", "pd_hi"))
  for (i in seq_len(nrow(x))) cat(sprintf("  %-6d %-8s %10s %10s %10s\n", x$grade[i], x$label[i],
                                          fmt_pct(x$pd_lo[i], 3), fmt_pct(x$pd_mid[i], 3), fmt_pct(x$pd_hi[i], 3)))
  invisible(x)
}

# -- alignment arithmetic --------------------------------------------------- #

#' Event ln(odds) of raw scores under an alignment
#' @keywords internal
#' @noRd
.pd_event_lnodds <- function(al, raw) {
  al$sign * (al$calibration$intercept + al$calibration$slope * as.double(raw))
}

#' New alignment from event ln(odds)* = a + b * event ln(odds)
#'
#' In the orientation of the alignment: I* = sign * a + b * I, S* = b * S.
#' @keywords internal
#' @noRd
.pd_align_new <- function(al, a, b, method) {
  I0 <- al$calibration$intercept; S0 <- al$calibration$slope
  cal <- list(method = paste0("pd_", method), intercept = al$sign * a + b * I0, slope = b * S0, r2 = NA_real_,
              n_bands = 0L, bands = NULL, note = sprintf("calibrated from I = %.6f, S = %.6f (%s)", I0, S0, al$calibration$method),
              before = list(intercept = I0, slope = S0, method = al$calibration$method), a = a, b = b)
  .scr_align_from(cal$intercept, cal$slope, al$base_score, al$base_odds, al$pdo, al$direction, cal)
}

#' PD of a production score: score -> raw through the scorecard alignment,
#' raw -> PD through the calibrated one
#' @keywords internal
#' @noRd
.pd_score_to_pd <- function(al_sc, al_cal, score) {
  raw <- (as.double(score) - al_sc$a) / al_sc$b
  stats::predict(al_cal, raw, type = "prob")
}

#' Production score at which the calibrated PD equals `pd`
#' @keywords internal
#' @noRd
.pd_pd_to_score <- function(al_sc, al_cal, pd) {
  ln_or <- al_cal$sign * stats::qlogis(as.double(pd))
  raw <- (ln_or - al_cal$calibration$intercept) / al_cal$calibration$slope
  al_sc$a + al_sc$b * raw
}

#' Implied AUC of a vector of PDs: P(PD_D > PD_ND) when the PDs are true
#' @keywords internal
#' @noRd
.pd_auc_implied <- function(p) {
  p <- sort(as.double(p))
  u <- rle(p)
  grp <- rep.int(seq_along(u$lengths), u$lengths)
  P <- as.double(rowsum(p, grp, reorder = FALSE))           # expected defaults per PD value
  Q <- u$lengths - P                                        # expected non-defaults
  cumQ <- cumsum(Q)
  sum(P * (cumQ - Q + Q / 2)) / (sum(P) * sum(Q))
}

#' Solve the intercept `a` so that mean(plogis(a + b * l)) equals `ct`
#' @keywords internal
#' @noRd
.pd_solve_a <- function(l, b, ct) {
  f <- function(a) mean(stats::plogis(a + b * l)) - ct
  # f increases in a; extend the bracket when b * l reaches beyond +-60
  stats::uniroot(f, c(-60, 60), tol = 1e-13, maxiter = 2000L, extendInt = "upX")$root
}

#' Solve (a, b) so that the mean PD equals `ct` and the implied AUC equals `auc_target`
#' @keywords internal
#' @noRd
.pd_solve_ab <- function(l, ct, auc_target) {
  g <- function(lb) { b <- exp(lb); a <- .pd_solve_a(l, b, ct); .pd_auc_implied(stats::plogis(a + b * l)) - auc_target }
  lo <- log(0.02); hi <- log(50)
  glo <- g(lo); ghi <- g(hi)
  note <- ""
  if (glo >= 0) { lb <- lo; note <- "AR target below the reachable range; slope at its lower bound" }
  else if (ghi <= 0) { lb <- hi; note <- "AR target above the reachable range; slope at its upper bound" }
  else lb <- stats::uniroot(g, c(lo, hi), tol = 1e-10, maxiter = 500L)$root
  b <- exp(lb)
  list(a = .pd_solve_a(l, b, ct), b = b, note = note)
}

# -- calibration ------------------------------------------------------------ #

#' Calibrate the alignment to a central tendency
#'
#' Re-anchors the probability of default of a scorecard to a long-run
#' average default rate (the central tendency, CT) without touching the
#' points: the result is a **new** alignment `(I*, S*)` such that
#' `predict(alignment, raw, type = "prob")` is the calibrated PD, while the
#' scorecard keeps its own alignment for the score. Four methods:
#'
#' \describe{
#'   \item{`"intercept"`}{The prior-correction shift of King and Zeng
#'     (2001), \eqn{\delta = \ln[\tau(1-\bar y) / ((1-\tau)\bar y)]}, added
#'     to the event ln(odds); `S` unchanged, so the rank order and every
#'     discrimination statistic are untouched. The closed form is exact on
#'     the odds; when the calibration sample is available the shift is
#'     refined by a one-dimensional root so that the mean PD equals the CT
#'     exactly (the closed form is reported as `shift_prior`).}
#'   \item{`"logodds_ab"`}{Tasche (2013): `ln(odds*) = a + b ln(odds)`, with
#'     `(a, b)` solving `mean(PD*) = CT` and implied accuracy ratio equal to
#'     `ar_target` (default: the accuracy ratio observed on the sample).
#'     With the observed accuracy ratio this is Tasche's quasi-moment
#'     matching (QMM) proper. The implied AUC is the probability that a
#'     default has a higher PD than a non-default when the PDs are true:
#'     each score carries weight `PD*` among the defaults and `1 - PD*`
#'     among the non-defaults, ties counted one half.}
#'   \item{`"qmm"`}{The outcome-free variant of the same two equations: the
#'     target accuracy ratio is the implied AR of the current PDs, so no
#'     outcome is needed and the implied discriminatory power of the
#'     uncalibrated curve is carried over to the new level.}
#'   \item{`"scaling"`}{`PD* = PD * CT / ybar`. The proportional rescaling is
#'     not a logit map, so the slope is the least-squares projection of
#'     `logit(PD*)` on the ln(odds) and the intercept is solved to the CT.}
#' }
#'
#' @param x An [scr_scorecard()] (uses the ln(odds) and outcome of
#'   `sample`), an [scr_align()] (pass `raw` and, for the two-parameter
#'   methods, `y`) or a numeric vector of event ln(odds) (aligned directly
#'   to the default 600/50/20 scale).
#' @param target The central tendency: a number in `(0, 1)` or an `scr_dr`
#'   from [scr_default_rate()] (its `lra$mean` is used). With `segment`, a
#'   named vector with one CT per segment.
#' @param sample_rate Event rate of the calibration sample; `NULL` uses the
#'   mean of `y`.
#' @param method `"intercept"`, `"logodds_ab"`, `"scaling"` or `"qmm"`;
#'   `NULL` uses `config$pd_calibration`.
#' @param ar_target Target accuracy ratio for `"logodds_ab"` and `"qmm"`.
#' @param segment Optional vector of segment labels, one per calibration
#'   row: one alignment per segment is fitted as well.
#' @param raw,y Raw ln(odds) and 0/1 outcome when `x` is not a scorecard.
#' @param sample Sample of the scorecard used for the calibration.
#'
#' @return An object of class `scr_pd_calibration`: `alignment` (the new
#'   `scr_align`), `alignment_before`, `method`, `ct`, `target_source`,
#'   `sample_rate`, `shift` (change of the event intercept), `shift_prior`
#'   (the closed-form King-Zeng shift), `slope_ratio` (`S* / S`),
#'   `mean_pd_before`, `mean_pd_after`, `ar_before`, `ar_after` (observed),
#'   `ar_implied_before`, `ar_implied_after`, `n`, `segments` (table and
#'   alignments when `segment` is given), `ledger`.
#'   Also `ar_target`, `note` and `sample`.
#'
#' @references
#' King, G. and Zeng, L. (2001). Logistic regression in rare events data.
#' *Political Analysis*, 9(2), 137-163.
#'
#' Tasche, D. (2013). The art of probability-of-default curve calibration.
#' *Journal of Credit Risk*, 9(4), 63-103.
#'
#' @family irb-pd
#' @examples
#' set.seed(1)
#' l <- stats::qlogis(0.12) + stats::rnorm(2000)
#' y <- stats::rbinom(2000, 1, stats::plogis(l))
#' cal <- scr_calibrate(l, target = 0.04, y = y)
#' cal
#' mean(predict(cal$alignment, l, type = "prob"))
#' scr_calibrate(l, target = 0.04, y = y, method = "logodds_ab", ar_target = 0.55)
#' @export
scr_calibrate <- function(x, target, sample_rate = NULL, method = NULL, ar_target = NULL, segment = NULL,
                          raw = NULL, y = NULL, sample = "holdout") {
  if (inherits(x, "scr_scorecard")) {
    al <- x$alignment
    s <- x$samples[[sample]]
    if (is.null(s)) stop("scr_calibrate(): sample '", sample, "' does not exist in the scorecard.", call. = FALSE)
    raw <- s$link; y <- s$y
    method <- method %||% x$config$pd_calibration
  } else if (inherits(x, "scr_align")) {
    al <- x
  } else if (is.numeric(x)) {
    al <- .scr_align_from(0, -1, 600, 50, 20, "higher_is_safer",
                          list(method = "direct", intercept = 0, slope = -1, r2 = NA_real_, n_bands = 0L, bands = NULL,
                               note = "numeric event ln(odds)"))
    raw <- as.double(x)
  } else stop("scr_calibrate() expects an scr_scorecard, an scr_align or a numeric vector of ln(odds).", call. = FALSE)
  method <- match.arg(method %||% "intercept", c("intercept", "logodds_ab", "scaling", "qmm"))
  if (!is.null(raw) && !is.null(y) && length(raw) != length(y)) stop("scr_calibrate(): `raw` and `y` must have the same length.", call. = FALSE)

  # -- central tendency ------------------------------------------------------ #
  target_source <- "numeric"
  if (inherits(target, "scr_dr")) {
    ct <- target$lra$mean; target_source <- sprintf("scr_dr long-run average (%d %sly cohorts)", target$lra$n_cohorts, target$by)
  } else ct <- target
  if (!is.numeric(ct) || anyNA(ct) || any(ct <= 0 | ct >= 1)) stop("scr_calibrate(): the central tendency must lie in (0, 1).", call. = FALSE)

  # -- segments: one alignment each, then the pooled one --------------------- #
  segments <- NULL
  if (!is.null(segment)) {
    if (is.null(raw) || length(segment) != length(raw)) stop("scr_calibrate(): `segment` must have one label per calibration row.", call. = FALSE)
    segment <- as.character(segment)
    lv <- sort(unique(segment))
    cts <- if (length(ct) > 1L || !is.null(names(ct))) {
      miss <- setdiff(lv, names(ct))
      if (length(miss)) stop("scr_calibrate(): no central tendency for segment(s) ", lst(miss), call. = FALSE)
      ct[lv]
    } else stats::setNames(rep(ct, length(lv)), lv)
    fits <- lapply(lv, function(g) {
      i <- segment == g
      scr_calibrate(al, target = unname(cts[[g]]), sample_rate = NULL, method = method, ar_target = ar_target,
                    raw = raw[i], y = if (is.null(y)) NULL else y[i])
    })
    names(fits) <- lv
    tab <- data.table::rbindlist(lapply(lv, function(g) {
      f <- fits[[g]]
      data.table::data.table(segment = g, n = f$n, sample_rate = f$sample_rate, ct = f$ct, shift = f$shift,
                             slope_ratio = f$slope_ratio, mean_pd_before = f$mean_pd_before, mean_pd_after = f$mean_pd_after)
    }))
    segments <- list(table = tab, alignments = lapply(fits, `[[`, "alignment"), fits = fits)
    ct <- sum(tab$ct * tab$n) / sum(tab$n)
  }
  ct <- unname(ct[1])

  # -- the calibration sample ---------------------------------------------- #
  if (!is.null(raw)) {
    ok <- is.finite(raw) & (if (is.null(y)) TRUE else !is.na(y))
    raw <- as.double(raw[ok]); if (!is.null(y)) y <- as.integer(y[ok])
    if (!length(raw)) stop("scr_calibrate(): no calibration row with a finite ln(odds) and an outcome.", call. = FALSE)
  }
  ybar <- sample_rate %||% (if (!is.null(y)) mean(y) else NULL)
  if (is.null(ybar)) stop("scr_calibrate(): give `sample_rate` or `y`.", call. = FALSE)
  .scr_num1(ybar, "sample_rate", lower = 0, upper = 1, open_lower = TRUE)
  if (ybar >= 1) stop("scr_calibrate(): `sample_rate` must be below 1.", call. = FALSE)
  shift_prior <- log(ct * (1 - ybar) / ((1 - ct) * ybar))
  l <- if (is.null(raw)) NULL else .pd_event_lnodds(al, raw)
  p0 <- if (is.null(l)) NULL else stats::plogis(l)
  auc_obs <- if (!is.null(y) && length(unique(y)) == 2L) .auc_ks(l, y)$auc else NA_real_
  auc_impl0 <- if (is.null(p0)) NA_real_ else .pd_auc_implied(p0)
  note <- ""

  a <- shift_prior; b <- 1
  if (identical(method, "intercept")) {
    if (!is.null(l)) a <- .pd_solve_a(l, 1, ct)
  } else {
    if (is.null(l)) stop("scr_calibrate(): method \"", method, "\" needs the calibration sample (`raw`).", call. = FALSE)
    if (identical(method, "scaling")) {
      k <- ct / ybar
      p1 <- pmin(pmax(k * p0, 1e-9), 1 - 1e-9)
      b <- unname(stats::coef(stats::lm(stats::qlogis(p1) ~ l))[2])
      if (!is.finite(b) || b <= 0) { b <- 1; note <- "scaling slope not identifiable; slope kept" }
      a <- .pd_solve_a(l, b, ct)
      note <- if (nzchar(note)) note else sprintf("proportional rescaling by %.4f projected on the ln(odds)", k)
    } else {
      ar_t <- ar_target %||% (if (identical(method, "qmm")) 2 * auc_impl0 - 1 else {
        if (is.na(auc_obs)) stop("scr_calibrate(): method \"logodds_ab\" needs `y` (both classes) or `ar_target`.", call. = FALSE)
        2 * auc_obs - 1 })
      .scr_num1(ar_t, "ar_target", lower = 0, upper = 1)
      sol <- .pd_solve_ab(l, ct, (ar_t + 1) / 2)
      a <- sol$a; b <- sol$b; note <- sol$note
    }
  }
  al_new <- .pd_align_new(al, a, b, method)
  p1 <- if (is.null(l)) NULL else stats::plogis(a + b * l)
  out <- list(
    alignment = al_new, alignment_before = al, method = method, ct = ct, target_source = target_source,
    sample_rate = ybar, shift = a, shift_prior = shift_prior, slope_ratio = b,
    mean_pd_before = if (is.null(p0)) NA_real_ else mean(p0),
    mean_pd_after  = if (is.null(p1)) NA_real_ else mean(p1),
    ar_before = if (is.na(auc_obs)) NA_real_ else 2 * auc_obs - 1,
    ar_after  = if (is.na(auc_obs) || is.null(p1)) NA_real_ else 2 * .auc_ks(p1, y)$auc - 1,
    ar_implied_before = if (is.na(auc_impl0)) NA_real_ else 2 * auc_impl0 - 1,
    ar_implied_after  = if (is.null(p1)) NA_real_ else 2 * .pd_auc_implied(p1) - 1,
    ar_target = if (method %in% c("logodds_ab", "qmm")) ar_t else NA_real_,
    n = if (is.null(l)) NA_integer_ else length(l), sample = if (inherits(x, "scr_scorecard")) sample else NA_character_,
    segments = segments, note = note,
    ledger = data.table::data.table(
      action = "pd_calibration", method = method,
      detail = sprintf("CT %s (%s) | sample rate %s | shift %+.6f (prior %+.6f) | slope x %.6f%s",
                       fmt_pct(ct, 3), target_source, fmt_pct(ybar, 3), a, shift_prior, b, if (nzchar(note)) paste0(" | ", note) else ""),
      date = format(Sys.Date())))
  structure(out, class = c("scr_pd_calibration", "list"))
}

#' @export
print.scr_pd_calibration <- function(x, ...) {
  cat(sprintf("<scr_pd_calibration> method %s | CT %s (%s) | sample rate %s%s\n", x$method, fmt_pct(x$ct, 3), x$target_source,
              fmt_pct(x$sample_rate, 3), if (is.na(x$n)) "" else sprintf(" | n %s", n_fmt(x$n))))
  cat(sprintf("  event ln(odds)* = %+.6f %+.6f * ln(odds)   [prior shift %+.6f]\n", x$shift, x$slope_ratio, x$shift_prior))
  cat(sprintf("  alignment: I %.6f -> %.6f | S %.6f -> %.6f\n", x$alignment_before$calibration$intercept, x$alignment$calibration$intercept,
              x$alignment_before$calibration$slope, x$alignment$calibration$slope))
  cat(sprintf("  mean PD %s -> %s | AR observed %s -> %s | AR implied %s -> %s\n", fmt_pct(x$mean_pd_before, 3), fmt_pct(x$mean_pd_after, 3),
              .fmt_num(x$ar_before), .fmt_num(x$ar_after), .fmt_num(x$ar_implied_before), .fmt_num(x$ar_implied_after)))
  if (nzchar(x$note)) cat("  note: ", x$note, "\n", sep = "")
  if (!is.null(x$segments)) {
    t <- x$segments$table
    cat("  segments:\n")
    for (i in seq_len(nrow(t))) cat(sprintf("    %-16s n %-6s CT %s shift %+.4f\n", t$segment[i], n_fmt(t$n[i]), fmt_pct(t$ct[i], 3), t$shift[i]))
  }
  invisible(x)
}

#' @keywords internal
#' @noRd
.fmt_num <- function(x, dig = 4) if (is.na(x)) "-" else sprintf(paste0("%.", dig, "f"), x)

# -- grades ------------------------------------------------------------------ #

#' Rating grades on the score
#'
#' Cuts the production score into grades whose PD is monotone. The grade
#' boundaries are score cut points, direction-aware: grade 1 is the safest
#' (the highest scores under `higher_is_safer`). Three constructions:
#' `"geometric"` builds a [scr_master_scale()] between percentiles 1
#' and 99 of the calibrated PD and converts its PD bounds into scores
#' through the calibrated alignment; `"quantile"` cuts equal-count score
#' bands (cut points moved half-way between neighbouring scores, so a
#' boundary never sits on an observed value); `"supplied"` grades by the PD
#' bands of a given master scale.
#'
#' Grades below `min_obligors` obligors or `min_defaults` defaults are
#' merged with the neighbour of closer default rate; the sequence of grade
#' PDs is then repaired by pool-adjacent-violators when `monotone = TRUE`,
#' and every merge is recorded in `repairs`. The grade PD (`pd_be`) is the
#' long-run average of the grade default rates when a default-rate series
#' by grade is given in `dr` (`pd_source = "lra"`), the sample default rate
#' of the grade otherwise, or the mean of the calibrated individual PDs
#' (`pd_source = "mean_pd"`). Concentration is reported as the Herfindahl
#' index, the coefficient of variation of the grade shares and the
#' Herfindahl-based `hi` index.
#'
#' @section Two-pass workflow with a default-rate series:
#'
#' The series must be keyed by the final grades of this same call. Run
#' [scr_grades()] once, grade the cohort panel with [predict.scr_grades()],
#' build the series with [scr_default_rate()] (`grade =`) and pass it as
#' `dr` in a second call with identical arguments (or in [scr_moc()] and
#' [scr_pd_validate()], which read it the same way).
#'
#' @param x An [scr_scorecard()].
#' @param calibration An [scr_calibrate()] object (or its alignment);
#'   `NULL` uses the scorecard's own alignment.
#' @param master_scale An [scr_master_scale()] for `method = "supplied"`
#'   (optional for `"geometric"`).
#' @param n_grades,method,min_obligors,min_defaults,pd_source `NULL` reads
#'   `pd_n_grades`, `pd_grade_method`, `pd_min_obligors`, `pd_min_defaults`
#'   and `pd_source` from the scorecard configuration.
#' @param monotone Repair non-monotone grade PDs by pooling.
#' @param sample Sample of the scorecard used to build the grades.
#' @param dr Optional `scr_dr` with a `grade` column keyed by the final
#'   grades (see the section above).
#'
#' @return An object of class `scr_grades`: `table` (`grade`, `label`,
#'   `score_lo`, `score_hi`, `pd_lo`, `pd_hi`, `n`, `share`, `defaults`,
#'   `dr`, `pd_mean`, `pd_be`, `merged_from`, and `n_series`, `t_series`
#'   when a series is given), `breaks` (ascending score cut points),
#'   `band_grade` (grade of every score band, ascending), `direction`,
#'   `method`, `pd_source`, `master_scale`, `alignment` (calibrated),
#'   `alignment_score` (the scorecard's), `concentration` (`hhi`, `cv`,
#'   `hi`, `k`), `repairs`, `ledger`, `moc` (empty, filled by
#'   [scr_moc()]), `dr` (the pooled series), `rows` (score, outcome and
#'   grade of the sample), `scorecard`, `sample`, `ct`, `sample_rate`.
#'   Also `calibration` (the `scr_pd_calibration` when one was given),
#'   `n_grades_requested`, `min_obligors`, `min_defaults`, `target` and
#'   `config`.
#'
#' @family irb-pd
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   use_lightgbm = FALSE, xgb_rounds = 40, n_boot = 10)
#' res <- scr_select(scr_demo, "default", config = cfg, drop = c("id", "churn"),
#'                   date_col = "ref_date")
#' sc <- scr_scorecard(res)
#' cal <- scr_calibrate(sc, target = 0.06)
#' gr <- scr_grades(sc, cal, n_grades = 7, min_defaults = 10)
#' gr
#' gr$table[, c("grade", "score_lo", "score_hi", "n", "dr", "pd_be")]
#' # grade a cohort panel with the score cut points (the panel score is a
#' # different scale; the demo only shows the mechanics)
#' head(predict(gr, score = scr_demo_panel$score))
#' @export
scr_grades <- function(x, calibration = NULL, master_scale = NULL, n_grades = NULL, method = NULL,
                       min_obligors = NULL, min_defaults = NULL, monotone = TRUE, pd_source = NULL,
                       sample = "holdout", dr = NULL) {
  check_scorecard(x, "scr_grades")
  cfg <- x$config
  old <- scr_verbose(isTRUE(cfg$verbose)); on.exit(scr_verbose(old), add = TRUE)
  n_grades <- as.integer(n_grades %||% cfg$pd_n_grades)
  method <- match.arg(method %||% cfg$pd_grade_method, c("geometric", "quantile", "supplied"))
  min_obligors <- as.integer(min_obligors %||% cfg$pd_min_obligors)
  min_defaults <- as.integer(min_defaults %||% cfg$pd_min_defaults)
  pd_source <- match.arg(pd_source %||% cfg$pd_source, c("lra", "mean_pd"))
  al_sc <- x$alignment
  al_cal <- if (is.null(calibration)) al_sc else if (inherits(calibration, "scr_pd_calibration")) calibration$alignment
            else if (inherits(calibration, "scr_align")) calibration
            else stop("scr_grades(): `calibration` must come from scr_calibrate() or scr_align().", call. = FALSE)
  s <- x$samples[[sample]]
  if (is.null(s)) stop("scr_grades(): sample '", sample, "' does not exist in the scorecard.", call. = FALSE)
  score <- as.double(s$score); y <- as.integer(s$y)
  dir <- x$direction
  safer <- identical(dir, "higher_is_safer")
  pd_i <- .pd_score_to_pd(al_sc, al_cal, score)
  ledger <- list()
  add <- function(action, detail) ledger[[length(ledger) + 1L]] <<- data.table::data.table(action = action, detail = detail, date = format(Sys.Date()))

  # -- score cut points ---------------------------------------------------- #
  ms <- NULL
  if (identical(method, "quantile")) {
    if (n_grades < 2L) stop("scr_grades(): `n_grades` must be at least 2.", call. = FALSE)
    breaks <- .pd_quantile_breaks(score, n_grades)
  } else {
    if (identical(method, "supplied")) {
      if (is.null(master_scale)) stop("scr_grades(): method \"supplied\" needs `master_scale`.", call. = FALSE)
      ms <- if (inherits(master_scale, "scr_master_scale")) master_scale else scr_master_scale(method = "supplied", grades = master_scale)
    } else {
      ms <- if (!is.null(master_scale)) {
        if (inherits(master_scale, "scr_master_scale")) master_scale else scr_master_scale(method = "supplied", grades = master_scale)
      } else {
        q <- stats::quantile(pd_i, c(0.01, 0.99), names = FALSE)
        scr_master_scale(pd_min = max(1e-6, q[1]), pd_max = min(0.999, max(q[2], q[1] * 1.5)), n_grades = n_grades)
      }
    }
    bounds <- ms$pd_hi[-nrow(ms)]
    sc_b <- .pd_pd_to_score(al_sc, al_cal, bounds)
    breaks <- sort(unique(sc_b))
  }
  k0 <- length(breaks) + 1L
  band <- findInterval(score, breaks, left.open = TRUE) + 1L
  grade0 <- if (safer) k0 - band + 1L else band
  # per original grade (ordered from the safest): counts, sums, score interval
  edges <- c(-Inf, breaks, Inf)
  f0 <- factor(grade0, levels = seq_len(k0))
  og <- data.table::data.table(grade0 = seq_len(k0), band = if (safer) k0 - seq_len(k0) + 1L else seq_len(k0),
                               n = tabulate(grade0, nbins = k0), d = as.integer(rowsum(y, f0)[, 1]),
                               sum_pd = as.double(rowsum(pd_i, f0)[, 1]))
  og[, score_lo := edges[band]]; og[, score_hi := edges[band + 1L]]
  if (!is.null(ms) && nrow(ms) == k0) og[, label0 := ms$label] else og[, label0 := as.character(grade0)]

  # -- merge below the minimum counts -------------------------------------- #
  groups <- as.list(seq_len(k0))
  repairs <- list()
  gstat <- function(g) list(n = sum(og$n[g]), d = sum(og$d[g]))
  repeat {
    if (length(groups) <= 1L) break
    st <- lapply(groups, gstat)
    nn <- vapply(st, `[[`, numeric(1), "n"); dd <- vapply(st, `[[`, numeric(1), "d")
    small <- which(nn == 0 | nn < min_obligors | dd < min_defaults)
    if (!length(small)) break
    i <- small[order(nn[small], dd[small])][1]
    m <- length(groups)
    cand <- c(if (i > 1L) i - 1L, if (i < m) i + 1L)
    rate <- ifelse(nn > 0, dd / pmax(nn, 1), NA_real_)
    j <- if (length(cand) == 1L) cand else {
      dist <- abs(rate[cand] - rate[i]); dist[is.na(dist)] <- Inf
      if (all(is.infinite(dist))) cand[which.min(nn[cand])] else cand[which.min(dist)]
    }
    repairs[[length(repairs) + 1L]] <- data.table::data.table(
      step = "min_counts", merged = .pd_label(og$label0, groups[[i]]), into = .pd_label(og$label0, groups[[j]]),
      reason = sprintf("n %d, defaults %d (minimum %d / %d)", as.integer(nn[i]), as.integer(dd[i]), min_obligors, min_defaults))
    a <- min(i, j); b <- max(i, j)
    groups[[a]] <- sort(c(groups[[a]], groups[[b]])); groups[[b]] <- NULL
  }

  # -- grade PD and the monotone repair ------------------------------------ #
  stat_fn <- function(g) {
    n <- sum(og$n[g])
    if (identical(pd_source, "mean_pd")) sum(og$sum_pd[g]) / max(1, n) else sum(og$d[g]) / max(1, n)
  }
  if (isTRUE(monotone)) {
    r <- .pd_pool_monotone(groups, stat_fn, og$label0, "monotone_sample")
    groups <- r$groups; repairs <- c(repairs, r$repairs)
  }
  K <- length(groups)
  series <- NULL; series_note <- NULL
  if (!is.null(dr)) {
    series <- .pd_series_match(dr, K, "scr_grades")
    lra_fn <- function(g) .pd_lra_of(series, g)
    if (isTRUE(monotone) && identical(pd_source, "lra")) {
      r <- .pd_pool_monotone(groups, lra_fn, og$label0, "monotone_lra", by_index = TRUE)
      if (length(r$repairs)) {
        # the series is keyed by the grades before this pooling: pool it too
        series <- .pd_series_pool(series, r$map)
        groups <- r$groups; repairs <- c(repairs, r$repairs)
        K <- length(groups)
      }
    }
    series_note <- sprintf("grade PD from the long-run average of %d cohorts (%s)", length(unique(series$cohort)), dr$by)
  } else if (identical(pd_source, "lra")) {
    series_note <- "no default-rate series given: the grade PD is the sample default rate of the grade (pass `dr` for a long-run average)"
  }
  add("grade_pd_source", series_note %||% "grade PD = mean of the calibrated individual PDs")

  # -- final table --------------------------------------------------------- #
  N <- length(score)
  tab <- data.table::rbindlist(lapply(seq_len(K), function(k) {
    g <- groups[[k]]
    n <- sum(og$n[g]); d <- sum(og$d[g])
    data.table::data.table(
      grade = k, label = .pd_label(og$label0, g),
      score_lo = min(og$score_lo[g]), score_hi = max(og$score_hi[g]),
      n = as.integer(n), share = n / N, defaults = as.integer(d), dr = if (n > 0) d / n else NA_real_,
      pd_mean = if (n > 0) sum(og$sum_pd[g]) / n else NA_real_,
      merged_from = if (length(g) > 1L) paste(og$label0[g], collapse = "+") else "")
  }))
  tab[, pd_be := if (identical(pd_source, "mean_pd")) pd_mean else if (!is.null(series)) vapply(seq_len(K), function(k) .pd_lra_of(series, k), numeric(1)) else dr]
  if (!is.null(series)) {
    agg <- series[, list(n_series = sum(.SD[["n"]]), t_series = .N), by = "grade"]
    tab[, n_series := agg$n_series[match(grade, agg$grade)]]
    tab[, t_series := agg$t_series[match(grade, agg$grade)]]
  }
  # implied PD at the score bounds (direction-aware)
  pd_at <- function(sc) ifelse(is.finite(sc), .pd_score_to_pd(al_sc, al_cal, ifelse(is.finite(sc), sc, 0)), NA_real_)
  lo_pd <- pd_at(if (safer) tab$score_hi else tab$score_lo); hi_pd <- pd_at(if (safer) tab$score_lo else tab$score_hi)
  tab[, pd_lo := ifelse(is.na(lo_pd), 0, lo_pd)]; tab[, pd_hi := ifelse(is.na(hi_pd), 1, hi_pd)]
  data.table::setcolorder(tab, c("grade", "label", "score_lo", "score_hi", "pd_lo", "pd_hi", "n", "share", "defaults", "dr", "pd_mean", "pd_be", "merged_from"))
  # final breaks and the grade of each ascending score band
  grade_of_band <- integer(k0)
  for (k in seq_len(K)) grade_of_band[og$band[groups[[k]]]] <- k
  keep <- which(diff(grade_of_band) != 0L)
  breaks_final <- breaks[keep]
  band_grade <- grade_of_band[c(1L, keep + 1L)]
  grade_row <- band_grade[findInterval(score, breaks_final, left.open = TRUE) + 1L]

  conc <- .pd_concentration(tab$n)
  rep_tab <- if (length(repairs)) data.table::rbindlist(repairs) else data.table::data.table(step = character(), merged = character(), into = character(), reason = character())
  add("grades", sprintf("%d grades (%s, from %d) on sample %s | min obligors %d, min defaults %d | %d merge(s) | HHI %.3f, CV %.3f, HI %s",
                        K, method, k0, sample, min_obligors, min_defaults, nrow(rep_tab), conc$hhi, conc$cv, .fmt_num(conc$hi, 3)))
  msg("  grades: %d (%s) | HHI %.3f | %d repair(s)", K, method, conc$hhi, nrow(rep_tab))
  if (K < 2L) warning("scr_grades(): a single grade remains after merging; lower the minimum counts or the number of grades.", call. = FALSE)
  structure(list(
    table = tab[], breaks = breaks_final, band_grade = band_grade, direction = dir, method = method, pd_source = pd_source,
    master_scale = ms, alignment = al_cal, alignment_score = al_sc, concentration = conc, repairs = rep_tab,
    ledger = data.table::rbindlist(ledger), moc = .pd_moc_empty(), dr = series,
    rows = data.table::data.table(score = score, y = y, grade = grade_row, pd = pd_i),
    scorecard = x, sample = sample, n_grades_requested = n_grades, min_obligors = min_obligors, min_defaults = min_defaults,
    ct = if (inherits(calibration, "scr_pd_calibration")) calibration$ct else NA_real_,
    sample_rate = mean(y), calibration = if (inherits(calibration, "scr_pd_calibration")) calibration else NULL,
    target = x$target, config = cfg
  ), class = c("scr_grades", "list"))
}

#' Equal-count score cut points moved half-way between neighbouring scores
#' @keywords internal
#' @noRd
.pd_quantile_breaks <- function(score, n_groups) {
  probs <- seq(0, 1, length.out = n_groups + 1L)[-c(1L, n_groups + 1L)]
  q <- stats::quantile(score, probs = probs, names = FALSE, type = 7)
  u <- sort(unique(score))
  b <- vapply(q, function(v) {
    i <- findInterval(v, u)                  # u[i] <= v < u[i + 1]
    if (i < 1L) return(u[1] - 0.5)
    if (i >= length(u)) return(u[length(u)] + 0.5)
    (u[i] + u[i + 1L]) / 2
  }, numeric(1))
  sort(unique(b))
}

#' @keywords internal
#' @noRd
.pd_label <- function(labels, g) paste(labels[g], collapse = "+")

#' Pool adjacent groups until `stat_fn` is non-decreasing along the groups
#'
#' `by_index = TRUE` passes the group index (1..K) instead of its members to
#' `stat_fn` and returns `map` (old index -> new index).
#' @keywords internal
#' @noRd
.pd_pool_monotone <- function(groups, stat_fn, labels, step, by_index = FALSE) {
  repairs <- list()
  idx <- lapply(seq_along(groups), function(k) k)     # current -> original indices
  repeat {
    K <- length(groups)
    if (K <= 1L) break
    v <- vapply(seq_len(K), function(k) if (by_index) stat_fn(idx[[k]]) else stat_fn(groups[[k]]), numeric(1))
    viol <- which(diff(v) < -1e-12)
    if (!length(viol)) break
    i <- viol[1]
    repairs[[length(repairs) + 1L]] <- data.table::data.table(
      step = step, merged = .pd_label(labels, groups[[i + 1L]]), into = .pd_label(labels, groups[[i]]),
      reason = sprintf("grade PD %s below the previous grade %s", fmt_pct(v[i + 1L], 3), fmt_pct(v[i], 3)))
    groups[[i]] <- sort(c(groups[[i]], groups[[i + 1L]])); groups[[i + 1L]] <- NULL
    idx[[i]] <- c(idx[[i]], idx[[i + 1L]]); idx[[i + 1L]] <- NULL
  }
  map <- integer(sum(lengths(idx)))
  for (k in seq_along(idx)) map[idx[[k]]] <- k
  list(groups = groups, repairs = repairs, map = map)
}

#' Match a default-rate series to the final grades
#' @keywords internal
#' @noRd
.pd_series_match <- function(dr, K, fn) {
  if (!inherits(dr, "scr_dr")) stop(fn, "(): `dr` must come from scr_default_rate().", call. = FALSE)
  if (!"grade" %in% names(dr$table)) stop(fn, "(): the default-rate series needs a `grade` column (scr_default_rate(..., grade = )).", call. = FALSE)
  t <- data.table::as.data.table(dr$table)
  g <- as.character(t[["grade"]])
  lv <- as.character(seq_len(K))
  bad <- setdiff(unique(g), lv)
  if (length(bad)) stop(fn, "(): grade label(s) ", lst(bad), " of the series do not match the final grades 1..", K,
                        ". Grade the panel with predict() on this object and rebuild the series.", call. = FALSE)
  data.table::data.table(cohort = t[["cohort"]], grade = as.integer(g), n = as.integer(t[["n"]]), defaults = as.integer(t[["defaults"]]))
}

#' Long-run average default rate of a set of grades in a series
#' @keywords internal
#' @noRd
.pd_lra_of <- function(series, g) {
  s <- series[series$grade %in% g]
  if (!nrow(s)) return(NA_real_)
  a <- s[, list(n = sum(.SD[["n"]]), d = sum(.SD[["defaults"]])), by = "cohort", .SDcols = c("n", "defaults")]
  a <- a[a$n > 0]
  if (!nrow(a)) return(NA_real_)
  mean(a$d / a$n)
}

#' Pool a series by an old -> new grade map
#' @keywords internal
#' @noRd
.pd_series_pool <- function(series, map) {
  s <- data.table::copy(series)
  s[, grade := map[grade]]
  s[, list(n = sum(.SD[["n"]]), defaults = sum(.SD[["defaults"]])), by = c("cohort", "grade"), .SDcols = c("n", "defaults")][order(cohort, grade)]
}

#' Concentration of the grade shares: HHI, CV and the Herfindahl-based index
#' @keywords internal
#' @noRd
.pd_concentration <- function(n) {
  K <- length(n); R <- n / sum(n)
  cv <- sqrt(K * sum((R - 1 / K)^2))
  list(hhi = sum(R^2), cv = cv, hi = if (K > 1L) 1 + log((cv^2 + 1) / K) / log(K) else NA_real_, k = K)
}

#' @keywords internal
#' @noRd
.pd_moc_empty <- function() {
  data.table::data.table(id = integer(), category = character(), method = character(), level = numeric(), grade = integer(),
                         pd_be = numeric(), value = numeric(), reason = character(), active = logical(), date = character())
}

#' Grade of a score under an scr_grades / scr_pd object
#' @keywords internal
#' @noRd
.pd_grade_of <- function(obj, score) {
  obj$band_grade[findInterval(as.double(score), obj$breaks, left.open = TRUE) + 1L]
}

#' Grade a score vector with the cut points of an scr_grades object
#'
#' @param object An [scr_grades()] object.
#' @param score Numeric production scores.
#' @param type `"grade"` (integer grade) or `"pd"` (calibrated individual PD).
#' @param ... Ignored.
#'
#' @return A vector of the length of `score`.
#'
#' @family irb-pd
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   use_lightgbm = FALSE, xgb_rounds = 40, n_boot = 10)
#' d <- scr_demo[, c("default", "ref_date", "ds_region", "ds_band", "vl_score_01",
#'                   "vl_score_02", "vl_score_05", "vl_score_10", "vl_hist_01")]
#' res <- scr_select(d, "default", config = cfg, date_col = "ref_date")
#' sc <- scr_scorecard(res)
#' gr <- scr_grades(sc, scr_calibrate(sc, target = 0.06), n_grades = 7, min_defaults = 10)
#' predict(gr, score = c(480, 560, 640))
#' predict(gr, score = c(480, 560, 640), type = "pd")
#' @export
predict.scr_grades <- function(object, score, type = c("grade", "pd"), ...) {
  type <- match.arg(type)
  if (identical(type, "grade")) .pd_grade_of(object, score) else .pd_score_to_pd(object$alignment_score, object$alignment, score)
}

#' @export
print.scr_grades <- function(x, ...) {
  t <- x$table; c0 <- x$concentration
  src <- if (identical(x$pd_source, "lra") && is.null(x$dr)) "lra (sample default rate; pass `dr` for the series)" else x$pd_source
  cat(sprintf("<scr_grades> target \"%s\" | %d grades (%s) on %s | PD source: %s | %s\n", x$target, nrow(t), x$method, x$sample,
              src, x$direction))
  cat(sprintf("  concentration: HHI %.3f | CV %.3f | HI %s | repairs %d%s\n", c0$hhi, c0$cv, .fmt_num(c0$hi, 3), nrow(x$repairs),
              if (!is.na(x$ct)) sprintf(" | calibrated to CT %s", fmt_pct(x$ct, 3)) else ""))
  cat(sprintf("  %-5s %-7s %9s %9s %6s %6s %5s %8s %8s %8s\n", "grade", "label", "score_lo", "score_hi", "n", "share", "def", "dr", "pd_mean", "pd_be"))
  for (i in seq_len(nrow(t))) cat(sprintf("  %-5d %-7s %9.2f %9.2f %6d %5.1f%% %5d %8s %8s %8s\n", t$grade[i], substr(t$label[i], 1, 7), t$score_lo[i], t$score_hi[i],
                                          t$n[i], 100 * t$share[i], t$defaults[i], fmt_pct(t$dr[i], 2), fmt_pct(t$pd_mean[i], 2), fmt_pct(t$pd_be[i], 2)))
  if (nrow(x$repairs)) for (i in seq_len(nrow(x$repairs))) cat(sprintf("  repair (%s): %s -> %s | %s\n", x$repairs$step[i], x$repairs$merged[i], x$repairs$into[i], x$repairs$reason[i]))
  m <- x$moc[x$moc$active]
  if (nrow(m)) {
    s <- m[, list(bp = 1e4 * mean(.SD[["value"]])), by = "category"]
    cat(sprintf("  margin of conservatism (active, mean over grades): %s\n", paste(sprintf("%s %.1f bp", s$category, s$bp), collapse = " | ")))
  }
  invisible(x)
}

# -- margin of conservatism -------------------------------------------------- #

#' Margin of conservatism, by category
#'
#' Appends entries to the MoC ledger of an [scr_grades()] object. Category
#' `"C"` (general estimation error) is quantified: `"ci_timeseries"` takes
#' the upper bound of a one-sided `level` interval of the long-run average
#' from the cohort series, \eqn{t_{q, T-1}\, sd(DR_t)/\sqrt{T}} per grade;
#' `"ci_binomial"` uses \eqn{z_q \sqrt{PD(1-PD)/n}} on the obligors (or
#' obligor-years when a series exists); `"bootstrap"` resamples the
#' outcomes of the sample within each grade (drawn as the resampled default
#' rate, Binomial(n, DR) / n, its exact distribution) and takes the `level`
#' quantile of the default rate above the estimate. Categories `"A"` (data and
#' methodological deficiencies) and `"B"` (changes in standards or
#' environment) are expert quantities: `value` (one number or one per
#' grade, in PD units) and a non-empty `reason` are mandatory. The ledger
#' is append-only: `A` and `B` entries accumulate, a new `C` supersedes the
#' previous one (kept with `active = FALSE`).
#'
#' @param x An [scr_grades()] object.
#' @param category `"A"`, `"B"` or `"C"`.
#' @param method For `"C"`: `"ci_timeseries"`, `"ci_binomial"` or
#'   `"bootstrap"`; `NULL` reads `config$pd_moc_method`.
#' @param level One-sided confidence level; `NULL` reads `config$pd_moc_level`.
#' @param value For `"A"`/`"B"`: the add-on in PD units, length 1 or one per grade.
#' @param reason Justification (mandatory for `"A"`/`"B"`).
#' @param dr Optional `scr_dr` by grade for `"ci_timeseries"`, keyed by the
#'   final grades of `x`; `NULL` uses the series already stored in `x$dr` by
#'   [scr_grades()].
#' @param n_boot,seed Bootstrap resamples and seed for `"bootstrap"`.
#'
#' @return The `scr_grades` object with the entries appended to `moc`.
#'
#' @family irb-pd
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   use_lightgbm = FALSE, xgb_rounds = 40, n_boot = 10)
#' res <- scr_select(scr_demo, "default", config = cfg, drop = c("id", "churn"),
#'                   date_col = "ref_date")
#' sc <- scr_scorecard(res)
#' gr <- scr_grades(sc, n_grades = 6, min_defaults = 10)
#' gr <- scr_moc(gr, "C", method = "ci_binomial")
#' gr <- scr_moc(gr, "A", value = 0.002, reason = "missing unlikeliness-to-pay trigger before 2024")
#' gr$moc
#' @export
scr_moc <- function(x, category = c("A", "B", "C"), method = NULL, level = NULL, value = NULL, reason = NULL,
                    dr = NULL, n_boot = 200L, seed = NULL) {
  if (!inherits(x, "scr_grades")) stop("scr_moc() expects an object from scr_grades().", call. = FALSE)
  category <- match.arg(category)
  cfg <- x$config
  level <- level %||% cfg$pd_moc_level
  .scr_num1(level, "level", lower = 0, upper = 1, open_lower = TRUE)
  t <- x$table; K <- nrow(t)
  if (category %in% c("A", "B")) {
    if (is.null(value) || is.null(reason) || !nzchar(trimws(paste(reason, collapse = "")))) {
      stop("scr_moc(): categories A and B need `value` (PD units) and a non-empty `reason`.", call. = FALSE)
    }
    value <- as.double(value)
    if (!length(value) %in% c(1L, K) || anyNA(value) || any(value < 0)) stop("scr_moc(): `value` must be non-negative, of length 1 or one per grade (", K, ").", call. = FALSE)
    v <- rep_len(value, K); method <- "manual"
  } else {
    method <- match.arg(method %||% cfg$pd_moc_method, c("ci_timeseries", "ci_binomial", "bootstrap"))
    if (!is.null(dr)) { x$dr <- .pd_series_match(dr, K, "scr_moc") }
    ser <- x$dr
    if (identical(method, "ci_timeseries")) {
      if (is.null(ser)) stop("scr_moc(): \"ci_timeseries\" needs a default-rate series by grade: pass `dr` (see ?scr_grades) or use \"ci_binomial\".", call. = FALSE)
      v <- vapply(seq_len(K), function(k) {
        s <- ser[ser$grade == k & ser$n > 0]
        if (nrow(s) < 2L) return(NA_real_)
        r <- s$defaults / s$n
        stats::qt(level, nrow(s) - 1L) * stats::sd(r) / sqrt(nrow(s))
      }, numeric(1))
      if (anyNA(v)) stop("scr_moc(): every grade needs at least two cohorts in the series for \"ci_timeseries\".", call. = FALSE)
    } else if (identical(method, "ci_binomial")) {
      n_use <- if (!is.null(ser)) vapply(seq_len(K), function(k) sum(ser$n[ser$grade == k]), numeric(1)) else as.double(t$n)
      v <- stats::qnorm(level) * sqrt(t$pd_be * (1 - t$pd_be) / pmax(n_use, 1))
    } else {
      seed <- seed %||% cfg$seed
      .scr_local_seed(seed)
      # the default rate of a with-replacement resample of the n_k outcomes of
      # grade k is exactly Binomial(n_k, DR_k) / n_k: drawn directly, O(K n_boot)
      n_boot <- as.integer(n_boot)
      rows <- x$rows
      n_k <- tabulate(rows$grade, nbins = K)
      d_k <- tabulate(rows$grade[rows$y %in% 1L], nbins = K)
      B <- vapply(seq_len(K), function(k) {
        if (n_k[k] > 0L) stats::rbinom(n_boot, n_k[k], d_k[k] / n_k[k]) / n_k[k] else rep(NA_real_, n_boot)
      }, numeric(n_boot))
      B <- matrix(B, nrow = n_boot)
      v <- pmax(0, apply(B, 2L, stats::quantile, probs = level, na.rm = TRUE, names = FALSE) - t$dr)
    }
    reason <- reason %||% sprintf("estimation error, %s at %.0f%% one-sided", method, 100 * level)
    # copy first: the ledger is shared with the object passed in, and := would edit it too
    if (any(x$moc$category == "C" & x$moc$active)) {
      x$moc <- data.table::copy(x$moc)
      x$moc[x$moc$category == "C", active := FALSE]
    }
  }
  new_id <- if (nrow(x$moc)) max(x$moc$id) + 1L else 1L
  entry <- data.table::data.table(id = new_id, category = category, method = method,
                                  level = if (category == "C") level else NA_real_, grade = t$grade, pd_be = t$pd_be,
                                  value = as.double(v), reason = paste(reason, collapse = " "), active = TRUE, date = format(Sys.Date()))
  x$moc <- data.table::rbindlist(list(x$moc, entry), use.names = TRUE)
  x$ledger <- data.table::rbindlist(list(x$ledger, data.table::data.table(
    action = sprintf("moc_%s", category), detail = sprintf("%s | mean %.1f bp | %s", method, 1e4 * mean(v), paste(reason, collapse = " ")),
    date = format(Sys.Date()))), use.names = TRUE)
  x
}

# -- the PD model ----------------------------------------------------------- #

#' The PD model: grades, margin of conservatism and the floor
#'
#' Assembles the final grade table: `pd_be` from [scr_grades()], the active
#' entries of the MoC ledger by category (`A` and `B` summed over their
#' entries, the latest `C`), `pd_moc = pd_be + A + B + C`, the PD floor of
#' the asset class under the framework of `params`, and
#' `pd_final = max(pd_moc, floor)`. With `philosophy = "pit"` the
#' through-the-cycle `pd_moc` is converted with the one-factor bridge of
#' [scr_pd_pit_ttc()] before the floor (`rho` and `z` required).
#'
#' @param grades An [scr_grades()] object, after [scr_moc()].
#' @param moc `NULL` uses `grades$moc`; otherwise a ledger in the same format.
#' @param params An [scr_irb_params()]; `NULL` uses `config$framework`.
#' @param asset_class Asset class of the floor; `NULL` uses `config$asset_class`.
#' @param philosophy `"ttc"` (default) or `"pit"`.
#' @param rho,z Asset correlation and systematic factor for `"pit"`.
#'
#' @return An object of class `scr_pd`: `table` (`grade`, `label`,
#'   `score_lo`, `score_hi`, `n`, `share`, `defaults`, `dr`, `pd_be`,
#'   `moc_a`, `moc_b`, `moc_c`, `pd_moc`, `pd_ttc`, `pd_pit`, `floor`,
#'   `pd_final`, `floor_applied`), `breaks`, `band_grade`, `direction`,
#'   `alignment`, `alignment_score`, `master_scale`, `asset_class`,
#'   `framework`, `floor`, `philosophy`, `rho`, `z`, `moc_ledger`,
#'   `calibration`, `concentration`, `portfolio` (weighted `pd_be`,
#'   `pd_moc`, `pd_final`, `moc_bp`, `share_at_floor`), `scorecard`,
#'   `ledger`, `model_card`.
#'   Also `params_modified`, `repairs`, `dr`, `rows`, `pd_source`,
#'   `grade_method`, `sample`, `target`, `config` and, after [scr_export()],
#'   `files`.
#'
#' @family irb-pd
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   use_lightgbm = FALSE, xgb_rounds = 40, n_boot = 10)
#' res <- scr_select(scr_demo, "default", config = cfg, drop = c("id", "churn"),
#'                   date_col = "ref_date")
#' sc <- scr_scorecard(res)
#' gr <- scr_moc(scr_grades(sc, n_grades = 6, min_defaults = 10), "C", method = "ci_binomial")
#' pd <- scr_pd(gr)
#' pd
#' pd$table[, c("grade", "pd_be", "moc_c", "pd_final", "floor_applied")]
#' head(predict(pd, score = c(480, 560, 640), type = "pd_final"))
#' head(scr_apply(pd, head(scr_demo, 5)))
#' cat(tail(scr_sql(pd), 8), sep = "\n")
#' @export
scr_pd <- function(grades, moc = NULL, params = NULL, asset_class = NULL, philosophy = c("ttc", "pit"), rho = NULL, z = NULL) {
  if (!inherits(grades, "scr_grades")) stop("scr_pd() expects an object from scr_grades().", call. = FALSE)
  philosophy <- match.arg(philosophy)
  cfg <- grades$config
  old <- scr_verbose(isTRUE(cfg$verbose)); on.exit(scr_verbose(old), add = TRUE)
  params <- .check_params(params %||% scr_irb_params(cfg$framework), "scr_pd")
  asset_class <- asset_class %||% cfg$asset_class
  floor <- .pd_floor_of(params, asset_class)
  ledger <- data.table::copy(grades$ledger)
  add <- function(action, detail) ledger <<- data.table::rbindlist(list(ledger, data.table::data.table(action = action, detail = detail, date = format(Sys.Date()))), use.names = TRUE)
  floor_note <- if (is.na(floor)) { floor <- 0; "no PD floor for this asset class in the framework" } else sprintf("PD floor %s (%s, %s)", fmt_pct(floor, 2), asset_class, params$framework)
  add("pd_floor", floor_note)
  if (isTRUE(params$modified)) add("params_modified", "the parameter tables were edited by the user")

  ml <- moc %||% grades$moc
  if (!is.null(moc) && !data.table::is.data.table(ml)) stop("scr_pd(): `moc` must be a ledger in the format of scr_moc().", call. = FALSE)
  t <- data.table::copy(grades$table)
  K <- nrow(t)
  act <- ml[ml$active %in% TRUE]
  sum_cat <- function(cat) vapply(t$grade, function(g) sum(act$value[act$category == cat & act$grade == g]), numeric(1))
  t[, moc_a := sum_cat("A")]; t[, moc_b := sum_cat("B")]; t[, moc_c := sum_cat("C")]
  if (!any(act$category == "C")) {
    add("moc_c_missing", "no estimation-error margin (category C) recorded: pd_moc = pd_be + A + B")
    msg("  note: no category C margin of conservatism recorded")
  }
  t[, pd_moc := pmin(1, pd_be + moc_a + moc_b + moc_c)]
  t[, pd_ttc := pd_moc]
  if (identical(philosophy, "pit")) {
    if (is.null(rho) || is.null(z)) stop("scr_pd(): philosophy \"pit\" needs `rho` and `z`.", call. = FALSE)
    .scr_num1(rho, "rho", lower = 0, upper = 1, open_lower = TRUE); .scr_num1(z, "z")
    # a grade without defaults and without a margin has pd_moc 0: the bridge keeps 0 and 1 fixed
    t[, pd_pit := .vasicek_pit(pd_moc, z, rho)]
    add("philosophy", sprintf("point-in-time bridge with rho %.3f and z %+.3f", rho, z))
  } else t[, pd_pit := NA_real_]
  base_pd <- if (identical(philosophy, "pit")) t$pd_pit else t$pd_ttc
  t[, floor := floor]
  t[, pd_final := pmin(1, pmax(base_pd, floor))]
  t[, floor_applied := base_pd < floor]
  data.table::setcolorder(t, c("grade", "label", "score_lo", "score_hi", "n", "share", "defaults", "dr", "pd_be", "moc_a", "moc_b", "moc_c",
                               "pd_moc", "pd_ttc", "pd_pit", "floor", "pd_final", "floor_applied"))
  w <- t$n / sum(t$n)
  portfolio <- list(n = sum(t$n), pd_be = sum(w * t$pd_be), pd_moc = sum(w * t$pd_moc), pd_final = sum(w * t$pd_final),
                    moc_bp = 1e4 * sum(w * (t$pd_moc - t$pd_be)), share_at_floor = sum(w[t$floor_applied]), n_grades = K,
                    philosophy = philosophy)
  add("pd_model", sprintf("%d grades | portfolio PD %s -> %s (MoC %.1f bp) -> %s | %s at the floor", K, fmt_pct(portfolio$pd_be, 3),
                          fmt_pct(portfolio$pd_moc, 3), portfolio$moc_bp, fmt_pct(portfolio$pd_final, 3), fmt_pct(portfolio$share_at_floor, 1)))
  msg("  PD model: %d grades | portfolio PD %s (MoC %.1f bp) | floor %s", K, fmt_pct(portfolio$pd_final, 3), portfolio$moc_bp, fmt_pct(floor, 2))
  sc <- grades$scorecard
  out <- structure(list(
    table = t[], breaks = grades$breaks, band_grade = grades$band_grade, direction = grades$direction,
    alignment = grades$alignment, alignment_score = grades$alignment_score, master_scale = grades$master_scale,
    asset_class = asset_class, framework = params$framework, params_modified = isTRUE(params$modified), floor = floor,
    philosophy = philosophy, rho = rho, z = z, moc_ledger = ml, calibration = grades$calibration,
    concentration = grades$concentration, repairs = grades$repairs, dr = grades$dr, rows = grades$rows, pd_source = grades$pd_source,
    grade_method = grades$method, sample = grades$sample, portfolio = portfolio, scorecard = sc, target = grades$target,
    ledger = ledger, config = cfg, model_card = NULL, files = NULL
  ), class = c("scr_pd", "list"))
  out$model_card <- .pd_model_card(out)
  out
}

#' @keywords internal
#' @noRd
.pd_model_card <- function(x) {
  cal <- x$calibration; sc <- x$scorecard; p <- x$portfolio
  act <- x$moc_ledger[x$moc_ledger$active %in% TRUE]
  list(
    package = sprintf("scorecraft %s", as.character(utils::packageVersion("scorecraft"))),
    fitted_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    target = x$target, model = "pd", framework = x$framework, params_modified = x$params_modified,
    asset_class = x$asset_class, pd_floor = x$floor, philosophy = x$philosophy,
    rho = x$rho %||% NA_real_, z = x$z %||% NA_real_,
    calibration_method = if (is.null(cal)) "none (scorecard alignment)" else cal$method,
    central_tendency = if (is.null(cal)) NA_real_ else cal$ct, sample_rate = if (is.null(cal)) x$rows[, mean(y)] else cal$sample_rate,
    calibration_shift = if (is.null(cal)) NA_real_ else cal$shift,
    grade_method = x$grade_method, n_grades = nrow(x$table), pd_source = x$pd_source, sample = x$sample,
    n_sample = p$n, n_repairs = nrow(x$repairs), hhi = x$concentration$hhi, cv = x$concentration$cv, hi = x$concentration$hi,
    moc_entries = nrow(act), moc_categories = paste(sort(unique(act$category)), collapse = ", "),
    moc_bp_weighted = p$moc_bp, portfolio_pd_be = p$pd_be, portfolio_pd_moc = p$pd_moc, portfolio_pd_final = p$pd_final,
    share_at_floor = p$share_at_floor, n_series_cohorts = if (is.null(x$dr)) 0L else length(unique(x$dr$cohort)),
    scorecard_features = paste(sc$features, collapse = ", "), scorecard_auc_holdout = sc$metrics[sc$metrics$sample == "holdout", ][["auc"]][1],
    direction = x$direction, base_score = sc$scale$base_score, base_odds = sc$scale$base_odds, pdo = sc$scale$pdo
  )
}

#' @export
print.scr_pd <- function(x, ...) {
  t <- x$table; p <- x$portfolio
  cat(sprintf("<scr_pd> target \"%s\" | %s | %s | floor %s | %d grades | %s\n", x$target, x$framework, x$asset_class, fmt_pct(x$floor, 2), nrow(t), x$philosophy))
  cat(sprintf("  calibration: %s%s | grades: %s, PD source %s | HHI %.3f\n",
              if (is.null(x$calibration)) "scorecard alignment" else x$calibration$method,
              if (is.null(x$calibration)) "" else sprintf(" to CT %s", fmt_pct(x$calibration$ct, 3)), x$grade_method, x$pd_source, x$concentration$hhi))
  cat(sprintf("  %-5s %9s %9s %6s %8s %7s %7s %7s %8s %8s %5s\n", "grade", "score_lo", "score_hi", "n", "pd_be", "moc_a", "moc_b", "moc_c", "pd_moc", "pd_final", "floor"))
  for (i in seq_len(nrow(t))) cat(sprintf("  %-5d %9.2f %9.2f %6d %8s %7.1f %7.1f %7.1f %8s %8s %5s\n", t$grade[i], t$score_lo[i], t$score_hi[i], t$n[i],
                                          fmt_pct(t$pd_be[i], 2), 1e4 * t$moc_a[i], 1e4 * t$moc_b[i], 1e4 * t$moc_c[i], fmt_pct(t$pd_moc[i], 2),
                                          fmt_pct(t$pd_final[i], 2), if (t$floor_applied[i]) "yes" else ""))
  cat(sprintf("  portfolio: PD_BE %s | MoC %.1f bp (A/B/C in bp) | PD_final %s | %s of obligors at the floor\n",
              fmt_pct(p$pd_be, 3), p$moc_bp, fmt_pct(p$pd_final, 3), fmt_pct(p$share_at_floor, 1)))
  invisible(x)
}

#' Predict grade and PD from an scr_pd object
#'
#' @param object An [scr_pd()] object.
#' @param newdata A table with the source columns of the scorecard, scored
#'   with [scr_apply()]; ignored when `score` is given.
#' @param score Production scores, as an alternative to `newdata`.
#' @param type `"grade"`, `"pd"` (calibrated individual PD), `"pd_final"`
#'   (grade PD after MoC and floor) or `"score"`.
#' @param ... Ignored.
#'
#' @return A vector of the length of the input.
#'
#' @family irb-pd
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   use_lightgbm = FALSE, xgb_rounds = 40, n_boot = 10)
#' d <- scr_demo[, c("default", "ref_date", "ds_region", "ds_band", "vl_score_01",
#'                   "vl_score_02", "vl_score_05", "vl_score_10", "vl_hist_01")]
#' res <- scr_select(d, "default", config = cfg, date_col = "ref_date")
#' sc <- scr_scorecard(res)
#' pd <- scr_pd(scr_grades(sc, n_grades = 6, min_defaults = 10))
#' predict(pd, score = c(480, 560, 640), type = "pd_final")
#' head(predict(pd, newdata = scr_demo[1:10, ], type = "grade"))
#' @export
predict.scr_pd <- function(object, newdata = NULL, score = NULL, type = c("grade", "pd", "pd_final", "score"), ...) {
  type <- match.arg(type)
  if (is.null(score)) {
    if (is.null(newdata)) stop("predict.scr_pd(): give `newdata` or `score`.", call. = FALSE)
    score <- scr_apply(object$scorecard, newdata)$score
  }
  score <- as.double(score)
  switch(type,
    score = score,
    grade = .pd_grade_of(object, score),
    pd = .pd_score_to_pd(object$alignment_score, object$alignment, score),
    pd_final = object$table$pd_final[.pd_grade_of(object, score)])
}

#' @rdname scr_apply
#' @export
scr_apply.scr_pd <- function(x, newdata, ...) {
  s <- scr_apply(x$scorecard, newdata)
  g <- .pd_grade_of(x, s$score)
  data.table::data.table(score = s$score, score_points = s$score_points, grade = g,
                         pd = .pd_score_to_pd(x$alignment_score, x$alignment, s$score),
                         pd_be = x$table$pd_be[g], pd_final = x$table$pd_final[g])
}

#' @rdname scr_sql
#' @export
scr_sql.scr_pd <- function(x, table = NULL, dialect = NULL, file = NULL, ...) {
  .scr_local_verbose(x)
  base <- scr_sql(x$scorecard, table = table, dialect = dialect)
  i_final <- which(base == "SELECT"); i_final <- i_final[length(i_final)]
  i_close <- max(which(base == ")" & seq_along(base) < i_final))
  body <- base[i_final:length(base)]
  body[length(body)] <- sub(";\\s*$", "", body[length(body)])
  br <- x$breaks; bg <- x$band_grade; K <- length(bg)
  whens <- function(vals) c(sprintf("WHEN score <= %s THEN %s", .sql_num(br), vals[-K]), sprintf("ELSE %s", vals[K]))
  pd_of <- .sql_num(x$table$pd_final[bg])
  out <- c(base[seq_len(i_close - 1L)], "),", "score_scr AS (", paste0("  ", body), ")", "",
           sprintf("-- Block 4: rating grade and final PD from the score cut points (%d grades, %s)", nrow(x$table), x$direction),
           "SELECT", "    s.*,",
           sprintf("    CASE %s END AS grade,", paste(whens(as.character(bg)), collapse = " ")),
           sprintf("    CASE %s END AS pd_final", paste(whens(pd_of), collapse = " ")),
           "FROM score_scr s;")
  .sql_out(.sql_lines(out), file)
}

# -- PIT / TTC bridge -------------------------------------------------------- #

#' One-factor bridge between point-in-time and through-the-cycle PD
#'
#' Vasicek's conditional default probability:
#' \deqn{PD_{PIT} = \Phi\left(\frac{\Phi^{-1}(PD_{TTC}) - \sqrt{\rho}\, z}{\sqrt{1-\rho}}\right),}
#' and its inverse for `to = "ttc"`. A positive `z` is a benign state (lower
#' PIT PD), a negative one a stressed state.
#'
#' @param pd Numeric PDs in `(0, 1)`.
#' @param z Systematic factor (standard normal scale).
#' @param rho Asset correlation in `(0, 1)`.
#' @param to `"pit"` (input is TTC) or `"ttc"` (input is PIT).
#'
#' @return A numeric vector of the length of `pd`.
#'
#' @references
#' Vasicek, O. (2002). The distribution of loan portfolio value. *Risk*,
#' 15(12), 160-162.
#'
#' @seealso [scr_pd_stress()], the same bridge with the systematic factor
#'   given as a quantile `q` rather than a value of `z`.
#' @family irb-pd
#' @examples
#' scr_pd_pit_ttc(c(0.01, 0.05), z = -2, rho = 0.15)
#' scr_pd_pit_ttc(scr_pd_pit_ttc(0.02, z = -1, rho = 0.1), z = -1, rho = 0.1, to = "ttc")
#' @export
scr_pd_pit_ttc <- function(pd, z, rho, to = c("pit", "ttc")) {
  to <- match.arg(to)
  .scr_num1(rho, "rho", lower = 0, upper = 1, open_lower = TRUE)
  if (rho >= 1) stop("`rho` must be below 1.", call. = FALSE)
  pd <- as.double(pd); z <- as.double(z)
  if (any(pd <= 0 | pd >= 1, na.rm = TRUE)) stop("`pd` must lie in (0, 1).", call. = FALSE)
  .pd_vasicek(pd, z, rho, to)
}

#' @keywords internal
#' @noRd
.pd_vasicek <- function(pd, z, rho, to = "pit") {
  # the forward bridge is shared with the capital module (.vasicek_pit in R/capital.R)
  if (identical(to, "pit")) .vasicek_pit(pd, z, rho)
  else stats::pnorm(stats::qnorm(pd) * sqrt(1 - rho) + sqrt(rho) * z)
}

# -- migration --------------------------------------------------------------- #

#' Migration matrix between two rating dates
#'
#' Counts `N_ij` of obligors in grade `i` at the first date and grade `j`
#' at the second, the row probabilities `p_ij`, the upper and lower matrix
#' weighted bandwidths
#' \deqn{MWB_{up} = \frac{\sum_{i<j} |i-j|\, N_i\, p_{ij}}{\sum_i \max(|i-K|, |i-1|)\, N_i \sum_{j>i} p_{ij}},}
#' (and the mirror image for downgrades), the `z` statistic of every
#' off-diagonal cell against its neighbour closer to the diagonal (a
#' significantly positive value means the probability does not decay away
#' from the diagonal) and the mobility summary. Values of `grade_t1`
#' outside `1..K` count as `default`, `NA` as `closed`; both stay out of
#' the bandwidths.
#'
#' @param grade_t0,grade_t1 Integer grades at the two dates, same length.
#' @param K Number of grades; `NULL` uses the largest grade observed.
#'
#' @return An object of class `scr_migration`: `matrix` (counts, `K` rows,
#'   `K + 2` columns), `p` (row probabilities), `n` (row totals),
#'   `mwb_upper`, `mwb_lower`, `z` (`K x K`), `n_significant` (cells with
#'   `z > 1.645`), `mobility` (`share_stable`, `share_up`, `share_down`,
#'   `mean_distance`, `share_default`, `share_closed`).
#'   Also `K`, the number of grades.
#'
#' @family irb-pd
#' @examples
#' set.seed(2)
#' g0 <- sample(1:5, 500, TRUE)
#' g1 <- pmin(5, pmax(1, g0 + sample(c(-1, 0, 0, 0, 1), 500, TRUE)))
#' g1[sample(500, 10)] <- NA
#' scr_migration(g0, g1, K = 5)
#' @export
scr_migration <- function(grade_t0, grade_t1, K = NULL) {
  if (length(grade_t0) != length(grade_t1)) stop("scr_migration(): the two grade vectors must have the same length.", call. = FALSE)
  g0 <- suppressWarnings(as.integer(as.character(grade_t0)))
  ok <- !is.na(g0)
  g0 <- g0[ok]; t1 <- grade_t1[ok]
  K <- as.integer(K %||% max(c(g0, suppressWarnings(as.integer(as.character(t1)))), na.rm = TRUE))
  if (any(g0 < 1L | g0 > K)) stop("scr_migration(): `grade_t0` must lie in 1..K.", call. = FALSE)
  t1c <- as.character(t1)
  g1 <- suppressWarnings(as.integer(t1c))
  col <- ifelse(is.na(t1c), K + 2L, ifelse(!is.na(g1) & g1 >= 1L & g1 <= K, g1, K + 1L))
  M <- matrix(0L, K, K + 2L, dimnames = list(as.character(seq_len(K)), c(as.character(seq_len(K)), "default", "closed")))
  M[] <- tabulate(g0 + (col - 1L) * K, nbins = K * (K + 2L))
  n_i <- rowSums(M)
  P <- M / pmax(n_i, 1)
  idx <- seq_len(K)
  D <- abs(outer(idx, idx, "-"))
  up <- outer(idx, idx, "<"); dn <- outer(idx, idx, ">")
  Pk <- P[, idx, drop = FALSE]
  norm_u <- sum(vapply(idx, function(i) max(abs(i - K), abs(i - 1)) * n_i[i] * sum(Pk[i, up[i, ]]), numeric(1)))
  norm_l <- sum(vapply(idx, function(i) max(abs(i - K), abs(i - 1)) * n_i[i] * sum(Pk[i, dn[i, ]]), numeric(1)))
  mwb_u <- if (norm_u > 0) sum((D * Pk * n_i)[up]) / norm_u else NA_real_
  mwb_l <- if (norm_l > 0) sum((D * Pk * n_i)[dn]) / norm_l else NA_real_
  # z of each off-diagonal cell against the neighbour closer to the diagonal
  Z <- matrix(NA_real_, K, K, dimnames = list(rownames(M), as.character(idx)))
  for (i in idx) for (j in idx) {
    if (j == i || n_i[i] == 0) next
    jn <- if (j > i) j - 1L else j + 1L
    p_far <- Pk[i, j]; p_near <- Pk[i, jn]
    v <- (p_near * (1 - p_near) + p_far * (1 - p_far) + 2 * p_near * p_far) / n_i[i]
    Z[i, j] <- if (v > 0) (p_far - p_near) / sqrt(v) else NA_real_
  }
  stay <- sum(diag(Pk * n_i)); in_grades <- sum(M[, idx])
  mobility <- list(share_stable = if (in_grades > 0) stay / in_grades else NA_real_,
                   share_up = if (in_grades > 0) sum((Pk * n_i)[up]) / in_grades else NA_real_,
                   share_down = if (in_grades > 0) sum((Pk * n_i)[dn]) / in_grades else NA_real_,
                   mean_distance = if (in_grades > 0) sum((D * Pk * n_i)) / in_grades else NA_real_,
                   share_default = sum(M[, K + 1L]) / max(1, sum(M)), share_closed = sum(M[, K + 2L]) / max(1, sum(M)))
  structure(list(matrix = M, p = P, n = n_i, K = K, mwb_upper = mwb_u, mwb_lower = mwb_l, z = Z,
                 n_significant = sum(Z > stats::qnorm(0.95), na.rm = TRUE), mobility = mobility),
            class = c("scr_migration", "list"))
}

#' @export
print.scr_migration <- function(x, ...) {
  m <- x$mobility
  cat(sprintf("<scr_migration> %d grades | %s obligors | stable %s | up %s | down %s | default %s | closed %s\n", x$K, n_fmt(sum(x$matrix)),
              fmt_pct(m$share_stable), fmt_pct(m$share_up), fmt_pct(m$share_down), fmt_pct(m$share_default), fmt_pct(m$share_closed)))
  cat(sprintf("  MWB upper %s | MWB lower %s | mean distance %s | %d cell(s) not decaying from the diagonal (z > 1.645)\n",
              .fmt_num(x$mwb_upper), .fmt_num(x$mwb_lower), .fmt_num(m$mean_distance, 3), x$n_significant))
  cn <- colnames(x$p)
  cat(sprintf("  %-5s", "from"), sprintf("%7s", cn), "\n")
  for (i in seq_len(x$K)) cat(sprintf("  %-5s", rownames(x$p)[i]), sprintf("%6.1f%%", 100 * x$p[i, ]), "\n")
  invisible(x)
}

# -- validation -------------------------------------------------------------- #

#' Cohort rows of a graded panel: population at every cohort start with the
#' 12-month outcome and the grade at the end of the window
#' @keywords internal
#' @noRd
.pd_cohorts <- function(dt, horizon, by) {
  dates <- sort(unique(dt$date))
  mth <- as.integer(format(dates, "%m"))
  starts <- switch(by, month = dates, quarter = dates[mth %in% c(1L, 4L, 7L, 10L)], year = dates[mth == 1L])
  starts <- starts[.add_months(starts, horizon) <= max(dates)]
  if (!length(starts)) stop("scr_pd_validate(): no cohort has a complete ", horizon, "-month window.", call. = FALSE)
  # every cohort at once (cohort-major, panel order within a cohort): one
  # rolling join for the outcome, one equi-join for the grade at the end
  pop <- dt[dt$default %in% 0L & !is.na(dt$grade) & as.double(dt$date) %in% as.double(starts)]
  pop <- pop[order(pop$date)]
  pop[, cohort := date]
  t1 <- .add_months_u(pop$cohort, horizon)
  pop[, y := .dr_outcome(dt[dt$default %in% 1L, c("id", "date"), with = FALSE], id, cohort, t1)]
  end <- dt[, list(id, t1 = date, g1 = grade)]
  q <- data.table::data.table(id = pop$id, t1 = t1)     # built outside: in `i`, t1 would be the column of `end`
  g1 <- end[q, on = c("id", "t1"), mult = "first", g1]
  pop[, grade_t1 := data.table::fifelse(y == 1L, "default", as.character(g1))]
  pop
}

#' Calibration tests on (N, D, PD): Jeffreys, binomial, normal
#' @keywords internal
#' @noRd
.pd_cal_tests <- function(N, D, PD, alpha) {
  N <- as.double(N); D <- as.double(D)
  ok <- N > 0 & is.finite(PD) & PD > 0 & PD < 1
  z <- p_j <- p_b <- crit <- rep(NA_real_, length(N))
  p_j[ok] <- stats::pbeta(PD[ok], D[ok] + 0.5, N[ok] - D[ok] + 0.5)
  p_b[ok] <- stats::pbinom(D[ok] - 1, N[ok], PD[ok], lower.tail = FALSE)
  crit[ok] <- stats::qbinom(1 - alpha, N[ok], PD[ok]) + 1
  z[ok] <- (D[ok] / N[ok] - PD[ok]) / sqrt(PD[ok] * (1 - PD[ok]) / N[ok])
  list(p_jeffreys = p_j, p_binomial = p_b, critical = crit, z = z, p_normal = stats::pnorm(z, lower.tail = FALSE))
}

#' Hosmer-Lemeshow statistic over the grades: chi2(K)
#'
#' The grade PDs are fixed before the validation sample is observed, so the
#' statistic is referred to a chi-square with K degrees of freedom (K - 2 is
#' the in-sample case, when the PDs are fitted on the same data): Hosmer and
#' Lemeshow (2000, sec. 5.2.2); BCBS Working Paper 14 (2005), sec. 5.
#' @keywords internal
#' @noRd
.pd_hl <- function(n, d, pd) {
  ok <- n > 0 & is.finite(pd) & pd > 0 & pd < 1
  n <- as.double(n[ok]); d <- as.double(d[ok]); pd <- pd[ok]
  chi2 <- sum((d - n * pd)^2 / (n * pd * (1 - pd)))
  df <- length(n)
  list(chi2 = chi2, df = df, p = if (df >= 1L) stats::pchisq(chi2, df, lower.tail = FALSE) else NA_real_)
}

#' Traffic light of a p-value
#' @keywords internal
#' @noRd
.pd_light <- function(p, lights) ifelse(is.na(p), NA_character_, ifelse(p <= lights[1], "red", ifelse(p <= lights[2], "amber", "green")))

#' DeLong standard error of an AUC, from mid-ranks
#'
#' The structural components of DeLong et al. (1988): for every event
#' V10_i is the share of non-events it outranks (ties count one half), for
#' every non-event V01_j the share of events that outrank it;
#' Var(AUC) = var(V10) / n1 + var(V01) / n0. This is the estimator of the
#' ECB instructions for reporting validation results (current AUC against
#' the initial one) and, unlike the Hanley-McNeil approximation, it holds
#' under heavy ties (a grade scale).
#' @keywords internal
#' @noRd
.pd_auc_se <- function(score, y, higher_is_event = TRUE) {
  ok <- !is.na(score) & !is.na(y)
  s <- as.double(score[ok]); y <- as.integer(y[ok])
  if (!higher_is_event) s <- -s
  e <- y == 1L; n1 <- sum(e); n0 <- sum(!e)
  if (n1 < 2L || n0 < 2L) return(NA_real_)
  r <- data.table::frank(s, ties.method = "average")
  v10 <- (r[e] - data.table::frank(s[e], ties.method = "average")) / n0
  v01 <- 1 - (r[!e] - data.table::frank(s[!e], ties.method = "average")) / n1
  sqrt(stats::var(v10) / n1 + stats::var(v01) / n0)
}

#' Validate a PD model on a cohort panel
#'
#' Runs the standard battery on a monthly panel with the default flag and
#' the grade (or the score) at every month: obligors non-defaulted at each
#' cohort start form the population, the outcome is a default within
#' `horizon` months, exactly as [scr_default_rate()] does.
#'
#' \describe{
#'   \item{Calibration}{Per grade (pooled over cohorts) and per cohort and
#'     grade: Jeffreys `p = F_Beta(PD; D + 1/2, N - D + 1/2)`, the binomial
#'     `P(X >= D)` with its critical count at `alpha`, the normal `z`, and
#'     the traffic light on the Jeffreys p-value. Portfolio: the same tests
#'     on the totals, Hosmer-Lemeshow over the grades (`K` degrees of
#'     freedom: the grade PDs are not fitted on the validation sample), the
#'     multi-period normal test over the cohort differences `DR_t - PD_t`
#'     (BCBS Working Paper 14, 2005) and the Brier score.}
#'   \item{Discrimination}{AUC, Gini and KS with a bootstrap interval
#'     ([scr_metrics()]) on the score when a `score` column exists,
#'     otherwise on the grade; the `S` statistic against `auc_init`
#'     (`(AUC_init - AUC_curr) / se`, with the DeLong standard error of
#'     the current AUC), `p = 1 - Phi(S)`.}
#'   \item{Stability}{PSI of the grade distribution against the development
#'     sample per cohort ([scr_psi()]); the migration matrix pooled over
#'     the cohorts whose end date is observed ([scr_migration()]); the
#'     concentration test on the coefficient of variation of the latest
#'     cohort against `cv_init`.}
#' }
#'
#' @param x An [scr_pd()] object.
#' @param newdata A `data.frame`/`data.table` panel, one row per `id` and month.
#' @param id,date,default Column names.
#' @param grade Column name of the grade at every month; `NULL` derives it
#'   from `score` with the cut points of `x`.
#' @param score Column name of the production score at every month, optional.
#' @param auc_init Development AUC; `NULL` uses the scorecard's hold-out AUC.
#' @param cv_init Development coefficient of variation; `NULL` uses the one of `x`.
#' @param tests Subset of the battery to run.
#' @param alpha Significance level of the binomial critical count.
#' @param lights Two p-value thresholds (red at or below the first, amber at
#'   or below the second, green above; the convention shared with the LGD
#'   and EAD validations); `NULL` reads `config$pd_lights`.
#' @param pd_column Grade PD tested: `"pd_final"` (default), `"pd_moc"` or `"pd_be"`.
#' @param horizon,by Cohort window in months and frequency (`NULL` reads `config$pd_dr_by`).
#' @param n_boot,seed Bootstrap resamples and seed of the discrimination interval.
#'
#' @return An object of class `scr_pd_validation`: `calibration` (per
#'   grade, pooled), `calibration_cohort` (per cohort and grade),
#'   `portfolio` (per cohort), `portfolio_tests` (list: `n`, `d`, `dr`,
#'   `pd`, `p_jeffreys`, `p_binomial`, `hl_chi2`, `hl_df`, `hl_p`,
#'   `multi_period_z`, `multi_period_p`, `brier`), `discrimination`,
#'   `stability` (`psi` table, `migration`, `concentration`), `summary`
#'   (one row per test with `statistic`, `p_value`, `light`), `light`
#'   (the worst light of the summary), `n_cohorts`, `alpha`, `lights`.
#'   `portfolio_tests` also carries `critical`, `z`, `p_normal`, `n_cohorts`
#'   and `pd_column`; the object also has `horizon`, `by`, `pd_column` and
#'   `target`.
#'
#' @family irb-pd
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   use_lightgbm = FALSE, xgb_rounds = 40, n_boot = 10)
#' res <- scr_select(scr_demo, "default", config = cfg, drop = c("id", "churn"),
#'                   date_col = "ref_date")
#' sc <- scr_scorecard(res)
#' pd <- scr_pd(scr_moc(scr_grades(sc, n_grades = 6, min_defaults = 10), "C", method = "ci_binomial"))
#' # the validation panel: default flag at every month plus the grade at the
#' # cohort start; here the behavioural score of the panel is graded with the
#' # cut points of the PD model
#' d <- scr_default(scr_demo_panel, "id", "ref_date", dpd = "dpd", config = cfg)
#' pnl <- merge(d$flags, scr_demo_panel[, c("id", "ref_date", "score")],
#'              by.x = c("id", "date"), by.y = c("id", "ref_date"))
#' pnl$grade <- predict(pd, score = pnl$score, type = "grade")
#' v <- scr_pd_validate(pd, pnl, id = "id", date = "date", default = "default",
#'                      grade = "grade", score = "score", by = "quarter")
#' v
#' v$summary
#' @export
scr_pd_validate <- function(x, newdata, id = "id", date = "date", default = "default", grade = NULL, score = NULL,
                            auc_init = NULL, cv_init = NULL,
                            tests = c("jeffreys", "binomial", "normal", "hl", "multi_period", "auc", "concentration", "psi", "migration"),
                            alpha = 0.05, lights = NULL, pd_column = c("pd_final", "pd_moc", "pd_be"), horizon = 12L, by = NULL,
                            n_boot = NULL, seed = NULL) {
  if (!inherits(x, "scr_pd")) stop("scr_pd_validate() expects an object from scr_pd().", call. = FALSE)
  tests <- match.arg(tests, several.ok = TRUE)
  pd_column <- match.arg(pd_column)
  cfg <- x$config
  lights <- lights %||% cfg$pd_lights
  by <- match.arg(by %||% cfg$pd_dr_by, c("month", "quarter", "year"))
  n_boot <- n_boot %||% cfg$n_boot; seed <- seed %||% cfg$seed
  .scr_num1(alpha, "alpha", lower = 0, upper = 1, open_lower = TRUE)
  raw <- data.table::as.data.table(newdata)
  need <- c(id, date, default, grade, score)
  miss <- setdiff(need, names(raw))
  if (length(miss)) stop("scr_pd_validate(): column(s) not found: ", lst(miss), call. = FALSE)
  if (is.null(grade) && is.null(score)) stop("scr_pd_validate(): give `grade` or `score` (a column of the panel).", call. = FALSE)
  dt <- data.table::data.table(id = as.character(raw[[id]]), date = as.Date(raw[[date]]), default = as.integer(raw[[default]]))
  dt[, score := if (is.null(score)) NA_real_ else as.double(raw[[score]])]
  dt[, grade := if (is.null(grade)) .pd_grade_of(x, dt$score) else as.integer(as.character(raw[[grade]]))]
  K <- nrow(x$table)
  if (any(!is.na(dt$grade) & (dt$grade < 1L | dt$grade > K))) stop("scr_pd_validate(): grades must lie in 1..", K, ".", call. = FALSE)
  rows <- .pd_cohorts(dt, as.integer(horizon), by)
  pd_g <- x$table[[pd_column]]

  # -- calibration ---------------------------------------------------------- #
  by_cg <- rows[, list(n = .N, d = sum(.SD[["y"]])), by = c("cohort", "grade"), .SDcols = "y"][order(cohort, grade)]
  by_g <- by_cg[, list(n = sum(.SD[["n"]]), d = sum(.SD[["d"]]), t = .N), by = "grade", .SDcols = c("n", "d")][order(grade)]
  full <- data.table::data.table(grade = seq_len(K))
  by_g <- merge(full, by_g, by = "grade", all.x = TRUE)
  by_g[is.na(n), `:=`(n = 0L, d = 0L, t = 0L)]
  by_g[, `:=`(dr = ifelse(n > 0, d / n, NA_real_), pd = pd_g[grade])]
  ct <- .pd_cal_tests(by_g$n, by_g$d, by_g$pd, alpha)
  by_g[, `:=`(p_jeffreys = ct$p_jeffreys, p_binomial = ct$p_binomial, critical = ct$critical, z = ct$z, p_normal = ct$p_normal)]
  # multi-period per grade
  mp <- by_cg[, list(mp_z = { r <- .SD[["d"]] / .SD[["n"]]; if (.N >= 2L && stats::sd(r) > 0) (mean(r) - pd_g[.BY[["grade"]]]) / (stats::sd(r) / sqrt(.N)) else NA_real_ }),
              by = "grade", .SDcols = c("n", "d")]
  by_g[, multi_period_z := mp$mp_z[match(grade, mp$grade)]]
  by_g[, multi_period_p := 1 - stats::pnorm(multi_period_z)]
  by_g[, light := .pd_light(p_jeffreys, lights)]
  by_cg[, `:=`(dr = d / n, pd = pd_g[grade])]
  ct2 <- .pd_cal_tests(by_cg$n, by_cg$d, by_cg$pd, alpha)
  by_cg[, `:=`(p_jeffreys = ct2$p_jeffreys, p_binomial = ct2$p_binomial, critical = ct2$critical, z = ct2$z, p_normal = ct2$p_normal)]
  by_cg[, light := .pd_light(p_jeffreys, lights)]
  # portfolio
  port <- by_cg[, list(n = sum(.SD[["n"]]), d = sum(.SD[["d"]]), pd = sum(.SD[["n"]] * .SD[["pd"]]) / sum(.SD[["n"]])), by = "cohort", .SDcols = c("n", "d", "pd")][order(cohort)]
  port[, dr := d / n]
  ct3 <- .pd_cal_tests(port$n, port$d, port$pd, alpha)
  port[, `:=`(p_jeffreys = ct3$p_jeffreys, p_binomial = ct3$p_binomial, critical = ct3$critical, z = ct3$z, p_normal = ct3$p_normal)]
  port[, light := .pd_light(p_jeffreys, lights)]
  N <- sum(by_g$n); Dn <- sum(by_g$d); PDp <- sum(by_g$n * by_g$pd) / N
  pt <- .pd_cal_tests(N, Dn, PDp, alpha)
  hl <- .pd_hl(by_g$n, by_g$d, by_g$pd)
  hl_chi2 <- hl$chi2; hl_df <- hl$df; hl_p <- hl$p
  T_c <- nrow(port)
  # normal test of BCBS WP 14: sum_t (DR_t - PD_t) / (sqrt(T) tau), tau the
  # standard deviation of the differences (the cohort PD moves with the mix)
  dd <- port$dr - port$pd
  mp_z <- if (T_c >= 2L && stats::sd(dd) > 0) mean(dd) / (stats::sd(dd) / sqrt(T_c)) else NA_real_
  brier <- sum(by_g$d * (1 - by_g$pd)^2 + (by_g$n - by_g$d) * by_g$pd^2, na.rm = TRUE) / N
  portfolio_tests <- list(n = N, d = Dn, dr = Dn / N, pd = PDp, p_jeffreys = pt$p_jeffreys, p_binomial = pt$p_binomial, critical = pt$critical,
                          z = pt$z, p_normal = pt$p_normal, hl_chi2 = hl_chi2, hl_df = hl_df, hl_p = hl_p,
                          multi_period_z = mp_z, multi_period_p = if (is.na(mp_z)) NA_real_ else 1 - stats::pnorm(mp_z), brier = brier,
                          n_cohorts = T_c, pd_column = pd_column)

  # -- discrimination ------------------------------------------------------- #
  disc <- NULL
  if ("auc" %in% tests) {
    use_score <- !all(is.na(rows$score))
    sv <- if (use_score) rows$score else as.double(rows$grade)
    hie <- if (use_score) identical(x$direction, "higher_is_riskier") else TRUE
    m <- scr_metrics(sv, rows$y, higher_is_event = hie, ci = TRUE, n_boot = n_boot, level = cfg$ci_level, seed = seed, nthread = cfg$nthread)
    auc_init <- auc_init %||% x$scorecard$metrics[x$scorecard$metrics$sample == "holdout", ][["auc"]][1]
    se <- if (is.na(m$auc)) NA_real_ else .pd_auc_se(sv, rows$y, higher_is_event = hie)
    S <- if (is.na(se) || se == 0) NA_real_ else (auc_init - m$auc) / se
    disc <- data.table::data.table(basis = if (use_score) "score" else "grade", n = m$n, events = m$events, auc = m$auc, auc_lo = m$auc_lo, auc_hi = m$auc_hi,
                                   gini = m$gini, gini_lo = m$gini_lo, gini_hi = m$gini_hi, ks = m$ks, ks_lo = m$ks_lo, ks_hi = m$ks_hi,
                                   auc_init = auc_init, se_auc = se, s_stat = S, p_value = if (is.na(S)) NA_real_ else 1 - stats::pnorm(S))
  }

  # -- stability ------------------------------------------------------------ #
  psi_tab <- NULL; mig <- NULL; conc <- NULL
  dev_grade <- x$rows$grade
  cohorts <- sort(unique(rows$cohort))
  if ("psi" %in% tests) {
    # one pass: the development grades converted once, the panel split once
    dev_chr <- as.character(dev_grade); lv_k <- as.character(seq_len(K))
    by_coh <- split(rows$grade, rows$cohort)
    psi_tab <- data.table::rbindlist(lapply(cohorts, function(c0) {
      g <- by_coh[[as.character(c0)]]
      r <- scr_psi(dev_chr, as.character(g), levels = lv_k, alpha = alpha)
      data.table::data.table(cohort = c0, n = length(g), psi = r$psi, flag_fixed = r$flag_fixed, critical = r$critical, flag_adjusted = r$flag_adjusted)
    }))
  }
  if ("migration" %in% tests) {
    obs <- rows[!is.na(rows$grade_t1) | rows$y == 1L]
    mig <- if (nrow(obs)) scr_migration(obs$grade, obs$grade_t1, K = K) else NULL
  }
  if ("concentration" %in% tests) {
    last <- rows[rows$cohort == max(cohorts)]
    cc <- .pd_concentration(tabulate(last$grade, nbins = K))
    cv0 <- cv_init %||% x$concentration$cv
    p_cv <- if (cc$cv > 0 && K > 1L) 1 - stats::pnorm(sqrt(K - 1) * (cc$cv - cv0) / sqrt(cc$cv^2 * (0.5 + cc$cv^2))) else NA_real_
    conc <- list(cohort = max(cohorts), cv = cc$cv, cv_init = cv0, hi = cc$hi, hhi = cc$hhi, p_value = p_cv,
                 shares = tabulate(last$grade, nbins = K) / nrow(last))
  }

  # -- summary --------------------------------------------------------------- #
  row <- function(test, level, stat, p) data.table::data.table(test = test, level = level, statistic = stat, p_value = p, light = .pd_light(p, lights))
  sm <- list()
  if ("jeffreys" %in% tests) { sm[[length(sm) + 1L]] <- row("jeffreys", "portfolio", portfolio_tests$dr, portfolio_tests$p_jeffreys)
    sm[[length(sm) + 1L]] <- row("jeffreys_grades_red", "grade", sum(by_g$light == "red", na.rm = TRUE),
                               if (all(is.na(by_g$p_jeffreys))) NA_real_ else min(by_g$p_jeffreys, na.rm = TRUE)) }
  if ("binomial" %in% tests) sm[[length(sm) + 1L]] <- row("binomial", "portfolio", portfolio_tests$critical, portfolio_tests$p_binomial)
  if ("normal" %in% tests) sm[[length(sm) + 1L]] <- row("normal", "portfolio", portfolio_tests$z, portfolio_tests$p_normal)
  if ("hl" %in% tests) {
    # two-sided by construction: a PD above the observed rate in every grade
    # (the conservative side, usual once the margin of conservatism is in)
    # cannot fail the calibration test
    hl_row <- row("hosmer_lemeshow", "portfolio", hl_chi2, hl_p)
    if (all(by_g$d <= by_g$n * by_g$pd + 1e-9, na.rm = TRUE) && !is.na(hl_p)) hl_row$light <- "green"
    sm[[length(sm) + 1L]] <- hl_row
  }
  if ("multi_period" %in% tests) sm[[length(sm) + 1L]] <- row("multi_period", "portfolio", mp_z, portfolio_tests$multi_period_p)
  if (!is.null(disc)) sm[[length(sm) + 1L]] <- row("auc_vs_initial", "portfolio", disc$s_stat, disc$p_value)
  if (!is.null(psi_tab)) { last_psi <- psi_tab[nrow(psi_tab)]
    sm[[length(sm) + 1L]] <- data.table::data.table(test = "psi_grades", level = "portfolio", statistic = last_psi$psi, p_value = NA_real_,
                                                    light = if (is.na(last_psi$flag_fixed)) NA_character_ else
                                                      switch(last_psi$flag_fixed, stable = "green", moderate = "amber", "red")) }
  if (!is.null(mig)) sm[[length(sm) + 1L]] <- data.table::data.table(test = "migration_mwb_upper", level = "portfolio", statistic = mig$mwb_upper, p_value = NA_real_, light = NA_character_)
  if (!is.null(conc)) sm[[length(sm) + 1L]] <- row("concentration_cv", "portfolio", conc$cv, conc$p_value)
  summary <- data.table::rbindlist(sm)
  worst <- if (any(summary$light == "red", na.rm = TRUE)) "red" else if (any(summary$light == "amber", na.rm = TRUE)) "amber" else "green"
  structure(list(calibration = by_g[], calibration_cohort = by_cg[], portfolio = port[], portfolio_tests = portfolio_tests,
                 discrimination = disc, stability = list(psi = psi_tab, migration = mig, concentration = conc),
                 summary = summary[], light = worst, n_cohorts = T_c, alpha = alpha, lights = lights, horizon = as.integer(horizon), by = by,
                 pd_column = pd_column, target = x$target),
            class = c("scr_pd_validation", "list"))
}

#' @export
print.scr_pd_validation <- function(x, ...) {
  pt <- x$portfolio_tests
  cat(sprintf("<scr_pd_validation> target \"%s\" | %d %sly cohorts, %d-month window | overall light: %s\n", x$target, x$n_cohorts, x$by, x$horizon, toupper(x$light)))
  cat(sprintf("  portfolio: N %s | D %s | DR %s vs %s %s | Jeffreys p %s | binomial p %s (critical %s) | HL chi2 %s (p %s) | multi-period z %s\n",
              n_fmt(pt$n), n_fmt(pt$d), fmt_pct(pt$dr, 2), x$pd_column, fmt_pct(pt$pd, 2), .fmt_num(pt$p_jeffreys), .fmt_num(pt$p_binomial),
              if (is.na(pt$critical)) "-" else format(pt$critical), .fmt_num(pt$hl_chi2, 2), .fmt_num(pt$hl_p), .fmt_num(pt$multi_period_z, 2)))
  g <- x$calibration
  cat(sprintf("  %-5s %7s %5s %8s %8s %9s %9s %-6s\n", "grade", "n", "d", "dr", "pd", "p_jeff", "p_binom", "light"))
  for (i in seq_len(nrow(g))) cat(sprintf("  %-5d %7d %5d %8s %8s %9s %9s %-6s\n", g$grade[i], g$n[i], g$d[i], fmt_pct(g$dr[i], 2), fmt_pct(g$pd[i], 2),
                                          .fmt_num(g$p_jeffreys[i]), .fmt_num(g$p_binomial[i]), g$light[i] %||% "-"))
  d <- x$discrimination
  if (!is.null(d)) cat(sprintf("  discrimination (%s): AUC %.4f [%.4f, %.4f] vs initial %.4f | S %s, p %s | KS %.4f\n", d$basis, d$auc, d$auc_lo, d$auc_hi,
                               d$auc_init, .fmt_num(d$s_stat, 2), .fmt_num(d$p_value), d$ks))
  s <- x$stability
  if (!is.null(s$psi)) { l <- s$psi[nrow(s$psi)]; cat(sprintf("  stability: grade PSI %.4f (%s, adjusted %s) at cohort %s", l$psi, l$flag_fixed, l$flag_adjusted, format(l$cohort))) }
  if (!is.null(s$migration)) cat(sprintf(" | MWB up %s / down %s", .fmt_num(s$migration$mwb_upper), .fmt_num(s$migration$mwb_lower)))
  if (!is.null(s$concentration)) cat(sprintf(" | CV %.3f vs %.3f (p %s)", s$concentration$cv, s$concentration$cv_init, .fmt_num(s$concentration$p_value)))
  cat("\n")
  invisible(x)
}

# -- export ------------------------------------------------------------------ #

#' @rdname scr_export
#' @export
scr_export.scr_pd <- function(x, dir, stamp = TRUE, validation = NULL, ...) {
  .scr_local_verbose(x)
  .need_openxlsx()
  out_dir <- .export_dir(dir, stamp)
  tag <- .file_tag(x$target)
  na_v <- data.frame(availability = "not_available", reason_code = "NO_VALIDATION_SUPPLIED", stringsAsFactors = FALSE)
  cal <- x$calibration
  cal_kv <- if (is.null(cal)) .kv_table(list(method = "none", note = "PD taken from the scorecard alignment")) else
    .kv_table(c(cal[c("method", "ct", "target_source", "sample_rate", "shift", "shift_prior", "slope_ratio", "mean_pd_before", "mean_pd_after",
                      "ar_before", "ar_after", "ar_implied_before", "ar_implied_after", "n")],
                list(intercept_before = cal$alignment_before$calibration$intercept, slope_before = cal$alignment_before$calibration$slope,
                     intercept_after = cal$alignment$calibration$intercept, slope_after = cal$alignment$calibration$slope,
                     a_after = cal$alignment$a, b_after = cal$alignment$b)))
  t <- x$table
  floors <- data.frame(framework = x$framework, asset_class = x$asset_class, floor = x$floor, params_modified = x$params_modified,
                       n_grades_floored = sum(t$floor_applied), share_obligors_at_floor = x$portfolio$share_at_floor, stringsAsFactors = FALSE)
  mig <- if (!is.null(validation$stability$migration)) {
    M <- validation$stability$migration$p
    cbind(data.frame(from = rownames(M), n = validation$stability$migration$n, stringsAsFactors = FALSE), as.data.frame(M, stringsAsFactors = FALSE))
  } else NULL
  stab <- if (is.null(validation)) NULL else {
    ps <- validation$stability$psi
    cc <- validation$stability$concentration
    rbind(if (!is.null(ps)) data.frame(item = paste0("psi_grades_", format(ps$cohort)), value = as.character(round(ps$psi, 6)), stringsAsFactors = FALSE),
          if (!is.null(cc)) .kv_table(cc[c("cohort", "cv", "cv_init", "hi", "hhi", "p_value")]),
          if (!is.null(validation$stability$migration)) .kv_table(c(validation$stability$migration[c("mwb_upper", "mwb_lower", "n_significant")],
                                                                    validation$stability$migration$mobility)))
  }
  sheets <- list(
    "PD_Grades"       = t,
    "Master_Scale"    = if (is.null(x$master_scale)) NULL else as.data.frame(x$master_scale),
    "Calibration"     = cal_kv,
    "MoC_Ledger"      = x$moc_ledger,
    "Floors"          = floors,
    "Grade_Series"    = x$dr,
    "Validation_Calibration"    = if (is.null(validation)) na_v else validation$calibration,
    "Validation_Cohorts"        = if (is.null(validation)) na_v else validation$calibration_cohort,
    "Validation_Discrimination" = if (is.null(validation) || is.null(validation$discrimination)) na_v else validation$discrimination,
    "Validation_Stability"      = stab %||% na_v,
    "Validation_Summary"        = if (is.null(validation)) na_v else validation$summary,
    "Migration"       = mig %||% na_v,
    "Decision_Ledger" = x$ledger,
    "Model_Card"      = .kv_table(x$model_card)
  )
  files <- list(pd = .scr_write_xlsx(sheets, file.path(out_dir, sprintf("pd_%s.xlsx", tag))),
                sql_pd = file.path(out_dir, sprintf("sql_pd_%s.sql", tag)))
  writeLines(scr_sql(x), files$sql_pd)
  for (f in files) msg("  %s", f)
  x$files <- files
  invisible(x)
}

# data.table column names used without quotes in this file
utils::globalVariables(c(
  "grade", "grade0", "label0", "band", "score_lo", "score_hi", "pd_lo", "pd_hi", "pd_be", "pd_mean", "n_series", "t_series",
  "active", "moc_a", "moc_b", "moc_c", "pd_moc", "pd_ttc", "pd_pit", "floor", "pd_final", "floor_applied",
  "cohort", "d", "dr", "pd", "p_jeffreys", "p_binomial", "critical", "z", "p_normal", "multi_period_z", "multi_period_p",
  "light", "grade_t1", ".BY", "g1", "y", "id", "date"
))
