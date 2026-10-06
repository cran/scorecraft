# ============================================================================ #
# align.R - Stage 5: alignment of a raw score to the declared scale
# ============================================================================ #
# Without a declared, auditable alignment, scorecards from distinct
# models/targets are not comparable: "700 points" only means the same thing
# if every raw score was taken to the SAME scale by the SAME procedure.
# Mechanism (the classic two-step alignment used in practitioner training):
#
#   ln(odds) = I + S * raw        (calibration regression by band)
#   score    = offset + factor * ln(odds) = a + b * raw
# ============================================================================ #

#' Stage 5: align a raw score to the declared scale
#'
#' Takes the raw score of **any** engine (the scorecard logit, the output of
#' a tree challenger, a legacy score) to the scale defined by `base_score`,
#' `base_odds` and `pdo`, recording `odds_orientation` on the object. This is
#' what makes two scorecards directly comparable. It runs automatically
#' inside [scr_scorecard()]; it is exposed to align other
#' scores to the same scale.
#'
#' @section Mechanism:
#'
#' With `method = "regression"` (default): the raw score is banded by
#' quantiles on the reference data, the empirical log-odds of every band is
#' computed with Laplace smoothing in the orientation `direction` implies,
#' and a regression weighted by band size fits
#' \deqn{\ln(\mathrm{odds}) = I + S \cdot \mathrm{raw}.}
#' This absorbs sample reweighting, miscalibration of the WOE fit and prior
#' shift. It then composes with the PDO map:
#' \deqn{\mathrm{factor} = \mathrm{pdo}/\ln 2,\quad
#'       \mathrm{offset} = \mathrm{base\_score} - \mathrm{factor}\cdot\ln(\mathrm{base\_odds}),}
#' \deqn{\mathrm{score} = \mathrm{offset} + \mathrm{factor}\,(I + S\cdot\mathrm{raw}) = a + b\cdot\mathrm{raw}.}
#'
#' With `method = "direct"`, the model is assumed calibrated: `I = 0` and
#' `S` is the sign of the direction (`-1` under `higher_is_safer`, `+1`
#' under `higher_is_riskier`), that is, `ln(odds)` is `raw` itself in the
#' right orientation.
#'
#' @section Odds orientation:
#'
#' `base_odds` is always expressed in the orientation `direction` implies:
#' non-event:event (`"safe:event"`) under `higher_is_safer`, event:non-event
#' (`"event:safe"`) under `higher_is_riskier`. The same word, `odds`,
#' changes meaning between the two, and that is the most common sign trap in
#' the literature; hence the object records `odds_orientation` explicitly.
#'
#' @param raw Raw score: an event logit (or any score on which a higher
#'   value means a higher probability of the event).
#' @param y 0/1 outcome vector (numeric or logical), same length as `raw`.
#' @param base_score Score at which the odds are `base_odds`.
#' @param base_odds Odds at `base_score`, positive, in the orientation of `direction`.
#' @param pdo Points that double the odds, positive.
#' @param direction `"higher_is_safer"` or `"higher_is_riskier"`.
#' @param method `"regression"` (default) or `"direct"`.
#' @param n_bands Bands of the calibration regression.
#' @param laplace Smoothing of the counts per band.
#' @param weights Optional non-negative weights per observation (sample
#'   reweighting), of the length of `raw`.
#'
#' @return An `scr_align` object with `base_score`, `base_odds`, `pdo`,
#'   `direction`, `odds_orientation`, `factor`, `offset`, `sign`,
#'   `calibration` (`method`, `intercept`, `slope`, `r2`, `n_bands`, `bands`)
#'   and the final coefficients `a` and `b` of `score = a + b * raw`. Use
#'   [predict.scr_align()] to apply it.
#'
#' @references
#' Siddiqi, N. (2006). *Credit Risk Scorecards*. Wiley, chapter 6.
#'
#' @family stages
#' @examples
#' set.seed(3)
#' y   <- stats::rbinom(4000, 1, 0.15)
#' raw <- stats::qlogis(0.15) + 1.3 * y + stats::rnorm(4000, sd = 1.2)
#' al  <- scr_align(raw, y, base_score = 600, base_odds = 50, pdo = 20)
#' al
#' head(predict(al, raw))
#' head(predict(al, raw, type = "prob"))
#'
#' # propensity: more points = more event, odds event:non-event
#' scr_align(raw, y, base_score = 500, base_odds = 1/9, pdo = 40,
#'           direction = "higher_is_riskier")
#' @export
scr_align <- function(raw, y, base_score = 600, base_odds = 50, pdo = 20,
                      direction = c("higher_is_safer", "higher_is_riskier"),
                      method = c("regression", "direct"), n_bands = 10L,
                      laplace = 0.5, weights = NULL) {
  direction <- match.arg(direction); method <- match.arg(method)
  .scr_num1(base_score, "base_score"); .scr_num1(base_odds, "base_odds", lower = 0, open_lower = TRUE)
  .scr_num1(pdo, "pdo", lower = 0, open_lower = TRUE)
  if (length(raw) != length(y)) stop("`raw` and `y` must have the same length.", call. = FALSE)
  if (!is.null(weights) && (!is.numeric(weights) || length(weights) != length(raw))) {
    stop("`weights` must be a numeric vector of the length of `raw`.", call. = FALSE)
  }
  y <- .scr_y01(y, "scr_align")
  ok <- is.finite(raw) & !is.na(y)
  if (!is.null(weights)) {
    if (any(weights[ok] < 0, na.rm = TRUE)) stop("`weights` must be non-negative.", call. = FALSE)
    ok <- ok & is.finite(weights)
  }
  if (!any(ok)) stop("scr_align(): no row with a finite `raw` and a non-missing `y`.", call. = FALSE)
  raw <- as.double(raw[ok]); y <- y[ok]
  w <- if (is.null(weights)) rep(1, length(raw)) else as.double(weights[ok])
  sgn <- if (identical(direction, "higher_is_safer")) -1 else 1

  cal <- list(method = method, intercept = 0, slope = sgn, r2 = NA_real_,
              n_bands = 0L, bands = NULL, note = "")
  if (identical(method, "regression")) {
    cal <- .align_regression(raw, y, w, sgn, n_bands, laplace)
    if (isTRUE(cal$fallback)) {
      warning("scr_align(): not enough bands for the calibration regression (", cal$note,
              "); using method = \"direct\".", call. = FALSE)
    } else if (sign(cal$slope) != sgn) {
      warning("scr_align(): the calibration slope (", sprintf("%.3f", cal$slope),
              ") has the opposite sign to the one expected for ", direction,
              " - the raw score ranks against the declared direction.", call. = FALSE)
    }
  }
  .scr_align_from(cal$intercept, cal$slope, base_score, base_odds, pdo, direction, cal)
}

#' Calibration regression by band: ln(odds) ~ raw
#' @keywords internal
#' @noRd
.align_regression <- function(raw, y, w, sgn, n_bands, laplace) {
  probs  <- seq(0, 1, length.out = n_bands + 1L)[-c(1L, n_bands + 1L)]
  breaks <- unique(c(-Inf, stats::quantile(raw, probs = probs, names = FALSE), Inf))
  band   <- cut(raw, breaks = breaks, include.lowest = TRUE)
  d <- data.table::data.table(band = band, raw = raw, y = y, w = w)[
    , .(n = sum(w), events = sum(w * y), raw_mean = stats::weighted.mean(raw, w)), by = band][order(band)]
  d[, ln_odds := if (sgn < 0) log((n - events + laplace) / (events + laplace))
                 else log((events + laplace) / (n - events + laplace))]
  use <- d[n > 0 & is.finite(ln_odds)]
  if (nrow(use) < 3L || stats::sd(use$raw_mean) == 0) {
    return(list(method = "direct", intercept = 0, slope = sgn, r2 = NA_real_,
                n_bands = nrow(use), bands = as.data.frame(d), fallback = TRUE,
                note = sprintf("%d usable band(s)", nrow(use))))
  }
  fit <- stats::lm(ln_odds ~ raw_mean, data = use, weights = n)
  cf  <- stats::coef(fit)
  d[, ln_odds_fit := cf[1] + cf[2] * raw_mean]
  list(method = "regression", intercept = unname(cf[1]), slope = unname(cf[2]),
       r2 = summary(fit)$adj.r.squared, n_bands = nrow(use), bands = as.data.frame(d),
       fallback = FALSE, note = "")
}

#' Build the alignment object from (I, S) and the scale
#' @keywords internal
#' @noRd
.scr_align_from <- function(intercept, slope, base_score, base_odds, pdo, direction, calibration = NULL) {
  factor <- pdo / log(2)
  offset <- base_score - factor * log(base_odds)
  structure(list(
    base_score = base_score, base_odds = base_odds, pdo = pdo, direction = direction,
    odds_orientation = if (identical(direction, "higher_is_safer")) "safe:event" else "event:safe",
    factor = factor, offset = offset,
    sign = if (identical(direction, "higher_is_safer")) -1 else 1,
    calibration = calibration %||% list(method = "direct", intercept = intercept, slope = slope,
                                        r2 = NA_real_, n_bands = 0L, bands = NULL, note = ""),
    a = offset + factor * intercept, b = factor * slope
  ), class = c("scr_align", "list"))
}

#' Apply an alignment to raw scores
#'
#' @param object An object from [scr_align()].
#' @param raw Raw scores on the same scale used in the fit.
#' @param type `"score"` (default) returns points; `"prob"` returns the
#'   calibrated event probability implied by the alignment.
#' @param ... Ignored.
#'
#' @return A numeric vector of the length of `raw`.
#'
#' @family production
#' @examples
#' set.seed(3)
#' y   <- stats::rbinom(2000, 1, 0.15)
#' raw <- stats::qlogis(0.15) + 1.3 * y + stats::rnorm(2000)
#' al  <- scr_align(raw, y)
#' head(predict(al, raw))
#' head(predict(al, raw, type = "prob"))
#' @export
predict.scr_align <- function(object, raw, type = c("score", "prob"), ...) {
  type <- match.arg(type)
  score <- object$a + object$b * as.double(raw)
  if (identical(type, "score")) return(score)
  .score_to_prob(object, score)
}

#' Event probability implied by an aligned score
#' @keywords internal
#' @noRd
.score_to_prob <- function(al, score) {
  ln_odds <- (score - al$offset) / al$factor
  if (al$sign < 0) 1 / (1 + exp(ln_odds)) else stats::plogis(ln_odds)
}

#' @export
print.scr_align <- function(x, ...) {
  cat(sprintf("<scr_align> %s points at odds %s:1 (%s), PDO %s | %s\n",
              format(x$base_score), format(x$base_odds), x$odds_orientation, format(x$pdo), x$direction))
  cat(sprintf("  factor = %.6f | offset = %.6f\n", x$factor, x$offset))
  cl <- x$calibration
  if (identical(cl$method, "regression")) {
    cat(sprintf("  calibration: ln(odds) = %.6f + %.6f * raw  (adj. R2 = %.4f, %d bands)\n",
                cl$intercept, cl$slope, cl$r2, cl$n_bands))
  } else {
    cat(sprintf("  calibration: direct (model assumed calibrated; I = %.3f, S = %+.0f)\n", cl$intercept, cl$slope))
  }
  cat(sprintf("  score = %.6f + %.6f * raw\n", x$a, x$b))
  invisible(x)
}
