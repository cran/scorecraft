# ============================================================================ #
# apply.R - score new data in R (the mirror of the production SQL)
# ============================================================================ #

#' Apply the WOE transformation or the scorecard to new data
#'
#' Materialises in R exactly what the production SQL does: the frozen Stage
#' 1 pre-processing (training median, special-population flags,
#' `"MISSING"`) followed by the frozen Stage 2 binning and, for a scorecard,
#' by the points. Nothing is refitted. The two paths, R and SQL, produce the
#' same numbers, and a test guarantees it.
#'
#' @param x An object from [scr_select()] (returns WOE/bin of the approved
#'   variables) or from [scr_scorecard()] (returns score and points).
#' @param newdata New table with the source columns of the requested
#'   variables. The target column is not needed.
#' @param ... Passed on to the methods.
#' @param features For `scr_result`: which variables to transform. Defaults to
#'   the approved ones.
#' @param what For `scr_result`: `"woe"`, `"bin"` or `"both"`. For
#'   `scr_scorecard`: `"score"`, `"points"`, `"woe"` or `"all"`. For
#'   `scr_lgd`: `"pool"` (pool and pool LGDs), `"lgd"` (adds the predicted
#'   LGD) or `"all"` (adds the cure probability and the severity). For
#'   `scr_ead`: `"all"` (default), `"ead"` (pool, measure, applied CCF,
#'   predicted EAD and the floor flag) or `"pool"` (pool and measure).
#'
#' @section Output columns:
#'
#' For `scr_result`: `<f>_woe` and/or `<f>_bin` per variable. For
#' `scr_scorecard`, `"score"` gives `link` (logit), `prob` (model
#' probability), `score` (exact, `a + b * logit`) and `score_points` (base
#' plus the whole points per bin); `"points"` gives `score`, `score_points`
#' and `<f>_points`; `"woe"` gives `link`, `score` and `<f>_woe`; `"all"`
#' gives everything.
#'
#' @section IRB models:
#'
#' `scr_pd` returns `score`, `score_points`, `grade`, `pd` (calibrated
#' individual PD), `pd_be` and `pd_final` of the grade. `scr_lgd` returns
#' `pool`, `lgd_lra`, `lgd_dt`, `lgd_final` and, with `what`, `p_cure`,
#' `severity` and `lgd_pred`. `scr_ead` returns `pool`, `measure`,
#' `utilisation`, `undrawn`, `ccf_applied`, `ead_model`, `ead_floor`,
#' `ead_predicted` and `ead_floor_binding`; the predicted EAD is never below
#' the drawn amount. [scr_capital()] reads `pd_final`, `lgd_final` and
#' `ead_predicted` from these outputs in its list form.
#'
#' @return A `data.table` with one row per row of `newdata`.
#'
#' @family production
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   xgb_rounds = 60, n_boot = 20)
#' res <- scr_select(scr_demo, "default", config = cfg, drop = "id",
#'                   date_col = "ref_date")
#' new <- head(scr_demo, 50)
#' str(scr_apply(res, new)[, 1:3])
#' sc <- scr_scorecard(res)
#' head(scr_apply(sc, new))
#' head(scr_apply(sc, new, what = "points"))
#' @export
scr_apply <- function(x, newdata, ...) UseMethod("scr_apply")

#' @rdname scr_apply
#' @export
scr_apply.scr_result <- function(x, newdata, features = scr_selected(x), what = c("woe", "bin", "both"), ...) {
  what <- match.arg(what)
  if (!length(features)) stop("no variable to apply.", call. = FALSE)
  missing_fit <- setdiff(features, names(x$fit$results))
  if (length(missing_fit)) stop("variable(s) without a fitted binning: ", lst(missing_fit), call. = FALSE)
  base <- .scr_preprocess(newdata, x$triage$ledger, features, x$config$special_values)
  apply_woe(x$fit, base, features, what)[]
}

#' @rdname scr_apply
#' @export
scr_apply.scr_scorecard <- function(x, newdata, what = c("score", "points", "woe", "all"), ...) {
  what <- match.arg(what)
  feats <- x$features
  base <- .scr_preprocess(newdata, x$ledger, feats, x$config$special_values)
  w <- apply_woe(x$fit, base, feats, "both")
  link <- .glm_link(x$coef, w, feats)
  al <- x$alignment
  out <- data.table::data.table(link = link, prob = stats::plogis(link), score = al$a + al$b * link)
  sp <- rep(x$base_points, nrow(w))
  pts_cols <- list()
  # a row in no fitted bin has WOE 0: its points are those of WOE 0, as in the SQL
  unb <- .sc_unbinned_points(x)
  for (f in feats) {
    p <- x$points[variable == f]
    pf <- p$points[match(w[[paste0(f, "_bin")]], p$bin)]
    pf[is.na(pf)] <- unb
    sp <- sp + pf
    pts_cols[[paste0(f, "_points")]] <- pf
  }
  out[, score_points := sp]
  # room for every <f>_points and <f>_woe column: set() cannot grow a table
  # past its allocated column slots (1024 by default)
  data.table::setalloccol(out, ncol(out) + 2L * length(feats) + getOption("datatable.alloccol", 1024L))
  if (what %in% c("points", "all")) for (nm in names(pts_cols)) data.table::set(out, j = nm, value = pts_cols[[nm]])
  if (what %in% c("woe", "all")) for (f in feats) data.table::set(out, j = paste0(f, "_woe"), value = w[[paste0(f, "_woe")]])
  if (identical(what, "points")) out <- out[, c("score", "score_points", names(pts_cols)), with = FALSE]
  if (identical(what, "woe")) out <- out[, c("link", "score", paste0(feats, "_woe")), with = FALSE]
  out[]
}

#' Reproduce in R the ledger pre-processing (the same as the SQL CTE)
#' @keywords internal
#' @noRd
.scr_preprocess <- function(newdata, ledger, features, sp) {
  dt  <- data.table::as.data.table(newdata)
  origin <- vapply(features, function(f) {
    s <- ledger[kind == "num_flag" & output == f, source]
    if (length(s)) s[1] else f
  }, character(1))
  missing <- setdiff(unique(origin), names(dt))
  if (length(missing)) stop("newdata lacks the source column(s): ", lst(missing), call. = FALSE)
  base <- data.table::data.table(.i = seq_len(nrow(dt)))
  data.table::setalloccol(base, length(features) + 1L + getOption("datatable.alloccol", 1024L))
  for (f in features) {
    src <- origin[[f]]
    if (identical(src, f)) {
      v <- dt[[f]]
      if (is.numeric(v)) {
        imp <- ledger[kind == "num_impute" & source == f, impute_value]
        if (length(imp) && is.finite(imp[1])) {
          bad <- is.na(v) | v %in% sp
          if (any(bad)) v[bad] <- imp[1]
        }
        data.table::set(base, j = f, value = as.double(v))
      } else {
        v <- as.character(v)
        if (nrow(ledger[kind == "cat_coalesce" & source == f]) && anyNA(v)) v[is.na(v)] <- "MISSING"
        data.table::set(base, j = f, value = v)
      }
    } else {
      data.table::set(base, j = f, value = .flag_levels(dt[[src]], sp))
    }
  }
  base[, ".i" := NULL]
  base[]
}

#' Reason codes: the variables that took the most points from each row
#'
#' For each row of `newdata`, the `k` variables whose contribution in points
#' fell furthest below the reference. The reference is the mean points of
#' the variable on the training population (`"mean"`, the Regulation B safe
#' harbour referenced to the average) or the maximum points of the variable
#' (`"max"`). Only applies to the additive scorecard; a tree challenger has
#' no reason codes.
#'
#' Under `higher_is_riskier` the shortfall is measured the other way round:
#' the reasons are the variables that added the most points.
#'
#' @param x An object from [scr_scorecard()].
#' @param newdata New table.
#' @param k Number of reasons per row.
#' @param reference `"mean"` (default) or `"max"`.
#'
#' @return A `data.table` with `reason_1` ... `reason_k` (variable names)
#'   and `shortfall_1` ... `shortfall_k` (points below the reference).
#'
#' @references
#' 12 CFR 1002.9 (Regulation B), official commentary to paragraph 9(b)(2).
#'
#' @family production
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   xgb_rounds = 60, n_boot = 20)
#' res <- scr_select(scr_demo, "default", config = cfg, drop = "id",
#'                   date_col = "ref_date")
#' sc <- scr_scorecard(res)
#' scr_reasons(sc, head(scr_demo, 5), k = 3)
#' @export
scr_reasons <- function(x, newdata, k = 4L, reference = c("mean", "max")) {
  if (!inherits(x, "scr_scorecard")) stop("scr_reasons() expects an object from scr_scorecard().", call. = FALSE)
  reference <- match.arg(reference)
  pts <- scr_apply(x, newdata, what = "points")
  feats <- x$features
  ref <- vapply(feats, function(f) {
    p <- x$points[variable == f]
    if (identical(reference, "max")) max(p$points) else sum(p$points * p$pct_train)
  }, numeric(1))
  # under higher_is_riskier more points is worse: the "reason" is the variable that ADDED most
  sgn <- if (identical(x$direction, "higher_is_safer")) 1 else -1
  M <- as.matrix(pts[, paste0(feats, "_points"), with = FALSE])
  short <- sweep(-M, 2L, -ref) * sgn   # (ref - points) * sgn
  k <- min(as.integer(k), length(feats))
  n <- nrow(M); nf <- length(feats)
  out <- data.table::data.table(.i = seq_len(n))
  data.table::setalloccol(out, 2L * k + 1L + getOption("datatable.alloccol", 1024L))
  if (k >= 1L && n > 0L) {
    # one global sort of the long (row, variable) table instead of one order()
    # and 2k set() calls per row: ties keep the order of the features, as
    # order() did row by row
    long <- data.table::data.table(i = rep(seq_len(n), times = nf), j = rep(seq_len(nf), each = n),
                                   s = as.vector(short))
    data.table::setorderv(long, c("i", "s", "j"), order = c(1L, -1L, 1L), na.last = TRUE)
    rk <- data.table::rowid(long$i)
    for (h in seq_len(k)) {
      at <- which(rk == h)   # sorted by row, so aligned with `out`
      data.table::set(out, j = paste0("reason_", h), value = feats[long$j[at]])
      data.table::set(out, j = paste0("shortfall_", h), value = long$s[at])
    }
  } else {
    for (j in seq_len(k)) { out[[paste0("reason_", j)]] <- NA_character_; out[[paste0("shortfall_", j)]] <- NA_real_ }
  }
  out[, .i := NULL]
  out[]
}
