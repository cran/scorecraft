# ============================================================================ #
# irb-binning.R - binning of drivers against a continuous bounded target
# ============================================================================ #
# The optimal-binning engine the scorecard uses accepts binary targets only.
# LGD and CCF are continuous and bounded, so this file provides a small
# supervised binner whose result has EXACTLY the shape of an `obwoe` object:
# `bin`, `cutpoints`, `woe`, `iv`, counts, `type`. The bin statistic of the
# target (its mean, or the logit of the mean) travels in the `woe` slot, so
# OptimalBinningWoE::obwoe_apply() and obwoe_sql() reproduce it in R and in
# every SQL dialect without a line of new code.
#
# Algorithm, numerics: candidate cut points at quantiles, greedy merge of
# the adjacent pair whose merge costs the least between-bin sum of squares,
# until the target number of bins, the minimum share and the minimum count
# hold; then pool-adjacent-violators on the bin means when monotonicity is
# required. Categoricals: levels ordered by mean, rare levels attached to
# the nearest neighbour, then the same merge on the ordered sequence.
# ============================================================================ #

#' Bin drivers against a continuous target (LGD, CCF)
#'
#' Supervised binning for a bounded continuous target, with the result in
#' the shape of an `obwoe` object, so that the [scr_apply()] and [scr_sql()]
#' machinery (`OptimalBinningWoE::obwoe_apply()` and
#' `OptimalBinningWoE::obwoe_sql()`) reproduces the bin statistic unchanged.
#' The `woe` slot of every bin carries the target mean of the bin (or its
#' logit with `scale = "logit"`); `iv` carries the bin's share of the
#' between-bin sum of squares, so `total_iv` is the eta-squared of the
#' driver, in `[0, 1]`.
#'
#' Numeric drivers must not contain missing values: run [scr_triage()] (or
#' impute) first, exactly as the scorecard pipeline does. Categorical
#' missing values become the level `"NA"`, as in the engine. When a
#' `holdout_idx` is given, the frozen bins are revalidated: the hold-out
#' bin means are recomputed, the PSI of the bin shares is reported with the
#' sample-size-adjusted critical value, a driver whose hold-out means
#' break the training order is flagged `UNSTABLE_HOLDOUT`, one whose bin
#' shares shift (fixed PSI flag `"shift"`, PSI at or above 0.25) is flagged
#' `PSI_ACTION`, and one with more than 1% of hold-out rows outside the
#' bins `UNBINNED_HOLDOUT`.
#'
#' @param data A `data.frame` or `data.table`.
#' @param target Column name of the continuous target.
#' @param features Column names of the drivers.
#' @param train_idx,holdout_idx Row indices; `NULL` uses every row for
#'   training and skips the revalidation.
#' @param min_bins,max_bins Target range of bins per driver.
#' @param min_share Minimum share of training rows per bin.
#' @param min_n Minimum number of training rows per bin.
#' @param monotone `"auto"` (direction from the Spearman sign),
#'   `"increasing"`, `"decreasing"` or `"none"`.
#' @param scale `"mean"` (bin mean in the `woe` slot) or `"logit"`.
#' @param nthread Parallel workers by driver, through the package backend.
#' @param alpha Alpha of the PSI critical value in the revalidation.
#'
#' @return An object of class `scr_cbins`: `fit` (the `obwoe`-shaped
#'   object), `summary` (one row per driver: `feature`, `type`, `n_bins`,
#'   `eta2`, `direction`, `converged`, and after revalidation `eta2_holdout`,
#'   `psi`, `psi_flag`, `holdout_ok`, `holdout_reason`), `holdout` (bin
#'   table per driver with train and hold-out means), `scale` and `target`.
#'   `summary` keeps the engine columns (`algorithm`, `total_iv`,
#'   `iterations`, `error`) and, after revalidation, `psi_critical`,
#'   `psi_flag_adjusted` and `pct_unbinned`.
#'
#' @family irb-ead
#' @examples
#' set.seed(1)
#' d <- data.frame(x = runif(600), g = sample(c("a", "b", "c", "d"), 600, TRUE))
#' d$y <- pmin(1, pmax(0, 0.2 + 0.6 * d$x + (d$g == "d") * 0.2 + rnorm(600, 0, 0.1)))
#' cb <- scr_bin_continuous(d, "y", c("x", "g"), train_idx = 1:400, holdout_idx = 401:600)
#' cb
#' cb$fit$results$x$bin
#' cb$fit$results$x$woe    # bin means of y
#' @export
scr_bin_continuous <- function(data, target, features, train_idx = NULL, holdout_idx = NULL,
                               min_bins = 2L, max_bins = 6L, min_share = 0.05, min_n = 30L,
                               monotone = c("auto", "increasing", "decreasing", "none"),
                               scale = c("mean", "logit"), nthread = 1L, alpha = 0.05) {
  monotone <- match.arg(monotone); scale <- match.arg(scale)
  dt <- data.table::as.data.table(data)
  miss <- setdiff(c(target, features), names(dt))
  if (length(miss)) stop("scr_bin_continuous(): column(s) not found: ", lst(miss), call. = FALSE)
  y <- as.double(dt[[target]])
  if (anyNA(y)) stop("scr_bin_continuous(): the target has missing values.", call. = FALSE)
  train_idx <- train_idx %||% seq_len(nrow(dt))
  fit <- .scr_bin_continuous(dt[train_idx], target, features, min_bins = min_bins, max_bins = max_bins,
                             min_share = min_share, min_n = min_n, monotone = monotone, scale = scale,
                             nthread = nthread)
  s <- data.table::as.data.table(fit$summary)
  holdout <- NULL
  if (!is.null(holdout_idx) && length(holdout_idx)) {
    ho <- .cbins_holdout(fit, dt[train_idx], dt[holdout_idx], target, alpha = alpha)
    holdout <- ho$bins
    s <- merge(s, ho$summary, by = "feature", all.x = TRUE, sort = FALSE)
  }
  structure(list(fit = fit, summary = s, holdout = holdout, scale = scale, target = target),
            class = c("scr_cbins", "list"))
}

#' @export
print.scr_cbins <- function(x, ...) {
  s <- x$summary
  cat(sprintf("<scr_cbins> %d driver(s) binned against '%s' (bin statistic: %s)\n", nrow(s), x$target, x$scale))
  for (i in seq_len(nrow(s))) {
    extra <- if ("holdout_ok" %in% names(s)) {
      sprintf(" | hold-out eta2 %.3f, PSI %.3f (%s)%s", s$eta2_holdout[i], s$psi[i], s$psi_flag[i],
              if (isTRUE(s$holdout_ok[i])) "" else paste0(" - ", s$holdout_reason[i]))
    } else ""
    cat(sprintf("  %-24s %-11s %d bins | eta2 %.3f | %s%s\n", s$feature[i], s$type[i], s$n_bins[i], s$eta2[i],
                s$direction[i], extra))
  }
  invisible(x)
}

#' The binner proper: an `obwoe`-shaped object for a continuous target
#' @keywords internal
#' @noRd
.scr_bin_continuous <- function(dt, target, features, min_bins = 2L, max_bins = 6L, min_share = 0.05,
                                min_n = 30L, monotone = "auto", scale = "mean", nthread = 1L) {
  y <- as.double(dt[[target]])
  entries <- .scr_lapply(features, function(f) {
    x <- dt[[f]]
    if (is.numeric(x)) {
      if (anyNA(x)) stop("driver '", f, "' has missing values; run scr_triage() or impute before binning.", call. = FALSE)
      .cbin_num(f, as.double(x), y, min_bins, max_bins, min_share, min_n, monotone, scale)
    } else {
      .cbin_cat(f, as.character(x), y, min_bins, max_bins, min_share, min_n, scale)
    }
  }, nthread = nthread)
  names(entries) <- features
  summary <- data.frame(
    feature = features,
    type = vapply(entries, `[[`, character(1), "type"),
    algorithm = "scr_continuous",
    n_bins = vapply(entries, function(e) length(e$bin), integer(1)),
    total_iv = vapply(entries, function(e) sum(e$iv), numeric(1)),
    converged = vapply(entries, `[[`, logical(1), "converged"),
    iterations = vapply(entries, `[[`, integer(1), "iterations"),
    error = FALSE, stringsAsFactors = FALSE)
  summary$eta2 <- summary$total_iv
  summary$direction <- vapply(entries, function(e) e$direction %||% "none", character(1))
  structure(list(results = entries, summary = summary, n_features = length(features), target = target,
                 algorithm = "scr_continuous", target_type = "continuous", scale = scale),
            class = "obwoe")
}

#' Bin statistic and share of the between-bin sum of squares
#' @keywords internal
#' @noRd
.cbin_stats <- function(idx, y, k, scale) {
  n_b <- tabulate(idx, nbins = k)
  s_b <- .cbin_sum_by(y, idx, k)
  m_b <- ifelse(n_b > 0, s_b / pmax(n_b, 1L), NA_real_)
  ss_tot <- sum((y - mean(y))^2)
  ss_b <- n_b * (m_b - mean(y))^2
  iv <- if (ss_tot > 0) ss_b / ss_tot else rep(0, k)
  w <- if (identical(scale, "logit")) stats::qlogis(pmin(pmax(m_b, 1e-4), 1 - 1e-4)) else m_b
  list(n = n_b, mean = m_b, woe = w, iv = ifelse(is.na(iv), 0, iv), sum = s_b)
}

#' Sum of `y` by bin index 1..k in one pass (0 for an empty bin)
#' @keywords internal
#' @noRd
.cbin_sum_by <- function(y, idx, k) {
  out <- numeric(k)
  ok <- !is.na(idx)
  if (!any(ok)) return(out)
  r <- rowsum(as.double(y[ok]), idx[ok], reorder = TRUE)
  out[as.integer(rownames(r))] <- r[, 1]
  out
}

#' Greedy merge of adjacent bins on a sequence of (n, sum) with a size floor
#'
#' Cost of merging bins a and b: n_a n_b / (n_a + n_b) (m_a - m_b)^2, the
#' between-bin sum of squares lost. Returns the group index of each input bin.
#' @keywords internal
#' @noRd
.cbin_merge <- function(n, s, max_bins, min_n, min_share, forced_pairs = NULL) {
  k <- length(n)
  groups <- as.list(seq_len(k))
  N <- sum(n)
  gn <- function(g) sum(n[g]); gs <- function(g) sum(s[g])
  it <- 0L
  repeat {
    it <- it + 1L
    m <- length(groups)
    if (m <= 1L) break
    nn <- vapply(groups, gn, numeric(1)); ss <- vapply(groups, gs, numeric(1)); mm <- ss / pmax(nn, 1)
    small <- which(nn < min_n | nn / N < min_share)
    if (length(small)) {
      i <- small[which.min(nn[small])]
      # attach the small group to the neighbour with the closer mean
      cand <- c(if (i > 1L) i - 1L, if (i < m) i + 1L)
      j <- cand[which.min(abs(mm[cand] - mm[i]))]
    } else if (m > max_bins) {
      cost <- vapply(seq_len(m - 1L), function(i) nn[i] * nn[i + 1L] / (nn[i] + nn[i + 1L]) * (mm[i] - mm[i + 1L])^2, numeric(1))
      i <- which.min(cost); j <- i + 1L
    } else break
    a <- min(i, j); b <- max(i, j)
    groups[[a]] <- c(groups[[a]], groups[[b]]); groups[[b]] <- NULL
  }
  g <- integer(k)
  for (gi in seq_along(groups)) g[groups[[gi]]] <- gi
  list(group = g, iterations = it)
}

#' Pool-adjacent-violators on bin means: merges adjacent bins that violate the order
#' @keywords internal
#' @noRd
.cbin_pava <- function(n, s, increasing) {
  groups <- as.list(seq_along(n))
  repeat {
    nn <- vapply(groups, function(g) sum(n[g]), numeric(1)); mm <- vapply(groups, function(g) sum(s[g]), numeric(1)) / pmax(nn, 1)
    if (length(groups) <= 1L) break
    d <- diff(mm)
    viol <- if (increasing) which(d < 0) else which(d > 0)
    if (!length(viol)) break
    i <- viol[1]
    groups[[i]] <- c(groups[[i]], groups[[i + 1L]]); groups[[i + 1L]] <- NULL
  }
  g <- integer(length(n))
  for (gi in seq_along(groups)) g[groups[[gi]]] <- gi
  g
}

#' Numeric driver
#' @keywords internal
#' @noRd
.cbin_num <- function(f, x, y, min_bins, max_bins, min_share, min_n, monotone, scale) {
  n <- length(x)
  ux <- unique(x)
  if (length(ux) < 2L || n < 2L * min_n) {
    return(.cbin_entry(f, "numerical", labels = .num_label(-Inf, Inf), cutpoints = numeric(), idx = rep(1L, n),
                       y = y, scale = scale, converged = FALSE, iterations = 0L, direction = "none"))
  }
  probs <- seq(0, 1, length.out = 22L)[-c(1L, 22L)]
  q <- unique(stats::quantile(x, probs = probs, type = 7, names = FALSE))
  q <- q[q > min(x) & q < max(x)]
  q <- round(q, 6)
  q <- unique(q)
  if (!length(q)) {
    return(.cbin_entry(f, "numerical", labels = .num_label(-Inf, Inf), cutpoints = numeric(), idx = rep(1L, n),
                       y = y, scale = scale, converged = FALSE, iterations = 0L, direction = "none"))
  }
  pre <- findInterval(x, q, left.open = TRUE) + 1L          # (q[i-1], q[i]] convention
  k <- length(q) + 1L
  st <- .cbin_stats(pre, y, k, scale)
  keep <- st$n > 0
  # drop empty prebins by collapsing their cut points
  if (!all(keep)) {
    q <- q[keep[-k] | FALSE]
    q <- q[seq_len(sum(keep) - 1L)]
    pre <- findInterval(x, q, left.open = TRUE) + 1L; k <- length(q) + 1L; st <- .cbin_stats(pre, y, k, scale)
  }
  mg <- .cbin_merge(st$n, st$sum, max_bins, min_n, min_share)
  g <- mg$group
  # monotone direction
  direction <- "none"
  if (!identical(monotone, "none")) {
    inc <- switch(monotone,
      auto = suppressWarnings(stats::cor(x, y, method = "spearman")) >= 0,
      increasing = TRUE, decreasing = FALSE)
    if (is.na(inc)) inc <- TRUE
    n_g <- vapply(seq_len(max(g)), function(i) sum(st$n[g == i]), numeric(1))
    s_g <- vapply(seq_len(max(g)), function(i) sum(st$sum[g == i]), numeric(1))
    g2 <- .cbin_pava(n_g, s_g, inc)
    g <- g2[g]
    direction <- if (inc) "increasing" else "decreasing"
  }
  # cut points = the prebin edges where the group changes
  edges <- q[which(diff(g) != 0)]
  idx <- findInterval(x, edges, left.open = TRUE) + 1L
  labels <- .num_labels(edges)
  .cbin_entry(f, "numerical", labels = labels, cutpoints = edges, idx = idx, y = y, scale = scale,
              converged = length(edges) + 1L >= min_bins, iterations = mg$iterations, direction = direction)
}

#' Categorical driver
#' @keywords internal
#' @noRd
.cbin_cat <- function(f, x, y, min_bins, max_bins, min_share, min_n, scale, sep = "%;%") {
  x[is.na(x)] <- "NA"
  if (any(grepl(sep, x, fixed = TRUE))) stop("driver '", f, "': a category contains the bin separator '", sep, "'.", call. = FALSE)
  lv <- sort(unique(x))
  li <- match(x, lv)
  n_l <- stats::setNames(as.double(tabulate(li, nbins = length(lv))), lv)
  m_l <- .cbin_sum_by(y, li, length(lv)) / n_l
  ord <- order(m_l, n_l)
  lv <- lv[ord]; n_l <- n_l[ord]; s_l <- (m_l * n_l)[ord]
  mg <- .cbin_merge(n_l, s_l, max_bins, min_n, min_share)
  g <- mg$group
  groups <- lapply(seq_len(max(g)), function(i) lv[g == i])
  labels <- vapply(groups, function(gr) paste(gr, collapse = sep), character(1))
  map <- stats::setNames(rep(seq_along(groups), lengths(groups)), unlist(groups))
  idx <- unname(map[x])
  .cbin_entry(f, "categorical", labels = labels, cutpoints = numeric(), idx = idx, y = y, scale = scale,
              converged = length(groups) >= min_bins, iterations = mg$iterations, direction = "ordered_by_mean")
}

#' Engine-style labels for a set of interior cut points
#' @keywords internal
#' @noRd
.num_labels <- function(edges) {
  e <- c(-Inf, edges, Inf)
  vapply(seq_len(length(e) - 1L), function(i) .num_label(e[i], e[i + 1L]), character(1))
}

#' Assemble an `obwoe`-shaped entry from bin memberships
#' @keywords internal
#' @noRd
.cbin_entry <- function(f, type, labels, cutpoints, idx, y, scale, converged, iterations, direction) {
  k <- length(labels)
  st <- .cbin_stats(idx, y, k, scale)
  list(id = as.numeric(seq_len(k)), bin = unname(labels), woe = unname(st$woe), iv = unname(st$iv),
       count = as.integer(st$n), count_pos = as.integer(round(st$sum)), count_neg = as.integer(st$n - round(st$sum)),
       cutpoints = as.double(cutpoints), converged = isTRUE(converged), iterations = as.integer(iterations),
       feature = f, type = type, algorithm = "scr_continuous",
       mean = unname(st$mean), direction = direction)
}

#' Hold-out revalidation of frozen continuous bins
#' @keywords internal
#' @noRd
.cbins_holdout <- function(fit, dt_train, dt_holdout, target, alpha = 0.05) {
  feats <- names(fit$results)
  y_ho <- as.double(dt_holdout[[target]])
  app_ho <- .cbins_apply_idx(fit, dt_holdout, feats)
  rows <- lapply(feats, function(f) {
    e <- fit$results[[f]]
    k <- length(e$bin)
    idx <- app_ho[[f]]
    ok <- !is.na(idx)
    st <- .cbin_stats(idx[ok], y_ho[ok], k, "mean")
    tr_share <- e$count / sum(e$count)
    ps <- if (k >= 2L) scr_psi(rep(seq_len(k), e$count), idx[ok], levels = seq_len(k), alpha = alpha)
          else list(psi = 0, flag_fixed = "stable", critical = NA_real_, flag_adjusted = "stable")
    order_ok <- if (k >= 2L && !identical(e$direction, "none")) {
      d <- diff(st$mean[st$n > 0])
      if (identical(e$direction, "decreasing")) all(d <= 1e-9) else all(d >= -1e-9)
    } else TRUE
    # scr_psi() flags "stable" / "moderate" / "shift": a shift (PSI >= 0.25) is the action level
    reason <- c(if (!order_ok) "UNSTABLE_HOLDOUT", if (identical(ps$flag_fixed, "shift")) "PSI_ACTION",
                if (mean(ok) < 0.99) "UNBINNED_HOLDOUT")
    list(
      bins = data.table::data.table(feature = f, bin = e$bin, n_train = e$count, mean_train = e$mean,
                                    share_train = tr_share, n_holdout = st$n, mean_holdout = st$mean,
                                    share_holdout = st$n / max(1L, sum(st$n))),
      summary = data.table::data.table(feature = f, eta2_holdout = sum(st$iv), psi = ps$psi, psi_flag = ps$flag_fixed,
                                       psi_critical = ps$critical, psi_flag_adjusted = ps$flag_adjusted,
                                       pct_unbinned = 1 - mean(ok),
                                       holdout_ok = !length(reason),
                                       holdout_reason = if (length(reason)) paste(reason, collapse = ";") else "OK"))
  })
  list(bins = data.table::rbindlist(lapply(rows, `[[`, "bins")),
       summary = data.table::rbindlist(lapply(rows, `[[`, "summary")))
}

#' Bin index (1..k) of every row for every driver of an obwoe-shaped fit
#'
#' Numerics follow the `(lo, hi]` convention of the engine; categoricals
#' match the `%;%` groups; `NA` for an unseen category or a missing numeric.
#' @keywords internal
#' @noRd
.cbins_apply_idx <- function(fit, dt, features = names(fit$results), sep = "%;%") {
  out <- lapply(features, function(f) {
    e <- fit$results[[f]]
    x <- dt[[f]]
    if (identical(e$type, "numerical")) {
      x <- as.double(x)
      i <- findInterval(x, e$cutpoints, left.open = TRUE) + 1L
      i[is.na(x)] <- NA_integer_
      i
    } else {
      x <- as.character(x); x[is.na(x)] <- "NA"
      groups <- strsplit(e$bin, sep, fixed = TRUE)
      map <- stats::setNames(rep(seq_along(groups), lengths(groups)), unlist(groups))
      unname(map[x])
    }
  })
  names(out) <- features
  data.table::as.data.table(out)
}

#' Bin statistic of every row: the `woe` slot looked up by bin index
#' @keywords internal
#' @noRd
.cbins_apply_value <- function(fit, dt, features = names(fit$results), unbinned = NA_real_) {
  idx <- .cbins_apply_idx(fit, dt, features)
  out <- lapply(features, function(f) {
    v <- fit$results[[f]]$woe[idx[[f]]]
    v[is.na(v)] <- unbinned
    v
  })
  names(out) <- paste0(features, "_woe")
  data.table::as.data.table(out)
}
