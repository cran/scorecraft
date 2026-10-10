# ============================================================================ #
# uplift.R - a score read on a treated and a control group
# ============================================================================ #
# One grouped pass counts the rows per arm and score cell. The bands, the
# Qini and uplift curves and their bootstrap are functions of those counts:
# a resample draws the cell counts of each arm and outcome from multinomial
# laws, at a cost proportional to the cells.
# ============================================================================ #

#' Uplift of a treatment along a score
#'
#' Reads a score on a treated and a control group (a campaign with a
#' hold-out): the difference of the event rates per score band with its
#' interval, the incremental events, the Qini and uplift curves with their
#' areas, and a check that the two groups have the same score distribution.
#'
#' @section Bands:
#'
#' The bands are cut on the score of both arms together, tie-safe as in
#' [scr_bands()], the event-richest end first (`band = 1`). Per band, with
#' \eqn{n_t}, \eqn{n_c} the rows with a known outcome and \eqn{r_t},
#' \eqn{r_c} the event rates of the treated and of the control:
#' \itemize{
#'   \item `uplift` = \eqn{r_t - r_c}, with the Newcombe hybrid score
#'     interval: from the Wilson intervals \eqn{[l_t, u_t]} and
#'     \eqn{[l_c, u_c]} of the two rates,
#'     \deqn{lo = r_t - r_c - \sqrt{(r_t - l_t)^2 + (u_c - r_c)^2}, \qquad
#'           hi = r_t - r_c + \sqrt{(u_t - r_t)^2 + (r_c - l_c)^2}.}
#'     Under weights the Wilson intervals use the Kish effective sizes.
#'   \item `incremental` = `uplift` \eqn{\times n_t}: the events of the
#'     treated that the treatment added; `cum_incremental` is its running
#'     sum from the first band (a band with an empty arm adds nothing) and
#'     `cum_pct_treated` the cumulative share of the treated.
#'   \item `type`, a convention of this package: `"persuadable"` when the
#'     interval lies above 0, `"negative"` when it lies below, `"no effect"`
#'     when it covers 0.
#' }
#'
#' @section Qini and AUUC:
#'
#' The score cells are accumulated from the event-rich end. With
#' \eqn{E_t(k)}, \eqn{E_c(k)} the cumulative events and \eqn{N_t(k)},
#' \eqn{N_c(k)} the cumulative rows of each arm after \eqn{k} cells:
#' \itemize{
#'   \item the Qini curve is
#'     \eqn{Q(k) = E_t(k) - E_c(k) N_t(k) / N_c(k)} (the second term is 0
#'     while no control row is in), read against \eqn{N_t(k)};
#'   \item `qini` is the area between the Qini curve and the straight line
#'     from the origin to its end point, by trapezoids from the origin,
#'     divided by the square of the treated volume:
#'     \deqn{qini = \frac{1}{N_t^2} \sum_k \frac{Q(k-1) + Q(k)}{2}
#'           (N_t(k) - N_t(k-1)) - \frac{Q(K)}{2 N_t}.}
#'     It is 0 for a score that orders at random and positive when the
#'     incremental events come first;
#'   \item the uplift curve is
#'     \eqn{U(k) = (E_t(k) / N_t(k) - E_c(k) / N_c(k)) (N_t(k) + N_c(k)) / N}
#'     (an arm with no row yet has rate 0), read against the cumulative
#'     share of all rows \eqn{(N_t(k) + N_c(k)) / N}; `auuc` is the area
#'     under it by trapezoids from the origin. A score that orders at random
#'     gives about half the overall uplift.
#' }
#' Rows with a missing outcome are left out of the rates and of the curves.
#'
#' `qini` is not Radcliffe's Q, which divides the same area by that of an
#' ideal curve; here the area is divided by \eqn{N_t^2} only, so it reads
#' as incremental events per treated row, averaged over the depth.
#'
#' The intervals of `qini` and `auuc` are percentile intervals of a
#' bootstrap on the counts, stratified by arm: each resample draws, for the
#' treated and for the control separately, the counts per score cell and
#' outcome from one multinomial law with the observed shares and the
#' unweighted size of the arm, the law of a row bootstrap within each arm.
#' The event totals of the arms vary from one resample to the next, so the
#' intervals carry the uncertainty of the overall uplift as well as that of
#' the ordering. Above `boot_cells` cells the resamples run on pooled cells
#' and are shifted to the point estimates, as in [scr_bands()]. A given
#' `seed` is local to the call.
#'
#' @section Randomization check:
#'
#' Under a randomized assignment the score has the same distribution in
#' both arms. `randomization` holds the PSI of the treated against the
#' control over the bands, with the n-adjusted critical value at
#' `1 - level` (see [scr_psi()]); when the PSI exceeds it, `flag` is `TRUE`
#' and `note` says so: the arms then differ along the score and the uplift
#' may reflect the assignment, not the treatment.
#'
#' @inheritParams scr_bands
#' @param x A `data.frame` with one row per case.
#' @param score,y,treat Column names of the score, of the 0/1 outcome (`NA`
#'   allowed) and of the 0/1 or logical treatment indicator (1 = treated).
#' @param n_bands Equal-share bands of the score when `cuts` is not given.
#' @param cuts Optional ascending cuts of the score, or an object from
#'   [scr_bands()] or [scr_tiers()], whose cuts, numbers and labels are then
#'   used, with its objective and direction.
#' @param level Confidence level of the intervals.
#' @param n_boot Bootstrap resamples of the Qini and AUUC intervals (`0`
#'   skips them).
#' @param objective `"propensity"` (default: the event is the outcome
#'   sought, and a higher score means a higher propensity) or `"risk"`.
#'
#' @return An object of class `c("scr_uplift", "list")`:
#'   \describe{
#'     \item{`table`}{One row per band, event-richest first: `band`,
#'       `label`, `score_lo`, `score_hi`, `n_t`, `n_c`, `events_t`,
#'       `events_c`, `rate_t`, `rate_c`, `uplift`, `uplift_lo`, `uplift_hi`,
#'       `incremental`, `cum_incremental`, `cum_pct_treated` and `type`.}
#'     \item{`summary`}{One row: `n_t`, `n_c`, `rate_t`, `rate_c`, `uplift`,
#'       `uplift_lo`, `uplift_hi`, `qini`, `qini_lo`, `qini_hi`, `auuc`,
#'       `auuc_lo`, `auuc_hi`, `psi`, `psi_critical` and `psi_flag`.}
#'     \item{`curve`}{The curves, at most about 1,000 rows spread evenly in
#'       the treated share (the areas use every cell): `cut` (the score
#'       boundary), `pct_treated`, `pct_all`, `n_t`, `n_c`, `qini`,
#'       `qini_random` (the straight line) and `uplift`.}
#'     \item{`randomization`}{A list: `psi`, `critical`, `flag` and `note`.}
#'     \item{`cuts`, `codes`, `labels`, `level`, `n_boot`, `objective`,
#'       `direction`, `score`, `target`, `treat`, `weighted`, `call`}{The
#'       cuts and the settings.}
#'   }
#'
#' @references
#' Newcombe, R. G. (1998). Interval estimation for the difference between
#' independent proportions: comparison of eleven methods. *Statistics in
#' Medicine*, 17(8), 873-890.
#'
#' Radcliffe, N. J. (2007). Using control groups to target on predicted
#' lift: building and assessing uplift models. *Direct Marketing Analytics
#' Journal*, 1, 14-21.
#'
#' Yurdakul, B. and Naranjo, J. (2020). Statistical properties of the
#' population stability index. *Journal of Risk Model Validation*, 14(4),
#' 89-100.
#'
#' @seealso [scr_bands()] for the bands, [scr_operating()] for the cut of a
#'   campaign under a budget.
#' @family score-studies
#' @examples
#' local({
#'   set.seed(1)
#'   n <- 8000
#'   x <- rnorm(n)
#'   treat <- rbinom(n, 1, 0.5)
#'   # the offer works on the customers with a high score only
#'   d <- data.frame(score = round(500 + 50 * x), treat = treat,
#'                   y = rbinom(n, 1, plogis(-1.5 + 0.5 * x + treat * pmax(x, 0))))
#'   up <- scr_uplift(d, n_bands = 5, n_boot = 50, seed = 1)
#'   print(up)
#'   up$summary[, c("uplift", "qini", "qini_lo", "qini_hi", "auuc")]
#' })
#' @export
scr_uplift <- function(x, ...) UseMethod("scr_uplift")

#' @rdname scr_uplift
#' @export
scr_uplift.data.frame <- function(x, score = "score", y = "y", treat = "treat", n_bands = 10L, cuts = NULL,
                                  level = 0.95, n_boot = 200L, seed = NULL, objective = "propensity",
                                  direction = NULL, weight = NULL, max_cells = 1e5, boot_cells = 1e4, ...) {
  fn <- "scr_uplift"
  .study_dots(list(...), fn)
  for (nm in c("score", "y", "treat", "weight")) {
    v <- get(nm)
    if (!is.null(v)) .study_chr1(v, nm, fn)
  }
  rd <- .study_one_direction(if (!missing(objective)) objective, direction, cuts, "propensity", fn)
  miss <- setdiff(c(score, y, treat, weight), names(x))
  if (length(miss)) stop(fn, "(): column(s) ", lst(miss), " not in `x`.", call. = FALSE)
  for (nm in c(score, weight)) {
    if (!is.numeric(x[[nm]])) stop(fn, "(): column '", nm, "' must be numeric.", call. = FALSE)
  }
  n_bands <- .study_whole(n_bands, "n_bands", fn, lower = 1)
  n_boot <- .study_whole(n_boot, "n_boot", fn)
  boot_cells <- .study_boot_cells(boot_cells, fn)
  level <- .study_level(level, fn)
  tr <- tryCatch(.scr_y01(x[[treat]], fn), error = function(e) NULL)
  if (is.null(tr) || anyNA(tr)) {
    stop(fn, "(): `treat` must be a 0/1 numeric or logical column without missing values.", call. = FALSE)
  }
  h <- .study_hist(x[[score]], x[[y]], w = if (!is.null(weight)) x[[weight]], by = list(sample = tr),
                   max_cells = max_cells, fn = fn)
  ht <- .study_cells(h, 1L); hc <- .study_cells(h, 0L)
  if (!nrow(ht) || !nrow(hc)) stop(fn, "(): both a treated and a control group are needed.", call. = FALSE)
  cells <- .study_collapse(h)
  side <- .study_side(rd$direction)
  bd <- .cross_bands(cuts, cells, n_bands, side, "cuts", fn)
  B <- length(bd$cuts) + 1L
  o <- if (side == "high") rev(seq_len(B)) else seq_len(B)
  a <- 1 - level

  # -- bands -------------------------------------------------------------------- #
  cols <- c("n", "w2", "n_y", "e", "w2_y")
  St <- lapply(.study_sum(ht, findInterval(ht$s, bd$cuts) + 1L, B, cols), `[`, o)
  Sc <- lapply(.study_sum(hc, findInterval(hc$s, bd$cuts) + 1L, B, cols), `[`, o)
  ut <- .uplift_diff(St$e, St$n_y, St$w2_y, Sc$e, Sc$n_y, Sc$w2_y, level)
  inc <- ut$uplift * St$n_y
  NT <- sum(St$n_y); NC <- sum(Sc$n_y)
  type <- ifelse(is.na(ut$lo), NA_character_, ifelse(ut$lo > 0, "persuadable", ifelse(ut$hi < 0, "negative", "no effect")))
  tab <- data.table::data.table(
    band = bd$codes[o], label = bd$labels[o], score_lo = c(-Inf, bd$cuts)[o], score_hi = c(bd$cuts, Inf)[o],
    n_t = St$n_y, n_c = Sc$n_y, events_t = St$e, events_c = Sc$e, rate_t = ut$rate_t, rate_c = ut$rate_c,
    uplift = ut$uplift, uplift_lo = ut$lo, uplift_hi = ut$hi, incremental = inc,
    cum_incremental = cumsum(ifelse(is.na(inc), 0, inc)),
    cum_pct_treated = if (NT > 0) cumsum(St$n_y) / NT else NA_real_, type = type)
  all_ <- .uplift_diff(sum(St$e), NT, sum(St$w2_y), sum(Sc$e), NC, sum(Sc$w2_y), level)

  # -- randomization: the score distribution of the treated against the control -- #
  ps <- .psi_counts(Sc$n / sum(Sc$n) * .study_kish(sum(Sc$n), sum(Sc$w2)),
                    St$n / sum(St$n) * .study_kish(sum(St$n), sum(St$w2)), bd$labels[o], a, c(0.10, 0.25))
  flag <- isTRUE(ps$psi > ps$critical)
  note <- if (flag) sprintf(paste0("the score distribution differs between the arms (PSI %s above the n-adjusted ",
                                   "critical value %s): the assignment may not be random along the score"),
                            .g3(ps$psi), .g3(ps$critical)) else NA_character_

  # -- curves from the cells, the event-rich end first --------------------------- #
  key <- if ("g" %in% names(cells)) "g" else "s"
  K <- nrow(cells)
  sel <- if (side == "high") rev(seq_len(K)) else seq_len(K)
  put <- function(hh, v) { z <- numeric(K); z[match(hh[[key]], cells[[key]])] <- v; z[sel] }
  te <- put(ht, ht$e); tn <- put(ht, pmax(ht$n_y - ht$e, 0))
  ce <- put(hc, hc$e); cn <- put(hc, pmax(hc$n_y - hc$e, 0))
  cv <- .uplift_curves(cbind(te), cbind(tn), cbind(ce), cbind(cn))
  qini <- auuc <- NA_real_
  q_ci <- a_ci <- c(NA_real_, NA_real_)
  curve <- NULL
  if (NT > 0 && NC > 0) {
    qini <- cv$qini; auuc <- cv$auuc
    bc <- .study_bcut(cells)
    cut <- if (side == "high") c(bc[rev(seq_len(K - 1L))], -Inf) else c(bc, Inf)
    keep <- if (K <= 1000L) seq_len(K) else sort(unique(c(.op_nearest(cv$pct_t[, 1], seq_len(1000L) / 1000), K)))
    curve <- data.table::data.table(
      cut = cut[keep], pct_treated = cv$pct_t[keep, 1], pct_all = cv$pct_all[keep, 1], n_t = cv$nt[keep, 1],
      n_c = cv$nc[keep, 1], qini = cv$q[keep, 1], qini_random = cv$q[K, 1] * cv$pct_t[keep, 1],
      uplift = cv$u[keep, 1])
    if (n_boot >= 2L) {
      # the unweighted rows with a known outcome of each arm: the sizes of the resamples
      bt <- .uplift_boot(te, tn, ce, cn, c(sum(ht$n_y_raw), sum(hc$n_y_raw)), n_boot, seed, boot_cells,
                         c(qini, auuc))
      if (!is.null(bt)) {
        q_ci <- stats::quantile(bt$qini, c(a / 2, 1 - a / 2), names = FALSE)
        a_ci <- stats::quantile(bt$auuc, c(a / 2, 1 - a / 2), names = FALSE)
      }
    }
  }
  summ <- data.table::data.table(
    n_t = NT, n_c = NC, rate_t = all_$rate_t, rate_c = all_$rate_c, uplift = all_$uplift, uplift_lo = all_$lo,
    uplift_hi = all_$hi, qini = qini, qini_lo = q_ci[1], qini_hi = q_ci[2], auuc = auuc, auuc_lo = a_ci[1],
    auuc_hi = a_ci[2], psi = ps$psi, psi_critical = ps$critical, psi_flag = flag)
  structure(list(table = tab, summary = summ, curve = curve,
                 randomization = list(psi = ps$psi, critical = ps$critical, flag = flag, note = note),
                 cuts = bd$cuts, codes = bd$codes, labels = bd$labels, bands = bd$source, level = level,
                 n_boot = if (all(is.na(q_ci))) 0L else n_boot, objective = rd$objective, direction = rd$direction,
                 score = score, target = y, treat = treat, weighted = !is.null(weight),
                 quantized = !is.null(attr(h, "edges")), call = sys.call()),
            class = c("scr_uplift", "list"))
}

#' Wilson score interval of a proportion (`x` events in `n`, both may be effective)
#' @keywords internal
#' @noRd
.uplift_wilson <- function(x, n, level) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  ok <- is.finite(n) & n > 0 & is.finite(x)
  lo <- hi <- rep(NA_real_, length(n))
  p <- pmin(pmax(x[ok] / n[ok], 0), 1); m <- n[ok]
  mid <- (p + z^2 / (2 * m)) / (1 + z^2 / m)
  half <- z * sqrt(p * (1 - p) / m + z^2 / (4 * m^2)) / (1 + z^2 / m)
  lo[ok] <- pmax(mid - half, 0); hi[ok] <- pmin(mid + half, 1)
  list(lo = lo, hi = hi)
}

#' Difference of two rates with the Newcombe hybrid score interval
#'
#' `e`, `n` and `w2` are the events, the volume with a known outcome and the
#' sum of squared weights of the treated (`t`) and of the control (`c`); the
#' Wilson intervals run on the Kish effective sizes. `NA` when an arm is
#' empty.
#' @keywords internal
#' @noRd
.uplift_diff <- function(et, nt, w2t, ec, nc, w2c, level) {
  rt <- ifelse(nt > 0, et / nt, NA_real_); rc <- ifelse(nc > 0, ec / nc, NA_real_)
  mt <- .study_kish(nt, w2t); mc <- .study_kish(nc, w2c)
  wt <- .uplift_wilson(rt * mt, mt, level); wc <- .uplift_wilson(rc * mc, mc, level)
  d <- rt - rc
  list(rate_t = rt, rate_c = rc, uplift = d, lo = d - sqrt((rt - wt$lo)^2 + (wc$hi - rc)^2),
       hi = d + sqrt((wt$hi - rt)^2 + (rc - wc$lo)^2))
}

#' Qini and uplift curves of every column of four count matrices
#'
#' Rows are the score cells in selection order (the event-rich end first),
#' columns are samples: the events and non-events of the treated (`te`,
#' `tn`) and of the control (`ce`, `cn`). Returns the curves and their
#' areas, `qini` and `auuc`, per column; the areas are `NA` when an arm is
#' empty.
#' @keywords internal
#' @noRd
.uplift_curves <- function(te, tn, ce, cn) {
  cs <- function(M) { r <- apply(M, 2L, cumsum); if (is.matrix(r)) r else matrix(r, nrow = nrow(M)) }
  K <- nrow(te)
  nt <- cs(te + tn); nc <- cs(ce + cn); et <- cs(te); ec <- cs(ce)
  NT <- nt[K, ]; NC <- nc[K, ]; N <- NT + NC
  # Qini: the control events rescaled to the treated volume; nothing while no control row is in
  q <- et - ifelse(nc > 0, ec * nt / nc, 0)
  u <- sweep((ifelse(nt > 0, et / nt, 0) - ifelse(nc > 0, ec / nc, 0)) * (nt + nc), 2L, N, "/")
  pct_t <- sweep(nt, 2L, NT, "/"); pct_all <- sweep(nt + nc, 2L, N, "/")
  # trapezoids from the origin
  lag0 <- function(M) rbind(0, M[-K, , drop = FALSE])
  area <- function(x, y) colSums((x - lag0(x)) * (y + lag0(y)) / 2)
  # the area less that under the straight line to the end point, per squared treated volume
  qini <- (area(nt, q) - q[K, ] * NT / 2) / NT^2
  auuc <- area(pct_all, u)
  bad <- !(NT > 0) | !(NC > 0)
  qini[bad] <- NA_real_; auuc[bad] <- NA_real_
  list(q = q, u = u, nt = nt, nc = nc, pct_t = pct_t, pct_all = pct_all, qini = qini, auuc = auuc)
}

#' Bootstrap of the Qini coefficient and of the AUUC on the counts
#'
#' Stratified by arm: each resample draws, per arm, the counts per cell and
#' outcome from one multinomial law over the 2K categories with the observed
#' shares and the unweighted size of the arm (`size`: treated, control),
#' rescaled to the weighted total of the arm. The event totals of the arms
#' are free, so the resamples carry the variance of the overall uplift.
#' Above `boot_cells` cells the resamples run on adjacent cells pooled into
#' `boot_cells` of equal share and are shifted by the difference between the
#' full and the pooled point estimates `point` (Qini, AUUC). The seed is
#' local to the call.
#' @keywords internal
#' @noRd
.uplift_boot <- function(te, tn, ce, cn, size, n_boot, seed, boot_cells, point) {
  cnt <- list(te, tn, ce, cn)
  size <- round(size)
  # an arm needs at least one row with a known outcome
  if (size[1] < 1 || size[2] < 1) return(NULL)
  .scr_local_seed(seed)
  shift <- c(0, 0)
  if (length(te) > boot_cells) {
    tot <- te + tn + ce + cn
    g <- pmin(pmax(ceiling(cumsum(tot) / sum(tot) * boot_cells), 1), boot_cells)
    cnt <- lapply(cnt, function(v) as.double(rowsum(v, g, reorder = TRUE)))
    pp <- .uplift_curves(cbind(cnt[[1]]), cbind(cnt[[2]]), cbind(cnt[[3]]), cbind(cnt[[4]]))
    shift <- point - c(pp$qini, pp$auuc)
  }
  K <- length(cnt[[1]])
  ev <- seq_len(K)
  # one arm: events in the first K rows, non-events in the last K
  draw <- function(e, ne, n, m) {
    tot <- sum(e) + sum(ne)
    stats::rmultinom(m, n, c(e, ne) / tot) * (tot / n)
  }
  chunk <- max(1L, min(as.integer(n_boot), as.integer(floor(5e5 / K))))
  qini <- auuc <- numeric(n_boot)
  done <- 0L
  while (done < n_boot) {
    m <- min(chunk, n_boot - done)
    Tm <- draw(cnt[[1]], cnt[[2]], size[1], m)
    Cm <- draw(cnt[[3]], cnt[[4]], size[2], m)
    r <- .uplift_curves(Tm[ev, , drop = FALSE], Tm[K + ev, , drop = FALSE], Cm[ev, , drop = FALSE],
                        Cm[K + ev, , drop = FALSE])
    qini[done + seq_len(m)] <- r$qini + shift[1]; auuc[done + seq_len(m)] <- r$auuc + shift[2]
    done <- done + m
  }
  list(qini = qini, auuc = auuc)
}

#' @export
print.scr_uplift <- function(x, ...) {
  s <- x$summary
  cat(sprintf("<scr_uplift> outcome \"%s\" | treatment \"%s\" | objective %s | %s\n", x$target, x$treat, x$objective,
              x$direction))
  cat(sprintf("  treated %s (rate %s) | control %s (rate %s) | uplift %s [%s, %s]\n", .study_n(s$n_t),
              .claims_pct(s$rate_t), .study_n(s$n_c), .claims_pct(s$rate_c), .mix_pp(s$uplift), .mix_pp(s$uplift_lo),
              .mix_pp(s$uplift_hi)))
  ci <- function(lo, hi) if (is.na(lo)) "" else sprintf(" [%s, %s]", .g3(lo), .g3(hi))
  cat(sprintf("  Qini %s%s | AUUC %s%s%s\n", if (is.na(s$qini)) "-" else .g3(s$qini), ci(s$qini_lo, s$qini_hi),
              if (is.na(s$auuc)) "-" else .g3(s$auuc), ci(s$auuc_lo, s$auuc_hi),
              if (x$n_boot > 0L) sprintf(" (%.0f%% bootstrap intervals, %d resamples)", 100 * x$level, x$n_boot) else ""))
  cat(sprintf("  randomization: PSI of the treated against the control %s (critical %s)\n", .study_f(s$psi, "%.4f"),
              .study_f(s$psi_critical, "%.4f")))
  if (isTRUE(x$randomization$flag)) cat(sprintf("  WARNING: %s\n", x$randomization$note))
  t <- x$table
  cat(sprintf("\nBands (event-richest first; %.0f%% Newcombe interval of the uplift)\n", 100 * x$level))
  cat(sprintf("  %4s %-22s %9s %9s %8s %8s %-34s %11s %8s  %s\n", "band", "score", "n_t", "n_c", "rate_t", "rate_c",
              "uplift [lo, hi]", "incremental", "cum", "type"))
  for (i in seq_len(nrow(t))) {
    cat(sprintf("  %4d %-22s %9s %9s %8s %8s %-34s %11s %8s  %s\n", t$band[i], substr(t$label[i], 1, 22),
                .study_n(t$n_t[i]), .study_n(t$n_c[i]), .claims_pct(t$rate_t[i]), .claims_pct(t$rate_c[i]),
                if (is.na(t$uplift[i])) "-" else sprintf("%s [%s, %s]", .mix_pp(t$uplift[i]), .mix_pp(t$uplift_lo[i]),
                                                         .mix_pp(t$uplift_hi[i])),
                .study_f(t$incremental[i], "%.1f"), .study_f(t$cum_incremental[i], "%.1f"),
                if (is.na(t$type[i])) "-" else t$type[i]))
  }
  invisible(x)
}

#' @rdname scr_export
#' @export
scr_export.scr_uplift <- function(x, dir, stamp = TRUE, ...) {
  .study_dots(list(...), "scr_export")
  .need_openxlsx()
  out_dir <- .export_dir(dir, stamp)
  tag <- .file_tag(x$target)
  settings <- .kv_table(list(score = x$score, target = x$target, treat = x$treat, objective = x$objective,
                             direction = x$direction, level = x$level, n_boot = x$n_boot, weighted = x$weighted,
                             randomization = x$randomization$note %||% NA_character_))
  cuts <- data.frame(cut = seq_along(x$cuts), score = x$cuts)
  sheets <- list(Summary = x$summary, Bands = x$table, Curve = x$curve, Cuts = cuts, Settings = settings)
  sheets <- lapply(sheets[!vapply(sheets, is.null, logical(1))], .study_sheet)
  files <- list(xlsx = .scr_write_xlsx(sheets, file.path(out_dir, sprintf("uplift_%s.xlsx", tag))))
  for (f in files) msg("  %s", f)
  x$files <- files
  invisible(x)
}
