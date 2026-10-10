# Score-study engine: count table, cuts, band table, count bootstrap and
# scr_bands(). Every statistic is checked against an independent formula.

study_df <- function(n = 3000, seed = 11) {
  set.seed(seed)
  x <- stats::rnorm(n)
  data.frame(score = round(600 + 40 * x), y = stats::rbinom(n, 1, stats::plogis(-1.8 - 1.1 * x)),
             smp = sample(c("dev", "oot"), n, TRUE), w = stats::runif(n, 0.2, 3),
             amt = stats::rexp(n, 1 / 500), stringsAsFactors = FALSE)
}

# -- the count table ------------------------------------------------------------ #

test_that(".study_hist equals a naive aggregation, with weights and missing outcomes", {
  set.seed(3)
  n <- 800
  s <- sample(1:60, n, TRUE) / 2
  y <- stats::rbinom(n, 1, 0.3)
  y[sample(n, 40)] <- NA
  s[sample(n, 15)] <- NA
  w <- stats::runif(n, 0, 3); w[sample(n, 25)] <- 0
  v <- stats::rexp(n)
  h <- .study_hist(s, y, w = w, value = v)
  ok <- !is.na(s) & w > 0
  u <- sort(unique(s[ok]))
  expect_equal(h$s, u)
  agg <- function(f) vapply(u, function(q) f(which(ok & s == q)), numeric(1))
  expect_equal(h$n, agg(function(i) sum(w[i])))
  expect_equal(h$n_y, agg(function(i) sum(w[i][!is.na(y[i])])))
  expect_equal(h$e, agg(function(i) sum((w * y)[i], na.rm = TRUE)))
  expect_equal(h$w2, agg(function(i) sum(w[i]^2)))
  expect_equal(h$w2_y, agg(function(i) sum(w[i][!is.na(y[i])]^2)))
  expect_equal(h$n_raw, agg(length))
  expect_equal(h$n_y_raw, agg(function(i) sum(!is.na(y[i]))))
  expect_equal(h$e_raw, agg(function(i) sum(y[i], na.rm = TRUE)))
  expect_equal(h$v, agg(function(i) sum((w * v)[i])))
  expect_equal(h$ve, agg(function(i) sum((w * v * y)[i], na.rm = TRUE)))
  # unweighted: the weighted columns are the counts
  h0 <- .study_hist(s, y)
  expect_equal(h0$n, h0$n_raw); expect_equal(h0$e, h0$e_raw); expect_equal(h0$w2_y, h0$n_y_raw)
  expect_equal(sum(h0$n), sum(!is.na(s)))
  # unit weights give the unweighted table
  h1 <- .study_hist(s, y, w = rep(1, n))
  expect_equal(h1$n, h0$n); expect_equal(h1$e, h0$e); expect_equal(h1$w2_y, h0$w2_y)
  expect_error(.study_hist(s, y, w = -w), "non-negative")
  expect_error(.study_hist(s, y + 2), "0/1")
})

test_that("the count table pools many distinct scores into cells that never straddle an edge", {
  set.seed(5)
  s <- stats::rnorm(4000); y <- stats::rbinom(4000, 1, 0.2)
  smp <- rep(c("a", "b"), 2000)
  h <- .study_hist(s, y, by = list(sample = smp), max_cells = 50, breaks = c(-Inf, 0.123456, Inf))
  edges <- attr(h, "edges")
  expect_false(is.null(edges))
  expect_true(0.123456 %in% edges)
  expect_lte(data.table::uniqueN(h$g), 50 + 1)
  expect_equal(sum(h$n), 4000); expect_equal(sum(h$e), sum(y))
  # no edge equals an observed score but the forced break, and every cell sits inside one bucket
  expect_false(any(setdiff(edges, 0.123456) %in% s))
  expect_identical(findInterval(h$s_lo, edges), findInterval(h$s_hi, edges))
  expect_identical(findInterval(h$s_lo, edges) + 1L, h$g)
  # every row lies in the cell of its bucket
  gi <- findInterval(s, edges) + 1L
  expect_equal(as.numeric(tapply(rep(1, 4000), list(gi), sum)),
               h[, list(n = sum(n)), keyby = "g"]$n)
  # the shares of the buckets are even (equal weighted share, up to one score)
  expect_lt(max(h[, list(n = sum(n)), keyby = "g"]$n), 4000 / 50 * 2)
  # the cuts of a pooled table are bucket edges
  cells <- .study_cells(h, "a")
  cs <- .study_cuts(cells, 10, "uniform", NULL, "high")
  expect_true(all(cs$cuts %in% edges))
})

test_that("pre-aggregated counts give the same table as the rows", {
  d <- study_df(1500)
  d <- d[!is.na(d$y), ]
  agg <- stats::aggregate(cbind(n = 1, events = d$y, value = d$amt, value_events = d$amt * d$y),
                          by = list(score = d$score, smp = d$smp), FUN = sum)
  bc <- scr_bands(agg, counts = TRUE, sample = "smp", value = "value", value_events = "value_events",
                  n_bands = 8, n_boot = 0)
  br <- scr_bands(d, sample = "smp", value = "amt", n_bands = 8, n_boot = 0)
  expect_equal(bc$cuts, br$cuts)
  expect_equal(bc$table, br$table)
  expect_equal(bc$summary, br$summary)
  expect_error(scr_bands(transform(agg, events = n + 1), counts = TRUE, n_boot = 0), "events <= n")
})

# -- cut points ------------------------------------------------------------------ #

test_that(".study_cuts never cuts on an observed score nor inside a tie", {
  d <- study_df(2000)
  h <- .study_hist(d$score, d$y)
  for (side in c("high", "low")) {
    cs <- .study_cuts(h, 10, "uniform", NULL, side)
    expect_false(any(cs$cuts %in% d$score))
    u <- sort(unique(d$score))
    # each cut sits strictly between two adjacent distinct scores
    for (ct in cs$cuts) {
      lo <- max(u[u < ct]); hi <- min(u[u > ct])
      expect_identical(which(u == lo) + 1L, which(u == hi))
    }
    expect_equal(cs$n_bands_effective, length(cs$cuts) + 1L)
  }
  # heavy ties: fewer bands than requested, and both counts reported
  tied <- .study_hist(rep(c(1, 2, 3), c(500, 30, 470)), stats::rbinom(1000, 1, 0.2))
  ct <- .study_cuts(tied, 10, "uniform", NULL, "high")
  expect_equal(ct$n_bands_requested, 10L)
  expect_lt(ct$n_bands_effective, 10L)
  expect_equal(ct$cuts, c(1.5, 2.5))
  # all tied: a single band
  one <- .study_hist(rep(7, 100), rep(0:1, 50))
  expect_identical(.study_cuts(one, 10, "uniform", NULL, "low")$n_bands_effective, 1L)
  expect_length(.study_cuts(one, 10, "uniform", NULL, "low")$cuts, 0)
})

test_that("tail spacing hits the requested shares within one cell, from the event-rich side", {
  set.seed(9)
  s <- round(stats::rnorm(20000) * 100)
  h <- .study_hist(s, stats::rbinom(20000, 1, 0.1))
  tp <- c(0.001, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50)
  for (side in c("high", "low")) {
    cs <- .study_cuts(h, 20, "tail", NULL, side)
    expect_equal(cs$n_bands_requested, length(tp) + 1L)
    share <- if (side == "high") vapply(cs$cuts, function(ct) mean(s >= ct), 1) else
      vapply(cs$cuts, function(ct) mean(s < ct), 1)
    cell <- max(h$n) / sum(h$n)
    for (t in tp) expect_lte(min(abs(share - t)), cell + 1e-12)
  }
})

# -- band table ------------------------------------------------------------------ #

test_that("the band table matches independent formulas", {
  d <- study_df(4000)
  b <- scr_bands(d, sample = "smp", n_bands = 10, n_boot = 0)
  for (nm in c("dev", "oot")) {
    t <- b$table[sample == nm]
    expect_equal(sum(t$pct), 1)
    expect_equal(t$cum_pct[nrow(t)], 1)
    expect_equal(t$capture[nrow(t)], 1)
    expect_equal(t$cum_nonevent_pct[nrow(t)], 1)
    expect_true(all(t$rate_lo <= t$rate & t$rate <= t$rate_hi))
    # Jeffreys bounds from qbeta
    expect_equal(t$rate_lo, stats::qbeta(0.025, t$events + 0.5, t$n - t$events + 0.5))
    expect_equal(t$rate_hi, stats::qbeta(0.975, t$events + 0.5, t$n - t$events + 0.5))
    bw <- .band_woe(t$events, t$n - t$events)
    expect_equal(t$woe, bw$log_odds)
    expect_equal(t$iv, (bw$pct_event - bw$pct_nonevent) * bw$log_odds)
    R <- sum(t$events) / sum(t$n)
    expect_equal(t$lift, t$rate / R)
    expect_equal(t$cum_lift, cumsum(t$events) / cumsum(t$n) / R)
    expect_equal(t$ks, abs(cumsum(t$events) / sum(t$events) - cumsum(t$n - t$events) / sum(t$n - t$events)))
    # rows: the event-richest band first (low scores under higher_is_safer)
    expect_true(all(diff(t$score_lo) > 0))
    expect_equal(t$odds, (t$n - t$events + 0.5) / (t$events + 0.5))
    # one-sided Fisher exact test of each band against the previous one
    for (i in 2:nrow(t)) {
      ft <- stats::fisher.test(matrix(c(t$events[i], t$n[i] - t$events[i], t$events[i - 1], t$n[i - 1] - t$events[i - 1]),
                                      2, byrow = TRUE), alternative = "greater")$p.value
      expect_equal(t$p_reversal[i], ft, tolerance = 1e-10)
    }
    expect_equal(t$p_reversal_adj, .study_holm(t$p_reversal))
    expect_true(is.na(t$p_reversal[1]))
    # counts by left-closed interval
    sc <- d$score[d$smp == nm]
    expect_equal(t$n, as.numeric(tabulate(findInterval(sc, b$cuts) + 1L, length(b$cuts) + 1L)))
  }
  # the summary PSI is the PSI of scr_psi() over the same bands (no score falls on a cut)
  ps <- scr_psi(d$score[d$smp == "dev"], d$score[d$smp == "oot"], breaks = c(-Inf, b$cuts, Inf))
  expect_equal(b$summary[sample == "oot", psi], ps$psi)
  expect_equal(sum(b$table[sample == "oot", psi]), ps$psi)
  expect_true(all(is.na(b$table[sample == "dev", psi])))
  expect_equal(b$summary$reversals, vapply(c("dev", "oot"), function(nm) sum(b$table[sample == nm, p_reversal_adj] < 0.05, na.rm = TRUE), 1L, USE.NAMES = FALSE))
})

test_that("weights act through the Kish effective size, value through the captured value", {
  d <- study_df(3000)
  b <- scr_bands(d, weight = "w", value = "amt", n_bands = 5, n_boot = 0)
  t <- b$table
  idx <- findInterval(d$score, b$cuts) + 1L
  B <- length(b$cuts) + 1L
  o <- seq_len(B)   # higher_is_safer: ascending intervals, event-richest first
  wsum <- as.numeric(tapply(d$w, factor(idx, levels = seq_len(B)), sum))[o]
  esum <- as.numeric(tapply(d$w * d$y, factor(idx, levels = seq_len(B)), sum))[o]
  w2 <- as.numeric(tapply(d$w^2, factor(idx, levels = seq_len(B)), sum))[o]
  expect_equal(t$n, wsum); expect_equal(t$events, esum)
  neff <- wsum^2 / w2
  x <- esum / wsum * neff
  expect_equal(t$rate_lo, stats::qbeta(0.025, x + 0.5, neff - x + 0.5))
  vs <- as.numeric(tapply(d$w * d$amt, factor(idx, levels = seq_len(B)), sum))[o]
  ves <- as.numeric(tapply(d$w * d$amt * d$y, factor(idx, levels = seq_len(B)), sum))[o]
  expect_equal(t$value, vs); expect_equal(t$value_events, ves)
  expect_equal(t$value_capture, cumsum(ves) / sum(ves))
  expect_equal(t$value_precision, cumsum(ves) / cumsum(vs))
})

# -- count bootstrap -------------------------------------------------------------- #

test_that("the count bootstrap has the point estimate and the law of the row bootstrap", {
  set.seed(21)
  n <- 700
  s <- round(stats::rnorm(n) * 4)
  y <- stats::rbinom(n, 1, stats::plogis(-1 + 0.3 * s))
  m <- scr_metrics(s, y, n_boot = 2000, seed = 4)
  idx <- data.table::frank(s, ties.method = "dense"); K <- max(idx)
  c1 <- tabulate(idx[y == 1], K); c0 <- tabulate(idx[y == 0], K)
  a <- .study_auc_boot(c1, c0, 2000, 0.95, seed = 4, keep = TRUE)
  expect_equal(a$auc, m$auc)
  expect_equal(a$ks, m$ks)
  expect_equal(a$auc, .auc_ks_counts(c1, c0)$auc)
  # the row bootstrap, stratified by outcome, as in scr_metrics()
  set.seed(8)
  i1 <- idx[y == 1]; i0 <- idx[y == 0]
  rb <- vapply(seq_len(2000), function(b) {
    .auc_ks_counts(tabulate(i1[sample.int(length(i1), replace = TRUE)], K),
                   tabulate(i0[sample.int(length(i0), replace = TRUE)], K))$auc
  }, numeric(1))
  se <- stats::sd(rb)
  expect_lt(abs(mean(a$boot_auc) - mean(rb)), 4 * se * sqrt(2 / 2000))
  expect_lt(abs(stats::sd(a$boot_auc) / se - 1), 0.08)
  expect_equal(a$auc_lo, m$auc_lo, tolerance = 0.01)
  expect_equal(a$auc_hi, m$auc_hi, tolerance = 0.01)
  expect_equal(a$gini_lo, 2 * a$auc_lo - 1)
  # reproducible with the seed, and the user's stream untouched
  b <- .study_auc_boot(c1, c0, 200, 0.95, seed = 4, keep = TRUE)
  expect_identical(b$boot_auc, .study_auc_boot(c1, c0, 200, 0.95, seed = 4, keep = TRUE)$boot_auc)
  set.seed(77); r1 <- stats::runif(3)
  set.seed(77); invisible(.study_auc_boot(c1, c0, 50, 0.95, seed = 1)); r2 <- stats::runif(3)
  expect_identical(r1, r2)
  set.seed(77); invisible(scr_bands(data.frame(score = s, y = y), n_boot = 30, seed = 2)); r3 <- stats::runif(3)
  expect_identical(r1, r3)
  # pooled resamples for many cells: same law, shifted to the full estimate
  p <- .study_auc_boot(c1, c0, 2000, 0.95, seed = 4, keep = TRUE, boot_cells = 10)
  expect_equal(p$auc, a$auc)
  expect_lt(abs(mean(p$boot_auc) - mean(a$boot_auc)), 4 * se * sqrt(2 / 2000))
  expect_lt(abs(stats::sd(p$boot_auc) / stats::sd(a$boot_auc) - 1), 0.1)
  # one class: nothing to estimate
  expect_true(is.na(.study_auc_boot(c1, 0 * c0, 100)$auc))
})

test_that("the bootstrap pooled above boot_cells stays within Monte Carlo error of the exact one", {
  skip_on_cran()   # the exact resamples over 2e5 score values take a few seconds
  set.seed(2024)
  n <- 2e5
  s <- stats::rnorm(n)
  y <- stats::rbinom(n, 1, stats::plogis(-1.5 + s))
  h <- .study_hist(s, y, max_cells = 1e6)
  expect_gt(nrow(h), 1.9e5)
  c1 <- h$e; c0 <- h$n_y - h$e
  ex <- .study_auc_boot(c1, c0, 200, 0.95, seed = 1, keep = TRUE, boot_cells = Inf)
  po <- .study_auc_boot(c1, c0, 200, 0.95, seed = 1, keep = TRUE)
  sig <- stats::sd(ex$boot_auc)
  expect_identical(po$auc, ex$auc)
  # The two runs draw independent resamples (different numbers of cells), so
  # their bounds differ by Monte Carlo error. Under normality, a 2.5% (or
  # 97.5%) percentile of B = 200 resamples has standard error
  # sqrt(p (1 - p) / B) / dnorm(1.96) = 0.189 sd; the difference of two
  # independent runs 0.267 sd. Four of those standard errors, 1.07 sd, is
  # the tolerance; the pooling bias itself is about 1e-7 here.
  tol <- 4 * sqrt(2 * 0.025 * 0.975 / 200) / stats::dnorm(stats::qnorm(0.975)) * sig
  expect_lt(abs(po$auc_lo - ex$auc_lo), tol)
  expect_lt(abs(po$auc_hi - ex$auc_hi), tol)
  # the spread: the sd of 200 resamples has a relative error of
  # 1 / sqrt(2 (B - 1)), so two runs differ by 0.071; four of those
  expect_lt(abs(stats::sd(po$boot_auc) / sig - 1), 4 * sqrt(2) / sqrt(2 * 199))
  # the centers: the mean difference has sd sqrt(2 / B) sd
  expect_lt(abs(mean(po$boot_auc) - mean(ex$boot_auc)), 4 * sqrt(2 / 200) * sig)
})

test_that("boot_cells reaches the bootstrap of scr_bands and scr_rag", {
  d <- study_df(2000)
  b0 <- scr_bands(d, n_bands = 5, n_boot = 50, seed = 3)
  expect_equal(b0$boot_cells, 1e4)
  # fewer score values than boot_cells: exact, the same as no pooling
  expect_identical(scr_bands(d, n_bands = 5, n_boot = 50, seed = 3, boot_cells = Inf)$summary, b0$summary)
  bp <- scr_bands(d, n_bands = 5, n_boot = 50, seed = 3, boot_cells = 5)
  expect_equal(bp$summary$auc, b0$summary$auc)
  expect_false(identical(bp$summary$auc_lo, b0$summary$auc_lo))
  expect_error(scr_bands(d, n_boot = 0, boot_cells = 1), "boot_cells")
  expect_error(scr_bands(d, n_boot = 0, boot_cells = 2.5), "boot_cells")
  r0 <- scr_rag(d, sample = "smp", n_boot = 50, seed = 3, boot_cells = Inf)
  expect_identical(r0$table, scr_rag(d, sample = "smp", n_boot = 50, seed = 3)$table)
  expect_false(identical(scr_rag(d, sample = "smp", n_boot = 50, seed = 3, boot_cells = 5)$table[metric == "gini_ratio", lo],
                         r0$table[metric == "gini_ratio", lo]))
})

test_that("the DeLong standard error from counts equals the row estimator", {
  set.seed(31)
  s <- round(stats::rnorm(500) * 3); y <- stats::rbinom(500, 1, stats::plogis(0.4 * s))
  idx <- data.table::frank(s, ties.method = "dense"); K <- max(idx)
  expect_equal(.study_delong_counts(tabulate(idx[y == 1], K), tabulate(idx[y == 0], K)), ref_auc_se(s, y), tolerance = 1e-12)
  expect_true(is.na(.study_delong_counts(c(1, 0), c(5, 5))))
})

# -- scr_bands -------------------------------------------------------------------- #

test_that("scr_bands on a scorecard: frozen on train, read on hold-out, config keys respected", {
  sc <- sc_demo()
  b <- scr_bands(sc, n_bands = 10, n_boot = 20, seed = 1)
  expect_s3_class(b, c("scr_study_bands", "scr_study"))
  expect_identical(unique(b$table$sample), c("train", "holdout"))
  expect_equal(b$summary$auc, sc$metrics$auc)
  expect_equal(b$summary$ks, sc$metrics$ks)
  expect_equal(b$summary$n, c(nrow(sc$samples$train), nrow(sc$samples$holdout)))
  # the same study from a data.frame of the scored samples
  d <- rbind(data.frame(smp = "train", sc$samples$train[, c("score", "y")]),
             data.frame(smp = "holdout", sc$samples$holdout[, c("score", "y")]))
  d$smp <- factor(d$smp, levels = c("train", "holdout"))
  bd <- scr_bands(d, sample = "smp", n_bands = 10, n_boot = 20, seed = 1)
  expect_equal(bd$cuts, b$cuts)
  expect_equal(bd$table[, -"sample"], b$table[, -"sample"])
  # config keys: study_bands and study_level when the argument is NULL
  sc2 <- sc; sc2$config$study_bands <- 4L; sc2$config$study_level <- 0.9
  b2 <- scr_bands(sc2, n_boot = 0)
  expect_equal(b2$n_bands_requested, 4L)
  expect_equal(b2$level, 0.9)
  # an older configuration without the keys falls back on the defaults
  sc3 <- sc; sc3$config$study_bands <- NULL
  expect_equal(scr_bands(sc3, n_boot = 0)$n_bands_requested, 20L)
  # breaks of the gains: left-closed here, so counts follow findInterval
  bg <- scr_bands(sc, breaks = sc$breaks, n_boot = 0)
  ho <- sc$samples$holdout$score
  expect_equal(bg$table[sample == "holdout", n],
               as.numeric(tabulate(findInterval(ho, sc$breaks[is.finite(sc$breaks)]) + 1L, length(sc$breaks) - 1L)))
  expect_output(print(b), "scr_study_bands")
  # the rate is printed with its interval, and the p-value column is the Holm-adjusted one
  expect_output(print(b), "rate \\[lo, hi\\]")
  expect_output(print(b), "p_rev_adj")
  hb <- b$table[sample == "holdout"]
  expect_output(print(b), sprintf("%.2f%% [%.2f%%, %.2f%%]", 100 * hb$rate[1], 100 * hb$rate_lo[1], 100 * hb$rate_hi[1]),
                fixed = TRUE)
  expect_output(print(b), sprintf("%.3f", hb$p_reversal_adj[2]), fixed = TRUE)
  expect_error(scr_bands(sc, sample = "oot"), "not in the scorecard")
  expect_error(scr_bands(sc, n_bandz = 3), "unused argument")
})

test_that("scr_bands under propensity orders from the high scores", {
  d <- study_df(2000)
  d$score <- -d$score
  b <- scr_bands(d, objective = "propensity", n_bands = 5, n_boot = 0)
  expect_identical(b$direction, "higher_is_riskier")
  t <- b$table
  expect_true(all(diff(t$score_lo) < 0))
  expect_gt(t$rate[1], t$rate[nrow(t)])
  expect_equal(t$odds, (t$events + 0.5) / (t$n - t$events + 0.5))
  expect_equal(b$summary$auc, scr_metrics(d$score, d$y, higher_is_event = TRUE, ci = FALSE)$auc)
  expect_identical(b$codes, rev(seq_len(length(b$cuts) + 1L)))
})

# -- edge cases ------------------------------------------------------------------- #

test_that("edge cases: single class, one band, all ties, zero weights, empty study, small n", {
  # single class: AUC and the shares undefined, the table still complete
  one_class <- data.frame(score = 1:100, y = 0)
  b <- scr_bands(one_class, n_bands = 4, n_boot = 10)
  expect_true(is.na(b$summary$auc)); expect_true(is.na(b$summary$iv))
  expect_equal(b$table$rate, rep(0, 4))
  expect_true(all(is.na(b$table$pct_event)))
  # one band
  d <- study_df(300)
  b1 <- scr_bands(d, n_bands = 1, n_boot = 0)
  expect_equal(nrow(b1$table), 1L)
  expect_identical(b1$table$label, "[-Inf, Inf)")
  expect_equal(b1$table$ks, 0)
  # all ties
  bt <- scr_bands(data.frame(score = rep(5, 40), y = rep(0:1, 20)), n_boot = 10, seed = 1)
  expect_equal(bt$summary$n_bands_effective, 1L)
  expect_equal(bt$summary$auc, 0.5)
  # zero weights leave the rows out
  zw <- data.frame(score = 1:100, y = rep(0:1, 50), w = rep(c(0, 1), each = 50))
  bz <- scr_bands(zw, weight = "w", n_bands = 5, n_boot = 0)
  expect_equal(sum(bz$table$n), 50)
  expect_true(all(bz$cuts > 50))
  expect_error(scr_bands(transform(zw, w = 0), weight = "w", n_boot = 0), "no scored row")
  # empty study sample: rows with a missing score only
  emp <- data.frame(score = c(1:50, rep(NA, 10)), y = rep(0:1, 30), smp = rep(c("a", "b"), c(50, 10)))
  be <- expect_silent(scr_bands(emp, sample = "smp", n_bands = 5, n_boot = 10))
  expect_equal(be$summary[sample == "b", n], 0)
  expect_true(is.na(be$summary[sample == "b", auc]))
  expect_equal(be$table[sample == "b", n], rep(0, 5))
  # missing outcomes count in the volume, not in the rates
  na_y <- data.frame(score = 1:20, y = c(rep(NA, 5), rep(0:1, 7), 1))
  bn <- scr_bands(na_y, n_bands = 2, n_boot = 0)
  expect_equal(sum(bn$table$n), 20)
  expect_equal(sum(bn$table$events), 8)
  expect_equal(bn$summary$rate, 8 / 15)
  # infinite scores are left out, like missing ones
  bi <- scr_bands(data.frame(score = c(1:10, Inf, -Inf, NaN), y = rep(0:1, length.out = 13)), n_bands = 2, n_boot = 0)
  expect_equal(sum(bi$table$n), 10)
  expect_true(all(is.finite(bi$cuts)))
  # small n
  bs <- scr_bands(data.frame(score = c(1, 2, 3), y = c(1, 0, 0)), n_bands = 3, n_boot = 5, seed = 1)
  expect_equal(bs$summary$auc, 1)
  expect_error(scr_bands(data.frame(s = 1), n_boot = 0), "not in `x`")
  expect_error(scr_bands(d, objective = "fraud"), "objective")
  expect_error(scr_bands(d, n_bands = 0), "n_bands")
  expect_error(scr_bands(d, reference = "dev"), "need a `sample` column")
})

# -- production ------------------------------------------------------------------- #

test_that("scr_apply and scr_sql assign the bands identically (duckdb)", {
  d <- study_df(2000)
  b <- scr_bands(d, sample = "smp", n_bands = 8, n_boot = 0)
  a <- scr_apply(b, d)
  expect_identical(a$tier, b$codes[findInterval(d$score, b$cuts) + 1L])
  expect_identical(a$tier_label, b$code_labels[findInterval(d$score, b$cuts) + 1L])
  # band numbers agree with the table: the band of every row has the row's label
  t <- b$table[sample == "dev"]
  expect_identical(t$label[match(a$tier, t$band)], a$tier_label)
  expect_false("tier" %in% names(d))   # the input is not modified
  expect_identical(scr_apply(b, c(500, NA))$tier, c(1L, NA))
  skip_if_not_installed("duckdb"); skip_if_not_installed("DBI")
  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  d$score[c(4, 9)] <- NA
  d$score[10] <- b$cuts[3]   # a score on a cut goes to the upper side in both
  DBI::dbWriteTable(con, "sc tab", d)
  sql <- scr_sql(b, table = "\"sc tab\"", dialect = "duckdb")
  got <- DBI::dbGetQuery(con, paste(sql, collapse = "\n"))
  exp <- scr_apply(b, d)
  expect_identical(as.integer(got$tier), exp$tier)
  expect_identical(got$tier_label, exp$tier_label)
})

test_that("scr_sql refuses an unknown dialect, as for a scorecard", {
  b <- scr_bands(study_df(500), n_bands = 4, n_boot = 0)
  e_study <- tryCatch(scr_sql(b, dialect = "foo"), error = conditionMessage)
  e_sc <- tryCatch(scr_sql(sc_demo(), dialect = "foo"), error = conditionMessage)
  expect_identical(e_study, e_sc)
  expect_match(e_study, "should be one of")
  for (dl in c("ansi", "postgres", "mysql", "mariadb", "sqlserver", "oracle", "spark", "hive", "databricks",
               "bigquery", "snowflake", "redshift", "duckdb", "sqlite")) {
    expect_true(any(grepl("CASE WHEN", scr_sql(b, dialect = dl))))
  }
})

test_that("scr_export writes one workbook with the study tables", {
  skip_if_not_installed("openxlsx")
  old <- scr_verbose(FALSE); on.exit(scr_verbose(old), add = TRUE)
  b <- scr_bands(study_df(1000), sample = "smp", n_bands = 5, n_boot = 0)
  out <- file.path(tempdir(), "scr-study-export")
  unlink(out, recursive = TRUE)
  ex <- scr_export(b, out, stamp = FALSE)
  expect_true(file.exists(ex$files$xlsx))
  expect_identical(openxlsx::getSheetNames(ex$files$xlsx), c("Summary", "Bands", "Cuts", "Settings"))
  bands <- openxlsx::read.xlsx(ex$files$xlsx, sheet = "Bands")
  expect_equal(nrow(bands), nrow(b$table))
  expect_equal(bands$n, b$table$n)
})

test_that("the score-study configuration keys are registered and validated", {
  k <- scr_config_keys(stage = 13)
  expect_setequal(k$key, c("study_bands", "study_level", "tier_min_pct", "tier_min_events", "tier_max_bins"))
  cfg <- scr_config(verbose = FALSE)
  expect_equal(c(cfg$study_bands, cfg$study_level, cfg$tier_min_pct, cfg$tier_min_events, cfg$tier_max_bins),
               c(20, 0.95, 0.05, 20, 100))
  expect_error(scr_config(study_level = 1), "study_level")
  expect_error(scr_config(study_bands = 0), "study_bands")
  expect_error(scr_config(tier_min_pct = 0.6), "tier_min_pct")
  expect_error(scr_config(tier_min_events = -1), "tier_min_events")
  expect_error(scr_config(tier_max_bins = 1000), "tier_max_bins")
  # a configuration saved before the keys existed is completed with the defaults
  old <- unclass(cfg)[setdiff(names(cfg), k$key)]
  expect_equal(.scr_validate_config(old)$tier_max_bins, 100L)
})

test_that("a numeric sample column is read in numeric order", {
  d <- study_df(2000, 15)
  d$k <- sample(c(9, 10, 100), nrow(d), TRUE)
  b <- scr_bands(d, sample = "k", n_bands = 5, n_boot = 0)
  # as text "10" would come before "9" and become the reference
  expect_identical(b$reference, "9"); expect_identical(b$samples, c("9", "10", "100"))
  expect_equal(b$cuts, scr_bands(d[d$k == 9, ], n_bands = 5, n_boot = 0)$cuts)
  # text and factors keep their order
  d$kc <- as.character(d$k)
  expect_identical(scr_bands(d, sample = "kc", n_bands = 5, n_boot = 0)$reference, "10")
  d$kf <- factor(d$k, levels = c(100, 9, 10))
  expect_identical(scr_bands(d, sample = "kf", n_bands = 5, n_boot = 0)$samples, c("100", "9", "10"))
  expect_identical(.study_levels(c(2.5, 10, 2.5, NA)), c("2.5", "10"))
  expect_identical(.study_levels(as.Date(c("2026-02-01", "2026-01-01"))), c("2026-01-01", "2026-02-01"))
  expect_null(.study_num_levels(list(c("a", "b")))); expect_null(.study_num_levels(list(factor(1:3))))
  expect_identical(.study_num_levels(list(c(10L, 9L), 100L)), c("9", "10", "100"))
})
