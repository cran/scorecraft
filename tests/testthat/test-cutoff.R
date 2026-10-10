test_that("the cut-off sweep freezes the cuts on train and behaves monotonically", {
  sc <- sc_demo()
  ct <- scr_cutoff(sc, n_cuts = 10)
  expect_s3_class(ct, "scr_cutoff")
  tb <- ct$table
  expect_identical(tb[sample == "train", cut], tb[sample == "holdout", cut])
  h <- tb[sample == "holdout"]
  expect_true(all(diff(h$pct_safe) <= 0))            # higher cut, fewer approved
  expect_true(all(diff(h$events_avoided_pct) >= 0))  # higher cut, more events avoided
  expect_true(all(h$ks_at_cut >= 0 & h$ks_at_cut <= 1))
  expect_equal(max(h$ks_at_cut), scr_score_metrics(sc)[sample == "holdout", ks], tolerance = 0.06)
  ct2 <- scr_cutoff(sc, cuts = c(520, 560))
  expect_equal(unique(ct2$table$cut), c(520, 560))
  expect_output(print(ct), "frozen on train")
})

test_that("the strategy table exposes the marginal expected profit and break-even", {
  sc <- sc_demo()
  st <- scr_strategy(sc, revenue_good = 1080, loss_bad = 4500)
  expect_equal(st$breakeven, 1080 / (1080 + 4500))
  d <- st$table
  expect_equal(d$ep_per_account, (1 - d$event_rate) * 1080 - d$event_rate * 4500)
  expect_equal(d$band_profit, d$n * d$ep_per_account)
  expect_true(all(d$decision[d$event_rate <= st$breakeven] == "approve"))
  expect_true(all(d$decision[d$event_rate > 1.25 * st$breakeven] == "decline"))
  expect_equal(d$cum_pct[nrow(d)], 1)
  expect_lt(d$event_rate[1], d$event_rate[nrow(d)])   # safest band first
  st2 <- scr_strategy(sc, decisions = rep("approve", nrow(d)))
  expect_true(all(st2$table$decision == "approve"))
  expect_error(scr_strategy(sc, decisions = "approve"), "one decision per band")
  expect_error(scr_strategy(sc, sample = "nope"), "does not exist")
  expect_output(print(st), "break-even")
})

test_that(".band_woe() gives exact shares and the event-oriented band WOE", {
  n  <- c(2675, 1531, 1710, 2011, 2231, 1901, 1782, 2176, 1511, 1912)
  ev <- c(148, 144, 219, 345, 522, 659, 638, 1098, 880, 1179)
  bw <- .band_woe(ev, n - ev)
  pinned <- c(-1.9903, -1.4178, -1.0708, -0.7273, -0.3387, 0.2135, 0.2633, 0.8657, 1.1799, 1.3226)
  expect_equal(round(bw$log_odds, 4), pinned)
  expect_lt(max(abs(bw$log_odds - pinned)), 1e-4)
  expect_equal(sum(bw$pct_event), 1)
  expect_equal(sum(bw$pct_nonevent), 1)
  expect_equal(bw$log_odds, log(bw$pct_event / bw$pct_nonevent))   # no zero: exact ratio
  expect_equal(bw$log_odds, log(bw$odds_event))
  # log_odds > 0 <=> band event rate above the overall rate <=> lift > 1
  expect_equal(sign(bw$log_odds), sign(ev / n - sum(ev) / sum(n)))
  D <- abs(cumsum(bw$pct_event) - cumsum(bw$pct_nonevent))
  expect_identical(which.max(D), 5L)
  expect_equal(round(max(D), 4), 0.4089)
})

test_that(".band_woe() smooths only when a band lacks a class and handles degenerate counts", {
  # a band without events: finite, smoothed ratio; the shares stay exact
  z <- .band_woe(c(0, 5, 10), c(10, 10, 10))
  expect_true(all(is.finite(z$log_odds)))
  expect_equal(z$pct_event, c(0, 5, 10) / 15)
  expect_equal(z$pct_nonevent, rep(1 / 3, 3))
  expect_equal(z$odds_event, ((c(0, 5, 10) + 0.5) / 16.5) / ((c(10, 10, 10) + 0.5) / 31.5))
  # a single class: that side and every ratio are NA
  s1 <- .band_woe(c(0, 0), c(3, 4))
  expect_true(all(is.na(s1$pct_event)))
  expect_true(all(is.na(s1$odds_event)) && all(is.na(s1$log_odds)))
  expect_equal(s1$pct_nonevent, c(3, 4) / 7)
  s0 <- .band_woe(c(0, 0), c(0, 0))
  expect_true(all(is.na(unlist(s0))))
  # an empty band carries no information and does not trigger the smoothing
  eb <- .band_woe(c(2, 0, 3), c(3, 0, 5))
  expect_true(is.na(eb$log_odds[2]))
  expect_equal(eb$log_odds[-2], log(c(2, 3) / 5 / (c(3, 5) / 8)))
  expect_equal(eb$pct_event, c(2, 0, 3) / 5)
  # a band with an unknown count is unknown on both sides and left out of the totals
  na <- .band_woe(c(2, NA, 3), c(3, 4, Inf))
  expect_true(all(is.na(c(na$pct_event[2:3], na$pct_nonevent[2:3], na$log_odds[2:3]))))
  expect_equal(na$pct_event[1], 1)
  expect_length(.band_woe(numeric(0), numeric(0))$log_odds, 0L)
})

test_that("the strategy table carries the event and non-event distributions", {
  sc <- sc_demo()
  st <- scr_strategy(sc, revenue_good = 1080, loss_bad = 4500)
  d <- st$table
  expect_identical(names(d), c("id", "band", "min_score", "max_score", "n", "pct", "events", "event_rate",
                               "pct_event", "pct_nonevent", "odds_event", "log_odds", "decision",
                               "ep_per_account", "band_profit", "cum_pct", "cum_event_rate", "cum_profit"))
  expect_equal(sum(d$pct_event), 1)
  expect_equal(sum(d$pct_nonevent), 1)
  expect_equal(sign(d$log_odds), sign(d$event_rate - sum(d$events) / sum(d$n)))
  expect_equal(d$log_odds, log(d$odds_event))
  # the same band WOE as the gains of the same sample and bands
  g <- scr_score_gains(sc, "holdout")
  expect_equal(d$log_odds, g$woe[match(d$band, g$band)])
  expect_identical(st$objective, "risk")
  expect_identical(st$rule, "breakeven")
  expect_output(print(st), "log_odds")
})

test_that("a fraud-like scorecard puts the low (safe) scores first and approves", {
  sf <- sc_fraud_demo()
  st <- scr_strategy(sf, revenue_good = 1080, loss_bad = 4500)
  d <- st$table
  expect_false(is.unsorted(d$min_score))
  expect_lt(d$event_rate[1], d$event_rate[nrow(d)])
  expect_true(all(d$decision %in% c("approve", "review", "decline")))
  expect_true(all(d$decision[d$event_rate <= st$breakeven] == "approve"))
  expect_equal(st$breakeven, 1080 / (1080 + 4500))
  expect_equal(d$ep_per_account, (1 - d$event_rate) * 1080 - d$event_rate * 4500)
  expect_equal(st$crossing$ks, max(scr_score_gains(sf, "holdout")$ks))
})

test_that("a propensity scorecard targets the most likely bands first", {
  sp <- sc_prop_demo()
  expect_identical(sp$config$objective, "propensity")
  st <- scr_strategy(sp, revenue_good = 100, loss_bad = 60)
  d <- st$table
  expect_identical(st$objective, "propensity")
  expect_false(is.unsorted(rev(d$min_score)))   # higher_is_riskier: high score first
  expect_gt(d$event_rate[1], d$event_rate[nrow(d)])
  expect_true(all(d$decision %in% c("target", "review", "skip")))
  expect_equal(d$ep_per_account, d$event_rate * 100 - (1 - d$event_rate) * 60)
  expect_equal(st$breakeven, 60 / (100 + 60))
  expect_true(all(d$decision[d$event_rate >= st$breakeven] == "target"))
  expect_true(all(d$decision[d$event_rate < st$breakeven] != "target"))
  expect_true(all(d$decision[1 - d$event_rate > 1.25 * 100 / 160] == "skip"))
  expect_equal(sum(d$pct_event), 1)
  expect_output(print(st), "target at or above")
})

# band edge of the crossing (the upper edge of the lower of the two rows
# around it) and the training scores on either side of it
crossing_edge <- function(st, breaks) {
  d <- st$table; k <- match(st$crossing$after_band, d$band)
  lower <- if (d$max_score[k] < d$min_score[k + 1L]) k else k + 1L
  min(breaks[breaks >= d$max_score[lower]])
}
frozen_cut <- function(x, edge) {
  tr <- x$samples$train$score
  lo <- max(tr[tr <= edge]); hi <- min(tr[tr > edge])
  list(lo = lo, hi = hi, cut = (lo + hi) / 2)
}

test_that("the crossing rule cuts at the KS boundary with contiguous decisions", {
  sc <- sc_demo()
  st <- scr_strategy(sc, rule = "crossing")
  d <- st$table; cr <- st$crossing
  k <- match(cr$after_band, d$band)
  expect_identical(d$decision, rep(c("approve", "decline"), c(k, nrow(d) - k)))
  expect_equal(cr$ks, max(scr_score_gains(sc, "holdout")$ks))
  expect_true(is.finite(cr$cut))
  # higher_is_safer: the approved bands sit above the band edge, the declined ones at or below it;
  # the cut is frozen on the training scores on either side of that edge
  e <- crossing_edge(st, sc$breaks); f <- frozen_cut(sc, e)
  expect_true(all(d$min_score[seq_len(k)] > e))
  expect_true(all(d$max_score[-seq_len(k)] <= e))
  expect_equal(cr$cut, f$cut)
  expect_true(f$lo <= e && e < f$hi && f$lo < cr$cut && cr$cut <= f$hi)
  expect_type(cr$single_crossing, "logical")
  # computed under either rule; decisions override the rule
  expect_identical(scr_strategy(sc)$crossing, cr)
  own <- scr_strategy(sc, rule = "crossing", decisions = rep("review", nrow(d)))
  expect_true(all(own$table$decision == "review"))
  expect_output(print(st), "cross at score")
  expect_error(scr_strategy(sc, rule = "nope"), "should be one of")
  expect_error(scr_strategy(sc, revenue_good = 0, loss_bad = 0), "cannot both be 0")

  sp <- sc_prop_demo()
  stp <- scr_strategy(sp, rule = "crossing")
  kp <- match(stp$crossing$after_band, stp$table$band)
  expect_identical(stp$table$decision, rep(c("target", "skip"), c(kp, nrow(stp$table) - kp)))
  expect_equal(stp$crossing$ks, max(scr_score_gains(sp, "holdout")$ks))
  # high score first: the targeted bands sit above the band edge
  ep <- crossing_edge(stp, sp$breaks)
  expect_true(all(stp$table$min_score[seq_len(kp)] > ep))
  expect_equal(stp$crossing$cut, frozen_cut(sp, ep)$cut)
})

test_that("the crossing handles an empty band, several sign changes and degenerate samples", {
  sc <- sc_demo()
  br <- c(-Inf, 10, 20, 30, Inf)
  # the band (10,20] is empty and no training score lies below the edge 10:
  # the cut is midway between the bands on either side, on the hold-out
  s2 <- sc
  s2$samples$holdout <- data.table::data.table(score = rep(c(40, 25, 5), each = 50),
                                               y = c(rep(0, 45), rep(1, 5), rep(0, 40), rep(1, 10), rep(0, 10), rep(1, 40)))
  st <- scr_strategy(s2, breaks = br, rule = "crossing")
  expect_equal(nrow(st$table), 3L)
  expect_equal(st$crossing$cut, 15)
  expect_identical(st$crossing$after_band, "(20,30]")
  expect_true(st$crossing$single_crossing)
  expect_identical(st$table$decision, c("approve", "approve", "decline"))
  # log_odds changes sign three times: still a contiguous cut, at the first maximum
  s3 <- sc
  s3$samples$holdout <- data.table::data.table(score = rep(c(35, 25, 15, 5), each = 50),
                                               y = rep(rep(c(0, 1), 4), c(45, 5, 25, 25, 44, 6, 25, 25)))
  st3 <- scr_strategy(s3, breaks = br, rule = "crossing")
  expect_false(st3$crossing$single_crossing)
  expect_identical(st3$table$decision, c("approve", "decline", "decline", "decline"))
  expect_equal(st3$crossing$cut, 30)
  expect_output(print(st3), "does not change sign exactly once")
  # a single band: no crossing
  one <- scr_strategy(sc, breaks = c(-Inf, Inf))
  expect_equal(nrow(one$table), 1L)
  expect_true(is.na(one$crossing$cut) && is.na(one$crossing$ks) && is.na(one$crossing$single_crossing))
  expect_error(scr_strategy(sc, breaks = c(-Inf, Inf), rule = "crossing"), "two bands")
  expect_output(print(one), "undefined")
  # a single class: shares of the missing side and every ratio are NA, the rule is refused
  s1 <- sc
  s1$samples$holdout <- data.table::copy(sc$samples$holdout)[, y := y * 0]
  st1 <- scr_strategy(s1)
  expect_true(all(is.na(st1$table$pct_event)) && all(is.na(st1$table$log_odds)))
  expect_equal(sum(st1$table$pct_nonevent), 1)
  expect_true(all(st1$table$decision == "approve"))
  expect_true(is.na(st1$crossing$ks))
  expect_error(scr_strategy(s1, rule = "crossing"), "both classes")
})

test_that("the crossing cut is frozen on train and scr_cutoff() at it splits as the crossing rule", {
  lab_good <- c("approve", "target"); lab_bad <- c("decline", "skip")
  # exact agreement on the sample of the strategy: the KS and both sides
  check <- function(x, st) {
    ct <- scr_cutoff(x, cuts = st$crossing$cut)
    h <- ct$table[sample == st$sample]
    expect_equal(h$ks_at_cut, st$crossing$ks, tolerance = 1e-12)
    n_good <- sum(st$table$n[st$table$decision %in% lab_good])
    n_bad <- sum(st$table$n[st$table$decision %in% lab_bad])
    expect_equal(n_good + n_bad, sum(st$table$n))
    # the safe side of scr_cutoff() is the good side under risk and the non-targeted side under propensity
    sides <- if (identical(st$objective, "propensity")) c(n_bad, n_good) else c(n_good, n_bad)
    expect_identical(c(h$n_safe, sum(st$table$n) - h$n_safe), sides)
    # the upper side is score >= cut
    s <- x$samples[[st$sample]]$score
    up <- if (st$table$min_score[1] > st$table$max_score[nrow(st$table)]) n_good else n_bad
    expect_identical(sum(s >= st$crossing$cut), as.integer(up))
  }
  # credit, fraud-like (risk, higher_is_riskier) and propensity
  for (x in list(sc_demo(), sc_fraud_demo(), sc_prop_demo())) {
    tr <- x$samples$train$score; s <- x$samples$holdout$score
    st_tr <- scr_strategy(x, sample = "train", rule = "crossing")
    check(x, st_tr)
    # the hold-out cut depends only on its band edge and the training scores
    st_ho <- scr_strategy(x, rule = "crossing")
    e <- crossing_edge(st_ho, x$breaks); f <- frozen_cut(x, e)
    expect_equal(st_ho$crossing$cut, f$cut)
    if (crossing_edge(st_tr, x$breaks) == e) expect_identical(st_ho$crossing$cut, st_tr$crossing$cut)
    # score >= cut reproduces the band split on train and on every score seen in training
    expect_identical(tr >= st_ho$crossing$cut, tr > e)
    off <- (s >= st_ho$crossing$cut) != (s > e)
    expect_false(any(off & s %in% tr))
    # on the hold-out a row changes side only strictly between the two training scores
    between <- s > f$lo & s < f$hi
    expect_true(all(between[off]))
    expect_lte(sum(off), sum(between))
    expect_identical(sum(s > e), sum(st_ho$table$n[st_ho$table$min_score > e]))
    if (!any(between)) check(x, st_ho)
  }
  # mass points on the band edges: a row at a right-closed edge belongs to the lower band
  s5 <- sc_demo()
  for (nm in c("train", "holdout")) {
    s5$samples[[nm]] <- data.table::data.table(score = rep(c(5, 10, 20, 30), each = 50),
                                               y = rep(rep(c(1, 0), 4), c(40, 10, 30, 20, 8, 42, 4, 46)))
  }
  st5 <- scr_strategy(s5, breaks = c(-Inf, 10, 20, Inf), rule = "crossing")
  expect_identical(st5$table$decision, c("approve", "approve", "decline"))
  expect_equal(st5$crossing$cut, 15)
  check(s5, st5)
  # a cut at the band edge (10) moves the rows scored 10 to the upper side
  expect_false(isTRUE(all.equal(scr_cutoff(s5, cuts = 10)$table[sample == "holdout"]$ks_at_cut, st5$crossing$ks)))
  # hold-out scores never seen in training, strictly between the training scores 10 and 20: the cut
  # stays at 15 (the hold-out midpoint would be 8.5) and those 50 rows fall on the other side
  s7 <- s5
  s7$samples$holdout <- data.table::copy(s5$samples$holdout)[score == 10, score := 12]
  st7 <- scr_strategy(s7, breaks = c(-Inf, 10, 20, Inf), rule = "crossing")
  expect_identical(st7$table$decision, c("approve", "approve", "decline"))
  expect_equal(st7$crossing$cut, 15)
  h7 <- scr_cutoff(s7, cuts = st7$crossing$cut)$table[sample == "holdout"]
  expect_equal(h7$n_safe, sum(st7$table$n[st7$table$decision == "approve"]) - 50L)
  # `breaks` as a number of intervals: the edges come from the sample, and so does the cut
  sc <- sc_demo()
  st6 <- scr_strategy(sc, breaks = 4, rule = "crossing")
  expect_true(is.finite(st6$crossing$cut))
  check(sc, st6)
})

test_that("under the crossing rule, scores outside the breaks get no decision", {
  sc <- sc_demo()
  s4 <- sc
  # the scores at 5 fall below the first break: a last row without band
  s4$samples$holdout <- data.table::data.table(score = rep(c(40, 25, 15, 5), each = 50),
                                               y = rep(rep(c(0, 1), 4), c(45, 5, 40, 10, 10, 40, 20, 30)))
  st <- scr_strategy(s4, breaks = c(10, 20, 30, Inf), rule = "crossing")
  d <- st$table
  expect_equal(nrow(d), 4L)
  expect_true(is.na(d$band[4]))
  expect_identical(d$decision, c("approve", "approve", "decline", NA_character_))
  expect_identical(st$crossing$after_band, "(20,30]")
  expect_equal(st$crossing$cut, 20)
  # the shares, hence the KS, are relative to the whole sample, outside row included
  expect_equal(sum(d$pct_event), 1)
  expect_equal(st$crossing$ks, abs(15 / 85 - 85 / 115))
  # the break-even rule still labels every row
  expect_false(anyNA(scr_strategy(s4, breaks = c(10, 20, 30, Inf))$table$decision))
})

test_that("a strategy object saved before the band WOE columns still prints", {
  st <- scr_strategy(sc_demo(), revenue_good = 1080, loss_bad = 4500)
  old <- st
  old$table <- old$table[, !c("pct_event", "pct_nonevent", "odds_event", "log_odds")]
  old$objective <- NULL; old$rule <- NULL; old$crossing <- NULL
  out <- capture.output(print(old))
  expect_false(any(grepl("log_odds", out)))
  expect_false(any(grepl("\\bNA\\b", out)))
  expect_true(any(grepl("break-even", out)))
  expect_true(any(grepl("objective risk | rule breakeven", out, fixed = TRUE)))
})

test_that("reject inference declares its scope and presents a sensitivity band", {
  sc <- sc_demo()
  rj <- scr_reject(sc)
  expect_s3_class(rj, "scr_reject")
  expect_equal(rj$scope$n_without_outcome, 0L)
  expect_match(rj$scope$statement, "WITH an observed outcome")
  expect_equal(sort(unique(rj$sensitivity$multiplier)), c(2, 4, 8))
  tot <- rj$sensitivity[band == "TOTAL"]
  expect_equal(nrow(tot), 3L)
  expect_equal(tot$rate_implied, tot$rate_dev)   # nothing without outcome: nothing changes
  expect_true(all(rj$coverage$coverage_flag %in% c("ok", "few_events", "no_outcome")))
  expect_output(print(rj), "sensitivity band")
})

test_that("reject inference with a through-the-door population raises the implied rate", {
  sc <- sc_demo()
  pop <- scr_demo
  acc <- seq_len(nrow(pop)) %in% res_demo()$split$holdout_idx
  rj <- scr_reject(sc, population = pop, accepted = acc, multipliers = c(2, 4))
  expect_equal(rj$scope$n_population, nrow(pop))
  expect_equal(rj$scope$n_without_outcome, sum(!acc))
  tot <- rj$sensitivity[band == "TOTAL"]
  expect_gt(tot[multiplier == 4, rate_implied], tot[multiplier == 2, rate_implied])
  expect_gt(tot[multiplier == 2, rate_implied], tot$rate_dev[1])
  expect_true(all(rj$coverage$coverage <= 1, na.rm = TRUE))
  expect_error(scr_reject(sc, population = pop, accepted = TRUE), "length")
})
