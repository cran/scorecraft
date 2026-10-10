test_that("the scorecard is additive: points sum to the score, exact and rounded", {
  sc <- sc_demo()
  expect_s3_class(sc, "scr_scorecard")
  s <- sc$samples$holdout
  expect_equal(s$score, sc$alignment$a + sc$alignment$b * s$link)
  a <- scr_apply(sc, scr_demo[sc$samples$holdout$y >= 0 & seq_len(nrow(scr_demo)) %in% res_demo()$split$holdout_idx, ], what = "points")
  expect_equal(a$score_points, sc$base_points + rowSums(a[, paste0(sc$features, "_points"), with = FALSE]))
  expect_equal(a$score_points, s$score_points)
  # rounded and exact score differ by less than the rounding error budget
  expect_lt(max(abs(a$score - a$score_points)), length(sc$features) * 0.5 + 0.5 + 1)
  expect_output(print(sc), "scr_scorecard")
})

test_that("the sign check leaves only positive coefficients and records removals", {
  sc <- sc_demo()
  b <- sc$coef[sc$features]
  expect_true(all(b > 0))
  expect_true(all(sc$sign_check$action %in% c("kept", "removed")))
  expect_setequal(sc$sign_check[action == "kept", variable], sc$features)
})

test_that("metrics report AUC above 0.5 in both directions, with a CI", {
  sc <- sc_demo()
  m <- scr_score_metrics(sc)
  expect_true(all(m$auc > 0.6))
  expect_true(all(m$auc_lo <= m$auc & m$auc <= m$auc_hi))
  expect_equal(m$direction, rep("higher_is_safer", 2))
  # higher_is_safer: a higher score means fewer events
  h <- sc$samples$holdout
  expect_lt(mean(h$score[h$y == 1]), mean(h$score[h$y == 0]))
  sp <- scr_scorecard(res_demo(), direction = "higher_is_riskier", n_boot = 10)
  hp <- sp$samples$holdout
  expect_gt(mean(hp$score[hp$y == 1]), mean(hp$score[hp$y == 0]))
  expect_true(all(scr_score_metrics(sp)$auc > 0.6))
  expect_equal(sp$odds_orientation, "event:safe")
  expect_equal(scr_score_metrics(sp)$auc, m$auc, tolerance = 1e-10)   # same ranking, mirrored scale
})

test_that("gains use bands frozen on train and run from the risky to the safe side", {
  sc <- sc_demo()
  g <- scr_score_gains(sc)
  expect_setequal(unique(g$sample), c("train", "holdout"))
  expect_identical(g[sample == "train", band], g[sample == "holdout", band])
  tr <- g[sample == "train"]
  expect_true(all(diff(tr$cum_pct) > 0))
  expect_gt(tr$event_rate[1], tr$event_rate[nrow(tr)])
  expect_equal(max(tr$ks), scr_score_metrics(sc)[sample == "train", ks], tolerance = 0.05)
  expect_equal(nrow(scr_score_gains(sc, "holdout")), sc$config$score_groups)
})

test_that("gains carry the band shares, the band WOE and odds in the orientation of the scale", {
  sc <- sc_demo()
  g <- scr_score_gains(sc)
  expect_identical(names(g), c("sample", "id", "band", "n", "pct", "events", "non_events", "event_rate",
                               "pct_event", "pct_nonevent", "woe", "min_score", "mean_score", "max_score",
                               "cum_pct", "cum_event_pct", "cum_nonevent_pct", "ks", "lift", "cum_lift",
                               "odds", "log_odds"))
  expect_equal(g[, sum(pct_event), by = sample]$V1, c(1, 1))
  expect_equal(g[, sum(pct_nonevent), by = sample]$V1, c(1, 1))
  expect_true(all(g$events > 0 & g$non_events > 0))
  expect_equal(g$woe, log(g$pct_event / g$pct_nonevent))
  expect_equal(sign(g$woe), sign(g$lift - 1))
  # higher_is_safer: non-events per event
  expect_equal(g$odds, (g$n - g$events + 0.5) / (g$events + 0.5))
  expect_true(all(g[, stats::cor(log_odds, mean_score), by = sample]$V1 > 0))
  # higher_is_riskier: events per non-event; log_odds still rises with the score
  gf <- scr_score_gains(sc_fraud_demo())
  expect_equal(gf$odds, (gf$events + 0.5) / (gf$n - gf$events + 0.5))
  expect_true(all(gf[, stats::cor(log_odds, mean_score), by = sample]$V1 > 0))
  expect_equal(gf$woe, log(gf$pct_event / gf$pct_nonevent))
  expect_equal(gf[, sum(pct_event), by = sample]$V1, c(1, 1))
})

test_that("the points table carries the event and non-event shares of every bin", {
  p <- sc_demo()$points
  i <- match("pos_rate", names(p))
  expect_identical(names(p)[i + 1:2], c("pct_event", "pct_nonevent"))
  tot <- p[, .(e = sum(pct_event), ne = sum(pct_nonevent)), by = variable]
  expect_equal(tot$e, rep(1, nrow(tot)))
  expect_equal(tot$ne, rep(1, nrow(tot)))
  expect_true(all(p$pct_event >= 0 & p$pct_nonevent >= 0))
})

test_that("stability, calibration and rank-order diagnostics are populated", {
  sc <- sc_demo()
  st <- sc$stability
  expect_true(is.finite(st$score$psi))
  expect_true(st$score$flag_adjusted %in% c("stable", "shift"))
  expect_setequal(st$variables$variable, sc$features)
  expect_true(all(is.finite(st$variables$points_shift)))
  cal <- sc$calibration$summary
  expect_true(cal$brier > 0 && cal$brier < 0.25)
  expect_equal(cal$slope, 1, tolerance = 0.4)
  expect_equal(nrow(sc$rank_order), sc$config$score_groups)
  expect_true(all(sc$rank_order$p_value[-1] >= 0 & sc$rank_order$p_value[-1] <= 1))
})

test_that("the rank-order p-value is the one-sided Fisher exact test against the previous band", {
  # riskiest band first; the last pair has no events
  g <- data.table::data.table(id = 1:6, band = letters[1:6], n = c(20L, 15L, 5L, 10L, 10L, 12L),
                              events = c(8L, 9L, 1L, 5L, 0L, 0L))
  g[, event_rate := events / n]
  ro <- .rank_order(g)
  expect_named(ro, c("id", "band", "n", "events", "event_rate", "prev_rate", "monotone", "p_value", "break_flag"))
  fisher_p <- vapply(2:6, function(i) {
    m <- matrix(c(g$events[i], g$events[i - 1L], g$n[i] - g$events[i], g$n[i - 1L] - g$events[i - 1L]), 2)
    stats::fisher.test(m, alternative = "greater")$p.value
  }, numeric(1))
  expect_equal(ro$p_value[-1], fisher_p)
  expect_true(is.na(ro$p_value[1]) && is.na(ro$monotone[1]) && is.na(ro$break_flag[1]))
  expect_equal(ro$p_value[6], 1)
  expect_identical(ro$monotone[-1], g$event_rate[-1] <= g$event_rate[-6])
  expect_identical(ro$break_flag[-1], !ro$monotone[-1] & ro$p_value[-1] < 0.05)
  # 5/10 after 1/5: a binomial test taking 1/5 as known flags a break
  # (p 0.033); with both rates estimated the rise is not significant
  expect_lt(stats::pbinom(4L, 10L, 0.2, lower.tail = FALSE), 0.05)
  expect_gt(ro$p_value[4], 0.05)
  expect_false(ro$break_flag[4])
  # a clear rise after a large band is still a break
  g2 <- data.table::data.table(id = 1:2, band = c("a", "b"), n = c(400L, 400L), events = c(40L, 80L))
  g2[, event_rate := events / n]
  expect_true(.rank_order(g2)$break_flag[2])
})

test_that("tied training scores that merge score bands are reported", {
  r <- res_demo()
  r$config$verbose <- TRUE
  # a single three-bin variable: three distinct scores, so at most three bands
  f <- names(which.min(vapply(scr_selected(r), function(v) length(r$fit$results[[v]]$bin), integer(1))))
  msgs <- capture_messages(s1 <- scr_scorecard(r, features = f))
  nb <- length(s1$breaks) - 1L
  expect_lt(nb, s1$config$score_groups)
  expect_true(any(grepl(sprintf("score bands: %d of %d requested (ties in the training score)", nb, s1$config$score_groups),
                        msgs, fixed = TRUE)))
  # no message when every requested band is kept
  r$config$score_groups <- 2L
  msgs2 <- capture_messages(s2 <- scr_scorecard(r, features = f))
  expect_equal(length(s2$breaks) - 1L, 2L)
  expect_false(any(grepl("score bands:", msgs2, fixed = TRUE)))
})

test_that("the challenger is aligned to the same scale and never pretends to be a scorecard", {
  sc <- scr_scorecard(res_demo(), challenger = "xgboost", n_boot = 10)
  ch <- sc$challenger
  expect_false(ch$supports_scorecard)
  expect_equal(ch$points, "NOT_APPLICABLE_ENGINE")
  expect_equal(ch$reason_codes, "NOT_APPLICABLE_ENGINE")
  expect_equal(ch$alignment$base_score, sc$scale$base_score)
  expect_equal(ch$alignment$odds_orientation, sc$odds_orientation)
  expect_gt(ch$metrics$auc, 0.6)
  expect_equal(nrow(ch$swapset), 3L)
  expect_equal(sc$model_card$challenger_supports_scorecard, FALSE)
  expect_output(print(sc), "supports_scorecard = FALSE")
})

test_that("distributed points spread the base over the characteristics", {
  sc <- scr_scorecard(res_demo(), points_style = "distributed", n_boot = 10)
  expect_equal(sc$base_points, 0)
  s <- sc$samples$train
  expect_lt(max(abs(s$score - s$score_points)), length(sc$features) * 0.5 + 1)
  expect_equal(mean(s$score), mean(sc_demo()$samples$train$score), tolerance = 1e-8)
})

test_that("scorecard overrides are validated and applied", {
  sc <- scr_scorecard(res_demo(), base_score = 700, base_odds = 30, pdo = 25, n_boot = 10)
  expect_equal(sc$scale$base_score, 700)
  expect_equal(sc$scale$factor, 25 / log(2))
  expect_error(scr_scorecard(res_demo(), features = "not_a_feature"), "without a fitted binning")
  expect_error(scr_scorecard(list()), "scr_select")
  expect_equal(scr_scorecard(res_demo(), align_method = "direct", n_boot = 10)$alignment$calibration$method, "direct")
})

test_that("reason codes point to the variables with the largest shortfall", {
  sc <- sc_demo()
  r <- scr_reasons(sc, head(scr_demo, 6), k = 3)
  expect_equal(nrow(r), 6L)
  expect_true(all(r$reason_1 %in% sc$features))
  expect_true(all(r$shortfall_1 >= r$shortfall_2 & r$shortfall_2 >= r$shortfall_3))
  expect_error(scr_reasons(res_demo(), scr_demo), "scr_scorecard")
})
