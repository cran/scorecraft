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
