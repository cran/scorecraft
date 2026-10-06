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
