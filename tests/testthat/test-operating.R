# Operating point: the curve against a manual cumulative computation on both
# sides, every constraint, the optimum and its shadow price on constructed
# cases, the daily capacity and the agreement with scr_strategy().

op_df <- function(n = 4000, seed = 51) {
  set.seed(seed)
  x <- stats::rnorm(n)
  data.frame(score = round(500 + 40 * x), y = stats::rbinom(n, 1, stats::plogis(-1.5 + 1.2 * x)),
             day = as.Date("2026-03-01") + sample(0:19, n, TRUE), w = stats::runif(n, 0.3, 2.5),
             amt = stats::rexp(n, 1 / 200), smp = sample(c("dev", "oot"), n, TRUE), stringsAsFactors = FALSE)
}

# Pool adjacent violators by repeated merging of the first violating pair,
# written independently of the package
pav_rates <- function(e, n) {
  groups <- as.list(seq_along(e))
  repeat {
    r <- vapply(groups, function(g) sum(e[g]) / sum(n[g]), 1)
    v <- which(diff(r) < 0)
    if (!length(v)) break
    i <- v[1]; groups[[i]] <- c(groups[[i]], groups[[i + 1L]]); groups[[i + 1L]] <- NULL
  }
  out <- numeric(length(e))
  for (g in groups) out[g] <- sum(e[g]) / sum(n[g])
  out
}

test_that("the curve equals a manual cumulative computation, side event from the high scores", {
  d <- op_df()
  op <- scr_operating(d, objective = "propensity", gain_event = 30, cost_select = 5, n_points = 1e6)
  expect_s3_class(op, "scr_operating")
  expect_identical(op$side, "event"); expect_true(op$select_high)
  cv <- op$curve
  u <- sort(unique(d$score), decreasing = TRUE)
  expect_equal(nrow(cv), length(u))
  # every cut sits midway between adjacent distinct scores; the last row selects everything
  expect_equal(cv$cut, c((u[-length(u)] + u[-1]) / 2, -Inf))
  sel_n <- vapply(cv$cut, function(ct) sum(d$score >= ct), 1)
  sel_e <- vapply(cv$cut, function(ct) sum(d$y[d$score >= ct]), 1)
  expect_equal(cv$n_sel, sel_n); expect_equal(cv$events_sel, sel_e)
  expect_equal(cv$depth, sel_n / nrow(d))
  expect_equal(cv$rate_sel, sel_e / sel_n)
  expect_equal(cv$rate_lo, ifelse(sel_e == 0, 0, stats::qbeta(0.025, sel_e + 0.5, sel_n - sel_e + 0.5)))
  expect_equal(cv$rate_hi, ifelse(sel_e == sel_n, 1, stats::qbeta(0.975, sel_e + 0.5, sel_n - sel_e + 0.5)))
  expect_equal(cv$capture, sel_e / sum(d$y))
  expect_equal(cv$lift, sel_e / sel_n / mean(d$y))
  expect_equal(cv$cost, 5 * sel_n)
  expect_equal(cv$value, 30 * sel_e - 5 * sel_n)
  # marginal rate: the smoothed rate of the score value just added (event-poor end first for the smoothing)
  ue <- as.numeric(tapply(d$y, factor(d$score, levels = rev(u)), sum))
  un <- as.numeric(tapply(d$y, factor(d$score, levels = rev(u)), length))
  expect_equal(cv$marginal_rate, rev(pav_rates(ue, un)))
  expect_true(all(diff(cv$marginal_rate) <= 1e-12))
  expect_true(all(cv$feasible))
  # no constraint: the best value
  expect_equal(op$optimum$value, max(cv$value))
  expect_identical(op$optimum$binding, "value")
})

test_that("the curve equals a manual computation on side safe and on side event from the low scores", {
  d <- op_df()
  # credit: higher_is_safer, approval from the high scores
  op <- scr_operating(d, objective = "risk", direction = "higher_is_safer", revenue_good = 100, loss_bad = 400,
                      cost_select = 2, n_points = 1e6)
  expect_identical(op$side, "safe"); expect_true(op$select_high)
  cv <- op$curve
  acc_n <- vapply(cv$cut, function(ct) sum(d$score >= ct), 1)
  acc_e <- vapply(cv$cut, function(ct) sum(d$y[d$score >= ct]), 1)
  expect_equal(cv$n_sel, acc_n); expect_equal(cv$events_sel, acc_e)
  expect_equal(cv$value, 100 * (acc_n - acc_e) - 400 * acc_e - 2 * acc_n)
  # side event under higher_is_safer: collections from the low scores, `score < cut`
  oe <- scr_operating(d, side = "event", objective = "risk", gain_event = 10, n_points = 1e6)
  expect_false(oe$select_high)
  u <- sort(unique(d$score))
  expect_equal(oe$curve$cut, c((u[-length(u)] + u[-1]) / 2, Inf))
  expect_equal(oe$curve$n_sel, vapply(oe$curve$cut, function(ct) sum(d$score < ct), 1))
  expect_equal(oe$curve$events_sel, vapply(oe$curve$cut, function(ct) sum(d$y[d$score < ct]), 1))
  # fraud: risk with higher_is_riskier defaults to side event from the high scores
  of <- scr_operating(d, objective = "risk", direction = "higher_is_riskier", n_points = 50)
  expect_identical(of$side, "event"); expect_true(of$select_high)
  expect_true(all(is.na(of$curve$value)))
  # without economics the optimum is the deepest feasible row
  expect_equal(of$optimum$depth, 1)
  expect_identical(of$optimum$binding, "end of the curve")
})

test_that("each constraint alone and combined picks the deepest feasible cut", {
  d <- op_df()
  base <- scr_operating(d, objective = "propensity", n_points = 1e6)
  cv <- base$curve
  last <- function(ok) max(which(ok))
  o1 <- scr_operating(d, objective = "propensity", max_n = 700, n_points = 1e6)
  expect_equal(o1$optimum$n_sel, cv$n_sel[last(cv$n_sel <= 700)])
  expect_identical(o1$optimum$binding, "max_n")
  o2 <- scr_operating(d, objective = "propensity", max_share = 0.25, n_points = 1e6)
  expect_equal(o2$optimum$depth, cv$depth[last(cv$depth <= 0.25)])
  o3 <- scr_operating(d, objective = "propensity", cost_select = 3, budget = 1500, n_points = 1e6)
  expect_equal(o3$optimum$n_sel, cv$n_sel[last(3 * cv$n_sel <= 1500)])
  expect_identical(o3$optimum$binding, "budget")
  o4 <- scr_operating(d, objective = "propensity", min_rate = 0.5, n_points = 1e6)
  expect_equal(o4$optimum$n_sel, cv$n_sel[last(cv$rate_sel >= 0.5)])
  expect_identical(o4$optimum$binding, "min_rate")
  expect_true(all(o4$curve$feasible == (cv$rate_sel >= 0.5)))
  # combined: the tightest wins, and every constraint is reported with its value at the optimum
  oc <- scr_operating(d, objective = "propensity", max_n = 700, max_share = 0.25, min_rate = 0.5, n_points = 1e6)
  expect_equal(oc$optimum$n_sel, min(o1$optimum$n_sel, o2$optimum$n_sel, o4$optimum$n_sel))
  expect_setequal(oc$constraints$constraint, c("max_n", "max_share", "min_rate"))
  expect_equal(oc$constraints[constraint == "max_n", at_optimum], oc$optimum$n_sel)
  expect_identical(oc$constraints$binding, oc$constraints$constraint %in% strsplit(oc$optimum$binding, ", ")[[1]])
  # max_rate on the safe side: the event rate among the accepted (here the low scores are safe)
  os <- scr_operating(d, objective = "risk", direction = "higher_is_riskier", side = "safe", max_rate = 0.12,
                      n_points = 1e6)
  expect_false(os$select_high)
  acc <- os$curve
  expect_equal(os$optimum$n_sel, acc$n_sel[last(acc$rate_sel <= 0.12)])
  # a limit met exactly is feasible
  ex <- scr_operating(d, objective = "propensity", max_n = cv$n_sel[10], n_points = 1e6)
  expect_equal(ex$optimum$n_sel, cv$n_sel[10])
})

test_that("optimum and shadow price on constructed cases with a known answer", {
  # propensity: rates 80%, 60%, 35%, 20%, 5% from the top score; gain 10, cost 4 per contact
  cnt <- data.frame(score = 5:1, n = 100, events = c(80, 60, 35, 20, 5))
  op <- scr_operating(cnt, counts = TRUE, objective = "propensity", gain_event = 10, cost_select = 4)
  expect_equal(op$curve$value, c(400, 600, 550, 350, 0))
  expect_equal(op$curve$marginal_rate, c(0.8, 0.6, 0.35, 0.2, 0.05))
  expect_equal(op$optimum$cut, 3.5); expect_equal(op$optimum$value, 600)
  expect_identical(op$optimum$binding, "value")
  expect_equal(op$optimum$shadow_price, 10 * 0.35 - 4)
  expect_equal(op$optimum$next_rate, 0.35)
  # a budget that allows a feasible but worse next row still binds: it keeps out the best row
  ob <- scr_operating(cnt, counts = TRUE, objective = "propensity", gain_event = 10, cost_select = 4, budget = 1000)
  expect_equal(ob$optimum$n_sel, 200)
  ob2 <- scr_operating(data.frame(score = 5:1, n = 100, events = c(80, 10, 90, 20, 5)), counts = TRUE,
                       objective = "propensity", gain_event = 10, cost_select = 4, budget = 1000)
  expect_equal(ob2$curve$value, c(400, 100, 600, 400, 50))
  expect_equal(ob2$optimum$n_sel, 100)
  expect_identical(ob2$optimum$binding, "budget")
  # a volume cap below the economic optimum binds; one more contact is worth 10 * 0.6 - 4
  om <- scr_operating(cnt, counts = TRUE, objective = "propensity", gain_event = 10, cost_select = 4, max_n = 150)
  expect_equal(om$optimum$cut, 4.5)
  expect_identical(om$optimum$binding, "max_n")
  expect_equal(om$optimum$shadow_price, 2)
  # ties in value: the smallest depth wins
  tie <- data.frame(score = 3:1, n = 100, events = c(80, 40, 10))
  ot <- scr_operating(tie, counts = TRUE, objective = "propensity", gain_event = 10, cost_select = 4)
  expect_equal(ot$curve$value, c(400, 400, 100))
  expect_equal(ot$optimum$n_sel, 100)
  # credit: approve from the safe end, revenue 100 per good, loss 1000 per bad
  cr <- data.frame(score = 4:1, n = 100, events = c(2, 5, 15, 40))
  oc <- scr_operating(cr, counts = TRUE, revenue_good = 100, loss_bad = 1000)
  expect_equal(oc$curve$value, c(7800, 12300, 5800, -28200))
  expect_equal(oc$optimum$cut, 2.5)
  expect_equal(oc$optimum$shadow_price, 100 * 0.85 - 1000 * 0.15)
  # a non-monotone score is smoothed for the margin: the PAV block rate of the next value
  nm <- data.frame(score = 4:1, n = 100, events = c(60, 20, 40, 10))
  on <- scr_operating(nm, counts = TRUE, objective = "propensity", gain_event = 10, cost_select = 3, max_n = 100)
  expect_equal(on$curve$marginal_rate, c(0.6, 0.3, 0.3, 0.1))
  expect_equal(on$optimum$shadow_price, 10 * 0.3 - 3)
  # every feasible cut loses: the optimum stays the best of them, with a note
  lose <- scr_operating(cnt, counts = TRUE, objective = "propensity", gain_event = 1, cost_select = 4)
  expect_equal(lose$optimum$n_sel, 100)
  expect_match(lose$message, "loses value")
})

test_that("the daily quantile equals a manual per-day computation", {
  d <- op_df()
  op <- scr_operating(d, objective = "propensity", date = "day", max_per_day = 25, day_quantile = 0.8, n_points = 60)
  expect_identical(op$n_days, 20L)
  days <- sort(unique(d$day))
  for (i in seq_len(nrow(op$curve))) {
    sel <- d$score >= op$curve$cut[i]
    vol <- vapply(days, function(dd) sum(sel & d$day == dd), 1)
    expect_equal(op$curve$day_q[i], stats::quantile(vol, 0.8, names = FALSE, type = 7))
    expect_equal(op$curve$pct_days_over[i], mean(vol > 25))
  }
  # the optimum is the deepest cut of every candidate whose daily quantile is within capacity
  u <- sort(unique(d$score), decreasing = TRUE)
  cuts <- c((u[-length(u)] + u[-1]) / 2, -Inf)
  q <- vapply(cuts, function(ct) {
    sel <- d$score >= ct
    stats::quantile(vapply(days, function(dd) sum(sel & d$day == dd), 1), 0.8, names = FALSE)
  }, 1)
  expect_true(all(diff(q) >= 0))
  expect_equal(op$optimum$cut, cuts[max(which(q <= 25))])
  expect_identical(op$optimum$binding, "max_per_day")
  expect_equal(op$constraints$at_optimum, q[max(which(q <= 25))])
  # weighted volumes per day, and rows without a date left out of the days only
  dw <- d; dw$day[1:30] <- NA
  ow <- scr_operating(dw, objective = "propensity", date = "day", weight = "w", n_points = 30)
  for (i in c(1, 10, nrow(ow$curve))) {
    sel <- dw$score >= ow$curve$cut[i] & !is.na(dw$day)
    vol <- vapply(days, function(dd) sum(dw$w[sel & dw$day %in% dd]), 1)
    expect_equal(ow$curve$day_q[i], stats::quantile(vol, 0.9, names = FALSE))
  }
  expect_equal(ow$curve$n_sel[nrow(ow$curve)], sum(dw$w))
  expect_true(all(is.na(ow$curve$pct_days_over)))
  expect_output(print(op), "day_q")
})

test_that("the credit side agrees with the cumulative profit of scr_strategy()", {
  sc <- sc_demo()
  st <- scr_strategy(sc, revenue_good = 1080, loss_bad = 4500)
  op <- scr_operating(sc, revenue_good = 1080, loss_bad = 4500, n_points = 1e6)
  expect_identical(op$side, "safe"); expect_identical(op$sample, "holdout")
  cv <- op$curve
  for (k in seq_len(nrow(st$table))) {
    j <- match(round(st$table$cum_pct[k] * nrow(sc$samples$holdout)), cv$n_sel)
    expect_false(is.na(j))
    expect_equal(cv$value[j], st$table$cum_profit[k])
  }
  # the scorecard dates give the days, and the configured level the intervals
  expect_identical(op$n_days, length(unique(sc$samples$holdout$date)))
  sc2 <- sc; sc2$config$study_level <- 0.8
  o8 <- scr_operating(sc2, n_points = 20)
  i <- nrow(o8$curve)
  x <- o8$curve$events_sel[i]; n <- o8$curve$n_sel[i]
  expect_equal(o8$curve$rate_lo[i], stats::qbeta(0.1, x + 0.5, n - x + 0.5))
  expect_error(scr_operating(sc, sample = "oot"), "not in the scorecard")
  expect_output(print(op), "scr_operating")
})

test_that("the value column replaces the gain, and weights act on volumes and intervals", {
  d <- op_df()
  ov <- scr_operating(d, objective = "propensity", value = "amt", cost_select = 4, n_points = 1e6)
  sel <- function(ct) d$score >= ct
  expect_equal(ov$curve$value, vapply(ov$curve$cut, function(ct) sum((d$amt * d$y)[sel(ct)]) - 4 * sum(sel(ct)), 1))
  expect_true(ov$economics)
  ow <- scr_operating(d, objective = "propensity", weight = "w", gain_event = 20, n_points = 1e6)
  ct <- ow$curve$cut[40]
  i <- d$score >= ct
  ny <- sum(d$w[i]); e <- sum(d$w[i] * d$y[i]); neff <- ny^2 / sum(d$w[i]^2)
  expect_equal(ow$curve$n_sel[40], ny); expect_equal(ow$curve$events_sel[40], e)
  expect_equal(ow$curve$rate_lo[40], stats::qbeta(0.025, e / ny * neff + 0.5, neff - e / ny * neff + 0.5))
  # zero weights leave the rows out
  dz <- d; dz$w[1:100] <- 0
  expect_equal(max(scr_operating(dz, objective = "propensity", weight = "w")$curve$n_sel), sum(dz$w))
  # a sample column: the curve on the study sample
  os <- scr_operating(d, objective = "propensity", sample = "smp", n_points = 1e6)
  expect_identical(os$sample, "oot")
  expect_equal(max(os$curve$n_sel), sum(d$smp == "oot"))
  expect_identical(scr_operating(d, objective = "propensity", sample = "smp", study = "dev")$sample, "dev")
})

test_that("thinning keeps about n_points rows and every row the constraints pick", {
  set.seed(3)
  d <- data.frame(score = stats::rnorm(20000), y = stats::rbinom(20000, 1, 0.2))
  op <- scr_operating(d, objective = "propensity", max_n = 3333, min_rate = 0.1, n_points = 50)
  full <- scr_operating(d, objective = "propensity", max_n = 3333, min_rate = 0.1, n_points = 1e6)
  expect_lt(nrow(op$curve), 70); expect_gt(nrow(op$curve), 40)
  expect_true(op$optimum$n_sel %in% op$curve$n_sel)
  expect_equal(op$optimum[, -c("binding", "shadow_price", "next_rate")],
               full$optimum[, -c("binding", "shadow_price", "next_rate")])
  expect_true(max(full$curve$n_sel[full$curve$n_sel <= 3333]) %in% op$curve$n_sel)
  expect_equal(op$curve$depth[nrow(op$curve)], 1)
  expect_true(all(op$curve$n_sel %in% full$curve$n_sel))
})

test_that("edge cases: infeasible constraints, single class, all ties, bad arguments", {
  d <- op_df(1000)
  expect_warning(op <- scr_operating(d, objective = "propensity", max_share = 1e-5), "no cut meets the constraints")
  expect_true(is.na(op$optimum$depth)); expect_true(is.na(op$optimum$binding))
  expect_match(op$message, "max_share")
  expect_output(print(op), "Optimum: NA")
  # single class: rates 0, the economics still defined
  one <- data.frame(score = 1:50, y = 0)
  o1 <- scr_operating(one, objective = "propensity", gain_event = 5, cost_select = 1)
  expect_true(all(o1$curve$rate_sel == 0)); expect_true(all(is.na(o1$curve$capture)))
  # no event, no lift: missing, not NaN
  expect_identical(o1$curve$lift, rep(NA_real_, nrow(o1$curve)))
  expect_equal(o1$optimum$n_sel, 1)
  # all ties: one row selecting everything
  ti <- scr_operating(data.frame(score = rep(3, 40), y = rep(0:1, 20)), objective = "propensity")
  expect_equal(nrow(ti$curve), 1L); expect_equal(ti$curve$cut, -Inf); expect_identical(ti$optimum$binding, "end of the curve")
  expect_error(scr_operating(d, objective = "propensity", revenue_good = 1), "does not apply|do\\(es\\) not apply")
  expect_error(scr_operating(d, gain_event = 1), "do\\(es\\) not apply")
  expect_error(scr_operating(d, objective = "propensity", budget = 10), "positive `cost_select`")
  expect_error(scr_operating(d, objective = "propensity", max_per_day = 10), "needs a date")
  expect_error(scr_operating(d, side = "both"), "`side`")
  expect_error(scr_operating(d, max_share = 1.5), "max_share")
  expect_error(scr_operating(d, study = "a"), "needs a `sample` column")
  expect_error(scr_operating(d, foo = 1), "unused argument")
  expect_error(scr_operating(data.frame(score = NA_real_, y = 1)), "no scored row")
})

test_that("scr_export writes the operating workbook", {
  skip_if_not_installed("openxlsx")
  old <- scr_verbose(FALSE); on.exit(scr_verbose(old), add = TRUE)
  op <- scr_operating(op_df(800), objective = "propensity", gain_event = 10, cost_select = 2, max_n = 200)
  out <- file.path(tempdir(), "scr-operating-export")
  unlink(out, recursive = TRUE)
  ex <- scr_export(op, out, stamp = FALSE)
  expect_identical(openxlsx::getSheetNames(ex$files$xlsx), c("Curve", "Optimum", "Constraints", "Settings"))
  expect_equal(nrow(openxlsx::read.xlsx(ex$files$xlsx, sheet = "Curve")), nrow(op$curve))
})

test_that("a constraint binds only when dropping it alone improves the optimum", {
  # a slack constraint next to a binding one: the budget stops the campaign, the share cap is far away
  ch <- sc_prop_demo()
  op <- scr_operating(ch, gain_event = 100, cost_select = 20, budget = 3000, max_share = 0.5)
  expect_identical(op$optimum$binding, "budget")
  expect_lt(op$optimum$depth, 0.5)
  expect_identical(op$constraints[constraint == "max_share", binding], FALSE)
  expect_identical(op$constraints[constraint == "budget", binding], TRUE)
  # dropping the budget alone improves the value; dropping the share cap alone changes nothing
  expect_gt(scr_operating(ch, gain_event = 100, cost_select = 20, max_share = 0.5)$optimum$value, op$optimum$value)
  expect_equal(scr_operating(ch, gain_event = 100, cost_select = 20, budget = 3000)$optimum$value, op$optimum$value)
  expect_output(print(op), "constraints: max_share 0.5, budget 3,000", fixed = TRUE)
  # a constraint that stops the curve at the row that is the best anyway does not bind
  cnt <- data.frame(score = 5:1, n = 100, events = c(90, 70, 50, 30, 10))
  ob <- scr_operating(cnt, counts = TRUE, objective = "propensity", gain_event = 10, cost_select = 4, budget = 1200)
  expect_equal(ob$curve$value, c(500, 800, 900, 800, 500))
  expect_equal(ob$optimum$n_sel, 300)
  expect_false(ob$curve$feasible[4])
  expect_identical(ob$optimum$binding, "value")
  expect_equal(ob$optimum$shadow_price, -1)
  expect_identical(ob$constraints$binding, FALSE)
  expect_identical(scr_operating(cnt, counts = TRUE, objective = "propensity", gain_event = 10,
                                 cost_select = 4)$optimum$n_sel, 300)
  # a tighter budget does bind, and one more contact is worth 10 * 0.5 - 4
  ot <- scr_operating(cnt, counts = TRUE, objective = "propensity", gain_event = 10, cost_select = 4, budget = 800)
  expect_identical(ot$optimum$binding, "budget")
  expect_equal(ot$optimum$shadow_price, 1)
  # two constraints that stop at the same row bind jointly: none improves the optimum alone
  oj <- scr_operating(cnt, counts = TRUE, objective = "propensity", gain_event = 10, cost_select = 4,
                      max_n = 150, max_share = 0.3)
  expect_equal(oj$optimum$n_sel, 100)
  expect_identical(oj$optimum$binding, "max_n, max_share")
  expect_identical(scr_operating(cnt, counts = TRUE, objective = "propensity", max_n = 150,
                                 max_share = 0.3)$optimum$binding, "max_n, max_share")
  # of two constraints at different rows only the tighter one binds, with or without economics
  o2 <- scr_operating(cnt, counts = TRUE, objective = "propensity", max_n = 150, max_share = 0.7)
  expect_identical(o2$optimum$binding, "max_n")
  expect_identical(o2$constraints$binding, c(TRUE, FALSE))
  # nothing binds without economics only when every row is selected
  expect_identical(scr_operating(cnt, counts = TRUE, objective = "propensity", max_n = 500)$optimum$binding,
                   "end of the curve")
})
