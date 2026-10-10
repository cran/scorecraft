# Maturity by band: Kaplan-Meier, Greenwood and the complementary log-log
# interval against survival::survfit() and a hand computation; the counts,
# the censoring patterns, the weights and the discrimination by horizon.

mat_df <- function(n = 2500, seed = 91) {
  set.seed(seed)
  x <- stats::rnorm(n)
  te <- ceiling(stats::rexp(n, 0.03 * exp(0.8 * x)))
  tc <- sample(1:40, n, TRUE)
  data.frame(score = round(500 - 40 * x), time = pmin(te, tc), event = as.integer(te <= tc),
             w = stats::runif(n, 0.2, 3), stringsAsFactors = FALSE)
}

# Kaplan-Meier with Greenwood's variance at one horizon, from the rows
km_manual <- function(time, event, h, w = rep(1, length(time)), level = 0.95) {
  tj <- sort(unique(time[event == 1 & time <= h]))
  surv <- 1; gw <- 0
  for (t in tj) {
    at <- time >= t
    n <- sum(w[at]); d <- sum(w[at & time == t & event == 1])
    neff <- n^2 / sum(w[at]^2)
    surv <- surv * (1 - d / n)
    gw <- gw + d / ((n - d) * neff)
  }
  z <- stats::qnorm(1 - (1 - level) / 2)
  sg <- sqrt(gw) / abs(log(surv))
  list(incidence = 1 - surv, se = surv * sqrt(gw), lo = 1 - surv^exp(-z * sg), hi = 1 - surv^exp(z * sg))
}

test_that("incidence, Greenwood error and interval equal survival::survfit()", {
  skip_if_not_installed("survival")
  d <- mat_df()
  hz <- c(3, 6, 12, 24, 36)
  m <- scr_maturity(d, horizons = hz, n_bands = 4, level = 0.9)
  expect_s3_class(m, "scr_maturity")
  f <- survival::survfit(survival::Surv(time, event) ~ 1, data = d, conf.type = "log-log", conf.int = 0.9)
  s <- summary(f, times = hz, extend = TRUE)
  al <- m$table[is.na(band)]
  expect_identical(al$label, rep("all", 5)); expect_equal(al$horizon, hz)
  expect_equal(al$incidence, 1 - s$surv); expect_equal(al$se, s$std.err)
  expect_equal(al$lo, 1 - s$upper); expect_equal(al$hi, 1 - s$lower)
  expect_equal(al$events, cumsum(s$n.event))
  # per band: band 1 is the event-richest, the low scores under higher_is_safer
  d$b <- findInterval(d$score, m$cuts) + 1L
  fb <- survival::survfit(survival::Surv(time, event) ~ b, data = d, conf.type = "log-log", conf.int = 0.9)
  sb <- summary(fb, times = hz, extend = TRUE)
  for (k in 1:4) {
    i <- sb$strata == paste0("b=", k)
    t <- m$table[band == k]
    expect_equal(t$incidence, 1 - sb$surv[i]); expect_equal(t$se, sb$std.err[i])
    expect_equal(t$lo, 1 - sb$upper[i]); expect_equal(t$hi, 1 - sb$lower[i])
    expect_equal(unique(t$n), sum(d$b == k))
  }
  expect_identical(m$table[horizon == 3 & !is.na(band), label], .study_labels(m$cuts))
  # weights: the weighted Kaplan-Meier estimate of survfit()
  mw <- scr_maturity(d, horizons = hz, n_bands = 4, weight = "w")
  fw <- survival::survfit(survival::Surv(time, event) ~ 1, data = d, weights = w)
  expect_equal(mw$table[is.na(band), incidence], 1 - summary(fw, times = hz, extend = TRUE)$surv)
})

test_that("a hand Kaplan-Meier gives the hand numbers", {
  # times 1 (event), 2 (censored), 3 (two events), 4 (censored), 5 (event)
  d <- data.frame(score = c(1, 2, 3, 4, 5, 6), time = c(1, 2, 3, 3, 4, 5), event = c(1, 0, 1, 1, 0, 1))
  m <- scr_maturity(d, horizons = c(0.5, 1, 3, 4.5, 5), n_bands = 1)
  t <- m$table[is.na(band)]
  # S(1) = 5/6, S(3) = 5/6 * 2/4, S(5) = 0
  expect_equal(t$incidence, c(0, 1 / 6, 1 - 5 / 12, 1 - 5 / 12, 1))
  gw1 <- 1 / (6 * 5); gw3 <- gw1 + 2 / (4 * 2)
  expect_equal(t$se[1:4], c(0, 5 / 6 * sqrt(gw1), 5 / 12 * sqrt(gw3), 5 / 12 * sqrt(gw3)))
  z <- stats::qnorm(0.975)
  sg <- sqrt(gw3) / abs(log(5 / 12))
  expect_equal(t$lo[3], 1 - (5 / 12)^exp(-z * sg)); expect_equal(t$hi[3], 1 - (5 / 12)^exp(z * sg))
  # no event yet: the interval has no upper end; every unit with the event: no interval
  expect_equal(t$lo[1], 0); expect_true(is.na(t$hi[1]))
  expect_true(is.na(t$se[5])); expect_true(is.na(t$lo[5])); expect_true(is.na(t$hi[5]))
  # events up to the horizon, censored before it, the others still followed
  expect_equal(t$events, c(0, 1, 3, 3, 4)); expect_equal(t$censored, c(0, 0, 1, 2, 2))
  expect_equal(t$at_risk, c(6, 5, 2, 1, 0))
  expect_equal(t$at_risk + t$events + t$censored, rep(6, 5))
  expect_equal(t$pct_of_final, t$incidence / 1)
  # the one band repeats the whole
  expect_equal(m$table[band == 1, incidence], t$incidence)
  expect_equal(km_manual(d$time, d$event, 3)$incidence, t$incidence[3])
})

test_that("the counts follow the censoring pattern", {
  d <- mat_df()
  hz <- c(5, 10, 20, 40)
  m <- scr_maturity(d, horizons = hz, n_bands = 3)
  d$b <- findInterval(d$score, m$cuts) + 1L
  for (k in 1:3) for (h in hz) {
    i <- d$b == k
    r <- m$table[band == k & horizon == h]
    expect_equal(r$events, sum(i & d$event == 1 & d$time <= h))
    expect_equal(r$censored, sum(i & d$event == 0 & d$time < h))
    expect_equal(r$at_risk, sum(i & (d$time > h | (d$time == h & d$event == 0))))
    expect_equal(r$n, sum(i))
    ref <- km_manual(d$time[i], d$event[i], h)
    expect_equal(r$incidence, ref$incidence); expect_equal(r$se, ref$se)
    expect_equal(r$lo, ref$lo); expect_equal(r$hi, ref$hi)
  }
  al <- m$table[is.na(band)]
  expect_equal(al$pct_of_final, al$incidence / al$incidence[4])
  # a unit censored at the horizon was followed that long; an event at the horizon is an event
  e <- data.frame(score = 1:4, time = c(10, 10, 12, 8), event = c(1, 0, 0, 0))
  r <- scr_maturity(e, horizons = 10, n_bands = 1)$table[is.na(band)]
  expect_equal(r$events, 1); expect_equal(r$censored, 1); expect_equal(r$at_risk, 2)
  expect_equal(r$incidence, 1 / 3)
  # without censoring the incidence is the share of events by the horizon
  nc <- data.frame(score = d$score, time = d$time, event = 1L)
  mn <- scr_maturity(nc, horizons = c(4, 15), n_bands = 2)
  expect_equal(mn$table[is.na(band), incidence], c(mean(nc$time <= 4), mean(nc$time <= 15)))
  expect_equal(mn$table[is.na(band), censored], c(0, 0))
  # censoring before the horizon lifts the estimate above the naive share
  al10 <- m$table[is.na(band) & horizon == 20]
  expect_gt(al10$incidence, mean(d$event == 1 & d$time <= 20))
  # a horizon before the first time: nothing has happened
  first <- scr_maturity(d, horizons = c(0.5, 20), n_bands = 2)$table[horizon == 0.5]
  expect_equal(first$incidence, rep(0, 3)); expect_equal(first$at_risk, first$n)
})

test_that("weights: Kish effective risk sets, equal weights and whole weights", {
  d <- mat_df(1500, 4)
  hz <- c(6, 12, 30)
  mw <- scr_maturity(d, horizons = hz, n_bands = 3, weight = "w", level = 0.9)
  d$b <- findInterval(d$score, mw$cuts) + 1L
  for (h in hz) {
    ref <- km_manual(d$time, d$event, h, d$w, level = 0.9)
    r <- mw$table[is.na(band) & horizon == h]
    expect_equal(r$incidence, ref$incidence); expect_equal(r$se, ref$se)
    expect_equal(r$lo, ref$lo); expect_equal(r$hi, ref$hi)
    i <- d$b == 2
    rb <- mw$table[band == 2 & horizon == h]
    expect_equal(rb$incidence, km_manual(d$time[i], d$event[i], h, d$w[i])$incidence)
    expect_equal(rb$se, km_manual(d$time[i], d$event[i], h, d$w[i])$se)
    expect_equal(rb$events, sum(d$w[i & d$event == 1 & d$time <= h]))
    expect_equal(rb$n, sum(d$w[i]))
  }
  # equal weights change nothing, the error included
  d$c5 <- 5
  m0 <- scr_maturity(d, horizons = hz, n_bands = 3)
  m5 <- scr_maturity(d, horizons = hz, n_bands = 3, weight = "c5")
  expect_equal(m5$table[, .(incidence, se, lo, hi, pct_of_final)], m0$table[, .(incidence, se, lo, hi, pct_of_final)])
  expect_equal(m5$table$n, 5 * m0$table$n)
  expect_equal(m5$discrimination$auc, m0$discrimination$auc)
  # whole weights give the incidence of the repeated rows
  d$k <- sample(1:3, nrow(d), TRUE)
  rp <- d[rep(seq_len(nrow(d)), d$k), ]
  mk <- scr_maturity(d, horizons = hz, cuts = m0$cuts, weight = "k")
  mr <- scr_maturity(rp, horizons = hz, cuts = m0$cuts)
  expect_equal(mk$table[, .(n, at_risk, events, censored, incidence)], mr$table[, .(n, at_risk, events, censored, incidence)])
  expect_equal(mk$discrimination, mr$discrimination)
  # a zero weight takes the row out
  d$w0 <- as.numeric(seq_len(nrow(d)) %% 2)
  mz <- scr_maturity(d, horizons = hz, cuts = m0$cuts, weight = "w0")
  expect_equal(mz$table[, .(n, incidence, se)], scr_maturity(d[d$w0 > 0, ], horizons = hz, cuts = m0$cuts)$table[, .(n, incidence, se)])
  expect_equal(mz$n_dropped, sum(d$w0 == 0))
})

test_that("the discrimination by horizon equals the row-level AUC on the complete windows", {
  d <- mat_df()
  hz <- c(3, 12, 24, 39)
  m <- scr_maturity(d, horizons = hz, n_bands = 4)
  for (h in hz) {
    full <- (d$event == 1 & d$time <= h) | d$time > h | (d$time == h & d$event == 0)
    yh <- as.integer(d$event == 1 & d$time <= h)[full]
    r <- m$discrimination[horizon == h]
    expect_equal(r$n, sum(full)); expect_equal(r$events, sum(yh))
    a <- .auc_ks(-d$score[full], yh)$auc
    expect_equal(r$auc, a); expect_equal(r$gini, 2 * a - 1)
    expect_equal(r$auc, scr_metrics(d$score[full], yh, higher_is_event = FALSE, ci = FALSE)$auc)
  }
  # weighted counts
  mw <- scr_maturity(d, horizons = 12, n_bands = 4, weight = "w")
  full <- (d$event == 1 & d$time <= 12) | d$time > 12 | (d$time == 12 & d$event == 0)
  yh <- as.integer(d$event == 1 & d$time <= 12)
  expect_equal(mw$discrimination$n, sum(d$w[full])); expect_equal(mw$discrimination$events, sum(d$w[full & yh == 1]))
  idx <- factor(data.table::frank(-d$score, ties.method = "dense"))
  c1 <- as.numeric(tapply(d$w * (full & yh == 1), idx, sum)); c0 <- as.numeric(tapply(d$w * (full & yh == 0), idx, sum))
  expect_equal(mw$discrimination$auc, .auc_ks_counts(c1, c0)$auc)
  # a score read the other way: the same curves and discrimination
  f <- d; f$score <- -f$score
  mf <- scr_maturity(f, horizons = hz, n_bands = 4, direction = "higher_is_riskier")
  expect_equal(mf$discrimination$auc, m$discrimination$auc)
  expect_equal(mf$table$incidence, m$table$incidence)
  expect_identical(mf$table[!is.na(band) & horizon == 3, band], 1:4)
})

test_that("bands: equal shares, explicit cuts and a score study", {
  d <- mat_df()
  m <- scr_maturity(d, horizons = c(6, 24), n_bands = 5)
  # tie-safe equal shares of the units, as scr_bands()
  expect_equal(m$cuts, scr_bands(data.frame(score = d$score, y = d$event), n_bands = 5, n_boot = 0)$cuts)
  # the event-richest band first: the incidence falls along the bands
  expect_true(all(diff(m$table[horizon == 24 & !is.na(band), incidence]) < 0))
  mc <- scr_maturity(d, horizons = c(6, 24), cuts = c(480, 520))
  expect_equal(mc$cuts, c(480, 520)); expect_identical(mc$bands, "cuts")
  expect_equal(mc$table[band == 1 & horizon == 6, n], sum(d$score < 480))
  # left-closed bands: a score on a cut is in the upper band
  expect_equal(mc$table[band == 3 & horizon == 6, n], sum(d$score >= 520))
  # a study brings its cuts, numbers, labels and direction
  tr <- scr_tiers(data.frame(score = d$score, y = as.integer(d$event == 1 & d$time <= 24)), n_tiers = 3)
  mt <- scr_maturity(d, horizons = c(6, 24), cuts = tr)
  expect_equal(mt$cuts, tr$cuts); expect_identical(mt$bands, "tiers")
  expect_identical(mt$table[horizon == 6 & !is.na(band), label], tr$table$label)
  expect_identical(mt$table[horizon == 6 & !is.na(band), band], tr$table$tier)
  expect_error(scr_maturity(d, horizons = 6, cuts = tr, direction = "higher_is_riskier"), "fitted with")
  # an explicit cut can leave a band empty
  me <- scr_maturity(d, horizons = 6, cuts = c(-1e6, 500))
  e <- me$table[band == 1]
  expect_equal(e$n, 0); expect_true(is.na(e$incidence)); expect_equal(e$at_risk, 0)
})

test_that("edge cases: a single class, tied scores, rows left out, bad input", {
  d <- mat_df(800, 7)
  # no event: incidence 0, no AUC
  z <- d; z$event <- 0L
  mz <- scr_maturity(z, horizons = c(5, 20), n_bands = 2)
  expect_equal(mz$table$incidence, rep(0, 6)); expect_equal(mz$table$lo, rep(0, 6))
  expect_true(all(is.na(mz$table$hi))); expect_true(all(is.na(mz$table$pct_of_final)))
  expect_true(all(is.na(mz$discrimination$auc))); expect_equal(mz$discrimination$events, c(0, 0))
  # every unit has the event: the curve reaches 1 and every window is complete
  a <- d; a$event <- 1L
  ma <- scr_maturity(a, horizons = c(5, 1000), n_bands = 2)
  expect_equal(ma$table[horizon == 1000, incidence], rep(1, 3))
  expect_equal(ma$discrimination$n, c(nrow(a), nrow(a)))
  expect_true(is.na(ma$discrimination$auc[2]))
  # every score tied: one band, AUC 0.5
  tie <- d; tie$score <- 3
  mt <- scr_maturity(tie, horizons = 10, n_bands = 5)
  expect_length(mt$cuts, 0L); expect_equal(nrow(mt$table), 2L); expect_equal(mt$discrimination$auc, 0.5)
  # rows without a score, a time or an event are left out
  na <- d; na$score[1:5] <- NA; na$time[6:9] <- NA; na$event[10:12] <- NA
  mn <- scr_maturity(na, horizons = 10, n_bands = 2)
  expect_equal(mn$n_dropped, 12L); expect_equal(mn$n_rows, nrow(d) - 12L)
  expect_equal(mn$table$incidence, scr_maturity(d[-(1:12), ], horizons = 10, n_bands = 2)$table$incidence)
  # the horizons are sorted and distinct
  expect_equal(scr_maturity(d, horizons = c(20, 5, 20), n_bands = 1)$horizons, c(5, 20))
  expect_error(scr_maturity(d), "horizons")
  expect_error(scr_maturity(d[0, ], horizons = 5), "no row has a score")
  # past the last follow-up of a band the curve is carried flat, with nobody at risk
  last <- max(d$time)
  fl <- scr_maturity(d, horizons = c(last, last + 50), n_bands = 2)$table
  expect_equal(fl[horizon == last + 50, incidence], fl[horizon == last, incidence])
  expect_equal(fl[horizon == last + 50, at_risk], rep(0, 3))
  expect_error(scr_maturity(d, horizons = c(5, -1)), "horizons")
  expect_error(scr_maturity(d, horizons = 5, time = "nope"), "not in `x`")
  expect_error(scr_maturity(transform(d, time = -time), horizons = 5), "non-negative")
  expect_error(scr_maturity(transform(d, event = 2), horizons = 5), "0/1")
  expect_error(scr_maturity(d, horizons = 5, foo = 1), "unused argument")
  expect_error(scr_maturity(transform(d, score = NA_real_), horizons = 5), "no row has a score")
  expect_error(scr_maturity(d, horizons = 5, weight = "w", level = 2), "level")
})

test_that("pooled score cells: the curves are exact, the discrimination close", {
  d <- mat_df(3000, 14)
  d$score <- d$score + stats::runif(nrow(d))
  ex <- scr_maturity(d, horizons = c(6, 24), n_bands = 4)
  po <- scr_maturity(d, horizons = c(6, 24), n_bands = 4, max_cells = 100)
  expect_true(po$quantized); expect_false(ex$quantized)
  # the bands are assigned on the rows: the incidence of each band is that of its rows
  d$b <- findInterval(d$score, po$cuts) + 1L
  for (k in 1:4) {
    expect_equal(po$table[band == k & horizon == 24, incidence], km_manual(d$time[d$b == k], d$event[d$b == k], 24)$incidence)
  }
  expect_equal(po$table[is.na(band), incidence], ex$table[is.na(band), incidence])
  expect_equal(po$discrimination$n, ex$discrimination$n)
  expect_equal(po$discrimination$auc, ex$discrimination$auc, tolerance = 5e-3)
})

test_that("print and export", {
  d <- mat_df(1000, 2)
  m <- scr_maturity(d, horizons = c(6, 12, 24), n_bands = 3)
  out <- utils::capture.output(print(m))
  expect_match(out[1], "^<scr_maturity> event \"event\" over \"time\"")
  expect_true(any(grepl("t=6 +t=12 +t=24", out)))
  expect_true(any(grepl("share of the final", out)))
  expect_true(any(grepl("^Discrimination at each horizon", out)))
  skip_if_not_installed("openxlsx")
  old <- scr_verbose(FALSE); on.exit(scr_verbose(old), add = TRUE)
  dir <- file.path(tempdir(), "scr-maturity-export")
  unlink(dir, recursive = TRUE)
  ex <- scr_export(m, dir, stamp = FALSE)
  expect_identical(basename(ex$files$xlsx), "maturity_event.xlsx")
  expect_identical(openxlsx::getSheetNames(ex$files$xlsx), c("Incidence", "Discrimination", "Cuts", "Settings"))
  expect_equal(nrow(openxlsx::read.xlsx(ex$files$xlsx, sheet = "Incidence")), 12L)
})
