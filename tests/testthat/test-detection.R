# Time to detection: constructed episodes with known answers, every
# threshold against a naive loop over the episodes, the loss accounting and
# the alert shares.

det_df <- function(n = 6000, n_ent = 500, n_hit = 60, seed = 121) {
  set.seed(seed)
  d <- data.frame(ent = sample(n_ent, n, TRUE), ts = stats::runif(n, 0, 30))
  hit <- sample(n_ent, n_hit)
  since <- stats::runif(n_ent, 5, 25)
  d$y <- as.integer(d$ent %in% hit & d$ts >= since[d$ent])
  d$score <- round(100 * stats::plogis(stats::rnorm(n, -1.5 + 1.5 * d$y)))
  d$amount <- round(stats::rexp(n, 1 / 60), 2)
  d
}

# one threshold, one episode at a time: sort the rows of the entity, start
# at its first event row, find the first alert
det_naive <- function(d, t, high = TRUE) {
  ev <- !is.na(d$y) & d$y == 1
  out <- lapply(unique(d$ent[ev]), function(e) {
    r <- d[d$ent == e, ]
    r <- r[order(r$ts), ]
    r <- r[which(r$y == 1)[1]:nrow(r), ]
    k <- which(if (high) r$score >= t else r$score < t)[1]
    isev <- !is.na(r$y) & r$y == 1
    tot <- sum(r$amount[isev])
    if (is.na(k)) return(data.frame(det = FALSE, evb = NA, tm = NA, loss = tot, tot = tot))
    b <- seq_len(k - 1L)
    data.frame(det = TRUE, evb = sum(isev[b]), tm = as.numeric(r$ts[k]) - as.numeric(r$ts[1]),
               loss = sum(r$amount[b][isev[b]]), tot = tot)
  })
  o <- do.call(rbind, out)
  list(episodes = nrow(o), detected = sum(o$det), median_events = stats::median(o$evb[o$det]),
       mean_events = mean(o$evb[o$det]), median_time = stats::median(o$tm[o$det]), loss_before = sum(o$loss),
       loss_total = sum(o$tot), alert_share = mean(if (high) d$score >= t else d$score < t))
}

toy <- function() {
  data.frame(
    ent = c("A", "A", "B", "B", "B", "B", "C", "C", "D", "D", "D"),
    ts = c(1, 2, 1, 2, 3, 4, 5, 6, 1, 2, 3),
    y = c(0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0),
    score = c(95, 90, 20, 30, 85, 40, 10, 15, 99, 50, 81),
    amount = c(10, 100, 50, 70, 10, 30, 200, 100, 5, 5, 5), stringsAsFactors = FALSE)
}

test_that("constructed episodes: detected at the first row, at a later row, never", {
  d <- toy()
  dt <- scr_detection(d, entity = "ent", time = "ts", amount = "amount", thresholds = c(80, 25, 12))
  expect_s3_class(dt, "scr_detection")
  t <- dt$table
  # the strictest threshold first
  expect_equal(t$threshold, c(80, 25, 12)); expect_equal(t$episodes, rep(3, 3))
  # 80: A at its first event row (the 95 before the start does not count), B at its third row, C never
  # 25: B at its second row; 12: A and B at their first rows, C at its second
  expect_equal(t$detected, c(2, 2, 3)); expect_equal(t$pct_detected, c(2, 2, 3) / 3)
  expect_equal(t$median_events_before, c(1, 0.5, 0)); expect_equal(t$mean_events_before, c(1, 0.5, 1 / 3))
  expect_equal(t$median_time, c(1, 0.5, 0))
  # loss before detection, the whole episode when never detected
  expect_equal(t$loss_total, rep(550, 3))
  expect_equal(t$loss_before, c(0 + 120 + 300, 0 + 50 + 300, 0 + 0 + 200))
  expect_equal(t$pct_loss_prevented, 1 - c(420, 350, 200) / 550)
  # the alert share is read on all rows, the entity without fraud included
  expect_equal(t$alert_share, c(5, 8, 10) / 11)
  expect_true(all(is.na(t$target_share)))
  # the episodes, with the first detection at the strictest threshold
  e <- dt$episodes
  expect_identical(e$entity, c("A", "B", "C")); expect_equal(e$start, c(2, 1, 5))
  expect_equal(e$rows, c(1, 4, 2)); expect_equal(e$events, c(1, 3, 2)); expect_equal(e$amount, c(100, 150, 300))
  expect_equal(e$peak_score, c(90, 85, 15))
  expect_identical(e$detected, c(TRUE, TRUE, FALSE))
  expect_equal(e$detect_time, c(2, 3, NA)); expect_equal(e$time_to_detect, c(0, 2, NA))
  expect_equal(e$events_before, c(0, 2, NA)); expect_equal(e$loss_before, c(0, 120, NA))
  expect_equal(dt$n_episode_rows, 7L); expect_equal(dt$n_rows, 11L)
  # the alert is `score >= threshold`: a score on the threshold alerts
  expect_equal(scr_detection(d, entity = "ent", time = "ts", thresholds = c(85, 85.01))$table$detected, c(1, 2))
  # from the low scores the alert is `score < threshold`
  f <- d; f$score <- 100 - f$score
  df <- scr_detection(f, entity = "ent", time = "ts", amount = "amount", thresholds = c(20.5, 75.5, 88.5),
                      direction = "higher_is_safer")
  tf <- df$table
  expect_equal(tf$threshold, c(20.5, 75.5, 88.5))
  expect_equal(tf$detected, t$detected); expect_equal(tf$loss_before, t$loss_before)
  expect_equal(tf$median_time, t$median_time); expect_equal(tf$alert_share, t$alert_share)
  expect_equal(df$episodes$peak_score, 100 - e$peak_score)
  # strict: a score equal to the threshold does not alert from the low side
  expect_equal(scr_detection(f, entity = "ent", time = "ts", thresholds = 15, direction = "higher_is_safer")$table$detected, 1)
  expect_equal(scr_detection(f, entity = "ent", time = "ts", thresholds = 15.01, direction = "higher_is_safer")$table$detected, 2)
})

test_that("every threshold equals a naive loop over the episodes", {
  d <- det_df()
  thr <- c(15, 30, 45.5, 60, 75, 90, 99, 101)
  dt <- scr_detection(d, entity = "ent", time = "ts", amount = "amount", thresholds = thr)
  t <- dt$table
  expect_equal(t$threshold, sort(thr, decreasing = TRUE))
  for (i in seq_len(nrow(t))) {
    ref <- det_naive(d, t$threshold[i])
    expect_equal(t$episodes[i], ref$episodes); expect_equal(t$detected[i], ref$detected)
    expect_equal(t$loss_before[i], ref$loss_before); expect_equal(t$loss_total[i], ref$loss_total)
    expect_equal(t$alert_share[i], ref$alert_share)
    if (ref$detected > 0) {
      expect_equal(t$median_events_before[i], ref$median_events); expect_equal(t$mean_events_before[i], ref$mean_events)
      expect_equal(t$median_time[i], ref$median_time)
    } else {
      expect_true(is.na(t$median_events_before[i])); expect_true(is.na(t$median_time[i]))
      expect_equal(t$pct_loss_prevented[i], 0)
    }
  }
  # a looser threshold never detects fewer episodes, nor later
  expect_true(all(diff(t$detected) >= 0)); expect_true(all(diff(t$loss_before) <= 1e-9))
  # the row order of the input does not matter
  sh <- d[sample(nrow(d)), ]
  expect_equal(scr_detection(sh, entity = "ent", time = "ts", amount = "amount", thresholds = thr)$table, t)
  # from the low scores
  f <- d; f$score <- 100 - f$score
  tf <- scr_detection(f, entity = "ent", time = "ts", amount = "amount", thresholds = 100 - thr,
                      direction = "higher_is_safer")$table
  for (i in seq_len(nrow(tf))) {
    ref <- det_naive(f, tf$threshold[i], high = FALSE)
    expect_equal(tf$detected[i], ref$detected); expect_equal(tf$loss_before[i], ref$loss_before)
    expect_equal(tf$alert_share[i], ref$alert_share)
    if (ref$detected > 0) expect_equal(tf$mean_events_before[i], ref$mean_events)
  }
  # the episodes: one per entity with an event, counted from its first event row
  e <- dt$episodes
  hit <- sort(unique(d$ent[d$y == 1]))
  expect_equal(e$entity, hit)
  one <- d[d$ent == hit[7], ]; one <- one[order(one$ts), ]; one <- one[which(one$y == 1)[1]:nrow(one), ]
  expect_equal(e$start[7], one$ts[1]); expect_equal(e$rows[7], nrow(one)); expect_equal(e$events[7], sum(one$y))
  expect_equal(e$amount[7], sum(one$amount[one$y == 1])); expect_equal(e$peak_score[7], max(one$score))
  expect_equal(sum(e$amount), t$loss_total[1]); expect_equal(sum(e$detected), t$detected[1])
})

test_that("the loss accounting", {
  d <- det_df(4000, 300, 40, 7)
  dt <- scr_detection(d, entity = "ent", time = "ts", amount = "amount", thresholds = c(40, 70))
  t <- dt$table
  ev <- d$y == 1
  expect_equal(t$loss_total, rep(sum(d$amount[ev]), 2))
  expect_equal(t$pct_loss_prevented, 1 - t$loss_before / t$loss_total)
  # at the strictest threshold: the loss before detection of the detected, all of it for the others
  e <- dt$episodes
  expect_equal(t$loss_before[1], sum(e$loss_before[e$detected]) + sum(e$amount[!e$detected]))
  expect_true(all(e$loss_before[e$detected] <= e$amount[e$detected] + 1e-9))
  # nothing is detected: everything is lost; everything alerts: nothing is lost
  none <- scr_detection(d, entity = "ent", time = "ts", amount = "amount", thresholds = 1000)$table
  expect_equal(none$detected, 0); expect_equal(none$loss_before, none$loss_total); expect_equal(none$pct_loss_prevented, 0)
  every <- scr_detection(d, entity = "ent", time = "ts", amount = "amount", thresholds = -5)$table
  expect_equal(every$pct_detected, 1); expect_equal(every$loss_before, 0); expect_equal(every$median_time, 0)
  expect_equal(every$mean_events_before, 0); expect_equal(every$alert_share, 1)
  # without an amount there is no loss
  na <- scr_detection(d, entity = "ent", time = "ts", thresholds = 70)
  expect_true(is.na(na$table$loss_before)); expect_true(is.na(na$table$pct_loss_prevented))
  expect_true(all(is.na(na$episodes$amount)))
  expect_equal(na$table$detected, t$detected[1])
  # a missing amount counts as 0
  ma <- d; ma$amount[which(ev)[1:3]] <- NA
  expect_equal(scr_detection(ma, entity = "ent", time = "ts", amount = "amount", thresholds = 70)$table$loss_total,
               sum(d$amount[ev]) - sum(d$amount[which(ev)[1:3]]))
})

test_that("alert shares become tie-safe thresholds on all rows", {
  d <- det_df()
  sh <- c(0.02, 0.10, 0.30)
  dt <- scr_detection(d, entity = "ent", time = "ts", alert_shares = sh)
  t <- dt$table
  expect_equal(t$target_share, sh)
  # each threshold sits between two distinct scores, at the boundary nearest to the share
  u <- sort(unique(d$score))
  act <- vapply(u, function(v) mean(d$score >= v), numeric(1))
  for (i in 1:3) {
    best <- which.min(abs(act - sh[i]))
    expect_equal(t$threshold[i], (u[best] + u[best - 1]) / 2)
    expect_equal(t$alert_share[i], act[best]); expect_equal(t$alert_share[i], mean(d$score >= t$threshold[i]))
    expect_equal(t$detected[i], det_naive(d, t$threshold[i])$detected)
  }
  # thresholds and shares together, duplicates kept once
  both <- scr_detection(d, entity = "ent", time = "ts", thresholds = c(50, t$threshold[2]), alert_shares = 0.10)$table
  expect_equal(nrow(both), 2L); expect_equal(sort(both$threshold), sort(c(50, t$threshold[2])))
  # a share of 1 alerts every row
  al <- scr_detection(d, entity = "ent", time = "ts", alert_shares = 1)$table
  expect_equal(al$threshold, -Inf); expect_equal(al$alert_share, 1); expect_equal(al$pct_detected, 1)
  # from the low scores
  f <- d; f$score <- 100 - f$score
  tf <- scr_detection(f, entity = "ent", time = "ts", alert_shares = sh, direction = "higher_is_safer")$table
  expect_equal(tf$threshold, 100 - t$threshold); expect_equal(tf$alert_share, t$alert_share)
  expect_equal(tf$detected, t$detected)
  # every score tied: the nearer of no row and every row
  tie <- d; tie$score <- 5
  tt <- scr_detection(tie, entity = "ent", time = "ts", alert_shares = c(0.2, 0.9))$table
  expect_equal(tt$alert_share, c(0, 1)); expect_equal(tt$pct_detected, c(0, 1))
})

test_that("time: dates, ties in the input order, rows before the start", {
  d <- det_df(3000, 200, 30, 9)
  # a date gives the time to detection in days
  dd <- d; dd$ts <- as.Date("2026-01-01") + floor(d$ts * 4)
  dn <- d; dn$ts <- floor(d$ts * 4)
  a <- scr_detection(dd, entity = "ent", time = "ts", thresholds = c(50, 80))
  b <- scr_detection(dn, entity = "ent", time = "ts", thresholds = c(50, 80))
  expect_equal(a$table, b$table)
  expect_s3_class(a$episodes$start, "Date"); expect_s3_class(a$episodes$detect_time, "Date")
  expect_equal(as.numeric(a$episodes$detect_time - a$episodes$start), a$episodes$time_to_detect)
  # a date-time gives seconds
  dp <- d; dp$ts <- as.POSIXct("2026-01-01", tz = "UTC") + round(d$ts * 3600)
  dh <- d; dh$ts <- round(d$ts * 3600)
  expect_equal(scr_detection(dp, entity = "ent", time = "ts", thresholds = 50)$table$median_time,
               scr_detection(dh, entity = "ent", time = "ts", thresholds = 50)$table$median_time)
  # rows at the same time keep their input order
  e1 <- data.frame(ent = 1, ts = c(1, 1), y = c(1, 0), score = c(10, 99))
  r1 <- scr_detection(e1, entity = "ent", time = "ts", thresholds = 50)
  expect_equal(r1$table$detected, 1); expect_equal(r1$episodes$events_before, 1); expect_equal(r1$episodes$rows, 2L)
  r2 <- scr_detection(e1[2:1, ], entity = "ent", time = "ts", thresholds = 50)
  expect_equal(r2$table$detected, 0); expect_equal(r2$episodes$rows, 1L)
  # character and factor entities
  dc <- d; dc$ent <- paste0("e", d$ent)
  expect_equal(scr_detection(dc, entity = "ent", time = "ts", thresholds = 50)$table$detected,
               scr_detection(d, entity = "ent", time = "ts", thresholds = 50)$table$detected)
  dc$ent <- factor(dc$ent)
  expect_equal(scr_detection(dc, entity = "ent", time = "ts", thresholds = 50)$table$detected,
               scr_detection(d, entity = "ent", time = "ts", thresholds = 50)$table$detected)
})

test_that("edge cases: no episode, one entity, rows left out, bad input", {
  d <- det_df(2000, 150, 20, 4)
  # no event row: no episode, the alert shares stand
  z <- d; z$y <- 0L
  dz <- scr_detection(z, entity = "ent", time = "ts", amount = "amount", thresholds = c(50, 80))
  expect_equal(dz$table$episodes, c(0, 0)); expect_equal(dz$table$detected, c(0, 0))
  expect_true(all(is.na(dz$table$pct_detected))); expect_true(all(is.na(dz$table$pct_loss_prevented)))
  expect_equal(dz$table$alert_share, c(mean(d$score >= 80), mean(d$score >= 50)))
  expect_null(dz$episodes); expect_equal(dz$n_episode_rows, 0L)
  expect_silent(utils::capture.output(print(dz)))
  # every row is one entity
  o <- d; o$ent <- 1
  do <- scr_detection(o, entity = "ent", time = "ts", thresholds = 60)
  ref <- det_naive(o, 60)
  expect_equal(do$table$episodes, 1); expect_equal(do$table$detected, ref$detected)
  expect_equal(do$table$mean_events_before, ref$mean_events)
  # rows without an entity, a time or a score are left out; a missing outcome is not an event
  na <- d; na$ent[1:4] <- NA; na$ts[5:7] <- NA; na$score[8:9] <- NA; na$y[10:30] <- NA
  dn <- scr_detection(na, entity = "ent", time = "ts", thresholds = 60)
  cl <- na[-(1:9), ]; cl$y[is.na(cl$y)] <- 0L
  expect_equal(dn$n_dropped, 9L)
  expect_equal(dn$table, scr_detection(cl, entity = "ent", time = "ts", thresholds = 60)$table)
  expect_error(scr_detection(d, time = "ts", thresholds = 50), "`entity` and `time` are needed")
  expect_error(scr_detection(d[0, ], entity = "ent", time = "ts", thresholds = 50), "no row has an entity")
  expect_error(scr_detection(d, entity = "ent", time = "ts"), "give `thresholds`")
  expect_error(scr_detection(d, entity = "ent", time = "nope", thresholds = 50), "not in `x`")
  expect_error(scr_detection(d, entity = "ent", time = "ts", thresholds = NA_real_), "thresholds")
  expect_error(scr_detection(d, entity = "ent", time = "ts", alert_shares = 1.5), "alert_shares")
  expect_error(scr_detection(d, entity = "ent", time = "ts", thresholds = 50, direction = "up"), "direction")
  expect_error(scr_detection(transform(d, ts = as.character(ts)), entity = "ent", time = "ts", thresholds = 50),
               "a number, a Date or a POSIXct")
  expect_error(scr_detection(transform(d, amount = -amount), entity = "ent", time = "ts", amount = "amount",
                             thresholds = 50), "non-negative")
  expect_error(scr_detection(d, entity = "ent", time = "ts", thresholds = 50, foo = 1), "unused argument")
  expect_error(scr_detection(transform(d, score = NA_real_), entity = "ent", time = "ts", thresholds = 50),
               "no row has an entity")
})

test_that("pooled score cells: the alert shares of the thresholds stay exact", {
  d <- det_df(5000, 400, 50, 17)
  d$score <- d$score + stats::runif(nrow(d))
  dt <- scr_detection(d, entity = "ent", time = "ts", thresholds = c(40.3, 77.7), alert_shares = 0.05, max_cells = 100)
  t <- dt$table
  # given thresholds are forced as cell edges; a share lands on an edge
  expect_equal(t$alert_share, vapply(t$threshold, function(v) mean(d$score >= v), numeric(1)))
  expect_equal(t[!is.na(target_share), alert_share], 0.05, tolerance = 0.21)
  for (i in seq_len(nrow(t))) expect_equal(t$detected[i], det_naive(d, t$threshold[i])$detected)
})

test_that("print and export", {
  d <- det_df(3000, 200, 30, 2)
  dt <- scr_detection(d, entity = "ent", time = "ts", amount = "amount", alert_shares = c(0.05, 0.2))
  out <- utils::capture.output(print(dt))
  expect_match(out[1], "^<scr_detection> entity \"ent\" \\| time \"ts\"")
  expect_true(any(grepl("episodes over", out))); expect_true(any(grepl("loss prevented", out)))
  skip_if_not_installed("openxlsx")
  old <- scr_verbose(FALSE); on.exit(scr_verbose(old), add = TRUE)
  dir <- file.path(tempdir(), "scr-detection-export")
  unlink(dir, recursive = TRUE)
  ex <- scr_export(dt, dir, stamp = FALSE)
  expect_identical(basename(ex$files$xlsx), "detection_y.xlsx")
  expect_identical(openxlsx::getSheetNames(ex$files$xlsx), c("Thresholds", "Episodes", "Settings"))
  expect_equal(nrow(openxlsx::read.xlsx(ex$files$xlsx, sheet = "Episodes")), nrow(dt$episodes))
  # no episode: the workbook has no episode sheet
  z <- d; z$y <- 0L
  ez <- scr_export(scr_detection(z, entity = "ent", time = "ts", thresholds = 50), dir, stamp = FALSE)
  expect_identical(openxlsx::getSheetNames(ez$files$xlsx), c("Thresholds", "Settings"))
})
