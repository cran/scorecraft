# Rules against the score: every set against logical operations on the
# rows, the retire flag, the value sums, the tie-safe alert share.

ovl_df <- function(n = 8000, seed = 111) {
  set.seed(seed)
  x <- stats::rnorm(n)
  y <- stats::rbinom(n, 1, stats::plogis(-4 + 1.5 * x))
  data.frame(score = round(100 * stats::plogis(x + stats::rnorm(n, sd = 0.5))), y = y,
             amount = round(stats::rexp(n, 1 / 80), 2),
             r_fast = as.integer(x > 1.5), r_dev = stats::rbinom(n, 1, ifelse(y == 1, 0.3, 0.02)),
             r_night = stats::runif(n) < 0.05, w = stats::runif(n, 0.2, 3), stringsAsFactors = FALSE)
}

jeff <- function(x, n, level = 0.95) {
  a <- (1 - level) / 2
  c(if (x <= 0) 0 else stats::qbeta(a, x + 0.5, n - x + 0.5), if (x >= n) 1 else stats::qbeta(1 - a, x + 0.5, n - x + 0.5))
}

test_that("every set equals the logical operation on the rows", {
  d <- ovl_df()
  rules <- c("r_fast", "r_dev", "r_night")
  ov <- scr_overlap(d, rules = rules, cut = 80, level = 0.9)
  expect_s3_class(ov, "scr_overlap")
  t <- ov$table
  expect_identical(t$rule, rules)
  a <- d$score >= 80
  for (k in seq_along(rules)) {
    r <- as.logical(d[[rules[k]]])
    sets <- list(rule = r, both = r & a, rule_only = r & !a, score_only = a & !r)
    for (s in names(sets)) {
      i <- sets[[s]]
      n <- sum(i); x <- sum(d$y[i])
      expect_equal(t[[paste0("n_", s)]][k], n, info = paste(rules[k], s))
      expect_equal(t[[paste0("events_", s)]][k], x, info = paste(rules[k], s))
      expect_equal(t[[paste0("precision_", s)]][k], x / n, info = paste(rules[k], s))
      ci <- jeff(x, n, 0.9)
      expect_equal(t[[paste0("precision_", s, "_lo")]][k], ci[1]); expect_equal(t[[paste0("precision_", s, "_hi")]][k], ci[2])
    }
    expect_equal(t$caught_by_score[k], sum(d$y[r & a]) / sum(d$y[r]))
    expect_identical(t$retire_candidate[k], sum(d$y[r & a]) / sum(d$y[r]) >= 0.95)
    # the rule splits into the part the score alerts and the rest
    expect_equal(t$n_both[k] + t$n_rule_only[k], t$n_rule[k])
    expect_equal(t$n_both[k] + t$n_score_only[k], sum(a))
  }
  # the summary: the score, any rule, either
  anyr <- d$r_fast == 1 | d$r_dev == 1 | d$r_night
  s <- ov$summary
  E <- sum(d$y)
  expect_equal(c(s$n, s$events), c(nrow(d), E))
  expect_equal(c(s$n_score, s$n_rules, s$n_any), c(sum(a), sum(anyr), sum(a | anyr)))
  expect_equal(c(s$share_score, s$share_rules, s$share_any), c(mean(a), mean(anyr), mean(a | anyr)))
  expect_equal(c(s$precision_score, s$precision_rules, s$precision_any), c(mean(d$y[a]), mean(d$y[anyr]), mean(d$y[a | anyr])))
  expect_equal(c(s$recall_score, s$recall_rules, s$recall_any), c(sum(d$y[a]), sum(d$y[anyr]), sum(d$y[a | anyr])) / E)
  expect_equal(s$incr_score, sum(d$y[a & !anyr]) / E); expect_equal(s$incr_rules, sum(d$y[anyr & !a]) / E)
  expect_identical(ov$n_patterns, nrow(unique(data.frame(a, d$r_fast, d$r_dev, d$r_night))))
  # under higher_is_safer the low scores are alerted: score < cut
  f <- d; f$score <- 100 - f$score
  of <- scr_overlap(f, rules = rules, cut = 20.5, direction = "higher_is_safer")
  expect_equal(of$summary$n_score, sum(f$score < 20.5)); expect_equal(of$summary$n_score, sum(d$score >= 80))
  expect_equal(of$table$n_both, t$n_both); expect_equal(of$table$caught_by_score, t$caught_by_score)
  # the direction follows the objective when left open
  expect_identical(scr_overlap(f, rules = rules, cut = 20.5, direction = NULL)$direction, "higher_is_safer")
})

test_that("the retire flag follows the share of the events of the rule caught by the score", {
  d <- ovl_df()
  a <- d$score >= 70
  # a rule inside the alerts of the score, one half inside, one without events
  d$r_in <- as.integer(a & d$y == 1 & seq_len(nrow(d)) %% 2 == 0)
  d$r_none <- as.integer(d$y == 0 & seq_len(nrow(d)) %% 50 == 0)
  d$r_off <- 0L
  ov <- scr_overlap(d, rules = c("r_in", "r_dev", "r_none", "r_off"), cut = 70)
  t <- ov$table
  expect_equal(t$caught_by_score[1], 1); expect_true(t$retire_candidate[1])
  expect_equal(t$n_rule_only[1], 0); expect_true(is.na(t$precision_rule_only[1]))
  cf <- sum(d$y[d$r_dev == 1 & a]) / sum(d$y[d$r_dev == 1])
  expect_lt(cf, 0.95)
  expect_equal(t$caught_by_score[2], cf); expect_false(t$retire_candidate[2])
  # the threshold decides: at or above `retire_at`
  expect_true(scr_overlap(d, rules = "r_dev", cut = 70, retire_at = cf)$table$retire_candidate)
  expect_false(scr_overlap(d, rules = "r_dev", cut = 70, retire_at = cf + 1e-6)$table$retire_candidate)
  # a rule without events, and one that never fires: nothing to judge
  expect_equal(t$events_rule[3], 0); expect_true(is.na(t$caught_by_score[3])); expect_true(is.na(t$retire_candidate[3]))
  expect_equal(t$precision_rule[3], 0); expect_equal(t$precision_rule_lo[3], 0)
  expect_equal(t$n_rule[4], 0); expect_true(is.na(t$precision_rule[4])); expect_true(is.na(t$retire_candidate[4]))
  expect_equal(t$n_score_only[4], sum(a))
})

test_that("the value sums are those of the events of each set", {
  d <- ovl_df()
  ov <- scr_overlap(d, rules = c("r_fast", "r_dev"), cut = 75, value = "amount")
  a <- d$score >= 75
  ve <- d$amount * d$y
  for (k in 1:2) {
    r <- d[[c("r_fast", "r_dev")[k]]] == 1
    expect_equal(ov$table$value_rule[k], sum(ve[r])); expect_equal(ov$table$value_both[k], sum(ve[r & a]))
    expect_equal(ov$table$value_rule_only[k], sum(ve[r & !a])); expect_equal(ov$table$value_score_only[k], sum(ve[a & !r]))
    expect_equal(ov$table$caught_by_score_value[k], sum(ve[r & a]) / sum(ve[r]))
  }
  anyr <- d$r_fast == 1 | d$r_dev == 1
  s <- ov$summary
  expect_equal(s$value_events, sum(ve))
  expect_equal(c(s$value_recall_score, s$value_recall_rules, s$value_recall_any),
               c(sum(ve[a]), sum(ve[anyr]), sum(ve[a | anyr])) / sum(ve))
  expect_equal(s$value_incr_score, sum(ve[a & !anyr]) / sum(ve)); expect_equal(s$value_incr_rules, sum(ve[anyr & !a]) / sum(ve))
  # without a value there is no value column
  expect_false(any(grepl("^value", names(scr_overlap(d, rules = "r_fast", cut = 75)$table))))
  # a missing amount counts as 0
  na <- d; na$amount[d$y == 1][1:5] <- NA
  expect_equal(scr_overlap(na, rules = "r_fast", cut = 75, value = "amount")$summary$value_events,
               sum(ve) - sum(d$amount[d$y == 1][1:5]))
})

test_that("the alert share gives a tie-safe cut and reports the share realized", {
  d <- ovl_df()
  ov <- scr_overlap(d, rules = "r_fast", alert_share = 0.05)
  # the boundary between two distinct scores nearest to the share
  u <- sort(unique(d$score))
  sh <- vapply(u, function(v) mean(d$score >= v), numeric(1))
  best <- which.min(abs(sh - 0.05))
  expect_equal(ov$cut, (u[best] + u[best - 1]) / 2)
  expect_equal(ov$summary$share_score, sh[best]); expect_equal(ov$summary$n_score, sum(d$score >= ov$cut))
  expect_identical(ov$alert_share, 0.05)
  # heavy ties: no tie is split, the share realized can sit far from the one asked
  t3 <- d; t3$score <- c(1, 2, 3)[1 + (d$score > 40) + (d$score > 70)]
  o3 <- scr_overlap(t3, rules = "r_fast", alert_share = 0.05)
  expect_equal(o3$cut, 2.5); expect_equal(o3$summary$share_score, mean(t3$score == 3))
  # every score tied: the nearer of nothing and everything
  t1 <- d; t1$score <- 9
  expect_equal(scr_overlap(t1, rules = "r_fast", alert_share = 0.3)$summary$n_score, 0)
  expect_equal(scr_overlap(t1, rules = "r_fast", alert_share = 0.8)$summary$n_score, nrow(d))
  expect_equal(scr_overlap(d, rules = "r_fast", alert_share = 1)$summary$share_score, 1)
  # from the low scores under higher_is_safer
  f <- d; f$score <- 100 - f$score
  of <- scr_overlap(f, rules = "r_fast", alert_share = 0.05, direction = "higher_is_safer")
  expect_equal(of$summary$n_score, ov$summary$n_score); expect_equal(of$cut, 100 - ov$cut)
  # the same cut as a tail band of scr_bands()
  expect_equal(ov$cut, scr_bands(d, spacing = "tail", tail_probs = 0.05, direction = "higher_is_riskier", n_boot = 0)$cuts)
})

test_that("weights, missing outcomes and missing flags", {
  d <- ovl_df()
  ow <- scr_overlap(d, rules = c("r_fast", "r_dev"), cut = 78, weight = "w", value = "amount")
  a <- d$score >= 78; r <- d$r_fast == 1
  expect_equal(ow$table$n_rule[1], sum(d$w[r])); expect_equal(ow$table$events_both[1], sum(d$w[r & a] * d$y[r & a]))
  p <- sum(d$w[r] * d$y[r]) / sum(d$w[r])
  expect_equal(ow$table$precision_rule[1], p)
  # the Jeffreys interval on the Kish effective size
  k <- sum(d$w[r])^2 / sum(d$w[r]^2)
  expect_equal(c(ow$table$precision_rule_lo[1], ow$table$precision_rule_hi[1]),
               c(stats::qbeta(0.025, p * k + 0.5, k - p * k + 0.5), stats::qbeta(0.975, p * k + 0.5, k - p * k + 0.5)))
  expect_equal(ow$table$value_both[1], sum((d$w * d$amount * d$y)[r & a]))
  expect_equal(ow$summary$recall_score, sum(d$w[a] * d$y[a]) / sum(d$w * d$y))
  expect_true(ow$weighted)
  # a zero weight and a missing score take the row out
  z <- d; z$w0 <- as.numeric(seq_len(nrow(d)) %% 2); z$score[c(2, 4, 6)] <- NA
  oz <- scr_overlap(z, rules = c("r_fast", "r_dev"), cut = 78, weight = "w0")
  keep <- z$w0 > 0 & !is.na(z$score)
  expect_equal(oz$table, scr_overlap(z[keep, ], rules = c("r_fast", "r_dev"), cut = 78)$table)
  expect_equal(oz$n_dropped, sum(!keep))
  # a missing outcome counts in the volume, not in the precision
  na <- d; na$y[seq(1, nrow(d), 9)] <- NA
  on <- scr_overlap(na, rules = "r_fast", cut = 78)
  expect_equal(on$table$n_rule, sum(r)); expect_equal(on$table$precision_rule, mean(na$y[r], na.rm = TRUE))
  # a missing flag is a rule that did not fire
  nf <- d; nf$r_fast[1:500] <- NA
  expect_equal(scr_overlap(nf, rules = "r_fast", cut = 78)$table$n_rule, sum(d$r_fast[-(1:500)]))
})

test_that("edge cases and bad input", {
  d <- ovl_df(2000, 5)
  # no event: no recall, the volumes stand
  z <- d; z$y <- 0L
  oz <- scr_overlap(z, rules = c("r_fast", "r_night"), cut = 80)
  expect_true(is.na(oz$summary$recall_score)); expect_true(is.na(oz$summary$incr_score))
  expect_equal(oz$summary$n_score, sum(d$score >= 80)); expect_equal(oz$table$precision_rule, c(0, 0))
  expect_true(all(is.na(oz$table$caught_by_score)))
  # a cut beyond the scores alerts nothing; one below alerts everything
  o0 <- scr_overlap(d, rules = "r_fast", cut = 1e6)
  expect_equal(o0$summary$n_score, 0); expect_equal(o0$table$n_both, 0); expect_equal(o0$table$caught_by_score, 0)
  o1 <- scr_overlap(d, rules = "r_fast", cut = -1)
  expect_equal(o1$table$caught_by_score, 1); expect_equal(o1$table$n_rule_only, 0)
  expect_equal(o1$summary$incr_rules, 0)
  expect_error(scr_overlap(d, cut = 80), "`rules` must be")
  expect_error(scr_overlap(d, rules = c("r_fast", "r_fast"), cut = 80), "`rules` must be")
  expect_error(scr_overlap(d, rules = "nope", cut = 80), "not in `x`")
  expect_error(scr_overlap(d, rules = "r_fast"), "one of them")
  expect_error(scr_overlap(d, rules = "r_fast", cut = 80, alert_share = 0.1), "one of them")
  expect_error(scr_overlap(d, rules = "r_fast", alert_share = 0), "^scr_overlap\\(\\): `alert_share` is outside")
  expect_error(scr_overlap(d, rules = "r_fast", alert_share = NA_real_), "^scr_overlap\\(\\): `alert_share` must be")
  expect_error(scr_overlap(d, rules = "amount", cut = 80), "0/1 or logical")
  expect_error(scr_overlap(d, rules = "r_fast", cut = 80, retire_at = 0), "^scr_overlap\\(\\): `retire_at` is outside")
  expect_error(scr_overlap(d[0, ], rules = "r_fast", cut = 80), "no row has a score")
  expect_error(scr_overlap(d, rules = "r_fast", cut = 80, foo = 1), "unused argument")
  expect_error(scr_overlap(transform(d, score = NA_real_), rules = "r_fast", cut = 80), "no row has a score")
})

test_that("pooled score cells: the alert cut sits on a cell edge and the sets stay exact", {
  d <- ovl_df(5000, 16)
  d$score <- d$score + stats::runif(nrow(d))
  ov <- scr_overlap(d, rules = c("r_fast", "r_dev"), alert_share = 0.1, max_cells = 100)
  a <- d$score >= ov$cut
  # within one cell of the share asked for
  expect_equal(ov$summary$share_score, 0.1, tolerance = 0.11)
  expect_equal(ov$summary$n_score, sum(a))
  expect_equal(ov$table$events_both, c(sum(d$y[a & d$r_fast == 1]), sum(d$y[a & d$r_dev == 1])))
})

test_that("print and export", {
  d <- ovl_df(3000, 2)
  ov <- scr_overlap(d, rules = c("r_fast", "r_dev"), alert_share = 0.05, value = "amount")
  out <- utils::capture.output(print(ov))
  expect_match(out[1], "^<scr_overlap> outcome \"y\" \\| score \"score\" \\(higher_is_riskier\\) \\| 2 rules")
  expect_true(any(grepl("incremental recall: the score over the rules", out)))
  expect_true(any(grepl("^  r_fast ", out)))
  skip_if_not_installed("openxlsx")
  old <- scr_verbose(FALSE); on.exit(scr_verbose(old), add = TRUE)
  dir <- file.path(tempdir(), "scr-overlap-export")
  unlink(dir, recursive = TRUE)
  ex <- scr_export(ov, dir, stamp = FALSE)
  expect_identical(basename(ex$files$xlsx), "overlap_y.xlsx")
  expect_identical(openxlsx::getSheetNames(ex$files$xlsx), c("Rules", "Summary", "Settings"))
  expect_equal(nrow(openxlsx::read.xlsx(ex$files$xlsx, sheet = "Rules")), 2L)
})
