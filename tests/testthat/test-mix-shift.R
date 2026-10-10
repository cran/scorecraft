# Mix and rate effects: the decomposition against a row-level computation,
# a hand example, the tests against base R, the PSI and the inputs.

mix_df <- function(n = 4000, seed = 71) {
  set.seed(seed)
  per <- sample(c("p1", "p2", "p3"), n, TRUE)
  x <- stats::rnorm(n, mean = c(p1 = 0, p2 = -0.3, p3 = 0.2)[per])
  data.frame(per = per, score = round(600 + 50 * x),
             y = stats::rbinom(n, 1, stats::plogis(-2 - x + 0.4 * (per == "p2"))),
             w = stats::runif(n, 0.2, 3), stringsAsFactors = FALSE)
}

# effects of one comparison from the rows, written independently of the package
mix_manual <- function(sb, yb, sc, yc, cuts, wb = rep(1, length(sb)), wc = rep(1, length(sc))) {
  B <- length(cuts) + 1L
  ib <- factor(findInterval(sb, cuts) + 1L, seq_len(B)); ic <- factor(findInterval(sc, cuts) + 1L, seq_len(B))
  nb <- as.numeric(tapply(wb, ib, sum)); nc <- as.numeric(tapply(wc, ic, sum))
  eb <- as.numeric(tapply(wb * yb, ib, sum)); ec <- as.numeric(tapply(wc * yc, ic, sum))
  nb[is.na(nb)] <- 0; nc[is.na(nc)] <- 0; eb[is.na(eb)] <- 0; ec[is.na(ec)] <- 0
  pb <- nb / sum(nb); pc <- nc / sum(nc)
  rb <- eb / nb; rc <- ec / nc
  list(pb = pb, pc = pc, rb = rb, rc = rc, mix = (pc - pb) * (rb + rc) / 2, rate = (rc - rb) * (pb + pc) / 2,
       nb = nb, nc = nc, eb = eb, ec = ec)
}

test_that("the effects equal a row-level computation and add up to the change of the rate", {
  d <- mix_df()
  ms <- scr_mix_shift(d, by = "per", n_bands = 6)
  expect_s3_class(ms, "scr_mix_shift")
  expect_identical(ms$base, "p1"); expect_identical(ms$groups, c("p2", "p3"))
  b <- d[d$per == "p1", ]
  for (g in c("p2", "p3")) {
    cc <- d[d$per == g, ]
    m <- mix_manual(b$score, b$y, cc$score, cc$y, ms$cuts)
    t <- ms$table[group == g]
    # higher_is_safer: band 1 is the lowest interval, the event-richest
    expect_identical(t$band, 1:6)
    expect_equal(t$pct_base, m$pb); expect_equal(t$pct_cmp, m$pc)
    expect_equal(t$rate_base, m$rb); expect_equal(t$rate_cmp, m$rc)
    expect_equal(t$mix_effect, m$mix); expect_equal(t$rate_effect, m$rate)
    expect_equal(t$total, m$mix + m$rate)
    expect_equal(t$n_base, m$nb); expect_equal(t$n_cmp, m$nc)
    s <- ms$summary[group == g]
    expect_equal(s$rate_base, mean(b$y)); expect_equal(s$rate_cmp, mean(cc$y))
    expect_equal(s$delta, mean(cc$y) - mean(b$y))
    # exact additivity, to rounding
    expect_equal(sum(t$mix_effect) + sum(t$rate_effect), s$delta, tolerance = 1e-13)
    expect_equal(s$mix_total + s$rate_total, s$delta, tolerance = 1e-13)
    expect_equal(s$share_mix, sum(m$mix) / s$delta)
    expect_equal(s$n_base, nrow(b)); expect_equal(s$n_cmp, nrow(cc))
  }
  # the bands are frozen on the base, tie-safe: the cuts of scr_bands() on the base rows
  expect_equal(ms$cuts, scr_bands(b, n_bands = 6, n_boot = 0)$cuts)
  expect_identical(ms$table[group == "p2", label], .study_labels(ms$cuts))
})

test_that("a hand example gives the hand numbers", {
  # score 1 is the risky band, score 2 the safe one
  cnt <- data.frame(smp = c("base", "base", "new", "new"), score = c(1, 2, 1, 2), n = c(100, 100, 150, 50),
                    events = c(30, 10, 60, 10))
  ms <- scr_mix_shift(cnt, by = "smp", counts = TRUE, breaks = 1.5)
  t <- ms$table
  expect_identical(t$band, 1:2)
  expect_equal(t$pct_base, c(0.5, 0.5)); expect_equal(t$pct_cmp, c(0.75, 0.25))
  expect_equal(t$rate_base, c(0.3, 0.1)); expect_equal(t$rate_cmp, c(0.4, 0.2))
  expect_equal(t$mix_effect, c(0.0875, -0.0375)); expect_equal(t$rate_effect, c(0.0625, 0.0375))
  expect_equal(t$total, c(0.15, 0))
  s <- ms$summary
  expect_equal(s$rate_base, 0.2); expect_equal(s$rate_cmp, 0.35); expect_equal(s$delta, 0.15)
  expect_equal(s$mix_total, 0.05); expect_equal(s$rate_total, 0.10); expect_equal(s$share_mix, 1 / 3)
  # the PSI of the shares and its n-adjusted critical value
  expect_equal(s$psi, (0.5 - 0.75) * log(0.5 / 0.75) + (0.5 - 0.25) * log(0.5 / 0.25))
  expect_equal(s$psi_critical, (1 / 200 + 1 / 200) * stats::qchisq(0.95, 1))
  # under higher_is_riskier the high score comes first
  r <- scr_mix_shift(cnt, by = "smp", counts = TRUE, breaks = 1.5, direction = "higher_is_riskier")
  expect_equal(r$table$rate_base, c(0.1, 0.3))
  expect_equal(r$summary$mix_total, 0.05)
})

test_that("a band empty in one sample takes the rate of the other", {
  cnt <- data.frame(smp = c("a", "a", "b", "b", "b"), score = c(1, 2, 1, 2, 3), n = c(100, 100, 80, 80, 40),
                    events = c(40, 10, 40, 16, 2))
  ms <- scr_mix_shift(cnt, by = "smp", counts = TRUE, breaks = c(1.5, 2.5))
  t <- ms$table
  # band 3 (score 3) does not exist in the base
  expect_equal(t$pct_base, c(0.5, 0.5, 0)); expect_true(is.na(t$rate_base[3]))
  expect_equal(t$rate_effect[3], 0)
  expect_equal(t$mix_effect[3], 0.2 * 0.05)
  expect_true(is.na(t$p_rate[3]))
  expect_equal(sum(t$total), ms$summary$delta, tolerance = 1e-13)
  expect_equal(ms$summary$delta, 58 / 200 - 50 / 200)
  # and the other way round: the comparison lacks a band of the base
  rv <- scr_mix_shift(cnt, by = "smp", counts = TRUE, breaks = c(1.5, 2.5), base = "b")
  tr <- rv$table
  expect_true(is.na(tr$rate_cmp[3])); expect_equal(tr$rate_effect[3], 0)
  expect_equal(tr$mix_effect[3], -0.2 * 0.05)
  expect_equal(sum(tr$total), -ms$summary$delta, tolerance = 1e-13)
  # the decomposition is antisymmetric in the two samples
  expect_equal(tr$mix_effect, -t$mix_effect); expect_equal(tr$rate_effect, -t$rate_effect)
  # the Holm adjustment leaves the missing test out
  expect_equal(t$p_rate_adj[1:2], stats::p.adjust(t$p_rate[1:2], "holm"))
})

test_that("the band tests are Fisher's exact test or the two-proportion z test", {
  # small expected counts: Fisher; large: the pooled z, equal to prop.test without correction
  e1 <- c(3, 40, 0, 7, 120); n1 <- c(20, 400, 15, 9, 1000)
  e2 <- c(9, 65, 0, 2, 90); n2 <- c(25, 420, 30, 8, 1100)
  p <- .mix_rate_test(e1, n1, e2, n2)
  ref <- vapply(seq_along(e1), function(i) {
    m <- matrix(c(e1[i], n1[i] - e1[i], e2[i], n2[i] - e2[i]), 2)
    N <- n1[i] + n2[i]; k <- e1[i] + e2[i]
    if (min(n1[i], n2[i]) * min(k, N - k) / N < 5) stats::fisher.test(m)$p.value else
      suppressWarnings(stats::prop.test(c(e1[i], e2[i]), c(n1[i], n2[i]), correct = FALSE)$p.value)
  }, numeric(1))
  expect_equal(p, ref)
  # both kinds occur, and a pair without events has p = 1
  expect_equal(p[3], 1)
  expect_true(is.na(.mix_rate_test(1, 0, 2, 10)))
  # in the study: the unweighted counts of each band, Holm across the bands
  d <- mix_df(1500, 5)
  ms <- scr_mix_shift(d, by = "per", n_bands = 8, compare = "p2", weight = "w")
  t <- ms$table
  b <- d[d$per == "p1", ]; cc <- d[d$per == "p2", ]
  m <- mix_manual(b$score, b$y, cc$score, cc$y, ms$cuts)
  expect_equal(t$p_rate, .mix_rate_test(m$eb, m$nb, m$ec, m$nc))
  expect_equal(t$p_rate_adj, stats::p.adjust(t$p_rate, "holm"))
  expect_equal(ms$summary$bands_changed, sum(t$p_rate_adj < 0.05))
})

test_that("by-groups: every group against the base, as the pairwise calls", {
  d <- mix_df()
  ms <- scr_mix_shift(d, by = "per", n_bands = 5, base = "p2")
  expect_identical(ms$groups, c("p1", "p3"))
  for (g in ms$groups) {
    one <- scr_mix_shift(d[d$per %in% c("p2", g), ], by = "per", n_bands = 5, base = "p2")
    expect_equal(ms$table[group == g], one$table)
    expect_equal(ms$summary[group == g], one$summary)
  }
  # `compare` selects, a factor keeps its level order, and rows without a group are left out
  expect_identical(scr_mix_shift(d, by = "per", compare = "p3")$groups, "p3")
  f <- d; f$per <- factor(f$per, levels = c("p3", "p1", "p2"))
  expect_identical(scr_mix_shift(f, by = "per")$base, "p3")
  na <- d; na$per[1:50] <- NA
  expect_equal(scr_mix_shift(na, by = "per", n_bands = 5)$summary,
               scr_mix_shift(d[-(1:50), ], by = "per", n_bands = 5)$summary)
  # a date column is a period
  dd <- d; dd$per <- as.Date(c(p1 = "2026-01-01", p2 = "2026-02-01", p3 = "2026-03-01")[d$per])
  r <- scr_mix_shift(dd, by = "per", n_bands = 5)
  expect_identical(r$base, "2026-01-01")
  expect_equal(r$summary$delta, scr_mix_shift(d, by = "per", n_bands = 5)$summary$delta)
})

test_that("the PSI is that of the band shares, with the n-adjusted critical value", {
  d <- mix_df()
  ms <- scr_mix_shift(d, by = "per", n_bands = 7, level = 0.9)
  b <- d[d$per == "p1", ]; cc <- d[d$per == "p3", ]
  m <- mix_manual(b$score, b$y, cc$score, cc$y, ms$cuts)
  s <- ms$summary[group == "p3"]
  expect_equal(s$psi, sum((m$pb - m$pc) * log(m$pb / m$pc)))
  expect_equal(s$psi_critical, (1 / nrow(b) + 1 / nrow(cc)) * stats::qchisq(0.9, 6))
  # the same index as scr_psi() on the rows, with the frozen cuts
  expect_equal(s$psi, scr_psi(b$score, cc$score, breaks = c(-Inf, ms$cuts, Inf))$psi)
})

test_that("weights: weighted shares and rates, and whole weights as repeated rows", {
  d <- mix_df(2500, 9)
  ms <- scr_mix_shift(d, by = "per", n_bands = 5, weight = "w", compare = "p2")
  b <- d[d$per == "p1", ]; cc <- d[d$per == "p2", ]
  m <- mix_manual(b$score, b$y, cc$score, cc$y, ms$cuts, b$w, cc$w)
  t <- ms$table
  expect_equal(t$mix_effect, m$mix); expect_equal(t$rate_effect, m$rate)
  expect_equal(ms$summary$delta, stats::weighted.mean(cc$y, cc$w) - stats::weighted.mean(b$y, b$w))
  expect_equal(sum(t$total), ms$summary$delta, tolerance = 1e-13)
  # the PSI runs on the Kish effective sizes
  kb <- sum(b$w)^2 / sum(b$w^2); kc <- sum(cc$w)^2 / sum(cc$w^2)
  expect_equal(ms$summary$psi_critical, (1 / kb + 1 / kc) * stats::qchisq(0.95, 4))
  # whole weights: the effects of the repeated rows
  d$k <- sample(1:3, nrow(d), TRUE)
  rep_d <- d[rep(seq_len(nrow(d)), d$k), ]
  br <- ms$cuts
  a <- scr_mix_shift(d, by = "per", weight = "k", breaks = br)
  r <- scr_mix_shift(rep_d, by = "per", breaks = br)
  expect_equal(a$table[, .(mix_effect, rate_effect, total)], r$table[, .(mix_effect, rate_effect, total)])
  # a zero weight takes the row out; a group left empty has nothing to decompose
  z <- d; z$w[z$per == "p3"] <- 0
  ze <- scr_mix_shift(z, by = "per", weight = "w", compare = "p3")
  expect_true(is.na(ze$summary$delta)); expect_equal(ze$summary$n_cmp, 0)
  expect_true(all(is.na(ze$table$total)))
  z2 <- d; z2$w0 <- as.numeric(seq_len(nrow(d)) %% 2)
  expect_equal(scr_mix_shift(z2, by = "per", weight = "w0", breaks = br)$table$total,
               scr_mix_shift(z2[z2$w0 > 0, ], by = "per", breaks = br)$table$total)
})

test_that("a score study and a scorecard are read as they are", {
  d <- mix_df()
  bd <- scr_bands(d, sample = "per", n_bands = 6, n_boot = 0)
  ms <- scr_mix_shift(bd)
  expect_equal(ms$table, scr_mix_shift(d, by = "per", n_bands = 6)$table)
  expect_identical(ms$base, "p1")
  # tiers: the tier numbers and labels, the event-richest tier first
  tr <- scr_tiers(d, sample = "per", n_tiers = 3)
  mt <- scr_mix_shift(tr, compare = "p2")
  expect_identical(mt$table$label, tr$table[sample == "p2", label])
  expect_identical(mt$table$band, tr$table[sample == "p2", tier])
  expect_equal(mt$table$rate_cmp, tr$table[sample == "p2", rate])
  expect_equal(sum(mt$table$total), mt$summary$delta, tolerance = 1e-13)
  expect_error(scr_mix_shift(scr_bands(d, n_bands = 4, n_boot = 0)), "nothing to compare")
  expect_error(scr_mix_shift(bd, base = "zz"), "not in the study")
  expect_error(scr_mix_shift(bd, n_bands = 3), "unused argument")

  sc <- sc_demo()
  m <- scr_mix_shift(sc, n_bands = 5)
  expect_identical(m$base, "train"); expect_identical(m$groups, "holdout")
  tr_ <- sc$samples$train; ho <- sc$samples$holdout
  mm <- mix_manual(tr_$score, tr_$y, ho$score, ho$y, m$cuts)
  expect_equal(m$table$mix_effect, mm$mix); expect_equal(m$table$rate_effect, mm$rate)
  expect_equal(m$summary$delta, mean(ho$y) - mean(tr_$y))
  expect_identical(m$target, "default")
  # by the stored dates: the scored rows of both samples, grouped by period
  md <- scr_mix_shift(sc, by = "date", n_bands = 5)
  al <- rbind(tr_, ho)
  dates <- sort(unique(as.character(al$date)))
  expect_identical(md$base, dates[1]); expect_identical(md$groups, dates[-1])
  b <- al[as.character(date) == dates[1]]; cc <- al[as.character(date) == dates[3]]
  mm <- mix_manual(b$score, b$y, cc$score, cc$y, md$cuts)
  t <- md$table[group == dates[3]]
  expect_equal(t$mix_effect, mm$mix); expect_equal(t$rate_effect, mm$rate)
  expect_equal(md$summary[group == dates[3], delta], mean(cc$y) - mean(b$y))
  expect_error(scr_mix_shift(sc, by = "date", base = "1999-01-01"), "not in the column")
  expect_error(scr_mix_shift(sc, by = "nope"), "not in the scored sample")
})

test_that("edge cases: one band, a single class, tied scores, missing outcomes, bad input", {
  d <- mix_df(1200, 3)
  # one band: no mix, the whole change is a rate effect
  one <- scr_mix_shift(d, by = "per", n_bands = 1, compare = "p2")
  expect_equal(nrow(one$table), 1L)
  expect_equal(one$table$mix_effect, 0); expect_equal(one$table$rate_effect, one$summary$delta)
  expect_true(is.na(one$summary$psi))
  # every score tied: one band whatever was asked
  tie <- d; tie$score <- 500
  expect_equal(nrow(scr_mix_shift(tie, by = "per", n_bands = 10, compare = "p3")$table), 1L)
  # a single class: no change, no share
  z <- d; z$y <- 0L
  r <- scr_mix_shift(z, by = "per", n_bands = 4)
  expect_equal(r$summary$delta, c(0, 0)); expect_true(all(is.na(r$summary$share_mix)))
  expect_equal(r$table$rate_effect, rep(0, 8)); expect_equal(r$table$p_rate, rep(1, 8))
  # rows with a missing outcome leave the shares and the rates
  na <- d; na$y[seq(1, nrow(d), 7)] <- NA
  expect_equal(scr_mix_shift(na, by = "per", breaks = c(560, 600, 640))$table,
               scr_mix_shift(na[!is.na(na$y), ], by = "per", breaks = c(560, 600, 640))$table)
  # a comparison without a known outcome: nothing to decompose
  nk <- d; nk$y[nk$per == "p3"] <- NA
  s <- scr_mix_shift(nk, by = "per", n_bands = 4)$summary
  expect_true(is.na(s$delta[2])); expect_false(is.na(s$delta[1]))
  expect_error(scr_mix_shift(d), "needs `by`")
  # no row, or no row with a group
  expect_error(scr_mix_shift(d[0, ], by = "per"), "`x` has no row with a value of 'per'")
  expect_error(scr_mix_shift(transform(d, per = NA_character_), by = "per"), "`x` has no row with a value of 'per'")
  expect_error(scr_mix_shift(d, by = "nope"), "not in `x`")
  expect_error(scr_mix_shift(d[d$per == "p1", ], by = "per"), "nothing to compare")
  expect_error(scr_mix_shift(d, by = "per", base = "zz"), "not in the column")
  expect_error(scr_mix_shift(d, by = "per", level = 1), "level")
  expect_error(scr_mix_shift(d, by = "per", n_bands = 0), "n_bands")
  expect_error(scr_mix_shift(d, by = "per", foo = 1), "unused argument")
  expect_error(scr_mix_shift(nk, by = "per", base = "p3"), "no row with a known outcome")
})

test_that("a data.table is filtered by its rows, whatever its columns are called", {
  d <- mix_df(2000, 31)
  d$per[1:60] <- NA
  ref <- scr_mix_shift(d, by = "per", n_bands = 5)
  # a column called `x`, and the group column itself called `by`: neither may be read as the filter
  dx <- data.table::as.data.table(d); dx[, x := seq_len(.N)]
  expect_equal(scr_mix_shift(dx, by = "per", n_bands = 5)$table, ref$table)
  db <- data.table::as.data.table(d); data.table::setnames(db, "per", "by")
  rb <- scr_mix_shift(db, by = "by", n_bands = 5)
  expect_equal(rb$table, ref$table); expect_identical(rb$by, "by")
  dk <- data.table::as.data.table(d); dk[, `:=`(keep = 0L, cols = "a", cn = 1)]
  expect_equal(scr_mix_shift(dk, by = "per", n_bands = 5, weight = "w")$summary,
               scr_mix_shift(d, by = "per", n_bands = 5, weight = "w")$summary)
  # the caller's table is left as it was
  expect_identical(names(dx), c(names(d), "x")); expect_equal(nrow(dx), nrow(d))
  # counts, with rows without a group
  cnt <- data.table::data.table(by = c("a", "a", "b", "b", NA), score = c(1, 2, 1, 2, 1), n = c(100, 100, 150, 50, 7),
                                events = c(30, 10, 60, 10, 1))
  expect_equal(scr_mix_shift(cnt, by = "by", counts = TRUE, breaks = 1.5)$summary$delta, 0.15)
})

test_that("numeric periods are read in numeric order", {
  d <- mix_df(3000, 32)
  d$k <- c(p1 = 9, p2 = 10, p3 = 100)[d$per]
  ms <- scr_mix_shift(d, by = "k", n_bands = 5)
  # as text "10" would come before "9" and become the base
  expect_identical(ms$base, "9"); expect_identical(ms$groups, c("10", "100"))
  expect_equal(ms$summary$delta, scr_mix_shift(d, by = "per", n_bands = 5)$summary$delta)
  expect_identical(scr_mix_shift(d, by = "k", base = 10)$groups, c("9", "100"))
  # a factor keeps its levels, text and dates the order of their labels
  f <- d; f$k <- factor(f$k, levels = c(100, 9, 10))
  expect_identical(scr_mix_shift(f, by = "k")$base, "100")
  ch <- d; ch$k <- as.character(ch$k)
  expect_identical(scr_mix_shift(ch, by = "k")$base, "10")
  # the scored samples of a scorecard grouped by a numeric column
  sc <- sc_demo()
  s2 <- sc
  s2$samples <- lapply(sc$samples, function(s) { s <- data.table::copy(s); s[, wk := 8 + as.integer(format(date, "%m")) * 2L]; s })
  m2 <- scr_mix_shift(s2, by = "wk", n_bands = 4)
  expect_identical(m2$base, "10"); expect_identical(m2$groups, as.character(c(12, 14, 16, 18, 20)))
  expect_equal(m2$summary$delta, scr_mix_shift(sc, by = "date", n_bands = 4)$summary$delta)
})

test_that("pooled score cells keep the decomposition exact", {
  d <- mix_df(3000, 12)
  d$score <- d$score + stats::runif(nrow(d))
  # more distinct scores than cells: the cuts sit on cell edges, so the bands hold whole cells
  ms <- scr_mix_shift(d, by = "per", n_bands = 5, max_cells = 40, compare = "p2")
  b <- d[d$per == "p1", ]; cc <- d[d$per == "p2", ]
  m <- mix_manual(b$score, b$y, cc$score, cc$y, ms$cuts)
  expect_equal(ms$table$mix_effect, m$mix); expect_equal(ms$table$rate_effect, m$rate)
  expect_equal(sum(ms$table$total), ms$summary$delta, tolerance = 1e-13)
  expect_length(ms$cuts, 4L)
})

test_that("print and export", {
  d <- mix_df(1500, 2)
  ms <- scr_mix_shift(d, by = "per", n_bands = 5)
  out <- utils::capture.output(print(ms))
  expect_match(out[1], "^<scr_mix_shift> target \"y\"")
  expect_true(any(grepl("rate [0-9.]+% -> [0-9.]+% \\([+-][0-9.]+ pp\\): mix [+-][0-9.]+ pp, rate [+-][0-9.]+ pp", out)))
  expect_equal(sum(grepl("^Largest band effects", out)), 2L)
  skip_if_not_installed("openxlsx")
  old <- scr_verbose(FALSE); on.exit(scr_verbose(old), add = TRUE)
  dir <- file.path(tempdir(), "scr-mix-export")
  unlink(dir, recursive = TRUE)
  ex <- scr_export(ms, dir, stamp = FALSE)
  expect_identical(basename(ex$files$xlsx), "mix_shift_y.xlsx")
  expect_identical(openxlsx::getSheetNames(ex$files$xlsx), c("Summary", "Bands", "Cuts", "Settings"))
  expect_equal(nrow(openxlsx::read.xlsx(ex$files$xlsx, sheet = "Bands")), 10L)
})
