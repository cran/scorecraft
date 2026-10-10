# Claims: exact binomial tests, Jeffreys bounds and Holm against base R,
# the verdict logic, score ranges, the floor, weights and the study inputs.

# Pool adjacent violators by repeated merging of the first violating pair,
# written independently of the package: the smoothed rate of every input
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

claims_df <- function(n = 5000, seed = 41) {
  set.seed(seed)
  x <- stats::rnorm(n)
  data.frame(score = round(500 + 50 * x), y = stats::rbinom(n, 1, stats::plogis(-0.5 + 1.2 * x)),
             smp = sample(c("dev", "oot"), n, TRUE), w = stats::runif(n, 0.2, 3), stringsAsFactors = FALSE)
}

test_that("p-values, bounds and Holm equal pbinom, qbeta and p.adjust on a manual filter", {
  d <- claims_df()
  cl <- data.frame(op = c(">=", ">=", "<=", "<=", ">="), rate = c(0.6, 0.85, 0.25, 0.08, 0.3),
                   score_lo = c(550, 560, NA, NA, 480), score_hi = c(NA, NA, 450, 430, 520))
  r <- scr_claims(d, cl, objective = "propensity", level = 0.9)
  t <- r$table
  expect_s3_class(r, "scr_claims")
  inr <- function(lo, hi) (is.na(lo) | d$score >= lo) & (is.na(hi) | d$score < hi)
  p <- pr <- bd <- numeric(5)
  for (k in 1:5) {
    i <- inr(cl$score_lo[k], cl$score_hi[k])
    x <- sum(d$y[i]); n <- sum(i)
    expect_equal(t$n[k], n); expect_equal(t$events[k], x); expect_equal(t$rate[k], x / n)
    if (cl$op[k] == ">=") {
      p[k] <- stats::pbinom(x - 1, n, cl$rate[k], lower.tail = FALSE)
      pr[k] <- stats::pbinom(x, n, cl$rate[k])
      bd[k] <- stats::qbeta(0.1, x + 0.5, n - x + 0.5)
    } else {
      p[k] <- stats::pbinom(x, n, cl$rate[k])
      pr[k] <- stats::pbinom(x - 1, n, cl$rate[k], lower.tail = FALSE)
      bd[k] <- stats::qbeta(0.9, x + 0.5, n - x + 0.5)
    }
    # the exact test of base R gives the same one-sided p-value
    alt <- if (cl$op[k] == ">=") "greater" else "less"
    expect_equal(t$p_value[k], stats::binom.test(x, n, cl$rate[k], alternative = alt)$p.value)
  }
  expect_equal(t$p_value, p); expect_equal(t$p_refute, pr); expect_equal(t$bound, bd)
  expect_equal(t$p_adj, stats::p.adjust(p, "holm"))
  expect_equal(t$p_refute_adj, stats::p.adjust(pr, "holm"))
  expect_identical(t$verdict, ifelse(t$p_adj < 0.1, "supported", ifelse(t$p_refute_adj < 0.1, "refuted", "not proven")))
  # no adjustment
  r0 <- scr_claims(d, cl, objective = "propensity", level = 0.9, adjust = "none")
  expect_equal(r0$table$p_adj, p)
  expect_equal(r0$table$p_refute_adj, pr)
  # the group text and the open ends
  expect_identical(t$group[1:3], c("score >= 550", "score >= 560", "score < 450"))
  expect_identical(t$group[5], "480 <= score < 520")
  expect_equal(t$score_lo[3], -Inf); expect_equal(t$score_hi[1], Inf)
})

test_that("the verdict: supported, refuted and not proven on constructed counts", {
  # one group of 1000 rows with 700 events at score 1, 1000 rows with 100 events at score 2
  cnt <- data.frame(score = c(1, 2), n = c(1000, 1000), events = c(700, 100))
  cl <- data.frame(op = c(">=", ">=", ">=", "<=", "<=", "<="), rate = c(0.65, 0.75, 0.69, 0.12, 0.08, 0.105),
                   score_lo = c(NA, NA, NA, 1.5, 1.5, 1.5), score_hi = c(1.5, 1.5, 1.5, NA, NA, NA))
  r <- scr_claims(cnt, cl, counts = TRUE, adjust = "none")
  expect_identical(r$table$verdict, c("supported", "refuted", "not proven", "supported", "refuted", "not proven"))
  # the bounds: lower for ">=", upper for "<=", from the Jeffreys posterior
  expect_equal(r$table$bound[1], stats::qbeta(0.05, 700.5, 300.5))
  expect_equal(r$table$bound[4], stats::qbeta(0.95, 100.5, 900.5))
  # the verdict of a claim and of its refutation can never both hold
  expect_false(any(r$table$p_adj < 0.05 & r$table$p_refute_adj < 0.05))
  # Holm can turn a marginal claim into "not proven"
  marg <- data.frame(op = ">=", rate = c(0.674, 0.674, 0.674, 0.674), score_hi = 1.5)
  p1 <- stats::pbinom(699, 1000, 0.674, lower.tail = FALSE)
  expect_lt(p1, 0.05); expect_gt(4 * p1, 0.05)
  expect_identical(unique(scr_claims(cnt, marg, counts = TRUE, adjust = "none")$table$verdict), "supported")
  expect_identical(unique(scr_claims(cnt, marg, counts = TRUE)$table$verdict), "not proven")
  # the statement names the sample, the rate, the bound and the verdict
  s <- r$table$statement[1]
  expect_match(s, "^On 'all' \\(n = 1,000\\), rows with score < 1.5 had an event rate of 70.0%")
  expect_match(s, "95% one-sided lower bound 67.6%\\); the claim 'rate >= 65%' is supported\\.$")
  expect_match(r$table$statement[4], "upper bound")
})

test_that("labels and numbers of bands and tiers select the same rows as the study table", {
  d <- claims_df()
  tr <- scr_tiers(d, sample = "smp", objective = "propensity", n_tiers = 3)
  cl <- data.frame(label = c("high", "1", "medium"), op = c(">=", "<=", ">="), rate = c(0.5, 0.2, 0.3),
                   name = c("high responds", NA, "medium"))
  r <- scr_claims(tr, cl)
  expect_identical(r$sample, "oot")
  t <- tr$table[sample == "oot"]
  i <- match(c("high", "low", "medium"), t$label)
  expect_equal(r$table$n, t$n[i]); expect_equal(r$table$events, t$events[i])
  expect_identical(r$table$group, c("high", "low", "medium"))
  expect_identical(r$table$name, c("high responds", "rate <= 20%", "medium"))
  expect_match(r$table$statement[1], "rows in tier 'high'")
  expect_match(r$table$statement[1], "the claim 'high responds' \\(rate >= 50%\\)")
  # another sample of the study, and the level of the study by default
  rd <- scr_claims(tr, cl, sample = "dev")
  expect_equal(rd$table$n, tr$table[sample == "dev"]$n[match(c("high", "low", "medium"), tr$table[sample == "dev"]$label)])
  expect_equal(rd$level, tr$level)
  # a band study: the band number and the interval label
  b <- scr_bands(d, objective = "propensity", n_bands = 5, n_boot = 0)
  rb <- scr_claims(b, data.frame(label = c(1, NA), op = ">=", rate = 0.5, score_lo = c(NA, b$table$score_lo[1])))
  expect_equal(rb$table$n[1], b$table$n[1])
  expect_equal(rb$table$n[2], b$table$n[1])
  rl <- scr_claims(b, data.frame(label = b$table$label[2], op = "<=", rate = 0.6))
  expect_equal(rl$table$events, b$table$events[2])
  # the numbered label of the tiers table selects the same tier
  expect_equal(scr_claims(tr, data.frame(label = "01.high", op = ">=", rate = 0.5))$table$n, r$table$n[1])
  expect_error(scr_claims(tr, data.frame(label = "huge", op = ">=", rate = 0.5)), "not a tier of the study")
  expect_error(scr_claims(tr, cl, sample = "zzz"), "not in the study")
  expect_error(scr_claims(tr, cl, foo = 1), "unused argument")
})

test_that("floor reads the weakest edge, found on the reference and tested on the study sample", {
  # propensity: the rate rises with the score; inside the group the low half is at 40%, the high half at 80%
  cnt <- data.frame(score = 1:10, n = 400, events = rep(c(160, 320), each = 5))
  cl <- data.frame(op = c(">=", "<="), rate = c(0.5, 0.7))
  av <- scr_claims(cnt, cl, counts = TRUE, objective = "propensity")
  fl <- scr_claims(cnt, cl, counts = TRUE, objective = "propensity", type = "floor")
  expect_equal(av$table$rate, c(0.6, 0.6))
  expect_identical(av$table$verdict, c("supported", "supported"))
  # ">=": the low-rate edge is the run of scores 1 to 5; "<=": the high-rate edge, scores 6 to 10
  expect_equal(fl$table$rate, c(0.4, 0.8))
  expect_equal(fl$table$n, c(2000, 2000))
  expect_equal(fl$table$score_hi[1], 5.5); expect_equal(fl$table$score_lo[2], 5.5)
  expect_identical(fl$table$verdict, c("refuted", "refuted"))
  expect_match(fl$table$statement[1], "low-rate edge \\(score < 5.5\\) of all rows")
  expect_match(fl$table$statement[2], "high-rate edge \\(score >= 5.5\\)")
  # under risk with higher_is_safer the event-poor end is the high scores
  rk <- data.frame(score = 10:1, n = 400, events = rep(c(160, 320), each = 5))
  fr <- scr_claims(rk, cl[1, ], counts = TRUE, objective = "risk", type = "floor")
  expect_equal(fr$table$rate, 0.4)
  expect_equal(fr$table$score_lo, 5.5)
  # a non-monotone group is smoothed first: the edge is the pooled block, not the single worst cell
  nm <- data.frame(score = 1:4, n = 1000, events = c(300, 200, 500, 600))
  fn <- scr_claims(nm, data.frame(op = ">=", rate = 0.2), counts = TRUE, objective = "propensity", type = "floor")
  expect_equal(fn$table$rate, 0.25)
  expect_equal(fn$table$n, 2000)
  # the edge comes from the reference: the study sample is read on the reference range
  d <- claims_df()
  b <- scr_bands(d, sample = "smp", objective = "propensity", n_bands = 4, n_boot = 0)
  f2 <- scr_claims(b, data.frame(op = ">=", rate = 0.5, score_lo = 520), type = "floor", floor_bins = 6)
  dev <- d[d$smp == "dev" & d$score >= 520, ]
  # the pre-bins: tie-safe equal shares of the group on the reference
  pc <- .study_cuts(.study_hist(dev$score, dev$y), 6, "uniform", NULL, "high")$cuts
  expect_false(any(pc %in% dev$score))
  pb <- findInterval(dev$score, pc) + 1L
  e <- as.numeric(tapply(dev$y, pb, sum)); n <- as.numeric(tapply(dev$y, pb, length))
  sm <- pav_rates(e, n)
  last <- max(which(sm <= sm[1]))
  hi <- if (last == length(n)) Inf else pc[last]
  i <- d$smp == "oot" & d$score >= 520 & d$score < hi
  expect_equal(f2$table$score_hi, hi)
  expect_equal(f2$table$n, sum(i))
  expect_equal(f2$table$events, sum(d$y[i]))
  # one pre-bin: the floor is the average
  expect_equal(scr_claims(b, data.frame(op = ">=", rate = 0.5, score_lo = 520), type = "floor", floor_bins = 1)$table$n,
               sum(d$smp == "oot" & d$score >= 520))
  expect_error(scr_claims(b, data.frame(op = ">=", rate = 0.5), floor_bins = 0), "floor_bins")
})

test_that("weights act through the Kish effective size, rounded only for the binomial", {
  d <- claims_df()
  cl <- data.frame(op = c(">=", "<="), rate = c(0.55, 0.15), score_lo = c(540, NA), score_hi = c(NA, 460))
  r <- scr_claims(d, cl, objective = "propensity", weight = "w")
  for (k in 1:2) {
    i <- (is.na(cl$score_lo[k]) | d$score >= cl$score_lo[k]) & (is.na(cl$score_hi[k]) | d$score < cl$score_hi[k])
    ny <- sum(d$w[i]); e <- sum(d$w[i] * d$y[i]); neff <- ny^2 / sum(d$w[i]^2)
    x <- e / ny * neff
    expect_equal(r$table$n[k], ny); expect_equal(r$table$events[k], e); expect_equal(r$table$n_eff[k], neff)
    if (cl$op[k] == ">=") {
      expect_equal(r$table$p_value[k], stats::pbinom(round(x) - 1, round(neff), 0.55, lower.tail = FALSE))
      expect_equal(r$table$bound[k], stats::qbeta(0.05, x + 0.5, neff - x + 0.5))
    } else {
      expect_equal(r$table$p_value[k], stats::pbinom(round(x), round(neff), 0.15))
      expect_equal(r$table$bound[k], stats::qbeta(0.95, x + 0.5, neff - x + 0.5))
    }
  }
  # the statement quotes the effective n of the test next to the weighted volume
  i <- d$score >= 540
  ny <- sum(d$w[i]); neff <- ny^2 / sum(d$w[i]^2)
  expect_match(r$table$statement[1],
               sprintf("On 'all' (effective n = %s of a weighted volume of %s), rows with score >= 540",
                       formatC(round(neff), big.mark = ",", format = "f", digits = 0),
                       formatC(round(ny), big.mark = ",", format = "f", digits = 0)), fixed = TRUE)
  expect_lt(neff, ny)
  expect_match(scr_claims(d, cl, objective = "propensity")$table$statement[1],
               sprintf("On 'all' (n = %s), rows", formatC(sum(i), big.mark = ",", format = "d")), fixed = TRUE)
  # unit weights give the unweighted result; zero weights leave the rows out
  expect_equal(scr_claims(transform(d, w = 1), cl, objective = "propensity", weight = "w")$table$p_value,
               scr_claims(d, cl, objective = "propensity")$table$p_value)
  dz <- d; dz$w[dz$score >= 540][1:50] <- 0
  rz <- scr_claims(dz, cl, objective = "propensity", weight = "w")
  expect_equal(rz$table$n[1], sum(dz$w[dz$score >= 540]))
})

test_that("score ranges stay exact when the scores are pooled into cells", {
  set.seed(7)
  n <- 20000
  x <- stats::rnorm(n)
  d <- data.frame(score = 500 + 50 * x, y = stats::rbinom(n, 1, stats::plogis(x)))
  cl <- data.frame(op = ">=", rate = 0.5, score_lo = 512.3456, score_hi = 600.25)
  r <- scr_claims(d, cl, objective = "propensity", max_cells = 50)
  i <- d$score >= 512.3456 & d$score < 600.25
  expect_equal(r$table$n, sum(i)); expect_equal(r$table$events, sum(d$y[i]))
  # a study pooled without the edges cannot answer a range that cuts a cell
  b <- scr_bands(d, objective = "propensity", max_cells = 50, n_boot = 0)
  expect_true(b$quantized)
  expect_error(scr_claims(b, cl), "fall inside a pooled cell")
})

test_that("scr_claims on a scorecard reads the hold-out and the configured level", {
  sc <- sc_demo()
  ho <- sc$samples$holdout
  cut <- stats::median(ho$score)
  cl <- data.frame(op = c("<=", ">="), rate = c(0.10, 0.20), score_lo = c(cut, NA), score_hi = c(NA, cut))
  r <- scr_claims(sc, cl)
  expect_identical(r$sample, "holdout"); expect_identical(r$reference, "train")
  expect_equal(r$table$n, c(sum(ho$score >= cut), sum(ho$score < cut)))
  expect_equal(r$table$events, c(sum(ho$y[ho$score >= cut]), sum(ho$y[ho$score < cut])))
  sc2 <- sc; sc2$config$study_level <- 0.9
  expect_equal(scr_claims(sc2, cl)$level, 0.9)
  expect_identical(scr_claims(sc, cl, sample = "train")$sample, "train")
  expect_output(print(r), "scr_claims")
  expect_output(print(r), "adjustment holm")
  expect_output(print(r), "On 'holdout'")
})

test_that("edge cases: single class, all ties, empty group, bad claims", {
  one <- data.frame(score = 1:200, y = 0)
  r <- scr_claims(one, data.frame(op = c(">=", "<="), rate = c(0.1, 0.1)))
  expect_identical(r$table$verdict, c("refuted", "supported"))
  expect_equal(r$table$bound, c(0, stats::qbeta(0.95, 0.5, 200.5)))
  tie <- data.frame(score = rep(5, 300), y = rep(c(1, 0, 0), 100))
  rt <- scr_claims(tie, data.frame(op = ">=", rate = 0.3), type = "floor")
  expect_equal(rt$table$rate, 1 / 3); expect_equal(rt$table$n, 300)
  # a floor edge that is the whole group says so instead of repeating the range
  expect_match(rt$table$statement, "rows at the low-rate edge (the whole group) of all rows", fixed = TRUE)
  flat <- data.frame(score = 1:4, n = 100, events = 50)
  rw <- scr_claims(flat, data.frame(op = ">=", rate = 0.4, score_lo = 1, score_hi = 5), counts = TRUE,
                   objective = "propensity", type = "floor")
  expect_match(rw$table$statement, "low-rate edge (the whole group) of rows with 1 <= score < 5", fixed = TRUE)
  expect_equal(rw$table$n, 400)
  # an empty group: no row, no test, "not proven"
  re <- scr_claims(claims_df(), data.frame(op = ">=", rate = 0.5, score_lo = 1e6))
  expect_equal(re$table$n, 0); expect_true(is.na(re$table$p_value)); expect_true(is.na(re$table$bound))
  expect_identical(re$table$verdict, "not proven")
  expect_match(re$table$statement, "had no row with a known outcome")
  # a floor on a group the reference does not reach
  rf <- scr_claims(claims_df(), data.frame(op = ">=", rate = 0.5, score_lo = 1e6), type = "floor")
  expect_identical(rf$table$verdict, "not proven")
  d <- claims_df(500)
  expect_error(scr_claims(d, data.frame(op = ">", rate = 0.5)), "`claims\\$op`")
  expect_error(scr_claims(d, data.frame(op = ">=", rate = 1)), "in \\(0, 1\\)")
  expect_error(scr_claims(d, data.frame(op = ">=")), "lacks the column")
  expect_error(scr_claims(d, data.frame(op = ">=", rate = 0.5, score_lo = 2, score_hi = 1)), "score_lo")
  expect_error(scr_claims(d, data.frame(op = ">=", rate = 0.5, label = 1, score_lo = 2)), "not both")
  expect_error(scr_claims(d, d[0, ]), "one row per claim")
  expect_error(scr_claims(d, data.frame(op = ">=", rate = 0.5), study = "x"), "need a `sample` column")
})

test_that("scr_export writes the claims workbook", {
  skip_if_not_installed("openxlsx")
  old <- scr_verbose(FALSE); on.exit(scr_verbose(old), add = TRUE)
  r <- scr_claims(claims_df(800), data.frame(op = ">=", rate = 0.5, score_lo = 520), objective = "propensity")
  out <- file.path(tempdir(), "scr-claims-export")
  unlink(out, recursive = TRUE)
  ex <- scr_export(r, out, stamp = FALSE)
  expect_identical(openxlsx::getSheetNames(ex$files$xlsx), c("Claims", "Settings"))
  expect_equal(nrow(openxlsx::read.xlsx(ex$files$xlsx, sheet = "Claims")), 1L)
})
