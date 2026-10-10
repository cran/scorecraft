# Uplift: the band numbers against manual filters, the Newcombe interval
# against its formula and a published value, the Qini and AUUC against a
# row-level computation, the bootstrap and the randomization check.

upl_df <- function(n = 5000, seed = 101, p_treat = 0.5) {
  set.seed(seed)
  x <- stats::rnorm(n)
  tr <- stats::rbinom(n, 1, p_treat)
  data.frame(score = round(500 + 50 * x), treat = tr,
             y = stats::rbinom(n, 1, stats::plogis(-1.5 + 0.5 * x + tr * pmax(x, 0))),
             w = stats::runif(n, 0.2, 3), stringsAsFactors = FALSE)
}

wilson <- function(x, n, level = 0.95) {
  z <- stats::qnorm(1 - (1 - level) / 2); p <- x / n
  mid <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  c(mid - half, mid + half)
}

newcombe <- function(x1, n1, x2, n2, level = 0.95) {
  w1 <- wilson(x1, n1, level); w2 <- wilson(x2, n2, level)
  p1 <- x1 / n1; p2 <- x2 / n2
  c(p1 - p2 - sqrt((p1 - w1[1])^2 + (w2[2] - p2)^2), p1 - p2 + sqrt((w1[2] - p1)^2 + (p2 - w2[1])^2))
}

# Qini coefficient and AUUC from the rows: one point per distinct score,
# every count taken again from the rows selected at that score
qini_rows <- function(score, y, treat, w = rep(1, length(score)), high = TRUE) {
  s <- if (high) score else -score
  u <- sort(unique(s), decreasing = TRUE)
  NT <- sum(w[treat == 1]); NC <- sum(w[treat == 0]); N <- NT + NC
  nt <- nc <- q <- up <- numeric(length(u))
  for (k in seq_along(u)) {
    i <- s >= u[k]
    nt[k] <- sum(w[i & treat == 1]); nc[k] <- sum(w[i & treat == 0])
    et <- sum(w[i & treat == 1] * y[i & treat == 1]); ec <- sum(w[i & treat == 0] * y[i & treat == 0])
    q[k] <- et - if (nc[k] > 0) ec * nt[k] / nc[k] else 0
    up[k] <- ((if (nt[k] > 0) et / nt[k] else 0) - (if (nc[k] > 0) ec / nc[k] else 0)) * (nt[k] + nc[k]) / N
  }
  trap <- function(x, v) sum(diff(c(0, x)) * (c(0, v[-length(v)]) + v) / 2)
  list(qini = trap(nt, q) / NT^2 - q[length(q)] / (2 * NT), auuc = trap((nt + nc) / N, up), q = q, nt = nt, nc = nc)
}

test_that("the band numbers equal manual filters and the Newcombe formula", {
  d <- upl_df()
  up <- scr_uplift(d, n_bands = 5, n_boot = 0, level = 0.9)
  expect_s3_class(up, "scr_uplift")
  t <- up$table
  # propensity: the high scores first
  expect_identical(t$band, 1:5); expect_identical(up$direction, "higher_is_riskier")
  expect_true(all(diff(t$score_lo) < 0))
  for (k in 1:5) {
    i <- d$score >= t$score_lo[k] & d$score < t$score_hi[k]
    xt <- sum(d$y[i & d$treat == 1]); nt <- sum(i & d$treat == 1)
    xc <- sum(d$y[i & d$treat == 0]); nc <- sum(i & d$treat == 0)
    expect_equal(t$n_t[k], nt); expect_equal(t$n_c[k], nc)
    expect_equal(t$events_t[k], xt); expect_equal(t$events_c[k], xc)
    expect_equal(t$rate_t[k], xt / nt); expect_equal(t$rate_c[k], xc / nc)
    expect_equal(t$uplift[k], xt / nt - xc / nc)
    ci <- newcombe(xt, nt, xc, nc, 0.9)
    expect_equal(t$uplift_lo[k], ci[1]); expect_equal(t$uplift_hi[k], ci[2])
    expect_equal(t$incremental[k], (xt / nt - xc / nc) * nt)
  }
  expect_equal(t$cum_incremental, cumsum(t$incremental))
  expect_equal(t$cum_pct_treated, cumsum(t$n_t) / sum(d$treat))
  expect_identical(t$type, ifelse(t$uplift_lo > 0, "persuadable", ifelse(t$uplift_hi < 0, "negative", "no effect")))
  # the whole sample
  s <- up$summary
  xt <- sum(d$y[d$treat == 1]); nt <- sum(d$treat); xc <- sum(d$y[d$treat == 0]); nc <- sum(d$treat == 0)
  expect_equal(s$uplift, xt / nt - xc / nc)
  expect_equal(c(s$uplift_lo, s$uplift_hi), newcombe(xt, nt, xc, nc, 0.9))
  expect_equal(c(s$n_t, s$n_c), c(nt, nc))
  # the bands are those of the pooled score, tie-safe
  expect_equal(up$cuts, scr_bands(d, n_bands = 5, n_boot = 0, objective = "propensity")$cuts)
  expect_identical(t$label, rev(.study_labels(up$cuts)))
})

test_that("the Newcombe interval reproduces the published example", {
  # Newcombe (1998), 56/70 against 48/80, method 10: 0.0524 to 0.3339
  r <- .uplift_diff(56, 70, 70, 48, 80, 80, 0.95)
  expect_equal(r$uplift, 0.2)
  expect_equal(round(c(r$lo, r$hi), 4), c(0.0524, 0.3339))
  expect_equal(c(r$lo, r$hi), newcombe(56, 70, 48, 80))
  # the Wilson interval is the score interval of prop.test() without correction
  w <- .uplift_wilson(c(56, 0, 80), c(70, 30, 80), 0.95)
  expect_equal(c(w$lo[1], w$hi[1]), as.numeric(stats::prop.test(56, 70, correct = FALSE)$conf.int))
  expect_equal(w$lo[2], 0); expect_equal(w$hi[3], 1)
  expect_true(is.na(.uplift_wilson(0, 0, 0.95)$lo))
})

test_that("the band type follows the interval", {
  # three score values: a clear gain, nothing, a clear loss
  cnt <- data.frame(score = rep(c(3, 2, 1), each = 4), treat = rep(c(1, 1, 0, 0), 3), y = rep(c(1, 0, 1, 0), 3),
                    k = c(300, 200, 100, 400,  250, 250, 250, 250,  100, 400, 300, 200))
  rows <- cnt[rep(seq_len(nrow(cnt)), cnt$k), ]
  up <- scr_uplift(rows, cuts = c(1.5, 2.5), n_boot = 0)
  expect_identical(up$table$type, c("persuadable", "no effect", "negative"))
  expect_equal(up$table$uplift, c(0.4, 0, -0.4))
  expect_equal(c(up$table$uplift_lo[1], up$table$uplift_hi[1]), newcombe(300, 500, 100, 500))
  expect_equal(up$table$incremental, c(200, 0, -200)); expect_equal(up$table$cum_incremental, c(200, 200, 0))
  # case weights: the same rates, the interval on the Kish effective sizes
  uw <- scr_uplift(cnt, weight = "k", cuts = c(1.5, 2.5), n_boot = 0)
  expect_equal(uw$table$uplift, c(0.4, 0, -0.4)); expect_equal(uw$table$incremental, c(200, 0, -200))
  kish <- function(w) sum(w)^2 / sum(w^2)
  ci <- newcombe(0.6 * kish(c(300, 200)), kish(c(300, 200)), 0.2 * kish(c(100, 400)), kish(c(100, 400)))
  expect_equal(c(uw$table$uplift_lo[1], uw$table$uplift_hi[1]), ci)
  # two effective rows per arm say nothing
  expect_identical(uw$table$type, rep("no effect", 3))
})

test_that("Qini and AUUC equal a row-level computation", {
  d <- upl_df()
  up <- scr_uplift(d, n_bands = 5, n_boot = 0)
  ref <- qini_rows(d$score, d$y, d$treat)
  expect_equal(up$summary$qini, ref$qini); expect_equal(up$summary$auuc, ref$auuc)
  # the curve: one row per score value here, ending at the overall incremental events
  cv <- up$curve
  expect_equal(nrow(cv), length(unique(d$score)))
  expect_equal(cv$qini, ref$q); expect_equal(cv$n_t, ref$nt); expect_equal(cv$n_c, ref$nc)
  K <- nrow(cv)
  expect_equal(cv$qini[K], sum(d$y[d$treat == 1]) - sum(d$y[d$treat == 0]) * sum(d$treat) / sum(d$treat == 0))
  expect_equal(cv$uplift[K], up$summary$uplift)
  expect_equal(cv$qini_random, cv$qini[K] * cv$pct_treated); expect_equal(cv$pct_treated[K], 1)
  # the rows selected at each cut are those of the curve
  expect_equal(cv$n_t[10], sum(d$treat == 1 & d$score >= cv$cut[10]))
  # weights, and a score that ranks the other way round
  uw <- scr_uplift(d, n_bands = 5, n_boot = 0, weight = "w")
  rw <- qini_rows(d$score, d$y, d$treat, d$w)
  expect_equal(uw$summary$qini, rw$qini); expect_equal(uw$summary$auuc, rw$auuc)
  f <- d; f$score <- -f$score
  uf <- scr_uplift(f, n_bands = 5, n_boot = 0, direction = "higher_is_safer")
  expect_equal(uf$summary$qini, ref$qini); expect_equal(uf$summary$auuc, ref$auuc)
  expect_equal(uf$table$uplift, up$table$uplift)
  # a score that orders at random: Qini near 0, AUUC near half the uplift; the true order: positive
  set.seed(3)
  rnd <- d; rnd$score <- stats::runif(nrow(d))
  ur <- scr_uplift(rnd, n_bands = 5, n_boot = 0)
  expect_lt(abs(ur$summary$qini), up$summary$qini / 4)
  expect_equal(ur$summary$auuc, ur$summary$uplift / 2, tolerance = 0.25)
  expect_gt(up$summary$qini, 0)
  # more than 1,000 score values: the curve is thinned, the areas are not
  cont <- d; cont$score <- d$score + stats::runif(nrow(d))
  uc <- scr_uplift(cont, n_bands = 5, n_boot = 0)
  expect_lte(nrow(uc$curve), 1001L); expect_equal(uc$curve$pct_treated[nrow(uc$curve)], 1)
  expect_equal(uc$summary$qini, qini_rows(cont$score, cont$y, cont$treat)$qini)
})

test_that("the bootstrap is reproducible, stratified by arm, and leaves the user's stream alone", {
  d <- upl_df(3000, 5)
  a <- scr_uplift(d, n_bands = 5, n_boot = 150, seed = 11)
  b <- scr_uplift(d, n_bands = 5, n_boot = 150, seed = 11)
  expect_identical(a$summary, b$summary)
  s <- a$summary
  expect_true(s$qini_lo < s$qini && s$qini < s$qini_hi); expect_true(s$auuc_lo < s$auuc && s$auuc < s$auuc_hi)
  expect_false(isTRUE(all.equal(s$qini_lo, scr_uplift(d, n_bands = 5, n_boot = 150, seed = 12)$summary$qini_lo)))
  expect_identical(a$n_boot, 150L)
  set.seed(77); r1 <- stats::runif(3)
  set.seed(77); invisible(scr_uplift(d, n_bands = 5, n_boot = 40, seed = 2)); r2 <- stats::runif(3)
  expect_identical(r1, r2)
  # no bootstrap: no interval, the point estimates unchanged
  z <- scr_uplift(d, n_bands = 5, n_boot = 0)
  expect_true(is.na(z$summary$qini_lo)); expect_equal(z$summary$qini, s$qini); expect_identical(z$n_boot, 0L)
  # the resampled intervals agree with a row bootstrap within each arm
  set.seed(9)
  grp <- split(seq_len(nrow(d)), d$treat)
  rb <- vapply(1:150, function(i) {
    j <- unlist(lapply(grp, function(g) g[sample.int(length(g), replace = TRUE)]))
    r <- qini_rows(d$score[j], d$y[j], d$treat[j])
    c(r$qini, r$auuc)
  }, numeric(2))
  expect_equal(c(s$qini_lo, s$qini_hi), as.numeric(stats::quantile(rb[1, ], c(0.025, 0.975))), tolerance = 0.15)
  expect_equal(c(s$auuc_lo, s$auuc_hi), as.numeric(stats::quantile(rb[2, ], c(0.025, 0.975))), tolerance = 0.15)
  # the event totals of the arms are free: the AUUC interval is about as wide as that of half the uplift
  expect_equal(s$auuc_hi - s$auuc_lo, (s$uplift_hi - s$uplift_lo) / 2, tolerance = 0.35)
  # pooled cells: the interval stays around the full point estimate
  p <- scr_uplift(d, n_bands = 5, n_boot = 150, seed = 11, boot_cells = 20)
  expect_equal(p$summary$qini, s$qini)
  expect_true(p$summary$qini_lo < s$qini && s$qini < p$summary$qini_hi)
  expect_equal(p$summary$qini_hi - p$summary$qini_lo, s$qini_hi - s$qini_lo, tolerance = 0.3)
})

test_that("the bootstrap intervals cover the population values at their nominal level", {
  skip_on_cran()
  gen <- function(n) {
    x <- stats::rnorm(n); tr <- stats::rbinom(n, 1, 0.5)
    data.frame(score = round(500 + 50 * x), treat = tr,
               y = stats::rbinom(n, 1, stats::plogis(-1.5 + 0.5 * x + tr * pmax(x, 0))))
  }
  # the population values, from a sample large enough to stand for them
  set.seed(2024)
  truth <- scr_uplift(gen(1e6), n_boot = 0)$summary
  set.seed(7)
  hit <- vapply(1:200, function(i) {
    s <- scr_uplift(gen(1500), n_boot = 150)$summary
    c(qini = s$qini_lo <= truth$qini && truth$qini <= s$qini_hi, auuc = s$auuc_lo <= truth$auuc && truth$auuc <= s$auuc_hi)
  }, logical(2))
  cover <- rowMeans(hit)
  # nominal 95%; 200 samples estimate a coverage to within about 0.03 (two standard errors), and the
  # percentile interval is a first-order one: 0.88 is the floor. A resample that fixed the event
  # totals of the arms covers the AUUC about 58% of the time
  expect_gte(cover[["auuc"]], 0.88); expect_gte(cover[["qini"]], 0.88)
  expect_lte(cover[["auuc"]], 0.995); expect_lte(cover[["qini"]], 0.995)
})

test_that("the randomization check flags arms that differ along the score", {
  d <- upl_df(6000, 8)
  up <- scr_uplift(d, n_bands = 8, n_boot = 0)
  # PSI of the treated against the control over the bands, with the n-adjusted critical value
  ib <- factor(findInterval(d$score, up$cuts) + 1L, seq_len(8))
  pt <- as.numeric(table(ib[d$treat == 1])) / sum(d$treat); pc <- as.numeric(table(ib[d$treat == 0])) / sum(d$treat == 0)
  expect_equal(up$randomization$psi, sum((pt - pc) * log(pt / pc)))
  expect_equal(up$randomization$critical, (1 / sum(d$treat) + 1 / sum(d$treat == 0)) * stats::qchisq(0.95, 7))
  expect_identical(up$randomization$flag, up$randomization$psi > up$randomization$critical)
  expect_false(up$randomization$flag); expect_true(is.na(up$randomization$note))
  expect_equal(up$summary$psi, up$randomization$psi); expect_false(up$summary$psi_flag)
  expect_false(any(grepl("WARNING", utils::capture.output(print(up)))))
  # the treatment went to the high scores: not a random assignment
  set.seed(4)
  b <- d; b$treat <- stats::rbinom(nrow(d), 1, stats::plogis((d$score - 500) / 60))
  ub <- scr_uplift(b, n_bands = 8, n_boot = 0)
  expect_true(ub$randomization$flag); expect_true(ub$summary$psi_flag)
  expect_match(ub$randomization$note, "the score distribution differs between the arms")
  expect_true(any(grepl("WARNING: the score distribution differs", utils::capture.output(print(ub)))))
})

test_that("cuts, a score study and whole weights", {
  d <- upl_df(3000, 6)
  uc <- scr_uplift(d, cuts = c(470, 530), n_boot = 0)
  expect_equal(uc$cuts, c(470, 530)); expect_equal(nrow(uc$table), 3L)
  expect_equal(uc$table$n_t[1], sum(d$treat == 1 & d$score >= 530))
  bd <- scr_bands(d, n_bands = 4, n_boot = 0, objective = "propensity")
  us <- scr_uplift(d, cuts = bd, n_boot = 0)
  expect_equal(us$cuts, bd$cuts); expect_identical(us$table$label, bd$table$label)
  expect_equal(us$table$n_t + us$table$n_c, bd$table$n)
  expect_error(scr_uplift(d, cuts = bd, objective = "risk"), "fitted with")
  # whole weights give the rates and the areas of the repeated rows
  d$k <- sample(1:3, nrow(d), TRUE)
  rp <- d[rep(seq_len(nrow(d)), d$k), ]
  a <- scr_uplift(d, cuts = c(470, 530), weight = "k", n_boot = 0)
  r <- scr_uplift(rp, cuts = c(470, 530), n_boot = 0)
  expect_equal(a$table[, .(n_t, n_c, rate_t, rate_c, uplift, incremental)], r$table[, .(n_t, n_c, rate_t, rate_c, uplift, incremental)])
  expect_equal(a$summary$qini, r$summary$qini); expect_equal(a$summary$auuc, r$summary$auuc)
  # a zero weight takes the row out
  d$w0 <- as.numeric(seq_len(nrow(d)) %% 2)
  expect_equal(scr_uplift(d, cuts = c(470, 530), weight = "w0", n_boot = 0)$table,
               scr_uplift(d[d$w0 > 0, ], cuts = c(470, 530), n_boot = 0)$table)
})

test_that("edge cases: a single class, one band, tied scores, an empty arm, missing outcomes", {
  d <- upl_df(1500, 3)
  # no event at all: no uplift, a flat curve
  z <- d; z$y <- 0L
  uz <- scr_uplift(z, n_bands = 3, n_boot = 20, seed = 1)
  expect_equal(uz$table$uplift, rep(0, 3)); expect_identical(uz$table$type, rep("no effect", 3))
  expect_equal(uz$summary$qini, 0); expect_equal(uz$summary$auuc, 0)
  expect_equal(c(uz$summary$qini_lo, uz$summary$qini_hi), c(0, 0))
  # one band is the whole sample
  one <- scr_uplift(d, n_bands = 1, n_boot = 0)
  expect_equal(nrow(one$table), 1L); expect_equal(one$table$uplift, one$summary$uplift)
  expect_true(is.na(one$randomization$psi)); expect_false(one$randomization$flag)
  # every score tied: one cell, the curve is the straight line
  tie <- d; tie$score <- 1
  ut <- scr_uplift(tie, n_bands = 5, n_boot = 20, seed = 1)
  expect_length(ut$cuts, 0L); expect_equal(nrow(ut$curve), 1L)
  expect_equal(ut$summary$qini, 0); expect_equal(ut$summary$auuc, ut$summary$uplift / 2)
  # one cell: the Qini curve is its own straight line in every resample, but the AUUC is half
  # the uplift, whose uncertainty the interval must carry
  ut2 <- scr_uplift(tie, n_bands = 5, n_boot = 300, seed = 1)$summary
  expect_equal(c(ut2$qini_lo, ut2$qini_hi), c(0, 0))
  expect_gt(ut2$auuc_hi - ut2$auuc_lo, 0.5 * (ut2$uplift_hi - ut2$uplift_lo) / 2)
  expect_equal(c(ut2$auuc_lo, ut2$auuc_hi), c(ut2$uplift_lo, ut2$uplift_hi) / 2, tolerance = 0.15)
  expect_true(ut2$auuc_lo < ut2$auuc && ut2$auuc < ut2$auuc_hi)
  # a band without control rows has no uplift and adds nothing
  e <- data.frame(score = c(rep(3, 40), rep(1, 80)), treat = c(rep(1, 40), rep(c(1, 0), each = 40)),
                  y = rep(c(1, 0), 60))
  ue <- scr_uplift(e, cuts = 2, n_boot = 0)
  expect_true(is.na(ue$table$uplift[1])); expect_true(is.na(ue$table$type[1]))
  expect_equal(ue$table$cum_incremental, c(0, 0)); expect_equal(ue$table$n_c, c(0, 40))
  # a missing outcome counts in the score distribution only
  na <- d; na$y[1:200] <- NA
  un <- scr_uplift(na, cuts = c(470, 530), n_boot = 0)
  uk <- scr_uplift(na[-(1:200), ], cuts = c(470, 530), n_boot = 0)
  expect_equal(un$table, uk$table); expect_equal(un$summary$qini, uk$summary$qini)
  expect_false(isTRUE(all.equal(un$randomization$psi, uk$randomization$psi)))
  # a logical treatment column
  l <- d; l$treat <- l$treat == 1
  expect_equal(scr_uplift(l, n_bands = 3, n_boot = 0)$table, scr_uplift(d, n_bands = 3, n_boot = 0)$table)
  expect_error(scr_uplift(d[d$treat == 1, ]), "both a treated and a control")
  expect_error(scr_uplift(transform(d, treat = ifelse(seq_along(treat) == 1, NA, treat))), "without missing values")
  expect_error(scr_uplift(transform(d, treat = treat + 1)), "0/1")
  expect_error(scr_uplift(d, treat = "nope"), "not in `x`")
  expect_error(scr_uplift(d, n_boot = -1), "n_boot")
  expect_error(scr_uplift(d, foo = 1), "unused argument")
  expect_error(scr_uplift(d[0, ]), "both a treated and a control")
  # an arm without a known outcome: rates and areas are missing, not an error
  nk <- d; nk$y[nk$treat == 0] <- NA
  uu <- scr_uplift(nk, n_bands = 3, n_boot = 10)
  expect_true(is.na(uu$summary$uplift)); expect_true(is.na(uu$summary$qini)); expect_null(uu$curve)
})

test_that("pooled score cells: exact bands, close areas", {
  d <- upl_df(4000, 15)
  d$score <- d$score + stats::runif(nrow(d))
  ex <- scr_uplift(d, n_bands = 4, n_boot = 0)
  po <- scr_uplift(d, n_bands = 4, n_boot = 0, max_cells = 150)
  expect_true(po$quantized); expect_false(ex$quantized)
  t <- po$table
  for (k in 1:4) {
    i <- d$score >= t$score_lo[k] & d$score < t$score_hi[k]
    expect_equal(t$n_t[k], sum(i & d$treat == 1)); expect_equal(t$events_c[k], sum(d$y[i & d$treat == 0]))
  }
  expect_equal(po$summary$qini, ex$summary$qini, tolerance = 0.02)
  expect_equal(po$summary$auuc, ex$summary$auuc, tolerance = 0.02)
  expect_lte(nrow(po$curve), 150L)
})

test_that("print and export", {
  d <- upl_df(2000, 2)
  up <- scr_uplift(d, n_bands = 4, n_boot = 30, seed = 1)
  out <- utils::capture.output(print(up))
  expect_match(out[1], "^<scr_uplift> outcome \"y\" \\| treatment \"treat\"")
  expect_true(any(grepl("Qini .* \\| AUUC .* \\(95% bootstrap intervals, 30 resamples\\)", out)))
  expect_true(any(grepl("randomization: PSI", out)))
  skip_if_not_installed("openxlsx")
  old <- scr_verbose(FALSE); on.exit(scr_verbose(old), add = TRUE)
  dir <- file.path(tempdir(), "scr-uplift-export")
  unlink(dir, recursive = TRUE)
  ex <- scr_export(up, dir, stamp = FALSE)
  expect_identical(basename(ex$files$xlsx), "uplift_y.xlsx")
  expect_identical(openxlsx::getSheetNames(ex$files$xlsx), c("Summary", "Bands", "Curve", "Cuts", "Settings"))
  expect_equal(nrow(openxlsx::read.xlsx(ex$files$xlsx, sheet = "Bands")), 4L)
})
