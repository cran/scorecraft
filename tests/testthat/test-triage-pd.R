# Regression tests of the PD / default / IRB-binning triage: each block pins
# one fixed behaviour against a brute-force or textbook reference.

test_that("the vectorised default engine equals the row state machine unit by unit", {
  set.seed(11)
  for (rep in 1:60) {
    nid <- sample(1:12, 1)
    ids <- rep(c("a", "B", "b", "A-1", "A1", "z z", "Z", "a_2", "10", "9", "é", "B.")[seq_len(nid)],
               times = sample(1:30, nid, TRUE))
    o <- order(ids, method = "radix")                # the C-locale order setorder() uses
    ids <- ids[o]; n <- length(ids)
    trig <- runif(n) < runif(1, 0, 0.5); utp <- trig & runif(n) < 0.3; restr <- runif(n) < 0.05
    P <- sample(0:4, 1); PR <- sample(c(P, P + 5L), 1)
    new <- .default_flags(ids, trig, utp, restr, P, PR)
    # reference: .default_run on every unit, placed back by row position
    ref <- list(integer(n), integer(n), character(n), integer(n), integer(n))
    for (u in unique(ids)) {
      ix <- which(ids == u)
      r <- .default_run(trig[ix], utp[ix], restr[ix], P, PR)
      ref[[1]][ix] <- r$default; ref[[2]][ix] <- r$ev; ref[[3]][ix] <- r$trigger
      ref[[4]][ix] <- r$months; ref[[5]][ix] <- r$cured
    }
    expect_equal(new, ref)
  }
  expect_equal(lengths(.default_flags(character(), logical(), logical(), logical(), 3L, 12L)), rep(0L, 5))
})

test_that("scr_default flags a unit the same whatever the other units' names (no split/collation misalignment)", {
  cfg <- cfg_test()
  pn <- data.table::as.data.table(scr_demo_panel)
  keep <- unique(pn$id)[1:40]
  pn <- pn[id %in% keep]
  # mixed case and punctuation: C-locale and locale collation disagree on these
  pn[, id2 := paste0(ifelse(match(id, keep) %% 2L == 0L, "b", "B"), "-", id)]
  d <- scr_default(pn, "id2", "ref_date", dpd = "dpd", restructured = "restructured", config = cfg)
  for (u in unique(pn$id2)[c(1, 2, 7, 20)]) {
    d1 <- scr_default(pn[id2 == u], "id2", "ref_date", dpd = "dpd", restructured = "restructured", config = cfg)
    expect_equal(d$flags[id == u, list(default, months_in_default, cured, trigger)],
                 d1$flags[, list(default, months_in_default, cured, trigger)])
  }
})

test_that("scr_default_rate equals the per-cohort brute force", {
  cfg <- cfg_test()
  d <- scr_default(scr_demo_panel, "id", "ref_date", dpd = "dpd", config = cfg)
  f <- data.table::copy(d$flags)
  f[, seg := substr(id, nchar(id), nchar(id))]
  dr <- scr_default_rate(f, segment = "seg", by = "quarter", config = cfg)
  starts <- sort(unique(dr$table$cohort))
  ref <- data.table::rbindlist(lapply(starts, function(t0) {
    t1 <- .add_months(t0, 12L)
    pop <- f[date == t0 & default == 0L]
    dd <- unique(f[default == 1L & date > t0 & date <= t1, id])
    pop[, list(n = .N, defaults = sum(id %in% dd)), by = "seg"][, cohort := t0]
  }))
  data.table::setorderv(ref, c("cohort", "seg"))
  expect_equal(dr$table$n, ref$n)
  expect_equal(dr$table$defaults, ref$defaults)
  expect_equal(.dr_outcome(data.table::data.table(id = c("a", "a"), date = as.Date(c("2020-03-01", "2021-06-01"))),
                           c("a", "a", "b"), as.Date(c("2020-01-01", "2020-03-01", "2020-01-01")),
                           as.Date(c("2021-01-01", "2021-03-01", "2021-01-01"))), c(1L, 0L, 0L))
})

test_that("Hosmer-Lemeshow on fixed grade PDs has K degrees of freedom", {
  hl <- .pd_hl(n = c(1000, 800, 500), d = c(10, 24, 40), pd = c(0.01, 0.03, 0.06))
  expect_equal(hl$chi2, 100 / 28.2)
  expect_equal(hl$df, 3L)
  expect_equal(hl$p, stats::pchisq(100 / 28.2, 3, lower.tail = FALSE))
  # grades without obligors are left out of the count
  expect_equal(.pd_hl(c(100, 0, 100), c(1, 0, 2), c(0.01, 0.5, 0.02))$df, 2L)
})

test_that("the implied AUC and the DeLong standard error match their brute-force definitions", {
  set.seed(3)
  p <- round(stats::runif(60, 0.01, 0.4), 2)          # ties on purpose
  bf <- 0
  for (i in seq_along(p)) for (j in seq_along(p)) bf <- bf + p[i] * (1 - p[j]) * ((p[i] > p[j]) + 0.5 * (p[i] == p[j]))
  expect_equal(.pd_auc_implied(p), bf / (sum(p) * sum(1 - p)))
  s <- sample(1:6, 80, TRUE); y <- stats::rbinom(80, 1, s / 8)
  xe <- s[y == 1]; xn <- s[y == 0]
  psi <- outer(xe, xn, function(a, b) (a > b) + 0.5 * (a == b))
  v10 <- rowMeans(psi); v01 <- colMeans(psi)
  expect_equal(.pd_auc_se(s, y), sqrt(stats::var(v10) / length(xe) + stats::var(v01) / length(xn)))
  expect_equal(.pd_auc_se(-s, y, higher_is_event = FALSE), .pd_auc_se(s, y))
  expect_true(is.na(.pd_auc_se(c(1, 2, 3), c(1, 0, 0))))
})

test_that("the intercept solver extends its bracket when b * ln(odds) is large", {
  l <- seq(-5, 5, length.out = 101)
  a <- .pd_solve_a(l, 50, 0.9)
  expect_equal(mean(stats::plogis(a + 50 * l)), 0.9, tolerance = 1e-9)
  expect_error(scr_calibrate(c(NA, Inf), target = 0.05, sample_rate = 0.1), "no calibration row")
})

test_that("the migration counts equal a cell-by-cell tabulation", {
  set.seed(5)
  g0 <- sample(1:4, 300, TRUE)
  g1 <- as.character(pmin(4, pmax(1, g0 + sample(-1:1, 300, TRUE))))
  g1[sample(300, 10)] <- NA; g1[sample(300, 5)] <- "default"
  m <- scr_migration(g0, g1, K = 4)
  col <- ifelse(is.na(g1), "closed", g1)
  ref <- table(factor(g0, 1:4), factor(col, c(1:4, "default", "closed")))
  expect_equal(unname(m$matrix), unname(matrix(as.integer(ref), 4)))
})

test_that("scr_moc does not edit the ledger of the object it was given", {
  gr <- scr_moc(pd_grades(), "C", method = "ci_binomial")
  before <- data.table::copy(gr$moc)
  gr2 <- scr_moc(gr, "C", method = "ci_binomial", level = 0.99)
  expect_equal(gr$moc, before)
  expect_true(all(gr$moc$active))
  expect_equal(sum(!gr2$moc$active), nrow(gr$table))
})

test_that("the bootstrap margin draws the resampled default rate as Binomial(n, DR) / n", {
  gr <- pd_grades(); K <- nrow(gr$table)
  g <- scr_moc(gr, "C", method = "bootstrap", n_boot = 50, seed = 9, level = 0.9)
  n_k <- tabulate(gr$rows$grade, K); d_k <- tabulate(gr$rows$grade[gr$rows$y == 1L], K)
  set.seed(9)
  B <- vapply(seq_len(K), function(k) stats::rbinom(50, n_k[k], d_k[k] / n_k[k]) / n_k[k], numeric(50))
  ref <- pmax(0, apply(B, 2L, stats::quantile, probs = 0.9, names = FALSE) - gr$table$dr)
  expect_equal(g$moc$value, ref)
  expect_true(all(g$moc$value >= 0))
})

test_that("the PIT bridge of scr_pd accepts a grade whose PD is zero", {
  gr <- pd_grades()
  gr$table <- data.table::copy(gr$table); gr$table$pd_be[1] <- 0
  pd <- scr_pd(gr, philosophy = "pit", rho = 0.12, z = -1)
  expect_equal(pd$table$pd_pit[1], 0)
  expect_equal(pd$table$pd_final[1], pd$floor)
  expect_equal(pd$table$pd_pit[-1], .vasicek_pit(pd$table$pd_moc[-1], -1, 0.12))
})

test_that("the validation uses the cohort differences in the multi-period test and matches the brute-force cohorts", {
  pd <- pd_model(); p <- pd_panel()
  p[, grade := predict(pd, score = score)]
  v <- scr_pd_validate(pd, p, grade = "grade", score = "score", by = "quarter", tests = c("jeffreys", "hl", "multi_period", "psi", "migration"))
  dd <- v$portfolio$dr - v$portfolio$pd
  expect_equal(v$portfolio_tests$multi_period_z, mean(dd) / (stats::sd(dd) / sqrt(length(dd))))
  expect_equal(v$portfolio_tests$hl_df, sum(v$calibration$n > 0))
  # cohort rows: population, outcome and grade at the end of the window, brute force
  dt <- data.table::data.table(id = as.character(p$id), date = as.Date(p$date), default = as.integer(p$default),
                               score = p$score, grade = as.integer(p$grade))
  rows <- .pd_cohorts(dt, 12L, "quarter")
  t0 <- sort(unique(rows$cohort))[2]; t1 <- .add_months(t0, 12L)
  pop <- dt[date == t0 & default == 0L & !is.na(grade)]
  dids <- unique(dt[default == 1L & date > t0 & date <= t1, id])
  r0 <- rows[cohort == t0]
  expect_equal(sort(r0$id), sort(pop$id))
  expect_equal(r0$y[order(r0$id)], as.integer(sort(pop$id) %in% dids))
  endg <- dt[date == t1][match(r0$id, id), grade]
  expect_equal(r0$grade_t1, ifelse(r0$y == 1L, "default", as.character(endg)))
})

test_that("the continuous binner flags a hold-out whose bin shares shift (PSI_ACTION) and its sums are exact", {
  set.seed(8)
  d <- data.table::data.table(x = c(stats::runif(1000), stats::runif(500, 0.8, 1)),
                              g = sample(c("a", "b", "c"), 1500, TRUE))
  d[, y := pmin(1, pmax(0, 0.1 + 0.6 * x + stats::rnorm(.N, 0, 0.05)))]
  cb <- scr_bin_continuous(d, "y", c("x", "g"), train_idx = 1:1000, holdout_idx = 1001:1500)
  sx <- cb$summary[feature == "x"]
  expect_equal(sx$psi_flag, "shift")
  expect_match(sx$holdout_reason, "PSI_ACTION")
  expect_false(sx$holdout_ok)
  # stable driver: no PSI_ACTION
  expect_false(grepl("PSI_ACTION", cb$summary[feature == "g"]$holdout_reason))
  # one-pass sums equal the per-bin sums
  idx <- c(3L, 1L, 3L, NA, 1L); yy <- c(0.2, 0.4, 0.1, 9, 0.3)
  expect_equal(.cbin_sum_by(yy, idx, 4L), c(0.7, 0, 0.3, 0))
  e <- cb$fit$results$g
  expect_equal(e$mean, unname(vapply(strsplit(e$bin, "%;%", fixed = TRUE), function(l) mean(d$y[1:1000][d$g[1:1000] %in% l]), numeric(1))))
})
