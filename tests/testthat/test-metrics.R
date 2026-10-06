test_that("scr_metrics recovers perfect, random and reversed separation", {
  y <- rep(0:1, each = 100)
  perfect <- scr_metrics(as.numeric(y), y, ci = FALSE)
  expect_equal(perfect$auc, 1)
  expect_equal(perfect$ks, 1)
  expect_equal(perfect$gini, 1)
  # a higher_is_safer score with perfect separation, declared as such
  expect_equal(scr_metrics(-as.numeric(y), y, higher_is_event = FALSE, ci = FALSE)$auc, 1)
  set.seed(1)
  rnd <- scr_metrics(stats::rnorm(4000), rep(0:1, 2000), ci = FALSE)
  expect_lt(abs(rnd$auc - 0.5), 0.05)
})

test_that("scr_metrics handles ties and large counts in double", {
  y <- rep(0:1, each = 60000)
  s <- rep(c(0, 1), each = 60000)   # perfect separation with massive ties
  m <- scr_metrics(s, y, ci = FALSE)
  expect_equal(m$auc, 1)
  expect_false(is.na(m$auc))
})

test_that("scr_metrics bootstrap CI contains the point estimate and is reproducible", {
  set.seed(2)
  y <- stats::rbinom(1000, 1, 0.3)
  s <- y + stats::rnorm(1000)
  m1 <- scr_metrics(s, y, n_boot = 30, seed = 7)
  m2 <- scr_metrics(s, y, n_boot = 30, seed = 7, nthread = 2L)
  expect_s3_class(m1, "scr_metrics")
  expect_equal(m1$n_boot, 30L)
  expect_true(m1$auc_lo <= m1$auc && m1$auc <= m1$auc_hi)
  expect_true(m1$ks_lo <= m1$ks && m1$ks <= m1$ks_hi)
  expect_equal(m1$gini_lo, 2 * m1$auc_lo - 1)
  expect_equal(unclass(m1), unclass(m2))   # same seeds drawn in the main process
  expect_output(print(m1), "bootstrap CI")
  expect_s3_class(as.data.frame(m1), "data.frame")
})

test_that("scr_metrics degrades to NA on a single class or empty input", {
  expect_true(is.na(scr_metrics(1:5, rep(1, 5))$auc))
  expect_true(is.na(scr_metrics(numeric(), integer())$auc))
  expect_error(scr_metrics(1:3, 1:2), "same length")
})

test_that("scr_iv is zero without groups and positive with signal, finite with a single-class group", {
  set.seed(3)
  y <- stats::rbinom(1000, 1, 0.3)
  expect_equal(scr_iv(rep("A", 1000), y), 0)
  g <- ifelse(stats::runif(1000) < 0.5 + 0.3 * y, "A", "B")
  expect_gt(scr_iv(g, y), 0.05)
  g2 <- c(rep("only_events", 20), rep("mixed", 980)); y2 <- c(rep(1L, 20), y[21:1000])
  expect_true(is.finite(scr_iv(g2, y2)))
  expect_true(is.infinite(scr_iv(g2, y2, laplace = 0)) || scr_iv(g2, y2, laplace = 0) > scr_iv(g2, y2))
})

test_that("scr_psi reports the fixed and the n-adjusted threshold side by side", {
  set.seed(4)
  base <- stats::rnorm(1000)
  same <- scr_psi(base, base)
  expect_equal(same$psi, 0)
  expect_equal(same$flag_fixed, "stable")
  expect_equal(same$flag_adjusted, "stable")
  # Yurdakul & Naranjo: n = m = 1000, B = 10 -> critical ~ 0.034 at 5%
  expect_equal(same$critical, (2 / 1000) * stats::qchisq(0.95, 9), tolerance = 1e-12)
  expect_equal(round(same$critical, 3), 0.034)
  shifted <- scr_psi(base, stats::rnorm(1000, mean = 0.6))
  expect_gt(shifted$psi, 0.05)
  expect_true(shifted$flag_adjusted == "shift")
  expect_equal(nrow(shifted$table), 10L)
  expect_output(print(shifted), "n-adjusted threshold")
})

test_that("scr_psi works on categorical vectors with frozen levels and empty bins", {
  a <- c(rep("x", 50), rep("y", 50))
  b <- c(rep("x", 90), rep("y", 10))
  p <- scr_psi(a, b, levels = c("x", "y", "z"))
  expect_equal(p$n_bins, 3L)
  expect_true(is.finite(p$psi))
  expect_true(is.na(scr_psi(character(), b)$psi))
})
