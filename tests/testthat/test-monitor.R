test_that("scr_monitor reports PSI, CSI with points shift and vintage per period", {
  sc <- sc_demo()
  mo <- scr_monitor(sc, scr_demo, date_col = "ref_date", target = "default", n_boot = 10)
  expect_s3_class(mo, "scr_monitor")
  expect_equal(nrow(mo$psi), 6L)
  expect_true(all(mo$psi$psi >= 0))
  expect_true(all(mo$psi$flag_adjusted %in% c("stable", "shift")))
  expect_equal(nrow(mo$csi), 6L * length(sc$features))
  # vl_late degrades only in the last period: its CSI there dominates
  vt <- mo$csi[variable == "vl_late"]
  expect_equal(vt$period[which.max(vt$csi)], "2026-06-01")
  expect_true(all(is.finite(mo$csi$points_shift)))
  expect_equal(nrow(mo$vintage), 6L)
  expect_true(all(mo$vintage$auc > 0.6))
  expect_true(all(c("psi_score_fixed_action", "threshold_source") %in% mo$plan$item))
  expect_output(print(mo), "performance by vintage")
})

test_that("scr_monitor works without a date or target and validates names", {
  sc <- sc_demo()
  mo <- scr_monitor(sc, head(scr_demo, 500))
  expect_equal(mo$periods, "all")
  expect_null(mo$vintage)
  expect_error(scr_monitor(sc, scr_demo, date_col = "nope"), "date_col")
  expect_error(scr_monitor(sc, scr_demo, target = "nope"), "target")
  expect_error(scr_monitor(res_demo(), scr_demo), "scr_scorecard")
})

test_that("the CSI of the monitor is the PSI of scr_psi() on the bin counts", {
  # reference: smoothing over every bin, k - 1 degrees of freedom
  csi_ref <- function(pt, cmp, alpha, thresholds) {
    n_new <- sum(cmp); n_tr <- sum(pt$count_train); k <- nrow(pt)
    sm <- if (any(pt$count_train == 0L) || any(cmp == 0L)) 0.5 else 0
    pb <- (pt$count_train + sm) / (n_tr + sm * k); pc <- (cmp + sm) / (n_new + sm * k)
    csi <- sum((pb - pc) * log(pb / pc))
    crit <- (1 / n_tr + 1 / max(1L, n_new)) * stats::qchisq(1 - alpha, df = max(1L, k - 1L))
    list(csi = csi, critical = crit,
         flag_fixed = if (csi < thresholds[1]) "stable" else if (csi < thresholds[2]) "moderate" else "shift",
         flag_adjusted = if (csi < crit) "stable" else "shift")
  }
  pt <- data.table::data.table(bin = c("a", "b", "c", "d"), count_train = c(120L, 300L, 80L, 500L),
                               points = c(10, 5, -3, 0))
  # every bin populated in both samples, and a bin empty in the new sample only: unchanged, bit for bit
  for (cmp in list(c(30L, 70L, 25L, 90L), c(30L, 0L, 25L, 90L))) {
    r <- .csi_row(pt, cmp, 0.05, c(0.10, 0.25)); o <- csi_ref(pt, cmp, 0.05, c(0.10, 0.25))
    expect_identical(r[c("csi", "critical", "flag_fixed", "flag_adjusted")], o)
    expect_identical(r$n, sum(cmp))
    expect_equal(r$points_shift, sum((cmp / sum(cmp) - pt$count_train / 1000) * pt$points))
  }
  # a bin empty in both samples: out of the smoothing and of the degrees of freedom
  pt0 <- data.table::copy(pt)[, count_train := c(120L, 0L, 80L, 500L)]
  for (cmp in list(c(30L, 0L, 25L, 90L), c(30L, 0L, 0L, 90L))) {
    r <- .csi_row(pt0, cmp, 0.05, c(0.10, 0.25))
    p <- scr_psi(rep(pt0$bin, pt0$count_train), rep(pt0$bin, cmp), levels = pt0$bin, alpha = 0.05)
    expect_identical(r$csi, p$psi)
    expect_identical(r$critical, p$critical)
    expect_identical(c(r$flag_fixed, r$flag_adjusted), c(p$flag_fixed, p$flag_adjusted))
    expect_equal(r$critical, (1 / 700 + 1 / sum(cmp)) * stats::qchisq(0.95, df = 2))
  }
  expect_false(identical(.csi_row(pt0, c(30L, 0L, 25L, 90L), 0.05)$csi, csi_ref(pt0, c(30L, 0L, 25L, 90L), 0.05, c(0.10, 0.25))$csi))
  # no row of the period in the bins: every figure is NA
  r0 <- .csi_row(pt, integer(4), 0.05)
  expect_identical(r0$n, 0L)
  expect_true(is.na(r0$csi) && is.na(r0$critical) && is.na(r0$flag_fixed) && is.na(r0$flag_adjusted) && is.na(r0$points_shift))
})
