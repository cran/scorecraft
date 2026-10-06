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
