test_that("scr_compare and scr_core summarise a list of results", {
  r1 <- res_demo()
  r2 <- scr_select(scr_demo, "churn", config = cfg_test(objective = "propensity"),
                   drop = c("id", "default"), date_col = "ref_date")
  cmp <- scr_compare(list(default = r1, churn = r2))
  expect_equal(nrow(cmp), 2L)
  expect_true(all(c("target", "approved", "auc", "auc_lo", "auc_hi", "relaxation") %in% names(cmp)))
  core <- scr_core(list(default = r1, churn = r2), min_targets = 2)
  expect_true(all(core$n_targets == 2))
  expect_true(all(core$feature %in% intersect(scr_selected(r1), scr_selected(r2))))
  expect_error(scr_compare(list(a = 1)), "no valid result")
  expect_error(scr_core(list(a = 1)), "no valid result")
})
