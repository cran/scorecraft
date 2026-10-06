test_that("parallel binning by column is identical to serial", {
  sp <- scr_split(scr_demo, "default", date_col = "ref_date", drop = c("id", "churn"))
  tr <- scr_triage(sp, cfg_test())
  b1 <- scr_bin(tr, cfg_test(nthread = 1L))
  b2 <- scr_bin(tr, cfg_test(nthread = 2L))
  expect_identical(names(b1$fit$results), names(b2$fit$results))
  for (f in names(b1$fit$results)) {
    expect_identical(b1$fit$results[[f]]$bin, b2$fit$results[[f]]$bin)
    expect_equal(b1$fit$results[[f]]$woe, b2$fit$results[[f]]$woe)
  }
  expect_equal(b1$holdout, b2$holdout)
  expect_identical(b1$pool, b2$pool)
  expect_output(print(b2), "binned")
})

test_that("screening, hold-out and pruning each leave a verdict with a reason", {
  bn <- scr_bin(scr_triage(scr_split(scr_demo, "default", date_col = "ref_date",
                                     drop = c("id", "churn")), cfg_test()), cfg_test())
  s <- bn$screen$summary
  expect_true(all(c("feature", "selected", "reason") %in% names(s)))
  expect_true(any(grepl("IV_BELOW_MIN", s$reason)))
  h <- bn$holdout
  expect_true(all(c("iv_holdout", "psi", "psi_critical", "psi_flag_adjusted", "holdout_ok", "holdout_reason") %in% names(h)))
  expect_true(all(h$psi_flag_adjusted[!is.na(h$psi)] %in% c("stable", "shift")))
  # vl_score_02 and vl_redundant are near copies: at most one survives
  expect_lte(sum(c("vl_score_02", "vl_redundant") %in% bn$prune$keep), 1L)
  expect_true(all(bn$pool %in% bn$prune$keep))
  expect_false(any(bn$derived_excluded %in% bn$pool))
})

test_that("a numeric-only algorithm falls back to jedi for categoricals with a message", {
  sp <- scr_split(head(scr_demo, 1500), "default", drop = c("id", "churn", "ref_date"))
  tr <- scr_triage(sp, cfg_test())
  expect_message(scr_bin(tr, cfg_test(algorithm = "ir", verbose = TRUE)), "not categorical")
  expect_equal(.algorithm_for(cfg_test(algorithm = "ir"), "categorical"), "jedi")
  expect_equal(.algorithm_for(cfg_test(algorithm = "ir"), "numerical"), "ir")
  expect_equal(.algorithm_for(cfg_test(algorithm = "ivb"), "numerical"), "jedi")
})

test_that(".merge_obwoe keeps the requested feature order", {
  sp <- scr_split(head(scr_demo, 1500), "default", drop = c("id", "churn", "ref_date"))
  tr <- scr_triage(sp, cfg_test())
  feats <- setdiff(names(tr$clean), "default")
  fit <- fit_binning(tr$clean[sp$train_idx], "default", feats, cfg_test(nthread = 2L))
  expect_identical(names(fit$results), feats[feats %in% names(fit$results)])
  expect_identical(fit$summary$feature, names(fit$results))
  expect_equal(fit$n_features, length(fit$results))
})
