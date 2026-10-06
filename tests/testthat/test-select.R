test_that("scr_select returns a coherent object", {
  res <- res_demo()
  expect_s3_class(res, "scr_result")
  expect_equal(res$target, "default")
  expect_output(print(res), "scr_result")
  expect_s3_class(summary(res), "scr_summary")
  expect_output(print(summary(res)), "Funnel")
  expect_s3_class(as.data.frame(res), "data.frame")
  pdf(NULL); on.exit(dev.off())
  expect_invisible(plot(res))
})

test_that("the funnel and the consensus never contradict each other", {
  res <- res_demo()
  sel <- scr_selected(res)
  f <- scr_funnel(res, cols = "all")
  expect_equal(nrow(f[approved == TRUE]), length(sel))
  expect_true(all(f[approved == TRUE, exit_stage] == "07.approved"))
  expect_setequal(f[approved == TRUE, feature], sel)
  expect_true(all(sel %in% res$prune$keep))
  expect_false(any(sel %in% res$screen$summary[!(selected %in% TRUE), feature]))
  expect_false(any(sel %in% res$holdout[!(holdout_ok %in% TRUE), feature]))
  # every input column is in the funnel, dropped ones included
  expect_true(all(setdiff(names(scr_demo), "default") %in% f$feature))
  expect_equal(f[feature == "id", exit_stage], "00.config")
  expect_equal(f[feature == "ref_date", exit_stage], "00.config")
})

test_that("pure noise is barred and genuine signal recovered", {
  sel <- scr_selected(res_demo())
  expect_lte(sum(noise_names() %in% sel), 1L)
  expect_gte(sum(signal_names() %in% sel), 4L)
  expect_lte(length(sel), res_demo()$config$target_max)
})

test_that("derived flags stay out of the deliverable and the SQL", {
  res <- res_demo()
  expect_false(any(grepl("__sp$", scr_selected(res))))
  expect_false(any(grepl("__sp", scr_sql(res))))
  f <- scr_funnel(res, cols = "all")
  out <- f[exit_stage == "05b.derived_excluded"]
  if (nrow(out)) expect_true(all(!is.na(out$iv_holdout)))
  res2 <- scr_select(scr_demo, "default", config = cfg_test(allow_derived_final = TRUE),
                     drop = c("id", "churn"), date_col = "ref_date")
  expect_equal(nrow(scr_funnel(res2, cols = "all")[exit_stage == "05b.derived_excluded"]), 0L)
})

test_that("objective does not touch the selection", {
  r_risk <- res_demo()
  r_prop <- scr_select(scr_demo, "default", config = cfg_test(objective = "propensity"),
                       drop = c("id", "churn"), date_col = "ref_date")
  expect_identical(scr_selected(r_risk), scr_selected(r_prop))
  expect_equal(r_risk$funnel$total_iv, r_prop$funnel$total_iv)
  expect_equal(r_prop$meta$score_direction, "higher_is_riskier")
})

test_that("scr_select does not modify the user's data when copy = TRUE", {
  d <- head(scr_demo, 1500)
  before <- vapply(d, function(x) class(x)[1], character(1))
  invisible(scr_select(d, "default", config = cfg_test(), drop = c("id", "churn"), date_col = "ref_date"))
  expect_identical(vapply(d, function(x) class(x)[1], character(1)), before)
})

test_that("accessors validate their input and the leakage audit prints", {
  res <- res_demo()
  expect_error(scr_selected(list()), "scr_select")
  expect_error(scr_funnel(res, cols = "nope"), "not in the funnel")
  expect_s3_class(scr_gains(res), "data.table")
  expect_true(all(scr_gains(res)$approved))
  expect_gt(nrow(scr_gains(res, only_selected = FALSE)), nrow(scr_gains(res)))
  lk <- scr_leakage(res)
  expect_s3_class(lk, "scr_leakage")
  expect_output(print(lk), "scr_leakage")
})

test_that("scr_select relaxes the consensus in named steps when the pool is short", {
  res <- scr_select(scr_demo, "default", config = cfg_test(target_min = 20, target_max = 25, min_votes = 2),
                    drop = c("id", "churn"), date_col = "ref_date")
  expect_true(res$consensus$meta$scarce || res$consensus$meta$relaxation != "none")
  expect_output(print(res), "Warnings")
})
