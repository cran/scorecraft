# Shared IRB core: parameter tables, the default engine, default rates and the
# continuous-target binner whose result has the shape of an obwoe object.

test_that("parameter presets are complete tables and edits are detected", {
  for (fw in c("bcb", "basel3_final", "crr3")) {
    p <- scr_irb_params(fw)
    expect_s3_class(p, "scr_irb_params")
    expect_setequal(names(p$pd_floor), c("asset_class", "floor"))
    expect_equal(p$pd_floor$floor[p$pd_floor$asset_class == "qrre_revolver"], 0.001)
    expect_true(is.na(p$pd_floor$floor[p$pd_floor$asset_class == "sovereign"]))
    expect_equal(p$ccf_sa$ccf, c(0.10, 0.40, 0.50, 1.00))
    expect_equal(p$correlation$retail_mortgage, 0.15)
    expect_equal(p$output_floor, 0.725)
    expect_output(print(p), fw)
    expect_false(.check_params(p, "t")$modified)
  }
  expect_equal(scr_irb_params("bcb")$lgd_firb$lgd[1], 0.75)
  expect_equal(scr_irb_params("basel3_final")$lgd_firb$lgd[1], 0.40)
  expect_equal(scr_irb_params("bcb")$correlation$sme$lo, 15)
  p <- scr_irb_params("bcb"); p$pd_floor$floor[7] <- 0.002
  expect_true(.check_params(p, "t")$modified)
  expect_error(.check_params(list(), "t"), "scr_irb_params")
  expect_equal(.pd_floor_of(scr_irb_params(), c("retail_other", "qrre_revolver")), c(5e-4, 1e-3))
  expect_error(.pd_floor_of(scr_irb_params(), "boats"), "asset_class")
})

test_that("the IRB configuration keys are registered, validated and documented", {
  cfg <- scr_config(verbose = FALSE)
  expect_equal(cfg$framework, "bcb")
  expect_equal(cfg$ccf_horizon_months, 12L)
  keys <- scr_config_keys()
  expect_true(all(c("default_days", "pd_calibration", "lgd_t_max", "ccf_u_star", "framework", "ecl_sicr_ratio") %in% keys$key))
  expect_setequal(setdiff(names(cfg), keys$key), character())
  expect_true(all(keys$stage[keys$key %in% c("pd_n_grades")] == 9))
  expect_error(scr_config(pd_calibration = "magic"), "pd_calibration")
  expect_error(scr_config(framework = "mars"), "framework")
  expect_error(scr_config(ccf_u_star = 2), "ccf_u_star")
  expect_error(scr_config(pd_lights = c(0.05, 0.01)), "pd_lights")
  expect_equal(nrow(scr_config_keys(stage = 12)), sum(keys$stage == 12))
})

test_that("the default engine applies the trigger, the materiality test and the probation", {
  cfg <- scr_config(verbose = FALSE)
  d <- scr_default(scr_demo_panel, "id", "ref_date", dpd = "dpd", arrears = "arrears", exposure = "exposure",
                   restructured = "restructured", config = cfg)
  expect_s3_class(d, "scr_default")
  expect_equal(nrow(d$flags), nrow(scr_demo_panel))
  expect_true(d$summary$n_events > 50)
  expect_output(print(d), "events")
  # every event starts on a trigger row and ends after the probation without triggers
  f <- merge(d$flags, data.table::as.data.table(scr_demo_panel)[, list(id, date = ref_date, dpd, arrears, exposure)],
             by = c("id", "date"))
  starts <- f[!is.na(event_id) & months_in_default == 1L]
  expect_true(all(starts$dpd >= 90L))
  expect_true(all(starts$arrears >= 100 & starts$arrears / starts$exposure > 0.01))
  # immaterial arrears never default (the obligors capped at 80 in the generator)
  imm <- f[arrears > 0 & arrears < 100 & dpd >= 90L]
  expect_true(nrow(imm) > 0)
  expect_true(all(imm[months_in_default == 1L, .N] == 0))
  # rows in default form contiguous runs; the run has at least `probation` quiet months at its end
  ev <- d$events
  expect_true(all(ev$months >= 1L))
  expect_equal(sum(f$default), sum(ev$months))
  # without materiality the same panel defaults more
  d0 <- scr_default(scr_demo_panel, "id", "ref_date", dpd = "dpd", config = cfg)
  expect_gt(d0$summary$n_events, d$summary$n_events)
  # utp alone works
  pu <- data.table::as.data.table(scr_demo_panel)[, list(id, ref_date, utp = dpd >= 120L)]
  du <- scr_default(pu, "id", "ref_date", utp = "utp", config = cfg)
  expect_true(all(du$events$trigger == "utp"))
  expect_error(scr_default(scr_demo_panel, "id", "ref_date", config = cfg), "dpd")
  expect_error(scr_default(scr_demo_panel, "id", "ref_date", dpd = "nope", config = cfg), "not found")
  # identical under every backend
  d_ser <- withr::with_options(list(scorecraft.parallel = "serial"), scr_default(scr_demo_panel, "id", "ref_date", dpd = "dpd", config = scr_config(verbose = FALSE, nthread = 2L)))
  d_psk <- withr::with_options(list(scorecraft.parallel = "psock"),  scr_default(scr_demo_panel, "id", "ref_date", dpd = "dpd", config = scr_config(verbose = FALSE, nthread = 2L)))
  expect_equal(d_ser$flags, d0$flags); expect_equal(d_psk$flags, d0$flags)
})

test_that("a longer probation keeps units in default longer and restructuring lengthens it", {
  cfg3 <- scr_config(verbose = FALSE, default_probation = 3L)
  cfg6 <- scr_config(verbose = FALSE, default_probation = 6L)
  d3 <- scr_default(scr_demo_panel, "id", "ref_date", dpd = "dpd", config = cfg3)
  d6 <- scr_default(scr_demo_panel, "id", "ref_date", dpd = "dpd", config = cfg6)
  expect_gt(sum(d6$flags$default), sum(d3$flags$default))
  expect_lte(d6$summary$n_events, d3$summary$n_events)
  # a hand-made unit: trigger for 2 months, then quiet: exits after 3 quiet months
  r <- .default_run(trig = c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE), utp = rep(FALSE, 8), restr = rep(FALSE, 8), probation = 3L, probation_restr = 12L)
  expect_equal(r$default, c(0, 1, 1, 1, 1, 1, 0, 0))
  expect_equal(r$cured, c(0, 0, 0, 0, 0, 1, 0, 0))
  expect_equal(r$trigger[2], "dpd")
  r2 <- .default_run(trig = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE), utp = rep(FALSE, 6), restr = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE), probation = 3L, probation_restr = 12L)
  expect_equal(sum(r2$default), 6L)   # restructured: probation of 12 never completes here
  r3 <- .default_run(trig = c(TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE), utp = rep(FALSE, 9), restr = rep(FALSE, 9), probation = 3L, probation_restr = 12L)
  expect_equal(max(r3$ev), 2L)        # re-default after the cure is a new event
})

test_that("the pulling effect defaults the sibling facilities of a defaulted obligor", {
  fac <- data.table::as.data.table(scr_demo_panel)[, list(id, ref_date, dpd, arrears, exposure)]
  sib <- data.table::copy(fac)[, id := paste0(id, "-b")][, dpd := 0L][, arrears := 0]
  both <- rbind(fac, sib)[, obligor := sub("-b$", "", id)]
  cfg <- scr_config(verbose = FALSE)
  dp <- scr_default(both, "id", "ref_date", dpd = "dpd", arrears = "arrears", exposure = "exposure", obligor = "obligor", config = cfg)
  expect_gt(dp$summary$n_pulled_rows, 0)
  expect_gt(sum(dp$flags[grepl("-b$", id), default]), 0)
  expect_true("pulling" %in% dp$events$trigger)
  df <- scr_default(both, "id", "ref_date", dpd = "dpd", arrears = "arrears", exposure = "exposure", obligor = "obligor",
                    config = scr_config(verbose = FALSE, default_level = "facility"))
  expect_equal(df$summary$n_pulled_rows, 0)
  expect_equal(sum(df$flags[grepl("-b$", id), default]), 0)
})

test_that("default rates by cohort, grade and exposure, with the long-run average and its benchmark", {
  cfg <- scr_config(verbose = FALSE)
  d <- scr_default(scr_demo_panel, "id", "ref_date", dpd = "dpd", config = cfg)
  dr <- scr_default_rate(d, by = "quarter", config = cfg)
  expect_s3_class(dr, "scr_dr")
  expect_equal(dr$horizon, 12L)
  # 36 months: cohorts with a full 12-month window are the first 24 months -> 8 quarters
  expect_equal(nrow(dr$table), 8L)
  expect_true(all(dr$table$dr >= 0 & dr$table$dr <= 1))
  expect_equal(dr$table$defaults / dr$table$n, dr$table$dr)
  expect_equal(dr$lra$mean, mean(dr$table$dr))
  expect_equal(dr$lra$benchmark, max(dr$lra$recent5_mean, dr$lra$mean))
  expect_false(dr$lra$flag_below_benchmark)
  expect_null(dr$lra$adjusted)
  expect_output(print(dr), "long-run average")
  # a unit counted at most once per cohort and only if not in default at the start
  t0 <- dr$table$cohort[1]
  pop <- d$flags[date == t0 & default == 0L]
  expect_equal(dr$table$n[1], nrow(pop))
  # by grade: rates decrease with the score
  fl <- merge(d$flags, data.table::as.data.table(scr_demo_panel)[, list(id, date = ref_date, score, exposure)], by = c("id", "date"))
  fl[, grade := cut(score, stats::quantile(score, 0:4 / 4), include.lowest = TRUE, labels = FALSE)]
  dr2 <- scr_default_rate(fl, id = "id", date = "date", default = "default", grade = "grade", exposure = "exposure", by = "year", config = cfg)
  expect_true(all(c("grade", "ead", "dr_weighted") %in% names(dr2$table)))
  g1 <- dr2$table[cohort == min(cohort)][order(grade)]
  expect_gt(g1$dr[1], g1$dr[4])
  expect_equal(sum(g1$n), dr2$portfolio[cohort == min(cohort), n])
  expect_error(scr_default_rate(d, grade = "x", config = cfg), "pass a table")
  expect_error(scr_default_rate(fl, id = "id", date = "date", default = "default", horizon = 48L, config = cfg), "complete")
  expect_equal(.add_months(as.Date("2024-01-31"), 1L), as.Date("2024-02-29"))
  expect_equal(.add_months(as.Date("2023-11-01"), 3L), as.Date("2024-02-01"))
})

test_that("the continuous binner returns an obwoe-shaped object that the engine reproduces in R and SQL", {
  set.seed(7); n <- 1500
  dd <- data.table::data.table(x = stats::runif(n), z = stats::rnorm(n),
                               g = sample(c("a", "b", "c", "d", "e"), n, TRUE, prob = c(.3, .3, .2, .15, .05)))
  dd[, y := pmin(1, pmax(0, 0.2 + 0.6 * x - 0.1 * z + (g == "d") * 0.2 + stats::rnorm(n, 0, 0.1)))]
  cb <- scr_bin_continuous(dd, "y", c("x", "z", "g"), train_idx = 1:1000, holdout_idx = 1001:1500, nthread = 2L)
  expect_s3_class(cb, "scr_cbins")
  expect_s3_class(cb$fit, "obwoe")
  e <- cb$fit$results$x
  expect_setequal(names(e)[1:13], c("id", "bin", "woe", "iv", "count", "count_pos", "count_neg", "cutpoints",
                                   "converged", "iterations", "feature", "type", "algorithm"))
  expect_equal(length(e$bin), length(e$cutpoints) + 1L)
  expect_true(all(diff(e$woe) > 0))              # increasing in x, after the monotone step
  expect_equal(sum(e$count), 1000L)
  expect_true(all(e$count >= 30L))
  expect_true(cb$summary$eta2[cb$summary$feature == "x"] > cb$summary$eta2[cb$summary$feature == "z"])
  expect_equal(cb$summary$direction[cb$summary$feature == "z"], "decreasing")
  expect_true(cb$summary$holdout_ok[cb$summary$feature == "x"])
  expect_output(print(cb), "eta2")
  # categorical: levels ordered by mean, 'd' alone with the highest mean, labels with the engine separator
  eg <- cb$fit$results$g
  expect_equal(eg$type, "categorical")
  expect_equal(eg$bin[which.max(eg$woe)], "d")
  # the engine reproduces the bin statistic
  app <- OptimalBinningWoE::obwoe_apply(as.data.frame(dd[1001:1500, list(x, z, g)]), cb$fit, keep_original = FALSE)
  mine <- .cbins_apply_value(cb$fit, dd[1001:1500])
  expect_equal(app$x_woe, mine$x_woe); expect_equal(app$z_woe, mine$z_woe); expect_equal(app$g_woe, mine$g_woe)
  # bin means are the means of y on train
  idx <- .cbins_apply_idx(cb$fit, dd[1:1000])
  expect_equal(as.vector(tapply(dd$y[1:1000], idx$x, mean)), e$woe, tolerance = 1e-12)
  # serial == parallel
  cb0 <- withr::with_options(list(scorecraft.parallel = "serial"),
    scr_bin_continuous(dd, "y", c("x", "z", "g"), train_idx = 1:1000, holdout_idx = 1001:1500, nthread = 2L))
  expect_equal(cb0$fit$results, cb$fit$results)
  # missing numeric values are refused, categorical NA becomes a level
  dd2 <- data.table::copy(dd); dd2$x[3] <- NA
  expect_error(scr_bin_continuous(dd2, "y", "x"), "missing values")
  dd3 <- data.table::copy(dd); dd3$g[1:80] <- NA
  cb3 <- scr_bin_continuous(dd3, "y", "g")
  expect_true(any(grepl("NA", cb3$fit$results$g$bin)))
  # logit scale
  cbl <- scr_bin_continuous(dd, "y", "x", scale = "logit")
  expect_equal(cbl$fit$results$x$woe, stats::qlogis(cbl$fit$results$x$mean), tolerance = 1e-12)
  # degenerate driver: one bin, not converged
  dd$k <- 1
  cbk <- scr_bin_continuous(dd, "y", "k")
  expect_equal(length(cbk$fit$results$k$bin), 1L); expect_false(cbk$fit$results$k$converged)
  skip_if_not_installed("duckdb"); skip_if_not_installed("DBI")
  con <- DBI::dbConnect(duckdb::duckdb()); on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  DBI::dbWriteTable(con, "t", as.data.frame(dd[1001:1500, list(x, z, g)]))
  sql <- paste(OptimalBinningWoE::obwoe_sql(cb$fit, table = "t", dialect = "duckdb", output = "woe"), collapse = "\n")
  r <- DBI::dbGetQuery(con, sql)
  expect_equal(r$x_woe, mine$x_woe, tolerance = 1e-12)
  expect_equal(r$g_woe, mine$g_woe, tolerance = 1e-12)
})
