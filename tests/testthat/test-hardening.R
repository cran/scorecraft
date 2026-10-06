# The three points that were only partially covered: PSOCK backend (Windows
# semantics) on every platform, the CSI timeline without scr_monitor(), and
# the monitoring plan as a contract that scr_monitor() actually reads.

test_that("the PSOCK backend gives results identical to serial and fork", {
  sp <- scr_split(head(scr_demo, 2000), "default", date_col = "ref_date", drop = c("id", "churn"))
  tr <- scr_triage(sp, cfg_test())
  feats <- setdiff(names(tr$clean), "default")
  run <- function(backend) withr::with_options(list(scorecraft.parallel = backend), {
    fit <- fit_binning(tr$clean[sp$train_idx], "default", feats, cfg_test(nthread = 2L))
    fit <- strip_failed_features(fit)
    binned <- fit$summary$feature
    w_tr <- apply_woe(fit, tr$clean[sp$train_idx], binned); w_ho <- apply_woe(fit, tr$clean[sp$holdout_idx], binned)
    ho <- holdout_check(w_tr, w_ho, tr$clean$default[sp$train_idx], tr$clean$default[sp$holdout_idx], binned, cfg_test(nthread = 2L))
    m <- scr_metrics(w_tr[[paste0(binned[1], "_woe")]], tr$clean$default[sp$train_idx], n_boot = 20, seed = 3, nthread = 2L)
    list(fit = fit, ho = ho, m = m)
  })
  expect_identical(.scr_backend(), if (.Platform$OS.type == "unix") "fork" else "psock")
  withr::with_options(list(scorecraft.parallel = "psock"), expect_identical(.scr_backend(), "psock"))
  a <- run("serial"); b <- run("fork"); c <- run("psock")
  for (f in names(a$fit$results)) {
    expect_identical(a$fit$results[[f]]$bin, c$fit$results[[f]]$bin)
    expect_equal(a$fit$results[[f]]$woe, c$fit$results[[f]]$woe)
    expect_equal(b$fit$results[[f]]$woe, c$fit$results[[f]]$woe)
  }
  expect_equal(a$ho, c$ho); expect_equal(b$ho, c$ho)
  expect_equal(unclass(a$m), unclass(c$m)); expect_equal(unclass(b$m), unclass(c$m))
  # a worker error is re-thrown with the failing item under PSOCK as well
  withr::with_options(list(scorecraft.parallel = "psock"),
    expect_error(.scr_lapply(1:3, function(i) if (i == 2) stop("boom") else i, nthread = 2L), "boom"))
  # the whole scorecard, PSOCK vs serial
  res <- res_demo()
  s1 <- withr::with_options(list(scorecraft.parallel = "serial"), scr_scorecard(res, n_boot = 10))
  s2 <- withr::with_options(list(scorecraft.parallel = "psock"), { cfg <- res$config; cfg$nthread <- 2L; r2 <- res; r2$config <- cfg; scr_scorecard(r2, n_boot = 10) })
  expect_equal(s1$points$points, s2$points$points)
  expect_equal(s1$metrics$auc, s2$metrics$auc)
  expect_equal(s1$metrics$auc_lo, s2$metrics$auc_lo)
  expect_equal(s1$stability$variables$csi, s2$stability$variables$csi)
})

test_that("triage is parallel by column and identical to the serial profile", {
  sp <- scr_split(scr_demo, "default", date_col = "ref_date", drop = c("id", "churn"))
  t1 <- withr::with_options(list(scorecraft.parallel = "serial"), scr_triage(sp, cfg_test(nthread = 2L)))
  t2 <- withr::with_options(list(scorecraft.parallel = "fork"),   scr_triage(sp, cfg_test(nthread = 2L)))
  t3 <- withr::with_options(list(scorecraft.parallel = "psock"),  scr_triage(sp, cfg_test(nthread = 2L)))
  expect_equal(t1$profile, t2$profile); expect_equal(t1$profile, t3$profile)
  expect_equal(t1$ledger, t2$ledger);   expect_equal(t1$ledger, t3$ledger)
  expect_identical(t1$keep, t2$keep);   expect_identical(t1$derived, t3$derived)
  expect_equal(t1$clean, t2$clean)
  # fork_only tasks stay serial under PSOCK, in the main process
  withr::with_options(list(scorecraft.parallel = "psock"),
    expect_identical(.scr_lapply(1:3, function(i) Sys.getpid(), nthread = 2L, fork_only = TRUE), as.list(rep(Sys.getpid(), 3))))
})

test_that("the fork backend caps the workers by the memory available", {
  expect_identical(withr::with_options(list(scorecraft.fork_mem_fraction = Inf), .scr_fork_cap(20L)), 20L)
  expect_identical(withr::with_options(list(scorecraft.fork_mem_fraction = "x"), .scr_fork_cap(4L)), 4L)
  skip_if_not(file.exists("/proc/meminfo") && file.exists("/proc/self/statm"), "Linux /proc only")
  expect_true(is.finite(.scr_mem_available()) && .scr_mem_available() > 0)
  expect_true(is.finite(.scr_rss()) && .scr_rss() > 0)
  k <- withr::with_options(list(scorecraft.fork_mem_fraction = 0.75), .scr_fork_cap(2L))
  expect_true(k >= 1L && k <= 2L)
  # a budget smaller than one resident copy leaves a single worker: the call runs in the parent
  tiny <- 0.5 * .scr_rss() / .scr_mem_available()
  withr::with_options(list(scorecraft.parallel = "fork", scorecraft.fork_mem_fraction = tiny), {
    expect_identical(.scr_fork_cap(8L), 1L)
    old <- scr_verbose(TRUE); on.exit(scr_verbose(old), add = TRUE)
    expect_message(.scr_fork_cap(8L), "capped at 1 of 8")
    expect_identical(.scr_lapply(1:3, function(i) Sys.getpid(), nthread = 2L), as.list(rep(Sys.getpid(), 3)))
  })
  # a generous budget keeps the requested workers
  withr::with_options(list(scorecraft.parallel = "fork", scorecraft.fork_mem_fraction = Inf),
    expect_false(any(unlist(.scr_lapply(1:2, function(i) Sys.getpid(), nthread = 2L)) == Sys.getpid())))
})

test_that("the CSI timeline by vintage comes from the stored hold-out bins and matches scr_monitor()", {
  sc <- sc_demo()
  expect_s3_class(sc$holdout_bins, "data.table")
  expect_setequal(names(sc$holdout_bins), sc$features)
  expect_equal(nrow(sc$holdout_bins), nrow(sc$samples$holdout))
  expect_true(all(vapply(sc$holdout_bins, function(v) all(!is.na(v)), logical(1))))
  tl <- .timelines(sc)
  periods <- sort(unique(as.character(sc$samples$holdout$date)))
  expect_equal(length(periods), 2L)
  expect_setequal(unique(tl$csi$period), periods)
  expect_equal(nrow(tl$csi), length(periods) * length(sc$features))
  expect_true(all(c("csi", "flag_fixed", "critical", "flag_adjusted", "points_shift") %in% names(tl$csi)))
  # the same numbers as scr_monitor() on the hold-out rows
  ho_idx <- res_demo()$split$holdout_idx
  mo <- scr_monitor(sc, scr_demo[ho_idx, ], date_col = "ref_date")
  a <- tl$csi[order(period, variable)]; b <- mo$csi[order(period, variable)]
  expect_equal(a$csi, b$csi, tolerance = 1e-12)
  expect_equal(a$points_shift, b$points_shift, tolerance = 1e-12)
  expect_identical(a$flag_adjusted, b$flag_adjusted)
  expect_equal(tl$psi$psi, mo$psi$psi, tolerance = 1e-12)
  # and the exported sheet is the timeline, not a single hold-out row
  skip_if_not_installed("openxlsx")
  out <- scr_export(sc, file.path(tempdir(), "scr-csi-tl"), stamp = FALSE)
  sheet <- openxlsx::read.xlsx(out$files$validation, sheet = "Stability_CSI_Timeline")
  expect_equal(nrow(sheet), length(periods) * length(sc$features))
  expect_setequal(unique(sheet$period), periods)
})

test_that("scr_monitor() reads the monitoring plan: from the scorecard, a table, or the strategy workbook", {
  sc <- sc_demo()
  plan <- sc$monitoring_plan
  expect_s3_class(plan, "scr_monitoring_plan")
  expect_true(all(c("psi_score_fixed_moderate", "min_events_per_period") %in% plan$item))
  mo0 <- scr_monitor(sc, scr_demo, date_col = "ref_date", target = "default", n_boot = 5)
  expect_equal(mo0$thresholds$psi, c(0.10, 0.25))
  expect_equal(mo0$thresholds$alpha, sc$config$psi_alpha)
  expect_true(all(mo0$psi$flag_fixed == "stable"))
  expect_identical(mo0$vintage$status, ifelse(mo0$vintage$events < 100, "insufficient", "ok"))
  # tighten the plan: the flags must follow the plan, not the configuration
  tight <- plan
  tight$value[tight$item == "psi_score_fixed_moderate"] <- "0.005"
  tight$value[tight$item == "psi_score_fixed_action"]   <- "0.02"
  tight$value[tight$item == "csi_variable_fixed_moderate"] <- "0.001"
  tight$value[tight$item == "csi_variable_fixed_action"]   <- "0.30"
  tight$value[tight$item == "psi_adjusted_alpha"] <- "0.20"
  tight$value[tight$item == "min_events_per_period"] <- "1000"
  mo1 <- scr_monitor(sc, scr_demo, date_col = "ref_date", target = "default", n_boot = 5, plan = tight)
  expect_equal(mo1$thresholds$psi, c(0.005, 0.02))
  expect_true(any(mo1$psi$flag_fixed != "stable"))
  expect_true(any(mo1$csi$flag_fixed == "moderate"))
  expect_true(all(mo1$vintage$status == "insufficient"))
  expect_true(all(mo1$psi$critical < mo0$psi$critical))   # alpha 0.20 lowers the critical value
  expect_output(print(mo1), "PSI 0.005/0.02")
  expect_output(print(mo1), "insufficient events")
  # an explicit alpha still wins over the plan
  mo2 <- scr_monitor(sc, scr_demo, date_col = "ref_date", alpha = 0.01, plan = tight)
  expect_equal(mo2$thresholds$alpha, 0.20); expect_true(all(mo2$psi$critical > mo1$psi$critical))
  expect_error(scr_monitor(sc, scr_demo, plan = data.frame(a = 1)), "item")
  bad <- plan; bad$value[bad$item == "psi_score_fixed_action"] <- "0.01"
  expect_error(scr_monitor(sc, scr_demo, plan = bad), "must exceed")
  bad2 <- plan; bad2$value[bad2$item == "psi_adjusted_alpha"] <- "x"
  expect_error(scr_monitor(sc, scr_demo, plan = bad2), "not numeric")
  # round trip through the strategy workbook written by scr_export()
  skip_if_not_installed("openxlsx")
  out <- file.path(tempdir(), "scr-plan"); unlink(out, recursive = TRUE)
  sc_t <- sc; sc_t$monitoring_plan <- tight
  ex <- scr_export(sc_t, out, stamp = FALSE)
  back <- openxlsx::read.xlsx(ex$files$strategy, sheet = "Monitoring_Plan")
  expect_equal(back$value[back$item == "psi_score_fixed_action"], "0.02")
  mo3 <- scr_monitor(sc, scr_demo, date_col = "ref_date", target = "default", n_boot = 5, plan = ex$files$strategy)
  expect_equal(mo3$thresholds, mo1$thresholds)
  expect_identical(mo3$psi$flag_fixed, mo1$psi$flag_fixed)
  expect_identical(mo3$csi$flag_fixed, mo1$csi$flag_fixed)
  expect_identical(mo3$vintage$status, mo1$vintage$status)
  # the plan written by export without a monitor object is the scorecard's own
  ex0 <- scr_export(sc, out, stamp = FALSE)
  back0 <- openxlsx::read.xlsx(ex0$files$strategy, sheet = "Monitoring_Plan")
  expect_equal(back0$value, plan$value)
})

test_that("workers never fail silently: errors, deaths, warnings, NULLs and data.tables round-trip on every backend", {
  for (b in c("fork", "psock")) withr::with_options(list(scorecraft.parallel = b, scorecraft.fork_mem_fraction = Inf), {
    if (b == "fork" && .Platform$OS.type != "unix") skip("fork is unix only")
    # a legitimate NULL is a value, not a dead worker
    expect_identical(.scr_lapply(1:3, function(i) NULL, nthread = 2L), vector("list", 3))
    # warnings raised inside a worker reach the parent, with the values intact
    expect_warning(v <- .scr_lapply(1:3, function(i) { if (i == 2L) warning("wobble ", i); i * 10L }, nthread = 2L), "wobble 2")
    expect_identical(v, list(10L, 20L, 30L))
    # errors are re-thrown with the failing item
    expect_error(.scr_lapply(1:3, function(i) if (i == 3L) stop("boom") else i, nthread = 2L), "item 3.*boom")
    # data.table threads: one per worker, never the parent's default
    expect_true(all(unlist(.scr_lapply(1:2, function(i) data.table::getDTthreads(), nthread = 2L)) == 1L))
    # a data.table returned by a worker accepts := without the selfref warning
    out <- .scr_lapply(1:2, function(i) data.table::data.table(a = i), nthread = 2L)
    expect_silent(out[[1]][, b := a + 1L])
    expect_identical(out[[1]]$b, 2L)
    # a worker killed by the system is reported, not returned as NULL
    expect_error(.scr_lapply(1:4, function(i) { if (i == 1L) tools::pskill(Sys.getpid(), tools::SIGKILL); i }, nthread = 2L),
                 "died|PSOCK")
  })
  # the real pipeline under a worker warning: nothing lost, nothing duplicated
  withr::with_options(list(scorecraft.parallel = "psock"), {
    r <- suppressWarnings(.scr_lapply(1:6, function(i) { if (i %% 2L == 0L) warning("even"); i }, nthread = 2L))
    expect_identical(unlist(r), 1:6)
    expect_warning(.scr_lapply(1:6, function(i) { if (i %% 2L == 0L) warning("even"); i }, nthread = 2L), "even")
  })
})
