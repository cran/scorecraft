# The eleven items the documentation audit flagged as needing a code change:
# each one with its argument validation and the flow it changes.

test_that("scr_iv() ignores NA for every group type and validates its arguments", {
  set.seed(1)
  y <- rbinom(400, 1, 0.3)
  g_chr <- sample(c("a", "b", "c"), 400, TRUE); g_chr[1:40] <- NA
  g_fac <- factor(g_chr)
  g_int <- as.integer(factor(g_chr))            # NA stays NA
  ref <- scr_iv(g_chr[!is.na(g_chr)], y[!is.na(g_chr)])
  expect_equal(scr_iv(g_chr, y), ref)
  expect_equal(scr_iv(g_fac, y), ref)
  expect_equal(scr_iv(g_int, y), ref)
  y2 <- y; y2[c(5, 9)] <- NA
  expect_equal(scr_iv(g_chr, y2), scr_iv(g_chr[!is.na(y2)], y[!is.na(y2)]))
  expect_equal(scr_iv(rep(NA_character_, 10), rbinom(10, 1, 0.5)), 0)
  expect_error(scr_iv(g_chr, y[-1]), "same length")
  expect_error(scr_iv(g_chr, y, laplace = -1), "laplace")
})

test_that("scr_classing_read() validates the separator and the spec carries it into the import", {
  res <- res_demo()
  lab <- scr_coarse_classing(res)
  f <- file.path(tempdir(), "spec-audit.csv")
  scr_classing_spec(lab, file = f)
  sp <- scr_classing_read(f)
  expect_identical(attr(sp, "sep"), "%;%")
  expect_error(scr_classing_read(f, sep = ""), "sep")
  expect_error(scr_classing_read(f, sep = c("a", "b")), "sep")
  expect_error(scr_classing_read(file.path(tempdir(), "does-not-exist.csv")), "existing")
  # a categorical bin with an empty category around the separator, and a category in two bins
  d <- utils::read.csv(f, stringsAsFactors = FALSE)
  cat_var <- unique(d$variable[d$type != "numeric"])[1]
  skip_if(is.na(cat_var), "no categorical variable in the lab")
  rows <- which(d$variable == cat_var)
  bad <- d; bad$categories[rows[1]] <- paste0(bad$categories[rows[1]], "%;%")
  f2 <- file.path(tempdir(), "spec-bad1.csv"); utils::write.csv(bad, f2, row.names = FALSE)
  expect_error(scr_classing_read(f2), "empty category")
  bad2 <- d; first_cat <- strsplit(bad2$categories[rows[1]], "%;%", fixed = TRUE)[[1]][1]
  bad2$categories[rows[2]] <- paste(bad2$categories[rows[2]], first_cat, sep = "%;%")
  f3 <- file.path(tempdir(), "spec-bad2.csv"); utils::write.csv(bad2, f3, row.names = FALSE)
  expect_error(scr_classing_read(f3), "two bins")
  # a spec read with another separator is refused by the lab
  sp2 <- sp; attr(sp2, "sep") <- "|"
  expect_error(scr_classing_import(lab, sp2), "separator")
  # the lab's own config$verbose governs the message
  lab_v <- lab; lab_v$result$config$verbose <- TRUE
  expect_message(scr_classing_import(lab_v, sp), "nothing to import")
  lab_q <- lab; lab_q$result$config$verbose <- FALSE
  expect_no_message(scr_classing_import(lab_q, sp))
})

test_that("TOO_MANY_BINS fires when a fit exceeds max_bins", {
  res <- res_demo()
  cfg <- res$config
  fit <- res$fit
  f <- intersect(names(fit$results), grep("^vl_score", names(res$data_clean), value = TRUE))[1]
  tr <- res$data_clean[res$split$train_idx]
  x <- tr[[f]]; y <- tr[[res$target]]
  br <- unique(stats::quantile(x, probs = seq(0.1, 0.9, by = 0.1), names = FALSE))
  fit$results[[f]] <- .manual_entry_num(f, x, y, breaks = br)
  s <- data.table::as.data.table(fit$summary); s[feature == f, n_bins := length(br) + 1L]; fit$summary <- as.data.frame(s)
  scr <- screen_features(fit, cfg)
  row <- scr$summary[feature == f]
  expect_gt(length(br) + 1L, cfg$max_bins)
  expect_true(grepl("TOO_MANY_BINS", row$reason))
  # the pipeline itself never produces more than max_bins, so the rule is silent there
  expect_false(any(grepl("TOO_MANY_BINS", res$screen$summary$reason)))
})

test_that("scr_psi() stores and prints its thresholds and validates them", {
  set.seed(2)
  a <- rnorm(2000); b <- rnorm(2000, 0.3)
  p <- scr_psi(a, b, thresholds = c(0.05, 0.20), alpha = 0.10)
  expect_equal(p$thresholds, c(0.05, 0.20)); expect_equal(p$alpha, 0.10)
  expect_output(print(p), "0.05/0.2")
  expect_error(scr_psi(a, b, thresholds = c(0.25, 0.10)), "thresholds")
  expect_error(scr_psi(a, b, thresholds = 0.1), "thresholds")
  expect_error(scr_psi(a, b, alpha = 1.5), "alpha")
  p0 <- scr_psi(numeric(), numeric())
  expect_true(is.na(p0$psi)); expect_equal(p0$thresholds, c(0.10, 0.25))
})

test_that("scr_default_rate() benchmarks an adjusted long-run average and keeps one mean", {
  cfg <- scr_config(verbose = FALSE)
  d <- scr_default(scr_demo_panel, "id", "ref_date", dpd = "dpd", config = cfg)
  dr <- scr_default_rate(d, config = cfg)
  expect_null(dr$lra$all_mean)
  expect_equal(dr$lra$benchmark, max(dr$lra$recent5_mean, dr$lra$mean))
  expect_false(dr$lra$flag_below_benchmark); expect_null(dr$lra$adjusted)
  lo <- scr_default_rate(d, lra_adjusted = dr$lra$benchmark / 2, config = cfg)
  expect_true(lo$lra$flag_below_benchmark)
  expect_output(print(lo), "BELOW the benchmark")
  hi <- scr_default_rate(d, lra_adjusted = dr$lra$benchmark + 0.01, config = cfg)
  expect_false(hi$lra$flag_below_benchmark)
  expect_output(print(hi), "at or above")
  expect_error(scr_default_rate(d, lra_adjusted = 1.5, config = cfg), "lra_adjusted")
  expect_error(scr_default_rate(d, horizon = 0, config = cfg), "horizon")
})

test_that("one asset_class key serves the PD floor, the LGD floor and the capital function", {
  cfg <- scr_config(verbose = FALSE)
  expect_equal(cfg$asset_class, "retail_other")
  expect_null(cfg$pd_asset_class); expect_null(cfg$capital_asset_class)
  expect_error(scr_config(asset_class = "boats"), "asset_class")
  expect_true("asset_class" %in% scr_config_keys(stage = 8)$key)
  pd <- pd_model()
  expect_equal(pd$asset_class, pd$config$asset_class %||% "retail_other")
  m <- scr_lgd_floor(lgd_demo(), asset_class = NULL)
  expect_equal(m$floors$asset_class %||% m$config$asset_class, "retail_other")
  port <- data.table::as.data.table(scr_demo_portfolio)[1:300]
  cq <- scr_capital(port, segment = "segment", config = scr_config(verbose = FALSE, asset_class = "qrre_revolver", capital_sensitivity = FALSE), keep_rows = TRUE)
  expect_true(all(cq$exposures$asset_class == "qrre_revolver"))
})

test_that("scr_lgd_downturn() always records a reason", {
  m <- lgd_demo()
  per <- data.frame(start = as.Date("2022-01-01"), end = as.Date("2023-12-31"))
  expect_error(scr_lgd_downturn(m, periods = per), "reason")
  expect_error(scr_lgd_downturn(m, periods = per, reason = ""), "reason")
  expect_error(scr_lgd_downturn(m, periods = per, reason = c("a", "b")), "reason")
  d <- scr_lgd_downturn(m, periods = per, reason = "stress window")
  expect_equal(d$downturn$reason, "stress window")
  expect_true(any(grepl("stress window", d$ledger$reason %||% d$ledger$detail %||% "")))
})

test_that("the traffic-light convention is the same in PD, LGD and EAD: red at or below, amber at or below", {
  l <- c(0.01, 0.05)
  expect_equal(.pd_light(c(0.005, 0.01, 0.03, 0.05, 0.051, NA), l), c("red", "red", "amber", "amber", "green", NA))
  v <- scr_pd_validate(pd_model(), pd_panel(), score = "score", tests = "jeffreys")
  expect_true(all(v$summary$light[!is.na(v$summary$p_value)] == .pd_light(v$summary$p_value[!is.na(v$summary$p_value)], v$lights)))
})

test_that("scr_apply.scr_ead() takes `what` like the PD and LGD methods", {
  e <- ead_model()
  nd <- head(scr_demo_ead, 20)
  a_all <- scr_apply(e, nd)
  a_ead <- scr_apply(e, nd, what = "ead")
  a_pool <- scr_apply(e, nd, what = "pool")
  expect_setequal(names(a_ead), c("pool", "measure", "ccf_applied", "ead_predicted", "ead_floor_binding"))
  expect_setequal(names(a_pool), c("pool", "measure"))
  expect_true(all(c("utilisation", "undrawn", "ead_model", "ead_floor") %in% names(a_all)))
  expect_equal(a_all$ead_predicted, a_ead$ead_predicted); expect_equal(a_all$pool, a_pool$pool)
  expect_error(scr_apply(e, nd, what = "points"), "arg")
})

test_that("verbosity reaches the database functions", {
  skip_if_not_installed("DBI"); skip_if_not_installed("RSQLite")
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:"); on.exit(DBI::dbDisconnect(con), add = TRUE)
  DBI::dbWriteTable(con, "t", head(scr_demo, 300))
  old <- scr_verbose(TRUE); on.exit(scr_verbose(old), add = TRUE)
  expect_message(scr_fetch(con, "t"), "SQL:")
  expect_silent(scr_fetch(con, "t", verbose = FALSE))
  expect_message(scr_fetch(con, "t", verbose = TRUE), "SQL:")
  expect_true(scr_verbose())                                   # restored after the call
  expect_error(scr_fetch(con, "t", verbose = "yes"), "verbose")
  scr_verbose(FALSE)
  expect_message(scr_fetch(con, "t", verbose = TRUE), "SQL:")
  expect_silent(scr_fetch(con, "t"))
  # scr_run() follows config$verbose and restores the state afterwards
  scr_verbose(TRUE)
  rs <- NULL
  expect_silent(rs <- scr_run(con, "t", targets = "default", config = cfg_test(), drop = c("id", "churn", "ref_date")))
  expect_true(scr_verbose())
  expect_s3_class(rs, "scr_runset")
})
