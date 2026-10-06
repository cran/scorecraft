# LGD module: workout engine, two-stage model, pools, downturn, floors, ELBE,
# production (R and SQL), validation and export.

test_that("the present value and the realised LGD of a hand-made workout are exact", {
  cfg <- lgd_cfg(lgd_discount_rate = 0.12, lgd_discount_add_on = 0)
  wo <- scr_workout(hand_defaults(), hand_cashflows(), config = cfg)
  expect_s3_class(wo, "scr_workout")
  r <- wo$rds; data.table::setkey(r, default_id)
  expect_equal(r["A", discount_rate], 0.12)
  # A: 600 at 6 months, 300 at 12, a cost of 50 at 1; monthly factor 1.01
  pv_rec <- 600 / 1.01^6 + 300 / 1.01^12
  expect_equal(r["A", pv_recovery], pv_rec, tolerance = 1e-12)
  expect_equal(r["A", pv_cost], 50 / 1.01, tolerance = 1e-12)
  expect_equal(r["A", lgd_raw], (1000 - pv_rec + 50 / 1.01) / 1000, tolerance = 1e-12)
  expect_equal(r["A", lgd_raw], 0.2180430, tolerance = 1e-6)
  expect_equal(r["A", months_in_default], 12L)
  # B: one recovery of 500 at 6 months on 2,000
  expect_equal(r["B", lgd_raw], (2000 - 500 / 1.01^6) / 2000, tolerance = 1e-12)
  # C: a cure with 100 paid at month 1 and the outstanding 900 as an artificial recovery at the cure (month 3)
  expect_true(r["C", is_cure])
  expect_equal(r["C", recovery_artificial], 900)
  expect_equal(r["C", lgd_raw], (1000 - 100 / 1.01 - 900 / 1.01^3) / 1000, tolerance = 1e-12)
  expect_equal(r["C", lgd_raw], 0.027458966, tolerance = 1e-6)
  expect_equal(r["C", months_in_default], 3L)
  # long-run averages: default-weighted is the plain mean, exposure-weighted uses the EAD
  lg <- r[c("A", "B", "C"), lgd_real]; e <- r[c("A", "B", "C"), ead]
  expect_equal(wo$summary$lra_default_weighted, mean(lg))
  expect_equal(wo$summary$lra_exposure_weighted, sum(lg * e) / sum(e))
  expect_lt(wo$summary$lra_default_weighted, wo$summary$lra_exposure_weighted)   # the large exposure carries the large loss
  expect_equal(wo$summary$n, 3L); expect_equal(wo$summary$cure_rate, 1 / 3)
  expect_equal(wo$funnel[rule == "INCOMPLETE_EXTRAPOLATED", n], 0L)
  expect_output(print(wo), "long-run average")
  # the rates table: the rate in force at the default date
  rt <- data.frame(date = as.Date(c("2023-06-01", "2024-01-01", "2024-06-01")), rate = c(0.05, 0.10, 0.20))
  wo2 <- scr_workout(hand_defaults(), hand_cashflows(), rates = rt, config = lgd_cfg(lgd_discount_add_on = 0.02))
  expect_equal(unique(wo2$rds$discount_rate), 0.12)
  expect_error(scr_workout(hand_defaults(), hand_cashflows(), config = lgd_cfg()), "lgd_discount_rate")
  # indirect costs by count: 100 to each of the three events
  wo3 <- scr_workout(hand_defaults(), hand_cashflows(), config = lgd_cfg(lgd_discount_rate = 0.12, lgd_discount_add_on = 0,
                                                                          lgd_cost_allocation = "count"), indirect_costs = 300)
  expect_equal(wo3$rds$cost_indirect, rep(100, 3))
  expect_equal(wo3$rds[default_id == "A", lgd_raw], r["A", lgd_raw] + 0.1, tolerance = 1e-12)
  expect_error(scr_workout(hand_defaults()[, -1], hand_cashflows(), config = cfg), "lacks column")
  d <- hand_defaults(); d$status[1] <- "gone"
  expect_error(scr_workout(d, hand_cashflows(), config = cfg), "status")
})

test_that("two defaults of one facility inside the window are one event; outside they are two", {
  d <- rbind(hand_defaults(), data.frame(default_id = "A2", facility_id = "f1", default_date = as.Date("2024-09-15"), ead = 900,
                                         product = "p", status = "closed", close_date = as.Date("2025-09-15"), stringsAsFactors = FALSE))
  d$status[1] <- "cured"; d$close_date[1] <- as.Date("2024-04-15")
  cf <- rbind(hand_cashflows()[hand_cashflows()$default_id != "A", ],
              data.frame(default_id = c("A", "A2"), date = as.Date(c("2024-02-15", "2025-03-15")), amount = c(100, 400),
                         type = "recovery", stringsAsFactors = FALSE))
  cfg <- lgd_cfg(lgd_discount_rate = 0.12, lgd_discount_add_on = 0)
  wo <- scr_workout(d, cf, config = cfg)           # window 9 months: the gap from 2024-04-15 to 2024-09-15 is 5
  expect_equal(nrow(wo$rds), 3L)
  expect_equal(wo$funnel[rule == "MULTIPLE_DEFAULT_MERGED", n], 1L)
  a <- wo$rds[default_id == "A"]
  expect_equal(a$merged_n, 2L); expect_equal(a$absorbed, "A2")
  expect_equal(a$status, "closed"); expect_false(a$is_cure)
  expect_equal(a$ead, 1000); expect_equal(a$months_in_default, 20L)
  expect_equal(a$recovery_artificial, 0)
  expect_equal(a$lgd_raw, (1000 - 100 / 1.01 - 400 / 1.01^14) / 1000, tolerance = 1e-12)
  wo3 <- scr_workout(d, cf, config = lgd_cfg(lgd_discount_rate = 0.12, lgd_discount_add_on = 0, lgd_cure_window = 3L))
  expect_equal(nrow(wo3$rds), 4L)
  expect_true(wo3$rds[default_id == "A", is_cure])
  expect_equal(wo3$funnel[rule == "MULTIPLE_DEFAULT_MERGED", n], 0L)
  expect_equal(.months_between(as.Date("2024-01-31"), as.Date("2024-02-29")), 0L)
  expect_equal(.months_between(as.Date("2024-01-15"), as.Date("2025-01-15")), 12L)
})

test_that("an incomplete workout receives exactly the remaining recovery of the closed profile", {
  d <- rbind(hand_defaults(), data.frame(default_id = "D", facility_id = "f4", default_date = as.Date("2024-07-15"), ead = 1000,
                                         product = "p", status = "open", close_date = as.Date(NA), stringsAsFactors = FALSE))
  cf <- rbind(hand_cashflows(), data.frame(default_id = "D", date = as.Date("2024-10-15"), amount = 100, type = "recovery",
                                           stringsAsFactors = FALSE))
  cfg <- lgd_cfg(lgd_discount_rate = 0.12, lgd_discount_add_on = 0)
  wo <- scr_workout(d, cf, config = cfg, obs_date = as.Date("2025-01-15"))
  r <- wo$rds[default_id == "D"]
  expect_true(r$is_incomplete); expect_equal(r$months_in_default, 6L)
  # profile of product p (closed A and B, EAD 3,000): 1,100 at month 6, 300 more at month 12
  prof <- wo$recovery_profile[product == "p"]
  expect_equal(prof[month == 6, cum_recovery], (1100 / 1.01^6) / 3000, tolerance = 1e-12)
  expect_equal(prof[month == 60, cum_recovery], (1100 / 1.01^6 + 300 / 1.01^12) / 3000, tolerance = 1e-12)
  expect_equal(r$recovery_extrapolated, 1000 * (300 / 1.01^12) / 3000, tolerance = 1e-12)
  expect_equal(wo$extrapolation$expected_further, r$recovery_extrapolated)
  expect_equal(r$lgd_raw, (1000 - 100 / 1.01^3 - 100 / 1.01^12) / 1000, tolerance = 1e-12)
  expect_equal(wo$funnel[rule == "INCOMPLETE_EXTRAPOLATED", n], 1L)
  # beyond t_max: closed with zero further recovery
  wo2 <- scr_workout(d, cf, config = lgd_cfg(lgd_discount_rate = 0.12, lgd_discount_add_on = 0, lgd_t_max = 6L), obs_date = as.Date("2025-01-15"))
  r2 <- wo2$rds[default_id == "D"]
  expect_false(r2$is_incomplete); expect_true(r2$closed_at_t_max); expect_equal(r2$recovery_extrapolated, 0)
  expect_equal(wo2$funnel[rule == "OPEN_BEYOND_T_MAX_CLOSED", n], 1L)
  # bounds: a negative raw LGD is floored, the raw value kept; the cap is optional
  d3 <- hand_defaults()[1, ]; cf3 <- data.frame(default_id = "A", date = as.Date("2024-02-15"), amount = 1500, type = "recovery")
  w3 <- scr_workout(d3, cf3, config = cfg)
  expect_lt(w3$rds$lgd_raw, 0); expect_equal(w3$rds$lgd_real, 0)
  expect_equal(w3$funnel[rule == "NEGATIVE_LGD_FLOORED", n], 1L)
  cf4 <- data.frame(default_id = "A", date = as.Date("2024-02-15"), amount = 1500, type = "direct_cost")
  w4 <- scr_workout(d3, cf4, config = lgd_cfg(lgd_discount_rate = 0.12, lgd_discount_add_on = 0, lgd_cap_at_one = TRUE))
  expect_gt(w4$rds$lgd_raw, 1); expect_equal(w4$rds$lgd_real, 1)
  expect_equal(w4$funnel[rule == "LGD_ABOVE_ONE_CAPPED", n], 1L)
})

test_that("the demo workout is consistent and reproducible", {
  wo <- wo_demo()
  expect_equal(nrow(wo$rds), 900L - wo$funnel[rule == "MULTIPLE_DEFAULT_MERGED", n])
  expect_equal(wo$funnel[rule == "MULTIPLE_DEFAULT_MERGED", n], 15L)
  expect_true(all(wo$rds$lgd_real >= 0))
  expect_true(all(wo$rds$lgd_real[wo$rds$is_cure] < 0.25))
  expect_gt(wo$summary$by_product[product == "unsecured", lra], wo$summary$by_product[product == "mortgage", lra])
  expect_true(all(wo$recovery_profile[, diff(cum_recovery) >= 0, by = product]$V1))
  expect_equal(nrow(wo$recovery_profile), 4L * (wo$config$lgd_t_max + 1L))
  expect_true(all(wo$rds$discount_rate > 0.05))
  expect_equal(wo$summary$lra_default_weighted, mean(wo$rds$lgd_real))
  wk <- scr_workout(scr_demo_lgd, scr_demo_lgd_cashflows, rates = scr_demo_rates, config = lgd_cfg(), keep_rows = TRUE)
  expect_true(is.data.frame(wk$cashflows)); expect_true(all(wk$cashflows$pv <= wk$cashflows$amount + 1e-9))
})

test_that("the two-stage model fits, scores every default, and the pools are monotone", {
  m <- lgd_demo()
  expect_s3_class(m, "scr_lgd")
  expect_equal(m$split$method, "cohort")
  expect_true(length(m$cure$features) >= 1L); expect_true(length(m$severity$features) >= 1L)
  expect_true(all(m$cure$sign_check[action == "kept", coef] > 0))
  expect_true(all(m$severity$sign_check[action == "kept", coef] > 0))
  expect_equal(m$severity$engine, "fractional_logit")
  expect_equal(nrow(m$scored), nrow(wo_demo()$rds))
  expect_true(all(m$scored$p_cure > 0 & m$scored$p_cure < 1))
  expect_true(all(m$scored$lgd_pred >= 0 & m$scored$lgd_pred <= 1))
  expect_gt(m$lgd_cure, 0); expect_lt(m$lgd_cure, 0.2)
  expect_setequal(m$metrics$sample, c("train", "holdout"))
  expect_true(all(m$metrics$gauc > 0.6)); expect_true(all(m$metrics$gauc_lo <= m$metrics$gauc & m$metrics$gauc <= m$metrics$gauc_hi))
  expect_true(all(m$metrics$lcr > 0.3)); expect_true(all(m$metrics$spearman > 0.4))
  p <- m$pools
  expect_true(all(diff(p$lra) > 0)); expect_true(all(diff(p$pred_hi[-nrow(p)]) > 0))
  expect_true(all(p$n >= m$config$lgd_min_defaults_bin))
  expect_equal(p$pred_lo[1], -Inf); expect_equal(p$pred_hi[nrow(p)], Inf)
  expect_false(any(p$pred_hi[-nrow(p)] %in% m$scored$lgd_pred))   # breaks sit between observed predictions
  expect_equal(p$lra_moc, p$lra + p$moc_c); expect_true(all(p$moc_c > 0))
  expect_equal(sum(p$n), m$split$n_train)
  expect_equal(as.integer(table(m$scored[sample == "train", pool])[as.character(p$pool)]), p$n)
  expect_true(all(m$scored$pool == .lgd_pool_of(m$scored$lgd_pred, p)))
  expect_equal(m$downturn$status, "provisional")
  expect_true(any(grepl("DOWNTURN_PENDING", m$ledger$reason)))
  expect_output(print(m), "gAUC")
  p4 <- scr_lgd_pools(m, n_pools = 4, min_defaults = 50)
  expect_true(nrow(p4) >= nrow(p)); expect_true(all(diff(p4$lra) > 0))
  # the metrics are the hand formulas
  s <- m$scored[sample == "holdout"]
  expect_equal(m$metrics[sample == "holdout", rmse], sqrt(mean((s$lgd_pred - s$lgd_real)^2)))
  expect_equal(m$metrics[sample == "holdout", somers_d], .lgd_somers(s$lgd_pred, s$lgd_real))
  # Somers' D and the loss capture ratio on tiny hand vectors
  expect_equal(.lgd_somers(c(1, 2, 3, 4), c(0.1, 0.2, 0.3, 0.4)), 1)
  expect_equal(.lgd_somers(c(4, 3, 2, 1), c(0.1, 0.2, 0.3, 0.4)), -1)
  expect_equal(.lgd_somers(c(1, 2, 3, 4), c(0.1, 0.1, 0.3, 0.3)), 1)     # ties on the realised value are excluded
  expect_equal(.lgd_lcr(c(0.9, 0.5, 0.1), c(0.9, 0.5, 0.1), c(1, 1, 1)), 1)
  expect_lt(.lgd_lcr(c(0.1, 0.5, 0.9), c(0.9, 0.5, 0.1), c(1, 1, 1)), 0)
  expect_error(scr_lgd(wo_demo(), drivers = "nope", config = lgd_cfg()), "not found")
  wo_na <- wo_demo(); wo_na$rds$ltv[1] <- NA
  expect_error(scr_lgd(wo_na, drivers = c("product", "ltv"), config = lgd_cfg()), "missing values")
})

test_that("the beta severity engine runs when betareg is installed and fails clearly otherwise", {
  skip_if_not_installed("betareg")
  mb <- scr_lgd(wo_demo(), drivers = c("product", "ltv", "prior_dpd_max"), config = lgd_cfg(lgd_severity = "beta"))
  expect_equal(mb$severity$engine, "beta")
  expect_true(is.finite(mb$severity$phi))
  expect_true(all(mb$scored$severity > 0 & mb$scored$severity < 1))
  expect_true(all(mb$severity$sign_check[action == "kept", coef] > 0))
})

test_that("downturn: observed impact per pool, the type-3 add-on, the reference value and the reasons", {
  m <- lgd_demo()
  per <- data.frame(start = as.Date("2022-01-01"), end = as.Date("2023-12-31"))
  d1 <- scr_lgd_downturn(m, periods = per, reason = "stress window")
  tb <- d1$downturn$table
  expect_equal(d1$downturn$status, "final")
  expect_true(all(tb$method_used %in% c("type1", "type3_fallback")))
  # observed: default-weighted realised LGD of the pool's training defaults inside the periods
  s <- m$scored[sample == "train" & default_date >= per$start & default_date <= per$end]
  obs <- s[, list(v = mean(lgd_real), n = .N), by = pool][order(pool)]
  expect_equal(tb$dt_observed[match(obs$pool, tb$pool)], obs$v)
  expect_equal(tb$n_downturn[match(obs$pool, tb$pool)], obs$n)
  expect_equal(tb$lgd_dt, pmin(1, pmax(tb$lra + tb$moc_c, tb$dt + tb$moc_c)))
  expect_true(all(tb$lgd_dt >= tb$lra + tb$moc_c - 1e-12))
  expect_equal(d1$pools$lgd_dt, tb$lgd_dt); expect_equal(d1$pools$lgd_final, tb$lgd_dt)
  # reference value: the mean of the two worst calendar years of the pool, training rows
  yr <- m$scored[sample == "train" & pool == 1, list(v = mean(lgd_real)), by = list(y = as.integer(format(default_date, "%Y")))]
  expect_equal(tb$reference_value[1], mean(sort(yr$v, decreasing = TRUE)[1:2]))
  d3 <- scr_lgd_downturn(m, method = "type3", reason = "downturn data too thin")
  expect_equal(d3$downturn$table$lgd_dt, pmin(1, m$pools$lra + 0.15 + m$pools$moc_c))
  expect_true(any(grepl("downturn data too thin", d3$ledger$reason)))
  d0 <- scr_lgd_downturn(m, method = "none", reason = "no downturn by decision")
  expect_equal(d0$downturn$table$lgd_dt, pmin(1, m$pools$lra + m$pools$moc_c))
  expect_error(scr_lgd_downturn(m, method = "type3"), "reason")
  expect_error(scr_lgd_downturn(m, periods = per), "reason")
  expect_error(scr_lgd_downturn(m, periods = per, reason = "  "), "reason")
  expect_error(scr_lgd_downturn(m, reason = "x"), "periods")
  expect_error(scr_lgd_downturn(m, periods = data.frame(start = as.Date("2023-01-01"), end = as.Date("2022-01-01")), reason = "x"), "start <= end")
})

test_that("floors: the blended input floor binds from below and edits of the parameters are recorded", {
  m <- lgd_final()
  f <- m$floors$table
  expect_equal(f$floor, rep(0.3 * 0.6 + 0.1 * 0.4, nrow(f)))
  expect_true(all(m$pools$lgd_final >= m$pools$floor)); expect_true(all(m$pools$lgd_final >= m$pools$lgd_dt))
  expect_equal(m$pools$lgd_final, pmax(m$pools$lgd_dt, m$pools$floor))
  expect_equal(f$binding, f$lgd_final > f$lgd_dt)
  expect_equal(m$floors$binding_share, sum(f$n[f$binding]) / sum(f$n))
  expect_false(m$floors$params_modified)
  # the mortgage floor applies to the whole exposure
  mm <- scr_lgd_floor(lgd_demo(), asset_class = "retail_mortgage", secured_share = 0.9)
  expect_equal(unique(mm$pools$floor), 0.05)
  # a modified parameter table that binds everywhere
  p <- scr_irb_params("bcb"); p$lgd_floor$unsecured[p$lgd_floor$asset_class == "retail_other"] <- 0.95
  mh <- scr_lgd_floor(lgd_demo(), params = p, asset_class = "retail_other")
  expect_true(mh$floors$params_modified); expect_equal(mh$floors$binding_share, 1)
  expect_true(all(mh$pools$lgd_final == 0.95))
  expect_error(scr_lgd_floor(lgd_demo(), asset_class = "boats"), "asset_class")
  expect_error(scr_lgd_floor(lgd_demo(), secured_share = 2), "secured_share")
})

test_that("ELBE and in-default LGD are consistent with the pools at tau = 0 and never below ELBE", {
  m <- lgd_final()
  e <- scr_elbe(m)
  expect_s3_class(e, "scr_elbe")
  expect_true(all(e$consistency$ok))
  expect_equal(e$consistency$elbe_0, m$pools$lra); expect_equal(e$consistency$lgd_in_default_0, m$pools$lgd_dt)
  expect_true(all(e$table$delta_ul >= 0)); expect_true(all(e$table$lgd_in_default >= e$table$elbe, na.rm = TRUE))
  expect_equal(nrow(e$table), nrow(m$pools) * length(m$config$lgd_elbe_grid))
  expect_true(all(e$table$recovered_share >= 0 & e$table$recovered_share <= 1 + 1e-12))
  expect_true(all(e$table[, diff(recovered_share) >= -1e-12, by = pool]$V1))
  e2 <- scr_elbe(m, grid = c(0, 12))
  expect_equal(e2$grid, c(0L, 12L))
  expect_output(print(e), "consistency")
  expect_error(scr_elbe(m, grid = -1), "non-negative")
})

test_that("scr_apply reproduces the scored sample and the SQL reproduces scr_apply on DuckDB and SQLite", {
  m <- lgd_final()
  ap <- scr_apply(m, wo_demo()$rds, what = "all")
  expect_equal(ap$lgd_pred, m$scored$lgd_pred); expect_equal(ap$pool, m$scored$pool)
  expect_equal(ap$lgd_final, m$pools$lgd_final[m$scored$pool])
  expect_named(scr_apply(m, head(scr_demo_lgd, 5)), c("pool", "lgd_lra", "lgd_dt", "lgd_final"))
  expect_named(scr_apply(m, head(scr_demo_lgd, 5), what = "lgd"), c("lgd_pred", "pool", "lgd_lra", "lgd_dt", "lgd_final"))
  expect_error(scr_apply(m, scr_demo_lgd[, c("product", "ltv")]), "driver")
  for (d in c("ansi", "databricks", "spark", "hive", "mysql", "mariadb", "sqlserver", "bigquery", "postgres", "oracle",
              "snowflake", "redshift", "duckdb", "sqlite")) {
    sq <- scr_sql(m, table = "t", dialect = d)
    txt <- paste(sq, collapse = "\n")
    expect_match(txt, "AS lgd_final"); expect_match(txt, "FROM t\n")
    for (f in union(m$cure$features, m$severity$features)) expect_match(txt, f, fixed = TRUE)
    if (d %in% c("sqlite", "sqlserver")) expect_false(grepl("GREATEST", txt)) else expect_match(txt, "GREATEST\\(lgd_dt, lgd_floor\\)")
  }
  tmp <- tempfile(fileext = ".sql"); expect_invisible(scr_sql(m, file = tmp)); expect_true(file.exists(tmp))
  skip_if_not_installed("DBI")
  new <- scr_demo_lgd
  exp <- scr_apply(m, new, what = "all")
  if (requireNamespace("duckdb", quietly = TRUE)) {
    con <- DBI::dbConnect(duckdb::duckdb(), config = list(threads = "1"))
    DBI::dbWriteTable(con, "lgd_t", new)
    got <- DBI::dbGetQuery(con, paste(scr_sql(m, table = "lgd_t", dialect = "duckdb"), collapse = "\n"))
    DBI::dbDisconnect(con, shutdown = TRUE)
    expect_identical(as.integer(got$pool), exp$pool)
    expect_equal(got$lgd_final, exp$lgd_final, tolerance = 1e-9)
    expect_equal(got$lgd_pred, exp$lgd_pred, tolerance = 1e-9)
    expect_equal(got$p_cure, exp$p_cure, tolerance = 1e-9)
  }
  if (requireNamespace("RSQLite", quietly = TRUE)) {
    con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
    d2 <- new; d2$default_date <- as.character(d2$default_date); d2$close_date <- as.character(d2$close_date)
    DBI::dbWriteTable(con, "lgd_t", d2)
    got <- DBI::dbGetQuery(con, paste(scr_sql(m, table = "lgd_t", dialect = "sqlite"), collapse = "\n"))
    DBI::dbDisconnect(con)
    expect_identical(as.integer(got$pool), exp$pool)
    expect_equal(got$lgd_final, exp$lgd_final, tolerance = 1e-9)
    expect_equal(got$lgd_pred, exp$lgd_pred, tolerance = 1e-9)
  }
})

test_that("the validation battery runs on the hold-out and on new data with traffic lights", {
  m <- lgd_final()
  v <- scr_lgd_validate(m)
  expect_s3_class(v, "scr_lgd_validation")
  expect_equal(v$sample, "holdout"); expect_equal(v$n, m$split$n_holdout)
  expect_equal(nrow(v$calibration), nrow(m$pools))
  expect_true(all(v$calibration$p >= 0 & v$calibration$p <= 1))
  expect_true(all(v$summary$light %in% c("green", "amber", "red", "grey")))
  expect_equal(nrow(v$summary), 9L)
  expect_equal(v$discrimination$gauc_init, m$metrics[sample == "train", gauc])
  expect_true(is.finite(v$portfolio$loss_shortfall))
  expect_equal(nrow(v$heterogeneity), nrow(m$pools) - 1L)
  expect_true(nrow(v$stability$drivers) == length(m$cure$features) + length(m$severity$features))
  expect_output(print(v), "calibration")
  # calibration by hand for the first pool
  s <- m$scored[sample == "holdout" & pool == 1]
  tt <- (mean(s$lgd_real) - m$pools$lra[1]) / (stats::sd(s$lgd_real) / sqrt(nrow(s)))
  expect_equal(v$calibration$t[1], tt); expect_equal(v$calibration$p[1], 1 - stats::pnorm(tt))
  v2 <- scr_lgd_validate(m, newdata = wo_demo())
  expect_equal(v2$sample, "newdata"); expect_equal(v2$n, nrow(wo_demo()$rds))
  expect_error(scr_lgd_validate(m, newdata = scr_demo_lgd), "lacks column")
})

test_that("results are identical under the serial and the PSOCK backend", {
  s <- lgd_demo()$scored[sample == "train"]
  m_ser <- withr::with_options(list(scorecraft.parallel = "serial"), .lgd_metrics(s$lgd_pred, s$lgd_real, s$ead, n_boot = 12L, seed = 7L, nthread = 2L))
  m_psk <- withr::with_options(list(scorecraft.parallel = "psock"),  .lgd_metrics(s$lgd_pred, s$lgd_real, s$ead, n_boot = 12L, seed = 7L, nthread = 2L))
  expect_equal(m_ser, m_psk)
  expect_true(is.finite(m_psk$gauc_lo))
  w_ser <- withr::with_options(list(scorecraft.parallel = "serial"), scr_workout(scr_demo_lgd, scr_demo_lgd_cashflows, rates = scr_demo_rates, config = lgd_cfg(nthread = 2L)))
  w_psk <- withr::with_options(list(scorecraft.parallel = "psock"),  scr_workout(scr_demo_lgd, scr_demo_lgd_cashflows, rates = scr_demo_rates, config = lgd_cfg(nthread = 2L)))
  expect_equal(w_ser$rds, w_psk$rds)
  skip_if_not(psock_has_dev_ns(), "the installed scorecraft lacks the IRB internals the PSOCK workers need")
  l_ser <- withr::with_options(list(scorecraft.parallel = "serial"), scr_lgd(w_ser, drivers = lgd_drivers(), config = lgd_cfg(nthread = 2L)))
  l_psk <- withr::with_options(list(scorecraft.parallel = "psock"),  scr_lgd(w_psk, drivers = lgd_drivers(), config = lgd_cfg(nthread = 2L)))
  expect_equal(l_ser$scored, l_psk$scored); expect_equal(l_ser$metrics, l_psk$metrics); expect_equal(l_ser$pools, l_psk$pools)
})

test_that("the LGD workbook round-trips through openxlsx", {
  skip_if_not_installed("openxlsx")
  m <- lgd_final()
  out <- file.path(tempdir(), "scr-lgd-export"); unlink(out, recursive = TRUE)
  r <- scr_export(m, out, stamp = FALSE)
  expect_setequal(names(r$files), c("xlsx", "sql"))
  expect_true(all(file.exists(unlist(r$files))))
  expect_equal(basename(r$files$xlsx), "lgd_model.xlsx")
  sheets <- openxlsx::getSheetNames(r$files$xlsx)
  expect_true(all(c("RDS_Funnel", "RDS_Summary", "Recovery_Profile", "Cure_Bins", "Cure_Coefficients", "Severity_Bins",
                    "Severity_Coefficients", "Sign_Check", "Pools", "Downturn", "Floors", "ELBE_Grid", "Validation",
                    "Model_Card", "Decision_Ledger") %in% sheets))
  pl <- openxlsx::read.xlsx(r$files$xlsx, sheet = "Pools")
  expect_equal(nrow(pl), nrow(m$pools)); expect_equal(pl$lgd_final, m$pools$lgd_final, tolerance = 1e-9)
  expect_equal(nrow(openxlsx::read.xlsx(r$files$xlsx, sheet = "ELBE_Grid")), nrow(scr_elbe(m)$table))
  expect_true(any(grepl("lgd_final", readLines(r$files$sql))))
  mc <- openxlsx::read.xlsx(r$files$xlsx, sheet = "Model_Card")
  expect_true(all(c("n_pools", "severity_engine", "floor_asset_class") %in% mc$item))
})
