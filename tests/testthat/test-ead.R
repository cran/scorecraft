# EAD/CCF module: the realised-CCF reference data set, the pools, the downturn
# and floor stack, application in R and SQL, validation and the workbook.

hand_args <- function(...) {
  utils::modifyList(list(facility_id = "facility_id", date_col = "ref_date", limit = "limit", drawn = "drawn",
                         defaulted = "defaulted", config = ead_cfg()), list(...))
}

test_that("the reference data set applies the identities, the measure rule and the named funnel rules", {
  h <- ead_hand()
  ed <- do.call(scr_ead_data, c(list(h), hand_args(keep_rows = TRUE)))
  expect_s3_class(ed, "scr_ead_data")
  r <- ed$rds; data.table::setkey(r, facility_id)
  # 10 events (A-H, L; I and J never flag themselves without the obligor); F and G excluded
  expect_equal(nrow(ed$rows), 9L)
  expect_setequal(r$facility_id, c("A", "B", "C", "D", "E", "H", "L"))
  expect_true(all(r$ref_date[r$facility_id != "E"] == as.Date("2023-01-01")))
  expect_true(all(r$default_date == as.Date("2024-01-01")))
  # A: CCF = (700 - 400) / 600 and the identity LF = u + CCF (1 - u)
  a <- r["A"]
  expect_equal(a$measure, "ulf"); expect_equal(a$ccf, 0.5); expect_equal(a$ccf_raw, 0.5)
  expect_equal(a$utilisation_ref + a$ccf * (1 - a$utilisation_ref), a$ead_realised / a$limit_ref)
  expect_equal(a$horizon_months, 12L); expect_false(a$fast_default); expect_equal(a$rule, "OK")
  # B: nothing undrawn -> limit factor, uncapped realised EAD above the limit
  b <- r["B"]
  expect_equal(b$measure, "lf"); expect_equal(b$ccf, 1.1); expect_equal(b$rule, "ZERO_UNDRAWN")
  # C: repaid before default -> raw -1.5 kept, floored at 0 and counted
  cc <- r["C"]
  expect_equal(cc$ccf_raw, -1.5); expect_equal(cc$ccf, 0); expect_equal(cc$rule, "NEGATIVE_CCF_FLOORED")
  # D: over the limit at the reference date -> LF = 1300 / 1000
  d <- r["D"]
  expect_equal(d$measure, "lf"); expect_equal(d$ccf, 1.3); expect_equal(d$rule, "OVER_LIMIT_AT_REF")
  # E: two months of history -> first snapshot, flagged fast default, kept
  e <- r["E"]
  expect_equal(e$ref_date, as.Date("2023-11-01")); expect_equal(e$horizon_months, 2L)
  expect_true(e$fast_default); expect_equal(e$rule, "FAST_DEFAULT"); expect_equal(e$ccf, 0.75)
  # H: utilisation 0.96 at or above u* -> LF by the measure rule
  hh <- r["H"]
  expect_equal(hh$measure, "lf"); expect_equal(hh$ccf, 1); expect_equal(hh$rule, "OK")
  # L: above one, kept uncapped and counted
  l <- r["L"]
  expect_equal(l$ccf, 7 / 6); expect_equal(l$rule, "CCF_ABOVE_ONE")
  # funnel
  f <- ed$funnel
  get <- function(rule) f$n[f$rule == rule]
  expect_equal(get("FAST_DEFAULT_EXCLUDED"), 1L); expect_equal(get("NOT_IN_SCOPE"), 1L)
  expect_equal(get("ZERO_UNDRAWN"), 1L); expect_equal(get("OVER_LIMIT_AT_REF"), 1L)
  expect_equal(get("NEGATIVE_CCF_FLOORED"), 1L); expect_equal(get("CCF_ABOVE_ONE"), 1L)
  expect_equal(get("FAST_DEFAULT"), 1L); expect_equal(get("OK"), 2L)
  expect_equal(f$action[f$rule == "ZERO_UNDRAWN"], "routed_to_lf")
  expect_equal(f$action[f$rule == "NOT_IN_SCOPE"], "excluded")
  expect_equal(sum(f$n), 9L); expect_equal(sum(f$share), 1)
  expect_equal(ed$rows$rule[ed$rows$facility_id == "G"], "NOT_IN_SCOPE")
  # averages: simple over the ULF rows (A, C, E, L), exposure-weighted
  expect_equal(ed$lra$simple, mean(c(0.5, 0, 0.75, 7 / 6)))
  # exposure-weighted: the floored per-row value weighted by the undrawn amount, so that the
  # floor applies to both averages (C contributes 0, not -300)
  expect_equal(ed$lra$exposure_weighted, (300 + 0 + 600 + 700) / (600 + 200 + 800 + 600))
  expect_equal(ed$lra$share_lf, 3 / 7)
  s <- ed$summary
  expect_equal(s[cohort == "ALL" & measure == "lf", n], 3L)
  expect_equal(s[cohort == "ALL" & measure == "ulf", ccf_simple], ed$lra$simple)
  expect_output(print(ed), "ZERO_UNDRAWN")
})

test_that("floor and cap of the realised value are explicit and logged; raw value always kept", {
  h <- ead_hand()
  ed <- do.call(scr_ead_data, c(list(h), hand_args(config = ead_cfg(ccf_floor_realised = NA, ccf_cap_realised = 1))))
  r <- ed$rds; data.table::setkey(r, facility_id)
  expect_equal(r["C"]$ccf, -1.5); expect_equal(r["C"]$rule, "OK")
  expect_equal(r["L"]$ccf, 1); expect_equal(r["L"]$ccf_raw, 7 / 6)
  expect_equal(ed$funnel$action[ed$funnel$rule == "CCF_ABOVE_ONE"], "capped")
  expect_false("NEGATIVE_CCF_FLOORED" %in% ed$funnel$rule)
  # a fixed measure: "lf" everywhere, "ulf" routes only the zero-undrawn rows
  ed_lf <- do.call(scr_ead_data, c(list(h), hand_args(config = ead_cfg(ccf_measure = "lf"))))
  expect_true(all(ed_lf$rds$measure == "lf"))
  expect_equal(ed_lf$rds[facility_id == "A", ccf], 0.7)
  ed_ulf <- do.call(scr_ead_data, c(list(h), hand_args(config = ead_cfg(ccf_measure = "ulf"))))
  expect_equal(ed_ulf$rds[facility_id == "H", measure], "ulf")
  expect_equal(ed_ulf$rds[facility_id == "B", measure], "lf")
  ed_eadf <- do.call(scr_ead_data, c(list(h), hand_args(config = ead_cfg(ccf_measure = "eadf"))))
  expect_equal(ed_eadf$rds[facility_id == "A", ccf], 700 / 400)
})

test_that("obligor-level defaults propagate, default dates can come as a table, and the horizons differ", {
  h <- ead_hand()
  ed <- do.call(scr_ead_data, c(list(h), hand_args(obligor_id = "obligor_id")))
  expect_true("I" %in% ed$rds$facility_id)
  i <- ed$rds[facility_id == "I"]
  expect_equal(i$ccf, 0); expect_equal(i$obligor_id, "oA"); expect_equal(i$default_date, as.Date("2024-01-01"))
  ed_f <- do.call(scr_ead_data, c(list(h), hand_args(obligor_id = "obligor_id", config = ead_cfg(default_level = "facility"))))
  expect_false("I" %in% ed_f$rds$facility_id)
  # default table
  dd <- data.frame(facility_id = "A", default_date = as.Date("2024-01-01"))
  ed_t <- do.call(scr_ead_data, c(list(h), hand_args(defaulted = NULL, default_date = dd)))
  expect_equal(nrow(ed_t$rds), 1L); expect_equal(ed_t$rds$ccf, 0.5)
  expect_error(do.call(scr_ead_data, c(list(h), hand_args(defaulted = NULL))), "default_date")
  expect_error(do.call(scr_ead_data, c(list(h), hand_args(limit = "nope"))), "not found")
  expect_error(do.call(scr_ead_data, c(list(h), hand_args(drivers = "ref_date"))), "reserved")
  # variable horizon: one row per month before default for A
  ed_v <- do.call(scr_ead_data, c(list(h[h$facility_id == "A", ]), hand_args(config = ead_cfg(ccf_horizon = "variable"))))
  expect_equal(nrow(ed_v$rds), 12L)
  expect_equal(sort(ed_v$rds$horizon_months), 1:12)
  expect_true(all(ed_v$rds$ccf == 0.5))
  # cohort horizon: the calendar cohort start (2023-01) is the reference of A
  ed_c <- do.call(scr_ead_data, c(list(h[h$facility_id == "A", ]), hand_args(config = ead_cfg(ccf_horizon = "cohort"))))
  expect_equal(ed_c$rds$ref_date, as.Date("2023-01-01"))
  # post-default drawings booked in the CCF: the maximum drawn amount over the event
  k <- ead_hand()[ead_hand()$facility_id == "A", ]
  k$drawn <- c(rep(400, 9), 700, 900, 800, 800); k$defaulted <- as.integer(seq_len(13) >= 10)
  ed_l <- do.call(scr_ead_data, c(list(k), hand_args()))
  ed_k <- do.call(scr_ead_data, c(list(k), hand_args(config = ead_cfg(post_default_drawings_in = "ccf"))))
  expect_equal(ed_l$rds$ead_realised, 700); expect_equal(ed_k$rds$ead_realised, 900)
  expect_equal(ed_l$rds$horizon_months, 9L); expect_true(ed_l$rds$fast_default)
})

test_that("the demo reference data set is pinned", {
  ed <- ead_demo()
  expect_equal(nrow(ed$rds), 197L)
  expect_equal(ed$meta$n_events, 197L)
  f <- ed$funnel
  expect_equal(f$n[match(c("FAST_DEFAULT", "OVER_LIMIT_AT_REF", "NEGATIVE_CCF_FLOORED", "CCF_ABOVE_ONE", "OK"), f$rule)],
               c(25L, 7L, 45L, 14L, 106L))
  expect_equal(ed$lra$n, 189L)
  expect_equal(ed$lra$simple, 0.325, tolerance = 2e-3)
  expect_equal(ed$lra$share_lf, 8 / 197)
  expect_true(all(ed$rds$ead_realised >= 0))
  expect_true(all(ed$rds[measure == "ulf", undrawn_ref] > 0))
  expect_true(all(ed$rds[measure == "ulf" & rule != "NEGATIVE_CCF_FLOORED", ccf == ccf_raw]))
  expect_true(all(c("product", "months_on_book", "dpd", "utilisation_ref") %in% names(ed$rds)))
})

test_that("pools are estimated on train by reference date, admitted by named rules and pinned", {
  m <- ead_model()
  expect_s3_class(m, "scr_ead")
  expect_equal(m$split$method, "out-of-time"); expect_equal(m$split$cutoff, "2023-11-01")
  expect_equal(m$split$n_train, 126L); expect_equal(m$split$n_holdout, 71L)
  expect_setequal(m$drivers$feature, c("utilisation_ref", "product", "months_on_book"))
  expect_true(all(m$drivers$reason %in% c("OK", "TOO_FEW_DEFAULTS", "NO_SEPARATION", "NOT_MONOTONIC", "UNSTABLE_HOLDOUT") |
                    grepl(";", m$drivers$reason)))
  expect_equal(m$survivors, "utilisation_ref")
  expect_equal(m$drivers[feature == "utilisation_ref", direction], "decreasing")
  p <- m$pools
  expect_equal(p$pool, c("P1", "P2", "LF"))
  expect_equal(p$n, c(98L, 23L, 5L))
  expect_equal(p$lra, c(0.2856519319, 0.4962956955, 0.8884166667), tolerance = 1e-9)
  # pools monotone in the long-run average; MoC, final and applied by construction
  expect_true(all(diff(p$lra[p$pool != "LF"]) > 0))
  expect_equal(p$moc_est, stats::qnorm(0.95) * p$se)
  expect_equal(p$ccf_final, pmax(p$lra, p$ccf_dt) + p$moc_est)
  expect_equal(p$ccf_floor[p$pool != "LF"], rep(0.5 * 0.40, 2)); expect_true(is.na(p$ccf_floor[p$pool == "LF"]))
  expect_equal(p$ccf_applied, c(0.3459067005, 0.6264385946, 1.0294951378), tolerance = 1e-9)
  expect_true(all(p$ccf_applied[p$pool != "LF"] >= 0.2))
  expect_equal(p$ccf_dt, p$lra); expect_true(all(p$downturn == "none"))
  # every training row in the main measure lands in a cell of the admitted driver
  expect_equal(sum(m$cells$n), 121L)
  expect_true(all(m$rds[measure == "lf", pool] == "LF"))
  expect_equal(sort(unique(m$rds$pool)), c("LF", "P1", "P2"))
  mt <- m$metrics
  expect_equal(mt$sample, c("train", "holdout"))
  expect_equal(mt[sample == "train", gauc], 0.5582217, tolerance = 1e-6)
  expect_equal(mt[sample == "holdout", rmse], 0.4034747, tolerance = 1e-6)
  expect_true(all(mt$gauc_lo <= mt$gauc & mt$gauc <= mt$gauc_hi))
  expect_true(all(is.finite(mt$cear)))
  expect_equal(mt$gauc, (mt$somers_d + 1) / 2)
  expect_equal(m$model_card$n_pools, 2L); expect_true(m$model_card$lf_pool)
  expect_equal(m$model_card$ccf_floor, 0.2)
  expect_true(any(m$ledger$step == "pools"))
  expect_output(print(m), "applied")
  expect_error(scr_ead(ead_demo(), drivers = "nope", config = ead_cfg()), "not in the reference")
  expect_error(scr_ead(ead_demo(), drivers = "product", config = ead_cfg(ccf_min_defaults = 500L)), "fewer than twice")
  # a modified parameter table is recorded
  pr <- scr_irb_params("bcb"); pr$ccf_floor_fraction <- 0.6
  m6 <- scr_ead(ead_demo(), drivers = "utilisation_ref", config = ead_cfg(), params = pr)
  expect_equal(m6$meta$floor, 0.24); expect_true(m6$params$modified)
  expect_true(any(grepl("params_modified", m6$ledger$detail)))
  # no driver admitted: a single pool with the long-run average
  m1 <- scr_ead(ead_demo(), drivers = "dpd", config = ead_cfg())
  expect_equal(m1$pools$pool, c("P1", "LF")); expect_equal(m1$survivors, character())
  expect_equal(m1$pools$lra[1], mean(m1$rds[sample == "train" & measure == "ulf", ccf]))
  expect_true(any(m1$ledger$reason == "NO_DRIVER_ADMITTED"))
})

test_that("Somers' D, gAUC and the cumulative EAD accuracy ratio are pinned on hand data", {
  expect_equal(.somers_d(c(1, 1, 2, 2), c(0.1, 0.2, 0.3, 0.4)), 4 / 6)
  expect_equal(.somers_d(c(1, 2), c(0.5, 0.2)), -1)
  expect_equal(.somers_d(c(1, 1, 1), c(0.5, 0.2, 0.1)), 0)
  expect_true(is.na(.somers_d(c(1, 2), c(0.5, 0.5))))
  g <- .ead_gauc(c(1, 1, 2, 2), c(0.1, 0.2, 0.3, 0.4), n_boot = 10, seed = 1)
  expect_equal(g$gauc, 5 / 6); expect_true(g$lo <= g$gauc && g$gauc <= g$hi)
  expect_equal(.ead_cear(c(4, 3, 2, 1), c(40, 30, 20, 10)), 1)
  expect_equal(.ead_cear(c(1, 1, 1, 1), c(40, 30, 20, 10)), 0)
  expect_true(is.na(.ead_cear(c(1, 2), c(0, 0))))
})

test_that("scr_apply predicts an EAD never below the drawn amount nor the standardised floor, with the LF branch", {
  m <- ead_model()
  a <- scr_apply(m, scr_demo_ead)
  expect_equal(nrow(a), nrow(scr_demo_ead))
  expect_setequal(names(a), c("pool", "measure", "utilisation", "undrawn", "ccf_applied", "ead_model", "ead_floor",
                              "ead_predicted", "ead_floor_binding"))
  undrawn <- pmax(scr_demo_ead$limit - scr_demo_ead$drawn, 0)
  expect_true(all(a$ead_predicted >= scr_demo_ead$drawn))
  expect_true(all(a$ead_predicted >= scr_demo_ead$drawn + 0.5 * 0.40 * undrawn - 1e-9))
  expect_equal(a$undrawn, undrawn)
  u <- scr_demo_ead$drawn / scr_demo_ead$limit
  expect_true(all(a$pool[u >= 0.95] == "LF")); expect_true(all(a$measure[u >= 0.95] == "lf"))
  expect_true(all(a$pool[u < 0.95] != "LF"))
  lf <- m$pools[pool == "LF", ccf_applied]
  i <- which(u >= 0.95)
  expect_equal(a$ead_predicted[i], pmax(scr_demo_ead$drawn[i], lf * scr_demo_ead$limit[i], scr_demo_ead$drawn[i] + 0.2 * undrawn[i]))
  j <- which(u < 0.95)
  expect_equal(a$ead_predicted[j], pmax(scr_demo_ead$drawn[j] + a$ccf_applied[j] * undrawn[j], scr_demo_ead$drawn[j] + 0.2 * undrawn[j]))
  expect_equal(a$ead_floor_binding, a$ead_floor >= a$ead_model & a$undrawn > 0)
  # hand rows: the LF branch, the floor branch and an over-limit row
  hand <- data.frame(limit = c(1000, 1000, 1000, 0), drawn = c(400, 960, 1300, 0))
  b <- scr_apply(m, hand)
  expect_equal(b$pool[2:4], c("LF", "LF", "LF"))
  expect_equal(b$ead_predicted[2], max(960, lf * 1000, 960 + 0.2 * 40))
  expect_equal(b$ead_predicted[3], 1300)
  expect_equal(b$ead_predicted[4], 0)
  expect_equal(b$ead_predicted[1], 400 + b$ccf_applied[1] * 600)
  expect_true(b$ccf_applied[1] %in% m$pools$ccf_applied)
  expect_error(scr_apply(m, data.frame(limit = 1)), "lacks")
  # with a binding floor the applied CCF is the floor
  pr <- scr_irb_params("bcb"); mf <- scr_ead(ead_demo(), drivers = "utilisation_ref", config = ead_cfg(ccf_sa_ccf = 1), params = pr)
  expect_equal(mf$pools[pool != "LF", ccf_applied], pmax(mf$pools[pool != "LF", ccf_final], 0.5))
  expect_equal(mf$pools[pool == "P1", ccf_applied], 0.5); expect_true(mf$pools[pool == "P1", floor_binding])
  af <- scr_apply(mf, data.frame(limit = 1000, drawn = 900))      # high utilisation below u*: the low-CCF pool
  expect_equal(af$pool, "P1"); expect_equal(af$ead_predicted, 900 + 0.5 * 100); expect_true(af$ead_floor_binding)
})

test_that("the downturn stack is type1 or type3 with a mandatory reason and a ledger row", {
  m <- ead_model()
  per <- data.frame(start = as.Date("2024-01-01"), end = as.Date("2024-12-01"))
  expect_error(scr_ead_downturn(m, per), "reason")
  expect_error(scr_ead_downturn(m, method = "type1", reason = "x"), "periods")
  d1 <- scr_ead_downturn(m, per, reason = "2024 as the stress year of the demo panel")
  expect_equal(d1$downturn$method, "type1")
  t <- d1$downturn$table
  # training rows only: the hold-out stays independent
  in_per <- m$rds[default_date >= per$start & default_date <= per$end]
  expect_true(any(in_per$sample == "holdout"))
  expect_equal(t$n_downturn, in_per[sample == "train", .N, by = pool][match(t$pool, pool), N])
  obs <- in_per[sample == "train", mean(ccf), by = pool]
  expect_equal(t$dt_observed, obs$V1[match(t$pool, obs$pool)])
  expect_equal(d1$pools$ccf_dt, pmax(d1$pools$lra, t$dt_observed))
  expect_equal(d1$pools$ccf_final, pmax(d1$pools$lra, d1$pools$ccf_dt) + d1$pools$moc_est)
  expect_true(all(d1$pools$ccf_applied >= m$pools$ccf_applied))
  expect_true(all(d1$pools$downturn == "type1"))
  expect_equal(d1$ledger[step == "downturn", reason], "2024 as the stress year of the demo panel")
  expect_equal(d1$model_card$downturn, "type1")
  d3 <- scr_ead_downturn(m, per, method = "type3", add_on = 0.1, reason = "add-on")
  expect_equal(d3$pools$ccf_dt, m$pools$lra + 0.1)
  expect_equal(d3$pools$ccf_applied[1], max(m$pools$lra[1] + 0.1 + m$pools$moc_est[1], 0.2))
  d0 <- scr_ead_downturn(d3, method = "none", reason = "reset")
  expect_equal(d0$pools$ccf_dt, m$pools$lra); expect_equal(d0$pools$ccf_applied, m$pools$ccf_applied)
  expect_equal(nrow(d0$ledger), nrow(m$ledger) + 2L)
  expect_output(print(d1), "type1")
})

test_that("the validation battery reports calibration, discrimination, back-test and stability with lights", {
  m <- ead_model()
  v <- scr_ead_validate(m)
  expect_s3_class(v, "scr_ead_validation")
  expect_equal(v$source, "holdout"); expect_equal(v$n, 71L)
  c <- v$calibration
  expect_equal(c$pool, c("P1", "P2", "LF", "TOTAL"))
  expect_equal(c$n, c(50L, 18L, 3L, 71L))
  tot <- c[pool == "TOTAL"]
  expect_equal(tot$t, -1.93344252, tolerance = 1e-6)
  expect_equal(tot$adequacy, 0.79054114, tolerance = 1e-6)
  expect_equal(tot$adequacy, sum(tot$ead_realised) / sum(tot$ead_predicted))
  expect_equal(tot$ead_realised, sum(m$rds[sample == "holdout", ead_realised]))
  expect_true(all(c$light_p %in% c("green", "amber", "red")))
  expect_true(all(c$light_adequacy %in% c("green", "amber", "red")))
  expect_equal(c$p, stats::pt(c$t, df = c$n_main - 1L, lower.tail = FALSE))
  expect_true(all(c$predicted[c$pool != "TOTAL"] == m$pools$ccf_applied))
  d <- v$discrimination
  expect_true(d$gauc_lo <= d$gauc && d$gauc <= d$gauc_hi)
  expect_equal(d$gauc_dev, m$metrics[sample == "train", gauc])
  expect_true(is.finite(d$cear))
  expect_equal(sort(unique(v$backtest$cohort)), c("2023", "2024"))
  expect_equal(sum(v$backtest$n), 71L)
  expect_equal(v$stability$item, c("pool", "utilisation_ref"))
  expect_true(all(v$stability$flag_fixed %in% c("stable", "moderate", "shift")))
  expect_equal(v$summary$test, c("calibration_t_total", "ead_adequacy_total", "gauc_vs_development", "pool_psi"))
  expect_true(all(grepl("convention|Yurdakul|bootstrap", v$summary$convention)))
  expect_output(print(v), "TOTAL")
  # the lights follow the thresholds: a red one when realised exceeds predicted strongly
  vv <- scr_ead_validate(m, lights = c(0.99, 0.995))
  expect_true(all(vv$calibration$light_p == "red"))
  expect_error(scr_ead_validate(m, lights = c(0.5, 0.1)), "increasing")
  # newdata as a reference data set built elsewhere
  v2 <- scr_ead_validate(m, newdata = ead_demo())
  expect_equal(v2$source, "newdata"); expect_equal(v2$n, 197L)
  expect_equal(v2$calibration[pool == "TOTAL", n], 197L)
  expect_error(scr_ead_validate(m, newdata = data.frame(x = 1)), "lacks")
})

test_that("R and SQL agree on pool, applied CCF and predicted EAD (DuckDB), in every dialect", {
  m <- ead_model()
  dialects <- c("ansi", "databricks", "spark", "hive", "mysql", "mariadb", "sqlserver", "bigquery", "postgres",
                "oracle", "snowflake", "redshift", "duckdb", "sqlite")
  for (dl in dialects) {
    txt <- paste(scr_sql(m, table = "t", dialect = dl), collapse = "\n")
    expect_match(txt, "WITH base_ead AS", fixed = TRUE)
    expect_match(txt, "AS ead_predicted", fixed = TRUE)
    expect_match(txt, "utilisation_idx", fixed = TRUE)
    if (dl == "sqlite") expect_match(txt, "MAX(drawn_amt", fixed = TRUE)
    else if (dl == "sqlserver") expect_match(txt, "VALUES", fixed = TRUE)
    else expect_match(txt, "GREATEST(drawn_amt", fixed = TRUE)
  }
  expect_error(scr_sql(m, dialect = "cobol"), "dialect")
  tmp <- tempfile(fileext = ".sql")
  expect_invisible(scr_sql(m, table = "t", file = tmp)); expect_true(file.exists(tmp))
  skip_if_not_installed("duckdb"); skip_if_not_installed("DBI")
  con <- DBI::dbConnect(duckdb::duckdb(), config = list(threads = "1"))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  DBI::dbWriteTable(con, "snap", scr_demo_ead)
  got <- DBI::dbGetQuery(con, paste(scr_sql(m, table = "snap", dialect = "duckdb"), collapse = "\n"))
  exp <- scr_apply(m, scr_demo_ead)
  expect_equal(nrow(got), nrow(exp))
  expect_identical(got$pool, exp$pool)
  expect_equal(got$ccf_applied, exp$ccf_applied, tolerance = 1e-9)
  expect_equal(got$ead_predicted, exp$ead_predicted, tolerance = 1e-9)
  expect_equal(got$undrawn, exp$undrawn)
  # a model with a downturn and the single-pool model also reproduce
  d1 <- scr_ead_downturn(m, data.frame(start = as.Date("2024-01-01"), end = as.Date("2024-12-01")), reason = "test")
  got2 <- DBI::dbGetQuery(con, paste(scr_sql(d1, table = "snap", dialect = "duckdb"), collapse = "\n"))
  expect_equal(got2$ead_predicted, scr_apply(d1, scr_demo_ead)$ead_predicted, tolerance = 1e-9)
  m1 <- scr_ead(ead_demo(), drivers = "dpd", config = ead_cfg())
  got1 <- DBI::dbGetQuery(con, paste(scr_sql(m1, table = "snap", dialect = "duckdb"), collapse = "\n"))
  expect_identical(got1$pool, scr_apply(m1, scr_demo_ead)$pool)
})

test_that("R and SQL agree on pool, applied CCF and predicted EAD (SQLite)", {
  skip_if_not_installed("RSQLite"); skip_if_not_installed("DBI")
  m <- ead_model()
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  d <- scr_demo_ead; d$ref_date <- as.character(d$ref_date)
  DBI::dbWriteTable(con, "snap", d)
  got <- DBI::dbGetQuery(con, paste(scr_sql(m, table = "snap", dialect = "sqlite"), collapse = "\n"))
  exp <- scr_apply(m, scr_demo_ead)
  expect_identical(got$pool, exp$pool)
  expect_equal(got$ccf_applied, exp$ccf_applied, tolerance = 1e-9)
  expect_equal(got$ead_predicted, exp$ead_predicted, tolerance = 1e-9)
})

test_that("the pools are identical under the serial and PSOCK backends", {
  skip_if_not_installed("withr")
  a <- withr::with_options(list(scorecraft.parallel = "serial"),
    scr_ead(ead_demo(), drivers = c("utilisation_ref", "product", "months_on_book"), config = ead_cfg(nthread = 2L)))
  b <- withr::with_options(list(scorecraft.parallel = "psock"),
    scr_ead(ead_demo(), drivers = c("utilisation_ref", "product", "months_on_book"), config = ead_cfg(nthread = 2L)))
  expect_equal(a$pools, b$pools); expect_equal(a$metrics, b$metrics); expect_equal(a$cells, b$cells)
  expect_equal(a$pools, ead_model()$pools)
})

test_that("the workbook round trip keeps every sheet and the SQL file", {
  skip_if_not_installed("openxlsx")
  m <- ead_model()
  out <- file.path(tempdir(), "scr-ead-export"); unlink(out, recursive = TRUE)
  x <- scr_export(m, out, stamp = FALSE)
  expect_setequal(names(x$files), c("xlsx", "sql"))
  expect_true(all(file.exists(unlist(x$files))))
  expect_equal(basename(x$files$xlsx), "ead_ccf.xlsx")
  sheets <- openxlsx::getSheetNames(x$files$xlsx)
  expect_true(all(c("RDS_Funnel", "RDS_Summary", "Driver_Bins", "Pools", "Downturn_MoC", "Holdout", "Calibration",
                    "Discrimination", "Backtest_Cohort", "Stability", "Model_Card", "Decision_Ledger") %in% sheets))
  pools <- openxlsx::read.xlsx(x$files$xlsx, sheet = "Pools")
  expect_equal(nrow(pools), nrow(m$pools)); expect_equal(pools$ccf_applied, m$pools$ccf_applied, tolerance = 1e-9)
  fun <- openxlsx::read.xlsx(x$files$xlsx, sheet = "RDS_Funnel")
  expect_equal(fun$n, m$funnel$n)
  mc <- openxlsx::read.xlsx(x$files$xlsx, sheet = "Model_Card")
  expect_true(all(c("horizon", "ccf_floor", "drivers_admitted") %in% mc$item))
  expect_true(any(grepl("ead_predicted", readLines(x$files$sql))))
  v <- scr_ead_validate(m)
  y <- scr_export(scr_ead_downturn(m, data.frame(start = as.Date("2024-01-01"), end = as.Date("2024-12-01")), reason = "t"),
                  out, stamp = FALSE, validation = v, tag = "dt")
  expect_equal(basename(y$files$xlsx), "ead_dt.xlsx")
  dt <- openxlsx::read.xlsx(y$files$xlsx, sheet = "Downturn_MoC")
  expect_true(all(c("dt_observed", "method", "periods") %in% names(dt)))
})
