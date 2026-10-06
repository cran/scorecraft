# Regression tests of the lgd-ead triage: one block per confirmed bug.

# -- LGD ---------------------------------------------------------------------- #

test_that("a cure with drawings after default returns them in the artificial recovery", {
  d <- data.frame(default_id = "C", facility_id = "f1", default_date = as.Date("2024-01-15"), ead = 1000,
                  product = "p", status = "cured", close_date = as.Date("2024-04-15"), stringsAsFactors = FALSE)
  cf <- data.frame(default_id = "C", date = as.Date(c("2024-02-15", "2024-02-15")), amount = c(200, 100),
                   type = c("drawing", "recovery"), stringsAsFactors = FALSE)
  wo <- scr_workout(d, cf, config = lgd_cfg(lgd_discount_rate = 0, lgd_discount_add_on = 0))
  r <- wo$rds
  # outstanding at the cure date = 1000 + 200 - 100; no discount, no cost: zero loss
  expect_equal(r$recovery_artificial, 1100)
  expect_equal(r$lgd_raw, 0)
  # a direct cost stays in the loss of the cure
  cf2 <- rbind(cf, data.frame(default_id = "C", date = as.Date("2024-03-15"), amount = 30, type = "direct_cost"))
  wo2 <- scr_workout(d, cf2, config = lgd_cfg(lgd_discount_rate = 0, lgd_discount_add_on = 0))
  expect_equal(wo2$rds$lgd_raw, 0.03)
})

test_that("the vectorised cash-flow sums and extrapolation match the per-event definitions", {
  wo <- wo_demo()
  cf <- scr_workout(scr_demo_lgd, scr_demo_lgd_cashflows, rates = scr_demo_rates, config = lgd_cfg(), keep_rows = TRUE)$cashflows
  s <- cf[type == "recovery", list(pv = sum(pv)), by = default_id]
  r <- wo$rds[match(s$default_id, wo$rds$default_id)]
  expect_equal(r$pv_recovery[!r$is_cure], s$pv[!r$is_cure])
  pa <- r$recovery_artificial / (1 + r$discount_rate / 12)^r$months_in_default
  expect_equal(r$pv_recovery[r$is_cure], s$pv[r$is_cure] + pa[r$is_cure])
  e <- wo$extrapolation
  expect_true(nrow(e) > 0)
  expect_equal(e$expected_further, wo$rds$ead[match(e$default_id, wo$rds$default_id)] * pmax(0, e$rho_t_max - e$rho_tau))
  p <- wo$recovery_profile
  i <- 1L
  expect_equal(e$rho_tau[i], p[product == e$profile_source[i] & month == e$months_in_default[i], cum_recovery])
})

test_that("the merge map scans only facilities with several defaults and keeps singles as roots", {
  d <- data.table::data.table(default_id = c("a1", "a2", "b1", "c1", "c2"), facility_id = c("a", "a", "b", "c", "c"),
                              default_date = as.Date(c("2020-01-01", "2020-06-01", "2020-01-01", "2020-01-01", "2022-01-01")),
                              close_date = as.Date(c("2020-03-01", NA, NA, "2020-02-01", NA)))
  m <- .lgd_merge_map(d, 9L)
  m <- m[order(default_id)]
  expect_equal(m$root_id, c("a1", "a1", "b1", "c1", "c2"))
})

test_that("a rates table with a missing rate is refused", {
  rt <- data.frame(date = as.Date(c("2023-01-01", "2024-01-01")), rate = c(0.05, NA))
  expect_error(scr_workout(hand_defaults(), hand_cashflows(), rates = rt, config = lgd_cfg()), "missing date or rate")
})

test_that("downturn and floor do not modify the pools of the caller's object", {
  m <- lgd_demo()
  p0 <- data.table::copy(m$pools)
  m2 <- scr_lgd_downturn(m, method = "type3", add_on = 0.25, reason = "test")
  expect_equal(m$pools, p0)
  p2 <- data.table::copy(m2$pools)
  invisible(scr_lgd_floor(m2, asset_class = "retail_other", secured_share = 0.4))
  expect_equal(m2$pools, p2)
})

test_that("the in-default LGD at tau = 0 equals the downturn LGD even with recoveries in month 0", {
  m <- lgd_final()
  m$workout$recovery_profile <- data.table::copy(m$workout$recovery_profile)[, cum_recovery := cum_recovery + 0.05]
  e <- scr_elbe(m)
  expect_true(all(e$consistency$ok))
  t0 <- e$table[months_since_default == 0]
  expect_true(all(t0$recovered_share == 0))
})

test_that("scr_lgd_pools refuses a non-positive number of pools", {
  expect_error(scr_lgd_pools(lgd_demo(), n_pools = 0), "positive integer")
})

test_that("scr_lgd_validate refuses newdata with a missing realised LGD", {
  m <- lgd_demo()
  nd <- data.table::copy(wo_demo()$rds)
  nd$lgd_real[1] <- NA
  expect_error(scr_lgd_validate(m, newdata = nd), "finite `lgd_real`")
})

# -- EAD ---------------------------------------------------------------------- #

ead_panel <- function() {
  s <- data.table::data.table(fid = rep(c("a", "b"), each = 15),
                              d = rep(seq(as.Date("2023-01-01"), by = "month", length.out = 15), 2), lim = 1000, dr = 500)
  s[fid == "a" & d >= as.Date("2024-02-01"), dr := 900]
  s[, dd := as.Date(NA)]; s[fid == "a", dd := as.Date("2024-02-01")]
  s[]
}

test_that("a default-date column is read in the rows' own order, not the caller's", {
  s <- ead_panel()
  set.seed(1); sh <- s[sample(.N)]
  cfg <- scr_config(verbose = FALSE)
  e1 <- scr_ead_data(s, facility_id = "fid", date_col = "d", limit = "lim", drawn = "dr", default_date = "dd", config = cfg)
  e2 <- scr_ead_data(sh, facility_id = "fid", date_col = "d", limit = "lim", drawn = "dr", default_date = "dd", config = cfg)
  expect_equal(nrow(e2$rds), 1L)
  expect_equal(e2$rds$facility_id, "a")
  expect_equal(e2$rds$ccf, e1$rds$ccf)
  expect_equal(e2$rds$ccf, 0.8)
})

test_that("two snapshots of one facility in the same month are refused", {
  s <- ead_panel()
  s <- rbind(s, data.table::copy(s[1])[, d := d + 10])
  expect_error(scr_ead_data(s, facility_id = "fid", date_col = "d", limit = "lim", drawn = "dr", default_date = "dd",
                            config = scr_config(verbose = FALSE)), "duplicated \\(facility, month\\)")
})

test_that("a driver cannot overwrite an internal panel column", {
  s <- ead_panel()[, date := d]
  expect_error(scr_ead_data(s, facility_id = "fid", date_col = "d", limit = "lim", drawn = "dr", default_date = "dd",
                            drivers = "date", config = scr_config(verbose = FALSE)), "internal columns")
})

test_that("the predicted EAD per row works on data.table 1.14 (no vector fcase default)", {
  pools <- data.table::data.table(pool = c("P1", "LF"), ccf_applied = c(0.5, 1.1))
  pr <- .ead_predict_rows(c("P1", "LF", "P1"), c("ulf", "lf", "eadf"), c(400, 1000, 100), c(1000, 1000, 1000), pools, 0.2)
  expect_equal(pr$ead_model, c(400 + 0.5 * 600, 1.1 * 1000, 0.5 * 100))
  expect_equal(pr$ead_predicted, c(700, 1100, 100 + 0.2 * 900))
})

test_that("with the limit factor as the main measure the rows use the cells, not an LF pool", {
  cfg <- ead_cfg(ccf_measure = "lf")
  ed <- scr_ead_data(scr_demo_ead, facility_id = "facility_id", date_col = "ref_date", limit = "limit", drawn = "drawn",
                     defaulted = "defaulted", drivers = c("product", "months_on_book"), config = cfg)
  expect_true(all(ed$rds$measure == "lf"))
  m <- scr_ead(ed, drivers = c("utilisation_ref", "product"), config = cfg)
  expect_false("LF" %in% m$pools$pool)
  expect_false(any(m$rds$pool == "LF"))
  a <- scr_apply(m, scr_demo_ead[1:50, ], what = "pool")
  expect_false(any(a$pool == "LF"))
})

test_that("scr_ead_validate without a hold-out leaves the model's rows untouched", {
  m <- ead_model()
  m$rds <- data.table::copy(m$rds)[, sample := "train"]
  nm <- names(m$rds)
  v <- scr_ead_validate(m)
  expect_equal(v$source, "train (no hold-out)")
  expect_equal(names(m$rds), nm)
})

test_that("the EAD downturn refuses a period that ends before it starts", {
  expect_error(scr_ead_downturn(ead_model(), periods = data.frame(start = as.Date("2024-12-01"), end = as.Date("2024-01-01")),
                                reason = "x"), "start <= end")
})

test_that("the vectorised reference rows keep the per-event rules at the edges", {
  s <- ead_panel()
  cfg <- scr_config(verbose = FALSE)
  # a default dated before the first snapshot: the default row itself, horizon 0, excluded
  s[fid == "a", dd := as.Date("2022-06-01")]
  e <- scr_ead_data(s, facility_id = "fid", date_col = "d", limit = "lim", drawn = "dr", default_date = "dd",
                    config = cfg, keep_rows = TRUE)
  expect_equal(e$rows$horizon_months, 0L)
  expect_equal(e$rows$rule, "FAST_DEFAULT_EXCLUDED")
  # post-default drawings in the CCF: the maximum over the flagged run, NA drawn amounts ignored
  k <- ead_panel()[, flag := as.integer(fid == "a" & d >= as.Date("2024-02-01"))]
  k[fid == "a" & d == as.Date("2024-03-01"), dr := 1200]
  k[fid == "a" & d == as.Date("2024-02-01"), dr := NA]
  ek <- scr_ead_data(k, facility_id = "fid", date_col = "d", limit = "lim", drawn = "dr", defaulted = "flag",
                     config = scr_config(verbose = FALSE, post_default_drawings_in = "ccf"))
  expect_equal(ek$rds$ead_realised, 1200)
  expect_equal(ek$rds$ref_date, as.Date("2023-02-01"))
})
