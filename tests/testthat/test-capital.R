# Expected loss, IRB risk weights, standardised comparison, portfolio capital,
# the Vasicek helpers, SQL, export and the accounting ECL.

test_that("expected loss is the primitive with the ELBE branch on defaulted rows", {
  expect_equal(scr_el(c(0.01, 0.02), 0.45, c(1000, 2000)), c(4.5, 18))
  expect_equal(scr_el(0.02, 0.45, 1000, defaulted = TRUE, elbe = 0.6), 600)
  expect_equal(scr_el(0.02, 0.45, 1000, defaulted = TRUE), 450)          # no ELBE: PD one, LGD
  expect_equal(scr_el(0.02, 0.45, 1000, defaulted = c(0, 1), elbe = c(NA, 0.3)), c(9, 300))
})

test_that("the risk-weight function reproduces the worked examples", {
  p <- scr_irb_params("bcb")
  r <- scr_irb_rw(0.01, 0.45, m = 2.5, asset_class = "corporate", params = p)
  expect_equal(round(r$r, 6), 0.192784)
  expect_equal(round(r$b, 6), 0.137486)
  expect_equal(r$ma, 1 / (1 - 1.5 * r$b), tolerance = 1e-12)   # M = 2.5: only the denominator remains
  expect_equal(round(r$k, 6), 0.073853)
  expect_equal(round(r$rw, 4), 0.9232)
  expect_equal(r$rwa, r$rw, tolerance = 1e-12)
  p2 <- p; p2$scaling_factor <- 1.06
  expect_equal(round(scr_irb_rw(0.01, 0.45, m = 2.5, asset_class = "corporate", params = p2)$rw, 4), 0.9786)
  expect_equal(round(scr_irb_rw(0.01, 0.20, asset_class = "retail_mortgage", params = p)$rw, 4), 0.2507)
  expect_equal(round(scr_irb_rw(0.01, 0.20, asset_class = "retail_mortgage", params = p)$k, 6), 0.020053)
  expect_equal(round(scr_irb_rw(0.02, 0.80, asset_class = "qrre_revolver", params = p)$rw, 4), 0.5142)
  expect_equal(round(scr_irb_rw(0.02, 0.80, asset_class = "qrre_revolver", params = p)$k, 6), 0.041135)
  o <- scr_irb_rw(0.02, 0.50, asset_class = "retail_other", params = p)
  expect_equal(round(o$r, 6), 0.094556)
  expect_equal(round(o$k, 6), 0.051544)
  expect_equal(round(o$rw, 4), 0.6443)
  # retail: no maturity adjustment, b not reported
  expect_true(is.na(o$b)); expect_equal(o$ma, 1)
  # K equals LGD times the stressed PD minus PD at 99.9 % (retail, no MA)
  expect_equal(o$k, 0.5 * (scr_pd_stress(0.02, o$r, 0.999) - 0.02), tolerance = 1e-12)
  # vectorised over classes with an EAD
  v <- scr_irb_rw(c(0.01, 0.02), c(0.20, 0.80), ead = c(100, 200), asset_class = c("retail_mortgage", "qrre_revolver"), params = p)
  expect_equal(nrow(v), 2L); expect_equal(v$rwa, v$rw * c(100, 200))
  expect_error(scr_irb_rw(0.01, 0.45, asset_class = "boats"), "asset_class")
  expect_error(scr_irb_rw(1.5, 0.45, asset_class = "corporate"), "pd")
  expect_error(scr_irb_rw(0.01, 0.45, asset_class = "corporate", params = list()), "scr_irb_params")
  expect_error(scr_irb_rw(0.01, 0.45, asset_class = "corporate", apply_floors = "ccf"), "apply_floors")
})

test_that("floors, maturity, defaulted rows, SME and FI adjustments behave as specified", {
  p <- scr_irb_params("bcb")
  f <- scr_irb_rw(1e-4, 0.5, asset_class = "retail_other", params = p)
  expect_equal(f$pd_used, 5e-4)
  fh <- attr(f, "floors_hit")
  expect_equal(fh$n[fh$floor == "pd_floor"], 1L)
  expect_equal(scr_irb_rw(1e-4, 0.5, asset_class = "qrre_revolver", params = p)$pd_used, 1e-3)
  expect_equal(scr_irb_rw(1e-4, 0.5, asset_class = "retail_other", params = p, apply_floors = FALSE)$pd_used, 1e-4)
  expect_equal(scr_irb_rw(1e-4, 0.5, asset_class = "retail_other", params = p, apply_floors = "lgd")$pd_used, 1e-4)
  # LGD floors: unsecured column; mortgages fall back to real estate; F-IRB has none
  l <- scr_irb_rw(0.01, 0.10, asset_class = c("retail_other", "qrre_revolver", "corporate", "retail_mortgage"), params = p)
  expect_equal(l$lgd_used, c(0.30, 0.50, 0.25, 0.10))
  expect_equal(attr(l, "floors_hit")$n[2], 3L)
  expect_equal(scr_irb_rw(0.01, 0.02, asset_class = "retail_mortgage", params = p)$lgd_used, 0.05)
  expect_equal(scr_irb_rw(0.01, 0.10, asset_class = "corporate", params = p, approach = "firb")$lgd_used, 0.10)
  expect_equal(scr_irb_rw(0.01, 0.05, asset_class = "corporate", params = p, collateral = "receivables")$lgd_used, 0.10)
  expect_equal(scr_irb_rw(0.01, 0.05, asset_class = "corporate", params = p, collateral = "receivables", secured_share = 0.5)$lgd_used, 0.175)
  expect_error(scr_irb_rw(0.01, 0.05, asset_class = "corporate", collateral = "gold"), "collateral")
  # maturity: default when missing, clipped to the range, fixed under F-IRB
  expect_equal(scr_irb_rw(0.01, 0.45, asset_class = "corporate", params = p)$m, 2.5)
  mm <- scr_irb_rw(0.01, 0.45, m = c(0.2, 3, 9), asset_class = "corporate", params = p)
  expect_equal(mm$m, c(1, 3, 5))
  expect_equal(attr(mm, "floors_hit")$n[3:4], c(1L, 1L))
  expect_true(mm$ma[1] < mm$ma[2] && mm$ma[2] < mm$ma[3])
  expect_equal(scr_irb_rw(0.01, 0.45, m = 4, asset_class = "corporate", params = p, approach = "firb")$m, 2.5)
  # defaulted rows
  d1 <- scr_irb_rw(0.05, 0.6, asset_class = "retail_other", defaulted = TRUE, elbe = 0.7, params = p)
  expect_equal(d1$k, 0); expect_equal(d1$pd_used, 1); expect_true(is.na(d1$r))
  expect_equal(scr_irb_rw(0.05, 0.6, asset_class = "retail_other", defaulted = TRUE, elbe = 0.5, params = p)$k, 0.1)
  expect_equal(scr_irb_rw(0.05, 0.6, asset_class = "retail_other", defaulted = TRUE, params = p)$k, 0)
  expect_equal(scr_irb_rw(0.05, 0.6, asset_class = "corporate", defaulted = 1, elbe = 0.5, params = p, approach = "firb")$k, 0)
  # SME adjustment lowers the correlation, more so for smaller firms; FI multiplier raises it
  rc <- scr_irb_rw(0.01, 0.45, asset_class = "corporate", params = p)$r
  rs <- scr_irb_rw(0.01, 0.45, asset_class = "corporate_sme", sales = c(15, 100, 300, 1000, NA), params = p)$r
  expect_true(all(rs[1:2] < rc))
  expect_equal(rs[1], rc - 0.04, tolerance = 1e-12)
  expect_equal(rs[3], rc, tolerance = 1e-12); expect_equal(rs[4], rc, tolerance = 1e-12)
  expect_equal(rs[5], rc, tolerance = 1e-12)   # missing sales: no firm-size relief (CRE31.9)
  expect_true(rs[1] < rs[2])
  rf <- scr_irb_rw(0.01, 0.45, asset_class = "bank", fi = TRUE, params = p)$r
  expect_equal(rf, 1.25 * rc, tolerance = 1e-12)
  expect_gt(scr_irb_rw(0.01, 0.45, asset_class = "bank", fi = TRUE, params = p)$rw, rc * 0 + 0.9232)
  # fi ignored on retail; hvcre uses its own ceiling
  expect_equal(scr_irb_rw(0.01, 0.45, asset_class = "retail_other", fi = TRUE, params = p)$r,
               scr_irb_rw(0.01, 0.45, asset_class = "retail_other", params = p)$r)
  expect_gt(scr_irb_rw(0.01, 0.45, asset_class = "hvcre", params = p)$r, rc)
  # basel3_final: SME bounds in EUR m
  rs2 <- scr_irb_rw(0.01, 0.45, asset_class = "corporate_sme", sales = c(5, 50), params = scr_irb_params("basel3_final"))$r
  expect_equal(rs2[1], rc - 0.04, tolerance = 1e-12); expect_equal(rs2[2], rc, tolerance = 1e-12)
})

test_that("the Vasicek helpers are consistent and bounded", {
  expect_equal(.vasicek_pit(0.02, 0, 0.15), stats::pnorm(stats::qnorm(0.02) / sqrt(0.85)))
  expect_equal(scr_pd_stress(0.02, 0.15, 0.5), .vasicek_pit(0.02, 0, 0.15))
  expect_true(scr_pd_stress(0.02, 0.15, 0.99) > scr_pd_stress(0.02, 0.15, 0.95))
  expect_true(scr_pd_stress(0.02, 0.15, 0.95) > 0.02)
  expect_equal(scr_pd_stress(c(0, 1), 0.15, 0.99), c(0, 1))
  expect_equal(.vasicek_pit(0.05, 1.5, 0.1), stats::pnorm((stats::qnorm(0.05) - sqrt(0.1) * 1.5) / sqrt(0.9)))
  expect_error(scr_pd_stress(0.02, 0.15, 1), "q")
  expect_error(.vasicek_pit(0.02, 0, 1), "rho")
})

test_that("standardised risk weights are looked up by class, band, rating and default status", {
  expect_equal(scr_sa_rw("retail_other"), 0.75)
  expect_equal(scr_sa_rw(c("qrre_revolver", "qrre_transactor")), c(0.75, 0.45))
  expect_equal(scr_sa_rw("retail_other", transactor = TRUE), 0.45)
  expect_equal(scr_sa_rw("retail_mortgage", ltv = c(0.3, 0.5, 0.55, 0.7, 0.85, 0.95, 1.2, NA)),
               c(0.20, 0.20, 0.25, 0.30, 0.40, 0.50, 0.70, 0.70))
  expect_equal(scr_sa_rw("corporate", rating = c("AAA", "AA-", "A+", "BBB", "BB", "B", "CCC", NA, "", "IG")),
               c(0.20, 0.20, 0.50, 0.75, 1.00, 1.50, 1.50, 1.00, 1.00, 0.65))
  expect_equal(scr_sa_rw(c("corporate_sme", "corporate"), sme = c(FALSE, TRUE)), c(0.85, 0.85))
  expect_equal(scr_sa_rw("corporate_sme", rating = "A"), 0.50)
  expect_equal(scr_sa_rw(c("bank", "sovereign")), c(1, 1))
  expect_equal(scr_sa_rw("retail_other", defaulted = TRUE, provision_ratio = c(0.1, 0.3, NA)), c(1.5, 1.0, 1.5))
  expect_equal(scr_sa_rw("retail_mortgage", ltv = 0.3, defaulted = TRUE), 1.0)
  expect_error(scr_sa_rw("boats"), "asset_class")
})

test_that("scr_capital aggregates the book, reconciles EL with provisions and measures the floors", {
  cap <- cap_demo()
  expect_s3_class(cap, "scr_capital")
  d <- scr_demo_portfolio
  ex <- cap$exposures
  expect_equal(nrow(ex), nrow(d))
  expect_equal(ex$id, d$id)
  t <- cap$totals
  expect_equal(t$n, 5000L)
  expect_equal(t$ead, sum(d$ead))
  expect_equal(t$rwa_irb, sum(ex$rwa))
  expect_equal(t$el, sum(ex$el))
  expect_equal(t$rwa_sa, sum(ex$rw_sa * ex$ead))
  expect_equal(sum(cap$segments$rwa_irb), t$rwa_irb)
  expect_equal(sum(cap$segments$el), t$el)
  expect_equal(sum(cap$segments$n), 5000L)
  expect_equal(t$rwa_reported, max(t$rwa_irb, 0.725 * t$rwa_sa))
  expect_equal(t$capital, 0.08 * t$rwa_reported)
  expect_equal(t$density, t$rwa_reported / t$ead)
  # exposure rows reproduce scr_irb_rw() directly
  i <- which(d$asset_class == "corporate")[1:20]
  r <- scr_irb_rw(d$pd[i], d$lgd[i], d$ead[i], m = d$m[i], asset_class = "corporate", defaulted = d$defaulted[i], elbe = d$elbe[i])
  expect_equal(ex$rwa[i], r$rwa)
  expect_equal(ex$k[i], r$k)
  # EL: floored PD on performing rows, ELBE on defaulted rows
  j <- which(d$defaulted == 1L)
  expect_equal(ex$el[j], d$elbe[j] * d$ead[j])
  expect_equal(ex$pd_used[j], rep(1, length(j)))
  expect_equal(ex$el[-j], ex$pd_used[-j] * ex$lgd_used[-j] * d$ead[-j])
  # provisions: shortfall / excess and the tier-2 cap
  expect_equal(t$provisions, sum(d$provision))
  expect_equal(t$shortfall, max(0, t$el - t$provisions)); expect_equal(t$excess, max(0, t$provisions - t$el))
  expect_equal(t$tier2_cap, 0.006 * t$rwa_irb)
  expect_equal(t$tier2_addback, min(t$excess, t$tier2_cap))
  expect_equal(cap$segments$shortfall_excess, cap$segments$provisions - cap$segments$el)
  # floors: the PD floor bites on the grades below it, its delta equals the no-floor bridge
  fl <- cap$floors
  expect_equal(fl$floor, c("pd_floor", "lgd_floor", "m_floor"))
  expect_equal(fl$n_hit[1], sum(d$pd < 5e-4 & d$defaulted == 0L & d$asset_class != "qrre_revolver") +
                              sum(d$pd < 1e-3 & d$defaulted == 0L & d$asset_class == "qrre_revolver"))
  expect_equal(fl$delta_rwa[1], t$rwa_irb - t$rwa_irb_no_floors)
  expect_gt(fl$delta_rwa[1], 0)
  expect_equal(fl$n_hit[2:3], c(0L, 0L))
  # concentration
  expect_equal(t$hhi, sum((d$ead / sum(d$ead))^2))
  expect_equal(t$n_eff, 1 / t$hhi)
  expect_equal(sum(cap$concentration$ead_share), 1)
  expect_equal(sum(cap$concentration$hhi_contribution), sum((cap$segments$ead / t$ead)^2))
  # sensitivity grid
  s <- cap$sensitivity
  expect_equal(s$shock, c("base", "pd_x1.10", "pd_x1.25", "pd_x1.50", "lgd_plus_5pp", "ead_x1.10", "no_floor", "r_x1.25",
                          "vasicek_q0.95", "vasicek_q0.99"))
  expect_equal(s$rwa[1], t$rwa_irb)
  expect_equal(s$delta_pct[s$shock == "ead_x1.10"], 0.10, tolerance = 1e-12)
  expect_equal(s$rwa[s$shock == "no_floor"], t$rwa_irb_no_floors)
  expect_true(all(diff(s$rwa[2:4]) > 0))
  expect_gt(s$rwa[s$shock == "vasicek_q0.99"], s$rwa[s$shock == "vasicek_q0.95"])
  expect_gt(s$rwa[s$shock == "r_x1.25"], t$rwa_irb)
  # pools: one per segment and grade; retail pools homogeneous, corporate ones not (maturity varies)
  expect_equal(nrow(cap$pools), data.table::uniqueN(d[, c("segment", "grade")]))
  expect_true(all(cap$pools$homogeneous[!cap$pools$segment %in% c("corporate_large", "corporate_sme")]))
  # ledger, model card, print
  expect_true(all(c("framework", "inputs", "floors", "output_floor", "provisions") %in% cap$ledger$action))
  expect_equal(cap$model_card$n_exposures, 5000L)
  expect_false(cap$model_card$params_modified)
  expect_output(print(cap), "output floor")
  expect_output(print(cap), "top segments")
  # no rows kept by default; constant asset class; no provisions
  c0 <- scr_capital(d[d$asset_class == "retail_other", ], config = cfg_test())
  expect_null(c0$exposures)
  expect_true(is.na(c0$totals$provisions)); expect_true(is.na(c0$totals$shortfall))
  expect_equal(c0$segments$segment, "portfolio")
  expect_output(print(c0), "portfolio")
  # switches
  c1 <- scr_capital(d[1:200, ], asset_class = "asset_class", config = cfg_test(capital_output_floor = FALSE, capital_sensitivity = FALSE))
  expect_null(c1$sensitivity); expect_true(is.na(c1$totals$rwa_sa)); expect_equal(c1$totals$rwa_reported, c1$totals$rwa_irb)
  expect_false(c1$totals$floor_binding)
  # errors
  expect_error(scr_capital(d, pd = "nope", config = cfg_test()), "not found")
  expect_error(scr_capital(d, asset_class = "boats", config = cfg_test()), "asset_class")
  expect_error(scr_capital(list(pd = NULL), config = cfg_test()), "data")
  expect_error(scr_capital(d, config = list()), "scr_config")
  # modified params are recorded
  p <- scr_irb_params("bcb"); p$scaling_factor <- 1.06
  cm <- scr_capital(d[1:300, ], asset_class = "asset_class", params = p, config = cfg_test())
  expect_true(cm$model_card$params_modified)
  expect_output(print(cm), "params modified")
})

test_that("the output floor binds when the IRB result is below 72.5 % of the standardised one", {
  cfg <- cfg_test()
  low <- data.frame(pd = 0.001, lgd = 0.10, ead = 1000, ltv = 0.4, segment = "mtg")
  cap <- scr_capital(low, segment = "segment", asset_class = "retail_mortgage", ltv = "ltv", config = cfg)
  t <- cap$totals
  expect_equal(t$rwa_sa, 200)                                    # 20 % band
  expect_lt(t$rwa_irb, 0.725 * t$rwa_sa)
  expect_true(t$floor_binding)
  expect_equal(t$rwa_reported, 0.725 * 200)
  expect_lt(t$headroom, 0)
  expect_output(print(cap), "BINDING")
  hi <- data.frame(pd = 0.10, lgd = 0.60, ead = 1000)
  ch <- scr_capital(hi, asset_class = "retail_other", config = cfg)
  expect_false(ch$totals$floor_binding)
  expect_equal(ch$totals$rwa_reported, ch$totals$rwa_irb)
  # a lower floor in an edited params object
  p <- scr_irb_params("bcb"); p$output_floor <- 0.5
  c2 <- scr_capital(low, asset_class = "retail_mortgage", ltv = "ltv", params = p, config = cfg)
  expect_equal(c2$totals$rwa_floor, 100)
})

test_that("the list form applies fitted models through scr_apply() and records the provenance", {
  registerS3method("scr_apply", "mock_pd_model", function(x, newdata, ...) data.table::data.table(pd_final = rep(x$pd, nrow(newdata))),
                   envir = asNamespace("scorecraft"))
  registerS3method("scr_apply", "mock_ead_model", function(x, newdata, ...) data.table::data.table(ead_predicted = newdata$ead * 1.1),
                   envir = asNamespace("scorecraft"))
  registerS3method("scr_apply", "mock_bad_model", function(x, newdata, ...) data.table::data.table(other = 1),
                   envir = asNamespace("scorecraft"))
  d <- scr_demo_portfolio[1:100, ]
  x <- list(pd = structure(list(pd = 0.03), class = "mock_pd_model"), ead = structure(list(), class = "mock_ead_model"), data = d)
  cap <- scr_capital(x, asset_class = "retail_other", config = cfg_test(), keep_rows = TRUE)
  expect_equal(cap$exposures$pd, rep(0.03, 100))
  expect_equal(cap$exposures$ead, d$ead * 1.1)
  expect_equal(cap$exposures$lgd, d$lgd)                            # lgd fell back to the column
  expect_equal(cap$model_card$pd_source, "scr_apply(mock_pd_model)$pd_final")
  expect_equal(cap$model_card$ead_source, "scr_apply(mock_ead_model)$ead_predicted")
  expect_equal(cap$model_card$lgd_source, "column")
  expect_true(grepl("mock_pd_model", cap$ledger$detail[cap$ledger$action == "inputs"]))
  expect_error(scr_capital(list(lgd = structure(list(), class = "mock_bad_model"), data = d), asset_class = "retail_other", config = cfg_test()), "lgd_final")
})

test_that("scr_capital is identical under the serial and PSOCK backends", {
  d <- scr_demo_portfolio[1:400, ]
  a <- withr::with_options(list(scorecraft.parallel = "serial"),
    do.call(scr_capital, c(list(d, config = cfg_test(nthread = 2L)), cap_args())))
  b <- withr::with_options(list(scorecraft.parallel = "psock"),
    do.call(scr_capital, c(list(d, config = cfg_test(nthread = 2L)), cap_args())))
  expect_equal(a$sensitivity, b$sensitivity)
  expect_equal(a$totals, b$totals)
  expect_equal(a$segments, b$segments)
})

test_that("the SQL reproduces R on DuckDB and SQLite, and every dialect is emitted", {
  cap <- cap_demo()
  d <- scr_demo_portfolio
  ex <- cap$exposures
  dialects <- c("ansi", "databricks", "spark", "hive", "mysql", "mariadb", "sqlserver", "bigquery", "postgres",
                "oracle", "snowflake", "redshift", "duckdb", "sqlite")
  for (dl in dialects) {
    s <- scr_sql(cap, table = "port", dialect = dl)
    expect_true(any(grepl("^WITH pool_params AS", s)), info = dl)
    expect_true(any(grepl("exposure_capital AS", s)), info = dl)
    expect_equal(sum(grepl("UNION ALL", s)), nrow(cap$pools) - 1L, info = dl)
    expect_false(any(grepl("qnorm|pnorm", s)), info = dl)
  }
  expect_true(any(grepl("FROM DUAL", scr_sql(cap, dialect = "oracle"))))
  expect_false(any(grepl("FROM DUAL", scr_sql(cap, dialect = "duckdb"))))
  agg <- scr_sql(cap, table = "port", dialect = "duckdb", level = "portfolio")
  expect_true(any(grepl("GROUP BY segment", agg)))
  f <- withr::local_tempfile(fileext = ".sql")
  expect_invisible(scr_sql(cap, file = f))
  expect_true(file.exists(f))
  # retail rows (homogeneous pools) match per exposure; every segment matches in aggregate
  check <- function(con, dialect) {
    DBI::dbWriteTable(con, "port", d)
    got <- DBI::dbGetQuery(con, paste(scr_sql(cap, table = "port", dialect = dialect), collapse = "\n"))
    got <- got[match(d$id, got$id), ]
    retail <- !d$segment %in% c("corporate_large", "corporate_sme")
    expect_equal(got$el[retail], ex$el[retail], tolerance = 1e-9)
    expect_equal(got$rwa[retail], ex$rwa[retail], tolerance = 1e-9)
    expect_equal(got$k[retail], ex$k[retail], tolerance = 1e-9)
    expect_equal(got$rw[retail], ex$rw[retail], tolerance = 1e-9)
    ag <- DBI::dbGetQuery(con, paste(scr_sql(cap, table = "port", dialect = dialect, level = "portfolio"), collapse = "\n"))
    ag <- ag[match(cap$segments$segment, ag$segment), ]
    expect_equal(ag$el, cap$segments$el, tolerance = 1e-9)
    expect_equal(ag$rwa, cap$segments$rwa_irb, tolerance = 1e-9)
    expect_equal(ag$n, cap$segments$n)
    expect_equal(sum(ag$rwa), cap$totals$rwa_irb, tolerance = 1e-9)
  }
  skip_if_not_installed("DBI")
  if (requireNamespace("duckdb", quietly = TRUE)) {
    con <- DBI::dbConnect(duckdb::duckdb(), config = list(threads = "1"))
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
    check(con, "duckdb")
  }
  skip_if_not_installed("RSQLite")
  con2 <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con2), add = TRUE)
  check(con2, "sqlite")
  # a single pool without segment or grade is a CROSS JOIN, also exact
  small <- d[d$asset_class == "retail_other", ][1:50, ]
  small$pd <- 0.02; small$lgd <- 0.5
  c1 <- scr_capital(small, asset_class = "retail_other", defaulted = "defaulted", elbe = "elbe", id = "id", config = cfg_test(), keep_rows = TRUE)
  s1 <- scr_sql(c1, table = "small", dialect = "sqlite")
  expect_true(any(grepl("CROSS JOIN", s1)))
  DBI::dbWriteTable(con2, "small", small)
  g1 <- DBI::dbGetQuery(con2, paste(s1, collapse = "\n"))
  g1 <- g1[match(small$id, g1$id), ]
  expect_equal(g1$rwa, c1$exposures$rwa, tolerance = 1e-9)
  expect_equal(g1$el, c1$exposures$el, tolerance = 1e-9)
})

test_that("the capital workbook round-trips through openxlsx", {
  skip_if_not_installed("openxlsx")
  cap <- cap_demo()
  out <- withr::local_tempdir()
  cap2 <- scr_export(cap, out, stamp = FALSE)
  expect_true(file.exists(cap2$files$xlsx)); expect_true(file.exists(cap2$files$sql))
  expect_equal(basename(cap2$files$xlsx), "capital_bcb.xlsx")
  sheets <- openxlsx::getSheetNames(cap2$files$xlsx)
  expect_true(all(c("Capital_Summary", "Capital_Config", "Segments_Reconciliation", "Pools", "Floors_Impact", "Output_Floor_Bridge",
                    "Sensitivity", "Concentration", "EL_vs_Provisions", "Exposures", "Model_Card", "Decision_Ledger") %in% sheets))
  seg <- openxlsx::read.xlsx(cap2$files$xlsx, sheet = "Segments_Reconciliation")
  expect_equal(nrow(seg), nrow(cap$segments))
  expect_equal(seg$rwa_irb[match(cap$segments$segment, seg$segment)], cap$segments$rwa_irb, tolerance = 1e-8)
  br <- openxlsx::read.xlsx(cap2$files$xlsx, sheet = "Output_Floor_Bridge")
  expect_equal(br$value[br$step == "rwa_reported"], cap$totals$rwa_reported, tolerance = 1e-8)
  ex <- openxlsx::read.xlsx(cap2$files$xlsx, sheet = "Exposures")
  expect_equal(nrow(ex), 5000L)
  # without rows kept the sheet carries an availability row
  c0 <- scr_capital(scr_demo_portfolio[1:50, ], asset_class = "asset_class", config = cfg_test())
  c0 <- scr_export(c0, file.path(out, "b"), stamp = FALSE)
  expect_false("Exposures" %in% openxlsx::getSheetNames(c0$files$xlsx))
  expect_true("Capital_Summary" %in% openxlsx::getSheetNames(c0$files$xlsx))
})

test_that("ECL of a flat hazard equals the closed form, with discounting, prepayment and stages", {
  cfg <- cfg_test(ecl_discount = "none")
  h <- 0.01; lgd <- 0.4; ead <- 1000
  e <- scr_ecl(h, lgd, ead, config = cfg, keep_rows = TRUE)
  expect_s3_class(e, "scr_ecl")
  expect_equal(e$totals$ecl_12m, lgd * ead * (1 - (1 - h)^12), tolerance = 1e-12)
  expect_equal(e$totals$ecl_life, e$totals$ecl_12m)                           # t_max defaults to the horizon
  expect_equal(e$exposures$stage, 1L)
  expect_equal(e$exposures$pd_12m, 1 - (1 - h)^12, tolerance = 1e-12)
  e36 <- scr_ecl(h, lgd, ead, t_max = 36L, config = cfg)
  expect_equal(e36$totals$ecl_12m, lgd * ead * (1 - (1 - h)^12), tolerance = 1e-12)
  expect_equal(e36$totals$ecl_life, lgd * ead * (1 - (1 - h)^36), tolerance = 1e-12)
  # a matrix with the same flat hazard gives the same answer; a hand path is exact
  em <- scr_ecl(matrix(h, 1, 36), lgd, ead, config = cfg)
  expect_equal(em$totals$ecl_life, e36$totals$ecl_life, tolerance = 1e-12)
  hp <- c(0.02, 0.03, 0.05)
  ep <- scr_ecl(matrix(hp, 1), matrix(c(0.5, 0.6, 0.7), 1), matrix(c(100, 90, 80), 1), config = cfg)
  expect_equal(ep$totals$ecl_life, 0.02 * 0.5 * 100 + 0.98 * 0.03 * 0.6 * 90 + 0.98 * 0.97 * 0.05 * 0.7 * 80, tolerance = 1e-12)
  # discounting lowers the loss by the monthly factors
  ed <- scr_ecl(matrix(hp, 1), 0.5, 100, eir = 0.12, config = cfg_test())
  s <- c(1, 0.98, 0.98 * 0.97)
  expect_equal(ed$totals$ecl_life, sum(s * hp * 0.5 * 100 * 1.12^(-(1:3) / 12)), tolerance = 1e-12)
  expect_lt(ed$totals$ecl_life, scr_ecl(matrix(hp, 1), 0.5, 100, config = cfg)$totals$ecl_life)
  # a vector is one flat hazard per exposure, a matrix one path per exposure
  ev <- scr_ecl(hp, 0.5, 100, config = cfg, keep_rows = TRUE)
  expect_equal(nrow(ev$exposures), 3L)
  expect_equal(ev$exposures$ecl_12m, 0.5 * 100 * (1 - (1 - hp)^12), tolerance = 1e-12)
  expect_equal(ed$discount, "eir")
  # prepayment shortens the survival
  pp <- scr_ecl(h, lgd, ead, prepay = 0.05, t_max = 24L, config = cfg)
  s_prev <- cumprod(c(1, rep(1 - h - 0.05, 23)))
  expect_equal(pp$totals$ecl_life, sum(s_prev * h * lgd * ead), tolerance = 1e-12)
  expect_lt(pp$totals$ecl_life, scr_ecl(h, lgd, ead, t_max = 24L, config = cfg)$totals$ecl_life)
  # stage rule: dpd and the PD deterioration ratio; stage 3 carries LGD x EAD
  st <- scr_ecl(rep(h, 5), lgd, ead, dpd = c(0, 30, 90, 0, 0), pd_orig = c(0.2, 0.2, 0.2, 0.05, 0.2), t_max = 24L,
                config = cfg, keep_rows = TRUE)
  expect_equal(st$exposures$stage, c(1L, 2L, 3L, 2L, 1L))
  expect_equal(st$exposures$ecl[1], st$exposures$ecl_12m[1])
  expect_equal(st$exposures$ecl[2], st$exposures$ecl_life[2])
  expect_equal(st$exposures$ecl[3], lgd * ead)
  expect_equal(st$stages$n, c(2L, 2L, 1L))
  expect_equal(st$stages$coverage, st$stages$ecl / st$stages$ead)
  expect_equal(st$totals$share_stage3, 0.2)
  expect_true(grepl("stage 3 if dpd >= 90", st$ledger$detail[st$ledger$action == "stage"]))
  sup <- scr_ecl(rep(h, 3), lgd, ead, stage = c(1, 2, 3), config = cfg)
  expect_equal(sup$stage_rule$source, "supplied")
  expect_equal(sup$stages$n, c(1L, 1L, 1L))
  s2 <- scr_ecl(h, lgd, ead, config = cfg_test(ecl_stage_dpd = c(15L, 60L)), dpd = 20)
  expect_equal(s2$stages$n[2], 1L)
  # scenarios: the weighted ECL is the weighted mean, weights normalised, z shocks move the hazard
  sc <- scr_ecl(h, lgd, ead, t_max = 12L, scenarios = list(base = list(), bad = list(pd_mult = 2), worse = list(lgd_add = 0.1, ead_mult = 1.1)),
                weights = c(2, 1, 1), config = cfg)
  expect_equal(sc$scenarios$weight, c(0.5, 0.25, 0.25))
  expect_equal(sc$scenarios$ecl[2], lgd * ead * (1 - (1 - 2 * h)^12), tolerance = 1e-12)
  expect_equal(sc$scenarios$ecl[3], 0.5 * 1100 * (1 - (1 - h)^12), tolerance = 1e-12)
  expect_equal(sc$totals$ecl, sum(sc$scenarios$weight * sc$scenarios$ecl), tolerance = 1e-12)
  expect_true(grepl("normalised", sc$ledger$reason[sc$ledger$action == "scenarios"]))
  z <- scr_ecl(h, lgd, ead, t_max = 12L, scenarios = list(good = list(z = 1), bad = list(z = -1)), rho = 0.15, config = cfg)
  expect_lt(z$scenarios$ecl[1], z$scenarios$ecl[2])
  expect_equal(z$scenarios$ecl[2], lgd * ead * (1 - (1 - .vasicek_pit(h, -1, 0.15))^12), tolerance = 1e-12)
  # segments and print
  d <- scr_demo_portfolio
  hd <- 1 - (1 - d$pd)^(1 / 12)
  e2 <- scr_ecl(hd, d$lgd, d$ead, eir = d$eir, dpd = d$dpd, pd_orig = d$pd_orig, t_max = 36L, segment = d$segment, id = d$id, config = cfg_test())
  expect_equal(sum(e2$segments$ecl), e2$totals$ecl)
  expect_equal(sum(e2$stages$n), 5000L)
  expect_null(e2$exposures)
  expect_output(print(e2), "stage 3")
  expect_output(print(sc), "scenarios: base 0.50")
  # errors
  expect_error(scr_ecl(1.2, lgd, ead, config = cfg), "pd_term")
  expect_error(scr_ecl(matrix(h, 2, 12), matrix(lgd, 2, 6), ead, config = cfg), "columns")
  expect_error(scr_ecl(h, lgd, ead, stage = 4, config = cfg), "stage")
  expect_error(scr_ecl(h, lgd, ead, scenarios = list(list(pd_mult = 2)), config = cfg), "named list")
  expect_error(scr_ecl(h, lgd, ead, scenarios = list(a = list(foo = 1)), config = cfg), "scenario key")
  expect_error(scr_ecl(h, lgd, ead, scenarios = list(a = list()), weights = c(1, 2), config = cfg), "weights")
  expect_error(scr_ecl(h, lgd, ead, prepay = 0.995, config = cfg), "prepay")
})

test_that("the foundation approach reads the supervisory LGD of the claim and retail can be non-granular", {
  bcb <- scr_irb_params("bcb"); b3 <- scr_irb_params("basel3_final")
  r1 <- scr_irb_rw(0.01, 0, m = 2.5, asset_class = "corporate", approach = "firb", claim = "senior_unsecured", params = bcb)
  expect_equal(r1$lgd_used, 0.75)
  r2 <- scr_irb_rw(0.01, 0, m = 2.5, asset_class = "corporate", approach = "firb", claim = "senior_unsecured_corporate", params = b3)
  expect_equal(r2$lgd_used, 0.40)
  expect_equal(r2$m, 2.5)
  r3 <- scr_irb_rw(0.01, 0.40, m = 2.5, asset_class = "corporate", approach = "airb", apply_floors = FALSE, params = b3)
  expect_equal(r2$k, r3$k, tolerance = 1e-12)
  expect_error(scr_irb_rw(0.01, 0, asset_class = "corporate", approach = "firb", claim = "boats"), "claim")
  r4 <- scr_irb_rw(0.01, 0.45, m = 2.5, asset_class = "corporate", approach = "airb", claim = "subordinated")
  expect_equal(r4$lgd_used, 0.45)
  expect_equal(scr_sa_rw("retail_other"), 0.75)
  expect_equal(scr_sa_rw("retail_other", granular = FALSE), 1.00)
  expect_equal(scr_sa_rw(c("retail_other", "retail_other"), transactor = c(TRUE, TRUE), granular = c(TRUE, FALSE)), c(0.45, 1.00))
  expect_equal(scr_sa_rw("retail_mortgage", ltv = 0.55, granular = FALSE), 0.25)
  port <- data.table::as.data.table(scr_demo_portfolio)[asset_class %in% c("corporate", "retail_other")][1:400]
  port[, claim := ifelse(asset_class == "corporate", "senior_unsecured", NA_character_)]
  cf <- cfg_test(capital_approach = "firb", capital_sensitivity = FALSE)
  cap_f <- scr_capital(port, segment = "segment", asset_class = "asset_class", m = "m", claim = "claim", config = cf, keep_rows = TRUE)
  expect_true(all(cap_f$exposures[asset_class == "corporate", lgd_used] == 0.75))
  cap_g <- scr_capital(port, segment = "segment", asset_class = "asset_class", m = "m", granular = FALSE,
                       config = cfg_test(capital_sensitivity = FALSE), keep_rows = TRUE)
  cap_0 <- scr_capital(port, segment = "segment", asset_class = "asset_class", m = "m",
                       config = cfg_test(capital_sensitivity = FALSE), keep_rows = TRUE)
  expect_gt(cap_g$totals$rwa_sa, cap_0$totals$rwa_sa)
  expect_equal(cap_g$totals$rwa_irb, cap_0$totals$rwa_irb)
})
