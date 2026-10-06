# PD module: master scale, calibration of the alignment, grades on the score,
# margin of conservatism, floors, migration, validation battery, PIT/TTC
# bridge, SQL equivalence, backends and the workbook.

test_that("the geometric master scale doubles the PD every grade and the supplied one is validated", {
  ms <- scr_master_scale(0.001, 0.064, n_grades = 7)
  expect_s3_class(ms, "scr_master_scale")
  expect_equal(attr(ms, "ratio"), 2)
  expect_equal(ms$pd_mid, 0.001 * 2^(0:6))
  expect_equal(ms$pd_hi[1:6], sqrt(ms$pd_mid[1:6] * ms$pd_mid[2:7]))
  expect_equal(ms$pd_lo[1], 0); expect_equal(ms$pd_hi[7], 1)
  expect_equal(ms$pd_lo[-1], ms$pd_hi[-7])
  expect_output(print(ms), "7 grades")
  s <- scr_master_scale(method = "supplied", grades = c(0.05, 0.001, 0.01))
  expect_equal(s$pd_mid, c(0.001, 0.01, 0.05))
  s2 <- scr_master_scale(method = "supplied", grades = data.frame(pd_lo = c(0, 0.02, 0.1), pd_hi = c(0.02, 0.1, 1), label = c("A", "B", "C")))
  expect_equal(s2$label, c("A", "B", "C"))
  expect_equal(s2$pd_mid[2], sqrt(0.02 * 0.1))
  expect_error(scr_master_scale(0.1, 0.01), "pd_max")
  expect_error(scr_master_scale(method = "supplied"), "grades")
  expect_error(scr_master_scale(method = "supplied", grades = data.frame(pd_lo = c(0, 0.5), pd_hi = c(0.2, 1))), "contiguous")
  expect_error(scr_master_scale(0.001, 0.1, labels = "a"), "labels")
})

test_that("the intercept shift moves the mean PD to the central tendency exactly and leaves the AUC unchanged", {
  sc <- sc_demo(); cal <- pd_cal()
  expect_s3_class(cal, "scr_pd_calibration")
  expect_s3_class(cal$alignment, "scr_align")
  ho <- sc$samples$holdout
  p <- predict(cal$alignment, ho$link, type = "prob")
  expect_equal(mean(p), 0.06, tolerance = 1e-10)
  expect_equal(cal$mean_pd_after, 0.06, tolerance = 1e-10)
  expect_equal(cal$slope_ratio, 1)
  expect_equal(cal$alignment$calibration$slope, sc$alignment$calibration$slope)
  # King-Zeng closed form is the prior shift, the exact root is close to it
  ybar <- mean(ho$y)
  expect_equal(cal$shift_prior, log(0.06 * (1 - ybar) / (0.94 * ybar)))
  expect_lt(abs(cal$shift - cal$shift_prior), 0.1)
  # the same rank order: the AUC of the calibrated PD equals the hold-out AUC of the score
  auc0 <- sc$metrics[sample == "holdout", auc]
  expect_equal(scr_metrics(p, ho$y, ci = FALSE)$auc, auc0, tolerance = 1e-12)
  expect_equal(cal$ar_after, cal$ar_before)
  # the scorecard itself is untouched and the score of the scorecard alignment is unchanged
  expect_equal(predict(sc$alignment, ho$link), ho$score)
  expect_output(print(cal), "mean PD")
  # an scr_dr as the target uses its long-run average
  dr <- pd_series()
  c2 <- scr_calibrate(sc, target = dr)
  expect_equal(c2$ct, dr$lra$mean)
  expect_match(c2$target_source, "long-run")
})

test_that("the two-parameter methods hit the CT and the accuracy-ratio target; scaling and segments work", {
  sc <- sc_demo(); ho <- sc$samples$holdout
  ab <- scr_calibrate(sc, target = 0.06, method = "logodds_ab", ar_target = 0.55)
  p <- predict(ab$alignment, ho$link, type = "prob")
  expect_equal(mean(p), 0.06, tolerance = 1e-10)
  expect_equal(ab$ar_implied_after, 0.55, tolerance = 1e-6)
  expect_equal(ab$ar_target, 0.55)
  expect_equal(ab$alignment$calibration$slope, ab$slope_ratio * sc$alignment$calibration$slope)
  # default target of logodds_ab is the observed AR; of qmm the implied AR of the current PDs
  ab0 <- scr_calibrate(sc, target = 0.06, method = "logodds_ab")
  expect_equal(ab0$ar_implied_after, ab0$ar_before, tolerance = 1e-6)
  q <- scr_calibrate(sc, target = 0.06, method = "qmm")
  expect_equal(q$ar_implied_after, q$ar_implied_before, tolerance = 1e-6)
  expect_equal(q$mean_pd_after, 0.06, tolerance = 1e-10)
  s <- scr_calibrate(sc, target = 0.06, method = "scaling")
  expect_equal(s$mean_pd_after, 0.06, tolerance = 1e-10)
  expect_match(s$note, "rescaling")
  expect_true(s$slope_ratio > 0 && s$slope_ratio < 1)
  # numeric ln(odds) and an alignment with only a sample rate (closed form)
  set.seed(1); l <- stats::qlogis(0.12) + stats::rnorm(500); y <- stats::rbinom(500, 1, stats::plogis(l))
  cn <- scr_calibrate(l, target = 0.04, y = y)
  expect_equal(mean(predict(cn$alignment, l, type = "prob")), 0.04, tolerance = 1e-10)
  ca <- scr_calibrate(sc$alignment, target = 0.04, sample_rate = 0.10)
  expect_equal(ca$shift, log(0.04 * 0.9 / (0.96 * 0.10)))
  expect_true(is.na(ca$mean_pd_after))
  # segments: one alignment each and the pooled one
  seg <- ifelse(ho$score > stats::median(ho$score), "hi", "lo")
  cs <- scr_calibrate(sc, target = c(hi = 0.03, lo = 0.10), segment = seg)
  expect_equal(nrow(cs$segments$table), 2L)
  expect_equal(mean(predict(cs$segments$alignments$hi, ho$link[seg == "hi"], type = "prob")), 0.03, tolerance = 1e-10)
  expect_equal(cs$ct, sum(cs$segments$table$ct * cs$segments$table$n) / nrow(ho))
  expect_output(print(cs), "segments")
  expect_error(scr_calibrate(sc, target = 1.5), "central tendency")
  expect_error(scr_calibrate(sc$alignment, target = 0.05, sample_rate = 0.1, method = "qmm"), "calibration sample")
  expect_error(scr_calibrate(sc, target = c(hi = 0.03), segment = seg), "segment")
  expect_error(scr_calibrate(list(), target = 0.05), "expects")
})

test_that("grades are monotone after repair, every grade satisfies the minimum counts and the cut points are direction-aware", {
  sc <- sc_demo(); gr <- pd_grades()
  expect_s3_class(gr, "scr_grades")
  t <- gr$table
  expect_true(all(diff(t$pd_be) > 0))
  expect_true(all(t$n >= 30L)); expect_true(all(t$defaults >= 10L))
  expect_equal(sum(t$n), nrow(sc$samples$holdout)); expect_equal(sum(t$share), 1)
  expect_true(nrow(gr$repairs) >= 1L)
  expect_true(any(nzchar(t$merged_from)))
  # grade 1 is the safest: the highest scores under higher_is_safer
  expect_equal(gr$direction, "higher_is_safer")
  expect_equal(t$score_hi[1], Inf); expect_equal(t$score_lo[nrow(t)], -Inf)
  expect_true(all(diff(t$score_lo) < 0))
  expect_equal(gr$band_grade, rev(seq_len(nrow(t))))
  expect_equal(length(gr$breaks), nrow(t) - 1L)
  expect_true(all(diff(gr$breaks) > 0))
  # implied PD bounds follow the master scale and the grade PD sits inside them
  expect_true(all(t$pd_lo < t$pd_hi))
  expect_true(all(diff(t$pd_lo) > 0))
  expect_equal(t$pd_lo[1], 0); expect_equal(t$pd_hi[nrow(t)], 1)
  # the rows' grade agrees with predict() and the table counts
  expect_equal(as.integer(tabulate(gr$rows$grade, nbins = nrow(t))), t$n)
  expect_identical(predict(gr, score = gr$rows$score), gr$rows$grade)
  expect_equal(predict(gr, score = gr$rows$score, type = "pd"), gr$rows$pd)
  # concentration
  R <- t$n / sum(t$n); K <- nrow(t)
  expect_equal(gr$concentration$hhi, sum(R^2))
  cv <- sqrt(K * sum((R - 1 / K)^2))
  expect_equal(gr$concentration$cv, cv)
  expect_equal(gr$concentration$hi, 1 + log((cv^2 + 1) / K) / log(K))
  expect_output(print(gr), "repair")
  # mean_pd source and monotone = FALSE
  gm <- scr_grades(sc, pd_cal(), n_grades = 7L, min_obligors = 30L, min_defaults = 10L, pd_source = "mean_pd")
  expect_equal(gm$table$pd_be, gm$table$pd_mean)
  expect_equal(gm$breaks, gr$breaks)
  expect_equal(mean(gm$rows$pd), 0.06, tolerance = 1e-10)
  # a supplied master scale keeps its labels through the merges
  ms <- scr_master_scale(0.005, 0.4, n_grades = 8)
  gs <- scr_grades(sc, pd_cal(), method = "supplied", master_scale = ms, min_defaults = 5)
  expect_true(all(gs$table$label %in% c(ms$label, gs$table$merged_from)))
  expect_true(all(diff(gs$table$pd_be) > 0))
  expect_error(scr_grades(sc, pd_cal(), method = "supplied"), "master_scale")
  expect_error(scr_grades(sc, calibration = list()), "calibration")
  expect_error(scr_grades(sc, sample = "nope"), "sample")
})

test_that("quantile grades have equal counts and cut points that never coincide with an observed score", {
  sc <- sc_demo()
  gq <- scr_grades(sc, pd_cal(), method = "quantile", n_grades = 5L, min_obligors = 1L, min_defaults = 1L, monotone = FALSE)
  expect_equal(nrow(gq$table), 5L)
  expect_equal(gq$table$n, rep(280L, 5))
  expect_false(any(gq$breaks %in% sc$samples$holdout$score))
  expect_false(any(gq$breaks %in% sc$samples$train$score))
  expect_null(gq$master_scale)
  b <- .pd_quantile_breaks(c(1, 2, 3, 4, 5, 6), 3)
  expect_equal(b, c(2.5, 4.5))
})

test_that("a default-rate series keyed by the final grades feeds the long-run grade PD and the time-series MoC", {
  sc <- sc_demo(); gr <- pd_grades(); dr <- pd_series()
  expect_true("grade" %in% names(dr$table))
  g2 <- scr_grades(sc, pd_cal(), n_grades = 7L, min_obligors = 30L, min_defaults = 10L, dr = dr)
  expect_equal(g2$breaks, gr$breaks)
  expect_true(all(c("n_series", "t_series") %in% names(g2$table)))
  # pd_be is the arithmetic mean of the cohort default rates of the grade
  s <- dr$table
  lra1 <- mean(s$defaults[s$grade == "1"] / s$n[s$grade == "1"])
  expect_equal(g2$table$pd_be[1], lra1)
  expect_true(all(diff(g2$table$pd_be) > 0))
  expect_equal(g2$table$t_series, rep(dr$lra$n_cohorts, nrow(g2$table)))
  expect_match(g2$ledger$detail[g2$ledger$action == "grade_pd_source"], "long-run")
  # a series with foreign labels is refused with guidance
  bad <- dr; bad$table$grade <- paste0("g", bad$table$grade)
  expect_error(scr_grades(sc, pd_cal(), n_grades = 7L, min_obligors = 30L, min_defaults = 10L, dr = bad), "do not match")
  expect_error(scr_moc(gr, "C", method = "ci_timeseries"), "needs a default-rate series")
  # category C by time series: t quantile of the cohort dispersion
  gc <- scr_moc(gr, "C", method = "ci_timeseries", dr = dr, level = 0.95)
  r1 <- s$defaults[s$grade == "1"] / s$n[s$grade == "1"]
  expect_equal(gc$moc$value[gc$moc$grade == 1L], stats::qt(0.95, length(r1) - 1) * stats::sd(r1) / sqrt(length(r1)))
  expect_false(is.null(gc$dr))
})

test_that("the MoC ledger is append-only: A/B need value and reason, C is superseded, A accumulates", {
  gr <- pd_grades(); K <- nrow(gr$table)
  expect_error(scr_moc(gr, "A", value = 0.001), "reason")
  expect_error(scr_moc(gr, "B", reason = "x"), "value")
  expect_error(scr_moc(gr, "A", value = -1, reason = "x"), "non-negative")
  expect_error(scr_moc(gr, "A", value = c(1, 2), reason = "x"), "length")
  g <- scr_moc(gr, "C", method = "ci_binomial", level = 0.95)
  expect_equal(g$moc$value, stats::qnorm(0.95) * sqrt(gr$table$pd_be * (1 - gr$table$pd_be) / gr$table$n))
  g <- scr_moc(g, "A", value = 0.001, reason = "missing trigger")
  g <- scr_moc(g, "A", value = seq_len(K) * 1e-4, reason = "external data")
  g <- scr_moc(g, "B", value = 0.002, reason = "underwriting change")
  g <- scr_moc(g, "C", method = "bootstrap", n_boot = 30, seed = 1)
  expect_equal(nrow(g$moc), 5L * K)
  expect_equal(sum(g$moc$category == "C" & !g$moc$active), K)
  expect_equal(g$moc$id[g$moc$category == "C" & g$moc$active][1], 5L)
  pd <- scr_pd(g)
  expect_equal(pd$table$moc_a, 0.001 + seq_len(K) * 1e-4)
  expect_equal(pd$table$moc_b, rep(0.002, K))
  expect_equal(pd$table$moc_c, g$moc$value[g$moc$active & g$moc$category == "C"])
  expect_equal(pd$table$pd_moc, pmin(1, pd$table$pd_be + pd$table$moc_a + pd$table$moc_b + pd$table$moc_c))
  expect_true(all(pd$table$moc_c >= 0))
  expect_output(print(g), "margin of conservatism")
})

test_that("scr_pd assembles the grade table with the floor of the framework and a model card", {
  pd <- pd_model()
  expect_s3_class(pd, "scr_pd")
  t <- pd$table
  expect_setequal(names(t), c("grade", "label", "score_lo", "score_hi", "pd_lo", "pd_hi", "n", "share", "defaults", "dr", "pd_mean", "pd_be",
                              "merged_from", "moc_a", "moc_b", "moc_c", "pd_moc", "pd_ttc", "pd_pit", "floor", "pd_final", "floor_applied"))
  expect_equal(pd$floor, 5e-4); expect_equal(pd$framework, "bcb"); expect_equal(pd$asset_class, "retail_other")
  expect_equal(t$pd_final, pmax(t$pd_moc, 5e-4))
  expect_false(any(t$floor_applied))
  expect_true(all(is.na(t$pd_pit)))
  expect_true(all(diff(t$pd_final) > 0))
  expect_equal(pd$portfolio$pd_be, sum(t$n * t$pd_be) / sum(t$n))
  expect_equal(pd$portfolio$moc_bp, 1e4 * sum(t$n * (t$pd_moc - t$pd_be)) / sum(t$n))
  expect_equal(pd$model_card$n_grades, nrow(t)); expect_equal(pd$model_card$calibration_method, "intercept")
  expect_equal(pd$model_card$central_tendency, 0.06)
  expect_true("pd_floor" %in% pd$ledger$action)
  expect_output(print(pd), "portfolio: PD_BE")
  # a binding floor from an edited parameter table is applied and recorded
  p <- scr_irb_params("bcb"); p$pd_floor$floor[p$pd_floor$asset_class == "retail_other"] <- 0.05
  pf <- scr_pd(pd_grades(), params = p)
  expect_true(pf$table$floor_applied[1]); expect_equal(pf$table$pd_final[1], 0.05)
  expect_true(pf$params_modified); expect_true("params_modified" %in% pf$ledger$action)
  expect_gt(pf$portfolio$share_at_floor, 0)
  expect_true("moc_c_missing" %in% pf$ledger$action)
  # sovereign has no floor; unknown asset classes fail
  ps <- scr_pd(pd_grades(), asset_class = "sovereign")
  expect_equal(ps$floor, 0)
  expect_error(scr_pd(pd_grades(), asset_class = "boats"), "asset_class")
  # point-in-time philosophy
  g <- scr_moc(pd_grades(), "C", method = "ci_binomial")
  expect_error(scr_pd(g, philosophy = "pit"), "rho")
  pp <- scr_pd(g, philosophy = "pit", rho = 0.1, z = -1)
  expect_equal(pp$table$pd_pit, scr_pd_pit_ttc(pp$table$pd_moc, z = -1, rho = 0.1))
  expect_equal(pp$table$pd_final, pmax(pp$table$pd_pit, 5e-4))
  expect_true(all(pp$table$pd_pit > pp$table$pd_ttc))
  expect_equal(pp$model_card$philosophy, "pit")
  expect_error(scr_pd(list()), "scr_grades")
})

test_that("predict, scr_apply and the SQL agree on grade and PD (DuckDB and SQLite)", {
  pd <- pd_model(); sc <- sc_demo()
  new <- head(scr_demo, 300)
  ap <- scr_apply(pd, new)
  expect_setequal(names(ap), c("score", "score_points", "grade", "pd", "pd_be", "pd_final"))
  expect_equal(ap$score, scr_apply(sc, new)$score)
  expect_identical(ap$grade, predict(pd, new))
  expect_identical(ap$grade, predict(pd, score = ap$score, type = "grade"))
  expect_equal(ap$pd_final, predict(pd, score = ap$score, type = "pd_final"))
  expect_equal(ap$pd, predict(pd, score = ap$score, type = "pd"))
  expect_equal(ap$pd_final, pd$table$pd_final[ap$grade])
  expect_equal(predict(pd, score = c(480, 640), type = "score"), c(480, 640))
  expect_error(predict(pd), "newdata")
  # boundaries: a score exactly on a cut point belongs to the lower band (score <= break)
  b <- pd$breaks
  expect_identical(predict(pd, score = b), pd$band_grade[seq_along(b)])
  expect_identical(predict(pd, score = b + 1e-9), pd$band_grade[seq_along(b) + 1L])
  sql <- scr_sql(pd)
  txt <- paste(sql, collapse = "\n")
  expect_match(txt, "score_scr AS"); expect_match(txt, "AS grade"); expect_match(txt, "AS pd_final")
  expect_match(paste(scr_sql(pd, table = "prd.t", dialect = "databricks"), collapse = "\n"), "FROM prd.t")
  tmp <- tempfile(fileext = ".sql"); expect_invisible(scr_sql(pd, file = tmp)); expect_true(file.exists(tmp))
  skip_if_not_installed("DBI")
  if (requireNamespace("duckdb", quietly = TRUE)) {
    con <- DBI::dbConnect(duckdb::duckdb(), config = list(threads = "1"))
    DBI::dbWriteTable(con, "scr_demo", scr_demo)
    got <- DBI::dbGetQuery(con, paste(scr_sql(pd, table = "scr_demo", dialect = "duckdb"), collapse = "\n"))
    DBI::dbDisconnect(con, shutdown = TRUE)
    exp <- scr_apply(pd, scr_demo)
    expect_equal(got$score, exp$score, tolerance = 1e-10)
    expect_identical(as.integer(got$grade), exp$grade)
    expect_equal(got$pd_final, exp$pd_final, tolerance = 1e-12)
  }
  if (requireNamespace("RSQLite", quietly = TRUE)) {
    con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
    d <- scr_demo; d$ref_date <- as.character(d$ref_date)
    DBI::dbWriteTable(con, "scr_demo", d)
    got <- DBI::dbGetQuery(con, paste(scr_sql(pd, table = "scr_demo", dialect = "sqlite"), collapse = "\n"))
    DBI::dbDisconnect(con)
    exp <- scr_apply(pd, scr_demo)
    expect_identical(as.integer(got$grade), exp$grade)
    expect_equal(got$pd_final, exp$pd_final, tolerance = 1e-12)
  }
})

test_that("the calibration tests are pinned: Jeffreys, binomial critical count, normal z, Hosmer-Lemeshow", {
  r <- .pd_cal_tests(N = 1000, D = 12, PD = 0.01, alpha = 0.05)
  expect_equal(r$p_jeffreys, 0.252381226151348, tolerance = 1e-12)
  expect_equal(r$critical, 16)
  expect_gt(r$p_binomial, 0.05)
  # P(X >= 16 | N = 1000, PD = 0.01): the first count at or below alpha
  expect_equal(.pd_cal_tests(1000, 16, 0.01, 0.05)$p_binomial, 0.0478705857579472, tolerance = 1e-12)
  expect_equal(r$z, (0.012 - 0.01) / sqrt(0.01 * 0.99 / 1000))
  expect_equal(r$p_normal, 1 - stats::pnorm(r$z))
  # a grade with more defaults than the critical count is rejected: D = 16 gives p <= alpha
  expect_lte(.pd_cal_tests(1000, 16, 0.01, 0.05)$p_binomial, 0.05)
  expect_gt(.pd_cal_tests(1000, 15, 0.01, 0.05)$p_binomial, 0.05)
  # Hosmer-Lemeshow on a hand table: only the third grade deviates, (40 - 30)^2 / (500 * 0.06 * 0.94);
  # the PDs are fixed before the validation sample is seen, so the reference is chi-square with K df
  hl <- .pd_hl(n = c(1000, 800, 500), d = c(10, 24, 40), pd = c(0.01, 0.03, 0.06))
  expect_equal(hl$chi2, 100 / 28.2)
  expect_equal(hl$df, 3L)
  expect_equal(hl$p, stats::pchisq(100 / 28.2, 3, lower.tail = FALSE))
  expect_equal(.pd_hl(c(100, 100), c(1, 2), c(0.01, 0.02))$df, 2L)
  expect_equal(.pd_light(c(0.005, 0.03, 0.2, NA), c(0.01, 0.05)), c("red", "amber", "green", NA))
})

test_that("the migration matrix and its bandwidths are pinned on a hand matrix", {
  N <- matrix(c(80, 15, 5, 10, 70, 20, 5, 15, 80), 3, byrow = TRUE)
  g0 <- rep(1:3, times = rowSums(N))
  g1 <- unlist(lapply(1:3, function(i) rep(1:3, times = N[i, ])))
  m <- scr_migration(g0, g1, K = 3)
  expect_s3_class(m, "scr_migration")
  expect_equal(unname(m$matrix[, 1:3]), N)
  expect_equal(unname(m$n), rowSums(N))
  expect_equal(m$mwb_upper, 0.75)
  expect_equal(m$mwb_lower, 0.70)
  expect_equal(m$mobility$share_stable, 230 / 300)
  expect_equal(m$mobility$share_up, 40 / 300); expect_equal(m$mobility$share_down, 30 / 300)
  expect_equal(m$mobility$mean_distance, (15 + 10 + 20 + 15 + 2 * 5 + 2 * 5) / 300)
  # z of the cell two steps away against its neighbour: p decays, so z is negative
  expect_lt(m$z[1, 3], 0)
  expect_equal(m$n_significant, 0L)
  expect_true(is.na(m$z[1, 1]))
  # default and closed columns stay out of the bandwidths
  m2 <- scr_migration(c(g0, 1, 2), c(g1, "D", NA), K = 3)
  expect_equal(m2$matrix[1, "default"], 1L); expect_equal(m2$matrix[2, "closed"], 1L)
  expect_equal(m2$mobility$share_default, 1 / 302)
  expect_equal(sum((m2$p * m2$n)[, 1:3]), 300)
  # an upgrade-only matrix has full upper bandwidth when every move is the longest one
  m3 <- scr_migration(c(1, 1, 3, 3), c(3, 3, 1, 1), K = 3)
  expect_equal(m3$mwb_upper, 1); expect_equal(m3$mwb_lower, 1)
  expect_output(print(m), "MWB upper 0.7500")
  expect_error(scr_migration(1:3, 1:2), "same length")
  expect_error(scr_migration(c(0, 1), c(1, 1), K = 2), "1..K")
})

test_that("the PIT/TTC bridge is Vasicek's formula, invertible and neutral at z = 0", {
  pd <- c(0.01, 0.05, 0.2)
  pit <- scr_pd_pit_ttc(pd, z = -2, rho = 0.15)
  expect_equal(pit, stats::pnorm((stats::qnorm(pd) + sqrt(0.15) * 2) / sqrt(0.85)))
  expect_equal(scr_pd_pit_ttc(pit, z = -2, rho = 0.15, to = "ttc"), pd)
  expect_equal(scr_pd_pit_ttc(pd, z = 0, rho = 0.15), stats::pnorm(stats::qnorm(pd) / sqrt(0.85)))
  expect_true(all(scr_pd_pit_ttc(pd, z = 2, rho = 0.15) < pd))
  # shared with the capital module
  expect_equal(scr_pd_pit_ttc(0.02, z = -1, rho = 0.1), .vasicek_pit(0.02, -1, 0.1))
  expect_error(scr_pd_pit_ttc(0.02, z = 0, rho = 1), "rho")
  expect_error(scr_pd_pit_ttc(1.2, z = 0, rho = 0.1), "pd")
})

test_that("the validation battery runs on the cohort panel with lights, discrimination and stability", {
  pd <- pd_model(); p <- pd_panel()
  p[, grade := predict(pd, score = score)]
  v <- scr_pd_validate(pd, p, id = "id", date = "date", default = "default", grade = "grade", score = "score", by = "quarter")
  expect_s3_class(v, "scr_pd_validation")
  expect_equal(v$n_cohorts, 8L)
  K <- nrow(pd$table)
  expect_equal(nrow(v$calibration), K)
  expect_true(all(c("p_jeffreys", "p_binomial", "critical", "z", "p_normal", "multi_period_z", "light") %in% names(v$calibration)))
  expect_true(all(v$calibration$light %in% c("red", "amber", "green")))
  # the pooled counts are the sums over cohorts and the cohort table matches scr_default_rate()
  dr <- scr_default_rate(p, "id", "date", "default", grade = "grade", by = "quarter", config = cfg_test())
  s <- dr$table[order(cohort, grade)]
  cc <- v$calibration_cohort[order(cohort, grade)]
  expect_equal(cc$n, s$n); expect_equal(cc$d, s$defaults)
  expect_equal(v$calibration$n, as.integer(tapply(s$n, as.integer(s$grade), sum)))
  pt <- v$portfolio_tests
  expect_equal(pt$n, sum(s$n)); expect_equal(pt$dr, sum(s$defaults) / sum(s$n))
  expect_equal(pt$p_jeffreys, stats::pbeta(pt$pd, pt$d + 0.5, pt$n - pt$d + 0.5))
  hl <- .pd_hl(v$calibration$n, v$calibration$d, v$calibration$pd)
  expect_equal(pt$hl_chi2, hl$chi2); expect_equal(pt$hl_df, K)
  dd <- v$portfolio$dr - v$portfolio$pd
  expect_equal(pt$multi_period_z, mean(dd) / (stats::sd(dd) / sqrt(8)))
  expect_equal(pt$brier, sum(v$calibration$d * (1 - v$calibration$pd)^2 + (v$calibration$n - v$calibration$d) * v$calibration$pd^2) / pt$n)
  # the tested PD is the final one by default
  expect_equal(v$calibration$pd, pd$table$pd_final)
  vb <- scr_pd_validate(pd, p, grade = "grade", pd_column = "pd_be", tests = c("jeffreys", "hl"))
  expect_equal(vb$calibration$pd, pd$table$pd_be)
  expect_null(vb$discrimination); expect_null(vb$stability$psi)
  expect_setequal(vb$summary$test, c("jeffreys", "jeffreys_grades_red", "hosmer_lemeshow"))
  # discrimination on the score with the S test against the hold-out AUC
  d <- v$discrimination
  expect_equal(d$basis, "score")
  expect_equal(d$auc_init, sc_demo()$metrics[sample == "holdout", auc])
  expect_equal(d$s_stat, (d$auc_init - d$auc) / d$se_auc)
  expect_equal(d$p_value, 1 - stats::pnorm(d$s_stat))
  expect_true(d$auc_lo <= d$auc && d$auc <= d$auc_hi)
  # without a score the discrimination is on the grade
  vg <- scr_pd_validate(pd, p, grade = "grade", tests = "auc")
  expect_equal(vg$discrimination$basis, "grade")
  expect_gt(vg$discrimination$auc, 0.6)
  # grading from the score column gives the same battery
  vs <- scr_pd_validate(pd, p, score = "score", by = "quarter")
  expect_equal(vs$calibration, v$calibration)
  # stability: PSI per cohort against the development grades, concentration test, migration
  expect_equal(nrow(v$stability$psi), 8L)
  expect_equal(v$stability$concentration$cv_init, pd$concentration$cv)
  expect_equal(length(v$stability$concentration$shares), K)
  expect_s3_class(v$stability$migration, "scr_migration")
  expect_equal(v$stability$migration$K, K)
  expect_true(all(v$summary$light[!is.na(v$summary$light)] %in% c("red", "amber", "green")))
  expect_true(v$light %in% c("red", "amber", "green"))
  expect_output(print(v), "overall light")
  expect_error(scr_pd_validate(pd, p), "grade")
  expect_error(scr_pd_validate(pd, p, grade = "nope"), "not found")
  expect_error(scr_pd_validate(pd, p, grade = "grade", horizon = 48L), "complete")
  expect_error(scr_pd_validate(list(), p), "scr_pd")
})

test_that("every PD function is identical under the serial and the PSOCK backend", {
  sc <- sc_demo(); p <- pd_panel(); dr <- pd_series()
  run <- function(backend) withr::with_options(list(scorecraft.parallel = backend), {
    s2 <- sc; s2$config$nthread <- 2L
    cal <- scr_calibrate(s2, target = 0.06, method = "logodds_ab")
    gr <- scr_grades(s2, cal, n_grades = 7L, min_obligors = 30L, min_defaults = 10L, dr = dr)
    gr <- scr_moc(gr, "C", method = "bootstrap", n_boot = 40, seed = 7)
    gr <- scr_moc(gr, "A", value = 0.001, reason = "x")
    pd <- scr_pd(gr)
    pp <- data.table::copy(p)[, grade := predict(pd, score = score)]
    v <- scr_pd_validate(pd, pp, grade = "grade", score = "score", n_boot = 20)
    list(cal = cal, gr = gr, pd = pd, v = v)
  })
  a <- run("serial"); b <- run("psock")
  expect_equal(a$cal$alignment, b$cal$alignment)
  expect_equal(a$gr$table, b$gr$table); expect_equal(a$gr$moc, b$gr$moc)
  expect_equal(a$pd$table, b$pd$table)
  expect_equal(a$v$calibration, b$v$calibration)
  expect_equal(a$v$discrimination, b$v$discrimination)
  expect_equal(a$v$summary, b$v$summary)
})

test_that("the workbook is written with every sheet and re-read", {
  skip_if_not_installed("openxlsx")
  pd <- pd_model(); p <- pd_panel()
  p[, grade := predict(pd, score = score)]
  v <- scr_pd_validate(pd, p, grade = "grade", score = "score")
  out <- file.path(tempdir(), "scr-pd-export"); unlink(out, recursive = TRUE)
  ex <- scr_export(pd, out, stamp = FALSE, validation = v)
  expect_setequal(names(ex$files), c("pd", "sql_pd"))
  expect_true(all(file.exists(unlist(ex$files))))
  sheets <- openxlsx::getSheetNames(ex$files$pd)
  expect_true(all(c("PD_Grades", "Master_Scale", "Calibration", "MoC_Ledger", "Floors", "Validation_Calibration",
                    "Validation_Discrimination", "Validation_Stability", "Migration", "Model_Card") %in% sheets))
  g <- openxlsx::read.xlsx(ex$files$pd, sheet = "PD_Grades")
  expect_equal(nrow(g), nrow(pd$table)); expect_equal(g$pd_final, pd$table$pd_final, tolerance = 1e-12)
  m <- openxlsx::read.xlsx(ex$files$pd, sheet = "MoC_Ledger")
  expect_equal(nrow(m), nrow(pd$moc_ledger))
  vc <- openxlsx::read.xlsx(ex$files$pd, sheet = "Validation_Calibration")
  expect_equal(vc$p_jeffreys, v$calibration$p_jeffreys, tolerance = 1e-12)
  mc <- openxlsx::read.xlsx(ex$files$pd, sheet = "Model_Card")
  expect_true("central_tendency" %in% mc$item)
  expect_true(any(grepl("AS pd_final", readLines(ex$files$sql_pd))))
  # without a validation the sheets carry an availability row, never a fabricated number
  ex0 <- scr_export(pd, out, stamp = FALSE)
  v0 <- openxlsx::read.xlsx(ex0$files$pd, sheet = "Validation_Calibration")
  expect_equal(v0$reason_code, "NO_VALIDATION_SUPPLIED")
})

test_that("the Hosmer-Lemeshow light is one-sided in practice and the grades print says where the PD came from", {
  pd <- pd_model()
  p  <- pd_panel()
  v_final <- scr_pd_validate(pd, p, score = "score", pd_column = "pd_final", tests = c("jeffreys", "hl"))
  hl <- v_final$summary[test == "hosmer_lemeshow"]
  g  <- v_final$calibration
  conservative <- all(g$d <= g$n * g$pd + 1e-9, na.rm = TRUE)
  if (conservative) expect_equal(hl$light, "green") else expect_equal(hl$light, .pd_light(hl$p_value, v_final$lights))
  v_be <- scr_pd_validate(pd, p, score = "score", pd_column = "pd_be", tests = c("hl"))
  expect_true(is.finite(v_be$summary[test == "hosmer_lemeshow", statistic]))
  gr <- pd_grades()
  expect_output(print(gr), if (is.null(gr$dr)) "sample default rate" else "PD source: lra")
})
