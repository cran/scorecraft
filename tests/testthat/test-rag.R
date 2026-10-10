# Red / amber / green lights: every rule on hand-made inputs, the roll-up,
# the "grey" rule, and the scorecard, data.frame and study paths.

rag_df <- function(n = 8000, seed = 25, shift = 0, slope = 1.1) {
  set.seed(seed)
  smp <- rep(c("dev", "new"), each = n / 2)
  x <- stats::rnorm(n) + ifelse(smp == "new", shift, 0)
  p <- stats::plogis(-1.8 - slope * x)
  data.frame(score = round(600 + 40 * x, 1), y = stats::rbinom(n, 1, p), prob = stats::plogis(-1.8 - 1.1 * x),
             smp = smp, seg = sample(c("A", "B"), n, TRUE), stringsAsFactors = FALSE)
}

plan_row <- function(m, objective = "risk") {
  pl <- data.table::as.data.table(scr_rag_plan(objective))
  pl[match(m, pl[["metric"]])]
}

test_that("each rule turns a value into the documented light", {
  # deviation convention (gini_ratio, higher better): lit on the upper bound,
  # amber or red only when the whole interval lies below a threshold
  r <- plan_row("gini_ratio")
  expect_identical(.rag_rule(r, 1.00, 0.96, 1.04)$light, "green")
  expect_identical(.rag_rule(r, 0.97, 0.93, 1.01)$light, "green")    # consistent with 0.95
  expect_identical(.rag_rule(r, 0.80, 0.70, 0.95)$light, "green")    # the upper bound reaches 0.95
  expect_identical(.rag_rule(r, 0.85, 0.80, 0.93)$light, "amber")    # below 0.95, not below 0.90
  expect_identical(.rag_rule(r, 0.80, 0.70, 0.89)$light, "red")
  expect_match(.rag_rule(r, 0.85, 0.80, 0.93)$reason, "upper bound 0.93 below 0.95, not below 0.9")
  expect_identical(.rag_rule(r, NA)$light, "grey")
  expect_identical(.rag_rule(plan_row("gini_ratio", "propensity"), 0.85, 0.81, 0.89)$light, "amber")
  # the mirror image when lower is better: lit on the lower bound
  lb <- data.table::copy(r)[, `:=`(higher_better = FALSE, green = 0.1, red = 0.2)]
  expect_identical(.rag_rule(lb, 0.15, 0.05, 0.25)$light, "green")
  expect_identical(.rag_rule(lb, 0.18, 0.12, 0.25)$light, "amber")
  expect_identical(.rag_rule(lb, 0.30, 0.22, 0.40)$light, "red")
  expect_match(.rag_rule(lb, 0.30, 0.22, 0.40)$reason, "lower bound 0.22 > 0.2")
  # effect and significance (score PSI)
  p <- plan_row("score_psi")
  expect_identical(.rag_rule(p, 0.30, bench = 0.01)$light, "red")
  expect_identical(.rag_rule(p, 0.30, bench = 0.50)$light, "green")    # large but not significant
  expect_identical(.rag_rule(p, 0.15, bench = 0.01)$light, "amber")
  expect_identical(.rag_rule(p, 0.05, bench = 0.001)$light, "green")   # significant but small
  expect_identical(.rag_rule(p, 0.30, bench = NA)$light, "green")
  # p-value lights
  q <- plan_row("auc_change_p")
  expect_identical(.rag_rule(q, 0.005)$light, "red")
  expect_identical(.rag_rule(q, 0.01)$light, "red")
  expect_identical(.rag_rule(q, 0.03)$light, "amber")
  expect_identical(.rag_rule(q, 0.20)$light, "green")
  # count rule
  k <- plan_row("rank_order")
  expect_identical(vapply(c(0, 1, 2, 5), function(v) .rag_rule(k, v)$light, ""), c("green", "amber", "red", "red"))
  f <- plan_row("woe_sign_flip")
  expect_identical(.rag_rule(f, 9)$light, "amber")   # never red
  # O/E under propensity: two-sided interval
  o <- plan_row("oe_ratio", "propensity")
  expect_identical(.rag_rule(o, 1.00, 0.95, 1.05)$light, "green")
  expect_identical(.rag_rule(o, 1.20, 1.12, 1.28)$light, "amber")
  expect_identical(.rag_rule(o, 1.40, 1.30, 1.50)$light, "red")
  expect_identical(.rag_rule(o, 0.70, 0.60, 0.79)$light, "red")
  expect_identical(.rag_rule(o, 1.30, 1.05, 1.55)$light, "green")    # consistent with the tolerance
  # O/E under risk: one-sided, lit on the lower bound (over-prediction is prudent)
  ok <- plan_row("oe_ratio")
  expect_identical(.rag_rule(ok, 0.70, 0.60, 0.79)$light, "green")
  expect_identical(.rag_rule(ok, 1.00, 0.95, 1.05)$light, "green")
  expect_identical(.rag_rule(ok, 1.20, 1.12, 1.28)$light, "amber")
  expect_identical(.rag_rule(ok, 1.40, 1.30, 1.50)$light, "red")
  # threshold rule on the value alone
  t <- plan_row("iv_ratio")
  expect_identical(vapply(c(0.9, 0.6, 0.4), function(v) .rag_rule(t, v)$light, ""), c("green", "amber", "red"))
  expect_identical(.rag_rule(plan_row("ks"), 0.4)$light, "none")
})

test_that("lights roll up worst-first, and 'grey' never becomes green", {
  expect_identical(.rag_worst(c("green", "amber", "red")), "red")
  expect_identical(.rag_worst(c("green", "grey")), "green")
  expect_identical(.rag_worst(c("grey", "grey")), "grey")
  expect_identical(.rag_worst(c("none", NA)), "grey")
  expect_identical(.rag_worst(character()), "grey")
  expect_identical(.rag_light(c(1, 0.92, 0.5, NA), green = 0.95, red = 0.9), c("green", "amber", "red", "grey"))
  expect_identical(.rag_light(0.2, green = 0.1, red = 0.25, higher_better = FALSE), "amber")
  expect_identical(.rag_light(c(0.9, 0.9, 0.8), c(0.85, 0.85, 0.7), c(0.96, 0.92, 0.85), green = 0.95, red = 0.9),
                   c("green", "amber", "red"))
  expect_identical(.rag_light(0.2, 0.05, 0.3, green = 0.1, red = 0.25, higher_better = FALSE), "green")
})

test_that("a stable sample is green; a shifted one raises stability, which caps the overall at amber", {
  stable <- scr_rag(rag_df(), prob = "prob", sample = "smp", n_boot = 100, seed = 1)
  expect_identical(stable$summary$overall, "green")
  expect_identical(stable$reference, "dev")
  gr <- stable$table[metric == "gini_ratio"]
  expect_identical(gr$light, if (gr$hi >= 0.95) "green" else if (gr$hi < 0.9) "red" else "amber")
  # the same model on a shifted population: large PSI, discrimination and calibration intact
  sh <- scr_rag(rag_df(shift = 0.8), prob = "prob", sample = "smp", n_boot = 100, seed = 1)
  s <- sh$summary
  expect_identical(s$stability, "red")
  expect_identical(s$overall, "amber")
  psi <- sh$table[metric == "score_psi"]
  expect_gt(psi$value, 0.25); expect_gt(psi$value, psi$benchmark)
  # a weaker model on the new sample: discrimination red, overall red
  weak <- scr_rag(rbind(rag_df()[1:4000, ], rag_df(slope = 0.3)[4001:8000, ]), prob = "prob", sample = "smp",
                  n_boot = 100, seed = 1)
  expect_identical(weak$summary$discrimination, "red")
  expect_identical(weak$summary$overall, "red")
})

test_that("the values follow independent formulas", {
  d <- rag_df()
  rg <- scr_rag(d, prob = "prob", sample = "smp", n_boot = 200, seed = 3)
  ref <- d[d$smp == "dev", ]; new <- d[d$smp == "new", ]
  m_r <- scr_metrics(ref$score, ref$y, higher_is_event = FALSE, ci = FALSE)
  m_n <- scr_metrics(new$score, new$y, higher_is_event = FALSE, ci = FALSE)
  t <- rg$table
  expect_equal(t[metric == "gini_ratio", value], m_n$gini / m_r$gini)
  expect_equal(t[metric == "gini_ratio", benchmark], m_r$gini)
  se <- sqrt(ref_auc_se(new$score, new$y, higher_is_event = FALSE)^2 + ref_auc_se(ref$score, ref$y, higher_is_event = FALSE)^2)
  expect_equal(t[metric == "auc_change_p", value], stats::pnorm((m_r$auc - m_n$auc) / se, lower.tail = FALSE))
  expect_equal(t[metric == "ks", value], m_n$ks)
  expect_identical(t[metric == "ks", light], "none")
  expect_equal(t[metric == "oe_ratio", value], sum(new$y) / sum(new$prob))
  expect_equal(t[metric == "oe_ratio", lo], stats::qbeta(0.025, sum(new$y) + 0.5, nrow(new) - sum(new$y) + 0.5) * nrow(new) / sum(new$prob))
  # under risk, O/E is lit on its lower bound
  oe <- t[metric == "oe_ratio"]
  expect_identical(oe$light, if (oe$lo <= 1.10) "green" else if (oe$lo > 1.25) "red" else "amber")
  # score PSI over the left-closed bands frozen on the reference
  ps <- scr_psi(findInterval(ref$score, rg$cuts), findInterval(new$score, rg$cuts), levels = 0:length(rg$cuts))
  expect_equal(t[metric == "score_psi", value], ps$psi)
  expect_equal(t[metric == "score_psi", benchmark], ps$critical)
  # band calibration: per-band rows without a light, then the count
  bc <- t[metric == "band_calibration"]
  expect_identical(nrow(bc), length(rg$cuts) + 2L)
  expect_true(all(bc[level != "score", light] == "none"))
  band <- findInterval(new$score, rg$cuts) + 1L
  pd <- tapply(new$prob, band, mean); dr <- tapply(new$y, band, sum); nn <- tapply(new$y, band, length)
  p1 <- stats::pbeta(pd, dr + 0.5, nn - dr + 0.5)
  expect_equal(bc[level == "score", value], sum(p1 <= 0.01))
  expect_equal(sort(bc[level != "score", benchmark]), sort(as.numeric(pd)))
})

test_that("too few events turn every light 'grey'; a single class too", {
  d <- rag_df(2000)
  g <- scr_rag(d, prob = "prob", sample = "smp", min_events = 5000, n_boot = 20, seed = 1)
  lit <- g$table[light != "none"]
  expect_true(all(lit$light == "grey"))
  expect_true(all(grepl("too few events", lit$reason)))
  expect_identical(unlist(g$summary[, c("discrimination", "calibration", "stability", "overall")], use.names = FALSE),
                   rep("grey", 4))
  expect_true(is.finite(g$table[metric == "gini_ratio", value]))   # values are still reported
  one <- d; one$y[one$smp == "new"] <- 0
  o <- scr_rag(one, sample = "smp", n_boot = 20, seed = 1)
  expect_true(is.na(o$table[metric == "auc_change_p", value]))
  expect_identical(o$summary$overall, "grey")
  # no expected probability: calibration is "grey", with the reason
  np <- scr_rag(d, sample = "smp", n_boot = 20, seed = 1)
  expect_identical(np$summary$calibration, "grey")
  expect_true(all(np$table[family == "calibration", reason] == "no expected probability"))
})

test_that("calibration is one-sided under risk: a conservative model is green, and red under propensity", {
  set.seed(41)
  n <- 6000
  x <- stats::rnorm(n)
  p <- stats::plogis(-1.8 - 1.1 * x)
  d <- data.frame(score = round(600 + 40 * x, 1), y = stats::rbinom(n, 1, p), cons = pmin(1, 1.7 * p), opt = p / 1.7)
  # conservative (O/E near 0.59): prudent under risk, so green
  rk <- scr_rag(d, prob = "cons", n_boot = 0)
  oe <- rk$table[metric == "oe_ratio"]
  expect_lt(oe$hi, 0.8)
  expect_identical(oe$light, "green")
  expect_identical(rk$table[metric == "band_calibration" & level == "score", light], "green")
  expect_identical(rk$summary$calibration, "green")
  # optimistic (O/E near 1.7): under-prediction is red under risk, lit on the lower bound
  ro <- scr_rag(d, prob = "opt", n_boot = 0)$table[metric == "oe_ratio"]
  expect_gt(ro$lo, 1.25)
  expect_identical(ro$light, "red")
  expect_match(ro$reason, "lower bound")
  # propensity: both directions count, so the conservative model is red
  # (the score is mirrored so that a higher score means more events)
  pp <- scr_rag(transform(d, score = -score), prob = "cons", objective = "propensity", n_boot = 0)
  expect_identical(pp$table[metric == "oe_ratio", light], "red")
  expect_identical(pp$summary$calibration, "red")
  # the plans carry the convention
  pr <- scr_rag_plan("risk"); pq <- scr_rag_plan("propensity")
  expect_identical(pr$rule[pr$metric == "oe_ratio"], "ci")
  expect_false(pr$higher_better[pr$metric == "oe_ratio"])
  expect_equal(c(pr$green[pr$metric == "oe_ratio"], pr$red[pr$metric == "oe_ratio"]), c(1.10, 1.25))
  expect_identical(pq$rule[pq$metric == "oe_ratio"], "interval")
})

test_that("the overall light says when calibration was not tested", {
  rg <- scr_rag(rag_df(), sample = "smp", n_boot = 100, seed = 1)
  expect_identical(rg$summary$discrimination, "green")
  expect_identical(rg$summary$calibration, "grey")
  expect_identical(rg$summary$overall, "green")
  expect_match(rg$summary$reason, "calibration not tested")
  expect_output(print(rg), "calibration not tested")
  # with an expected probability nothing needs saying
  expect_identical(scr_rag(rag_df(), prob = "prob", sample = "smp", n_boot = 100, seed = 1)$summary$reason, "")
  # raised by stability; too few events
  sh <- scr_rag(rag_df(shift = 0.8), prob = "prob", sample = "smp", n_boot = 100, seed = 1)
  expect_match(sh$summary$reason, "raised to amber by stability")
  few <- scr_rag(rag_df(2000), prob = "prob", sample = "smp", min_events = 5000, n_boot = 0)
  expect_identical(few$summary$reason, "too few events")
})

test_that("the two-sample S-test keeps its size; a fixed reference AUC would not", {
  # H0: the reference and the study samples come from the same model. With
  # R replications the rejection rate at 5% has standard error
  # sqrt(0.05 * 0.95 / R) = 0.0049 for R = 2000; three of them give [0.035, 0.065].
  set.seed(7)
  R <- 2000; n <- 500
  cells <- function() {
    x <- round(stats::rnorm(n) * 4)
    y <- stats::rbinom(n, 1, stats::plogis(-1 + 0.25 * x))
    idx <- x + 21L   # scores -20..20 on fixed cells
    e <- tabulate(idx[y == 1L], 41L); m <- tabulate(idx, 41L)
    list(e = e, n_y = m, e_raw = e, n_y_raw = m)
  }
  sims <- replicate(R, {
    hs <- cells(); hr <- cells()
    st <- .rag_auc_change(hs, hr, "higher_is_riskier")
    # the one-sample form: the reference AUC taken as fixed
    c1 <- hs$e; c0 <- hs$n_y - hs$e
    se1 <- .study_delong_counts(c1, c0)
    c(two = st$p < 0.05, one = stats::pnorm((st$auc_ref - st$auc_study) / se1, lower.tail = FALSE) < 0.05)
  })
  size <- rowMeans(sims)
  expect_gt(size[["two"]], 0.035); expect_lt(size[["two"]], 0.065)
  expect_gt(size[["one"]], 0.08)   # about 12% in theory with equal sample sizes
  # se is the two-sample one
  hs <- cells(); hr <- cells()
  st <- .rag_auc_change(hs, hr, "higher_is_riskier")
  expect_equal(st$se, sqrt(.study_delong_counts(hs$e, hs$n_y - hs$e)^2 + .study_delong_counts(hr$e, hr$n_y - hr$e)^2))
})

test_that("by groups: one set of lights per group from one table, same group of the reference", {
  d <- rag_df()
  rg <- scr_rag(d, prob = "prob", sample = "smp", by = "seg", n_boot = 50, seed = 2)
  expect_identical(rg$summary$group, c("A", "B"))
  expect_identical(unique(rg$table$group), c("A", "B"))
  # the segment of the new sample is read against the same segment of the reference
  a_ref <- d[d$smp == "dev" & d$seg == "A", ]; a_new <- d[d$smp == "new" & d$seg == "A", ]
  g_r <- scr_metrics(a_ref$score, a_ref$y, higher_is_event = FALSE, ci = FALSE)$gini
  g_n <- scr_metrics(a_new$score, a_new$y, higher_is_event = FALSE, ci = FALSE)$gini
  expect_equal(rg$table[group == "A" & metric == "gini_ratio", value], g_n / g_r)
  # without a sample column every group is read against the whole data
  whole <- scr_rag(d, prob = "prob", by = "seg", n_boot = 50, seed = 2)
  g_all <- scr_metrics(d$score, d$y, higher_is_event = FALSE, ci = FALSE)$gini
  a <- d[d$seg == "A", ]
  expect_equal(whole$table[group == "A" & metric == "gini_ratio", value],
               scr_metrics(a$score, a$y, higher_is_event = FALSE, ci = FALSE)$gini / g_all)
  # a period absent from the reference is read against the whole reference
  d$per <- ifelse(d$smp == "dev", "2025", ifelse(seq_len(nrow(d)) %% 2 == 0, "2026-1", "2026-2"))
  pr <- scr_rag(d, prob = "prob", sample = "smp", by = "per", n_boot = 50, seed = 2)
  expect_identical(pr$summary$group, c("2026-1", "2026-2"))
  ref <- d[d$smp == "dev", ]; p1 <- d[d$per == "2026-1", ]
  expect_equal(pr$table[group == "2026-1" & metric == "gini_ratio", value],
               scr_metrics(p1$score, p1$y, higher_is_event = FALSE, ci = FALSE)$gini /
                 scr_metrics(ref$score, ref$y, higher_is_event = FALSE, ci = FALSE)$gini)
})

test_that("scr_rag on a scorecard covers the variables and reuses the scorecard's numbers", {
  sc <- sc_demo()
  rg <- scr_rag(sc, n_boot = 50, seed = 1)
  expect_s3_class(rg, "scr_rag")
  expect_true(all(c("discrimination", "calibration", "stability", "variables") %in% rg$table$family))
  t <- rg$table
  m <- sc$metrics
  expect_equal(t[metric == "gini_ratio", value], m[sample == "holdout", gini] / m[sample == "train", gini])
  ho <- sc$samples$holdout
  expect_equal(t[metric == "oe_ratio", value], sum(ho$y) / sum(.score_to_prob(sc$alignment, ho$score)))
  # CSI per variable equals the stored CSI of the scorecard
  csi <- t[metric == "csi"]
  expect_setequal(csi$level, sc$features)
  sv <- sc$stability$variables
  expect_equal(csi$value[match(sv$variable, csi$level)], sv$csi)
  # IV ratio with the same smoothing on both samples
  f <- sc$features[1]
  r <- sc$fit$results[[f]]
  idx <- sc$holdout_bins[[f]]
  ev <- tabulate(idx[ho$y == 1], length(r$bin)); cnt <- tabulate(idx, length(r$bin))
  iv <- function(bw) sum((bw$pct_event - bw$pct_nonevent) * bw$log_odds)
  expect_equal(t[metric == "iv_ratio" & level == f, value],
               iv(.band_woe(ev, cnt - ev)) / iv(.band_woe(r$count_pos, r$count_neg)))
  expect_output(print(rg), "scr_rag")
  # by period of the scored samples
  rb <- scr_rag(sc, by = "date", n_boot = 20, seed = 1)
  expect_gt(nrow(rb$summary), 1L)
  expect_true(all(rb$summary$group %in% as.character(unique(ho$date))))
  # variables need train against holdout
  rs <- scr_rag(sc, sample = "train", reference = "holdout", n_boot = 20, seed = 1)
  expect_true(all(rs$table[family == "variables", light] == "grey"))
})

test_that("scr_rag on a score study reuses its cuts and its count table", {
  d <- rag_df()
  b <- scr_bands(d, sample = "smp", n_bands = 10, n_boot = 0)
  rg <- scr_rag(b, n_boot = 50, seed = 1)
  expect_equal(rg$cuts, b$cuts)
  expect_equal(rg$table[metric == "score_psi", value], b$summary[sample == "new", psi])
  expect_equal(rg$table[metric == "rank_order", value], b$summary[sample == "new", reversals])
  expect_identical(rg$summary$calibration, "grey")
  # a study of a scorecard keeps the expected probability
  sc <- sc_demo()
  rs <- scr_rag(scr_tiers(sc, n_tiers = 3), n_boot = 20, seed = 1)
  expect_false(identical(rs$summary$calibration, "grey"))
  expect_error(scr_rag(b, sample = "zzz"), "not in the study")
})

test_that("an edited plan drives the lights; a malformed plan is refused", {
  d <- rag_df()
  pl <- scr_rag_plan()
  pl$green[pl$metric == "score_psi"] <- 0
  pl$red[pl$metric == "score_psi"] <- 0
  pl <- pl[pl$metric != "ks", ]
  rg <- scr_rag(d, prob = "prob", sample = "smp", plan = pl, n_boot = 20, seed = 1)
  expect_false("ks" %in% rg$table$metric)
  psi <- rg$table[metric == "score_psi"]
  expect_identical(psi$light, if (psi$value > psi$benchmark) "red" else "green")
  # the interval lights follow the edited thresholds too: a Gini ratio policy
  # above any reachable value turns red, a lenient one green
  gp <- scr_rag_plan()
  gp$green[1] <- 2; gp$red[1] <- 1.9
  hi <- scr_rag(d, prob = "prob", sample = "smp", plan = gp, n_boot = 50, seed = 1)$table[metric == "gini_ratio"]
  expect_lt(hi$hi, 1.9)
  expect_identical(hi$light, "red")
  gp$green[1] <- 0.5; gp$red[1] <- 0.4
  expect_identical(scr_rag(d, prob = "prob", sample = "smp", plan = gp, n_boot = 50, seed = 1)$table[metric == "gini_ratio", light],
                   "green")
  bad <- scr_rag_plan(); bad$rule[1] <- "magic"
  expect_error(scr_rag(d, sample = "smp", plan = bad), "rule")
  bad <- scr_rag_plan(); bad$metric[1] <- "auc_ratio"
  expect_error(scr_rag(d, sample = "smp", plan = bad), "unknown check")
  expect_error(scr_rag(d, sample = "smp", plan = scr_rag_plan()[, -2]), "lacks")
  expect_identical(scr_rag_plan("propensity")$green[1], 0.90)
})

test_that("scr_export writes the lights, alone or with a study", {
  skip_if_not_installed("openxlsx")
  old <- scr_verbose(FALSE); on.exit(scr_verbose(old), add = TRUE)
  d <- rag_df(2000)
  b <- scr_bands(d, sample = "smp", n_bands = 5, n_boot = 0)
  rg <- scr_rag(d, prob = "prob", sample = "smp", n_boot = 20, seed = 1)
  out <- file.path(tempdir(), "scr-rag-export")
  unlink(out, recursive = TRUE)
  ex <- scr_export(rg, out, stamp = FALSE)
  expect_identical(openxlsx::getSheetNames(ex$files$xlsx), c("Lights", "Light_Summary", "Light_Plan", "Settings"))
  eb <- scr_export(b, out, stamp = FALSE, rag = rg)
  expect_true(all(c("Bands", "Lights", "Light_Summary") %in% openxlsx::getSheetNames(eb$files$xlsx)))
  expect_error(scr_export(b, out, stamp = FALSE, rag = b), "scr_rag")
})

test_that("the groups of a numeric `by` are listed in numeric order", {
  d <- rag_df(6000)
  set.seed(3)
  d$k <- sample(c(9, 10, 100), nrow(d), TRUE)
  rg <- scr_rag(d, prob = "prob", sample = "smp", by = "k", n_boot = 0)
  expect_identical(rg$summary$group, c("9", "10", "100"))
  expect_identical(unique(rg$table$group), c("9", "10", "100"))
  # text labels keep their order; the checks of each group are the same
  d$kc <- as.character(d$k)
  rc <- scr_rag(d, prob = "prob", sample = "smp", by = "kc", n_boot = 0)
  expect_identical(rc$summary$group, c("10", "100", "9"))
  for (g in c("9", "10", "100")) {
    expect_equal(rg$table[group == g, value], rc$table[group == g, value])
    expect_identical(rg$summary[group == g, overall], rc$summary[group == g, overall])
  }
  # without a sample column, and for the dates of a scorecard, the order is unchanged
  expect_identical(scr_rag(d, prob = "prob", by = "k", n_boot = 0)$summary$group, c("9", "10", "100"))
  rb <- scr_rag(sc_demo(), by = "date", n_boot = 0)
  expect_identical(rb$summary$group, sort(rb$summary$group))
})
