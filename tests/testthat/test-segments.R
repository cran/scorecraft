# One score on many segments: the AUC and its DeLong error against the
# row-level estimators, the indirect standardization and the slope against
# manual computations and glm(), the chi-square test and the actions.

seg_df <- function(n = 5000, seed = 81) {
  set.seed(seed)
  seg <- sample(c("app", "store", "web"), n, TRUE, c(0.5, 0.3, 0.2))
  x <- stats::rnorm(n)
  data.frame(seg = seg, score = round(600 + 50 * x),
             y = stats::rbinom(n, 1, stats::plogis(-2 - x + 0.5 * (seg == "web"))),
             w = stats::runif(n, 0.2, 3), per = sample(c("h1", "h2"), n, TRUE), stringsAsFactors = FALSE)
}

test_that("AUC, KS and the DeLong error equal the row-level estimators", {
  d <- seg_df()
  sg <- scr_segments(d, segment = "seg", level = 0.9)
  expect_s3_class(sg, "scr_segments")
  t <- sg$table
  expect_identical(t$segment, c("app", "store", "web")); expect_identical(unique(t$group), "all")
  z <- stats::qnorm(0.95)
  for (s in t$segment) {
    r <- d[d$seg == s, ]
    row <- t[segment == s]
    # higher_is_safer: the events sit at the low scores
    m <- .auc_ks(-r$score, r$y)
    expect_equal(row$auc, m$auc); expect_equal(row$ks, m$ks); expect_equal(row$gini, 2 * m$auc - 1)
    se <- ref_auc_se(r$score, r$y, higher_is_event = FALSE)
    expect_equal(row$auc_se, se)
    expect_equal(row$auc_lo, m$auc - z * se); expect_equal(row$auc_hi, m$auc + z * se)
    expect_equal(row$n, nrow(r)); expect_equal(row$events, sum(r$y)); expect_equal(row$rate, mean(r$y))
  }
  # the pooled rows
  p <- sg$pooled
  expect_equal(p$auc, .auc_ks(-d$score, d$y)$auc)
  expect_equal(p$auc_se, ref_auc_se(d$score, d$y, higher_is_event = FALSE))
  expect_equal(p$n, nrow(d)); expect_equal(p$rate, mean(d$y))
  # a score read the other way has the mirrored AUC and the same error
  f <- d; f$score <- -f$score
  tf <- scr_segments(f, segment = "seg", direction = "higher_is_riskier")$table
  expect_equal(tf$auc, t$auc); expect_equal(tf$auc_se, t$auc_se)
})

test_that("the expected events, the O/E ratio and the offset equal the manual standardization", {
  d <- seg_df()
  sg <- scr_segments(d, segment = "seg", n_bands = 8)
  cuts <- sg$cuts
  B <- length(cuts) + 1L
  ib <- factor(findInterval(d$score, cuts) + 1L, seq_len(B))
  Rb <- as.numeric(tapply(d$y, ib, mean))
  for (s in sg$table$segment) {
    i <- d$seg == s
    nsb <- as.numeric(table(ib[i]))
    expd <- sum(nsb * Rb)
    row <- sg$table[segment == s]
    expect_equal(row$expected, expd)
    x <- sum(d$y[i]); n <- sum(i)
    expect_equal(row$oe_ratio, x / expd)
    expect_equal(row$oe_lo, stats::qbeta(0.025, x + 0.5, n - x + 0.5) * n / expd)
    expect_equal(row$oe_hi, stats::qbeta(0.975, x + 0.5, n - x + 0.5) * n / expd)
    expect_equal(row$offset, stats::qlogis(x / n) - stats::qlogis(expd / n))
    # the PSI of the segment against the pooled rows, over the same bands
    ps <- scr_psi(d$score, d$score[i], breaks = c(-Inf, cuts, Inf))
    expect_equal(row$psi, ps$psi)
    # the segment is part of the pooled rows: 1/n_s - 1/N, not 1/n_s + 1/N
    expect_equal(row$psi_critical, (1 / n - 1 / nrow(d)) * stats::qchisq(0.95, B - 1))
  }
  # the expected events of the segments add up to the events
  expect_equal(sum(sg$table$expected), sum(d$y))
  # the pooled bands are the equal-share bands of scr_bands() on all rows
  expect_equal(cuts, scr_bands(d, n_bands = 8, n_boot = 0)$cuts)
})

test_that("the slope equals glm() on the rows, with and without weights", {
  d <- seg_df()
  sg <- scr_segments(d, segment = "seg")
  b_all <- unname(stats::coef(stats::glm(y ~ score, stats::binomial(), d))[2])
  expect_equal(sg$pooled$slope, b_all, tolerance = 1e-7)
  for (s in sg$table$segment) {
    b <- unname(stats::coef(stats::glm(y ~ score, stats::binomial(), d[d$seg == s, ]))[2])
    expect_equal(sg$table[segment == s, slope], b, tolerance = 1e-7)
    expect_equal(sg$table[segment == s, slope_ratio], b / b_all, tolerance = 1e-7)
  }
  sw <- scr_segments(d, segment = "seg", weight = "w")
  bw_all <- unname(stats::coef(suppressWarnings(stats::glm(y ~ score, stats::binomial(), d, weights = w)))[2])
  bw <- unname(stats::coef(suppressWarnings(stats::glm(y ~ score, stats::binomial(), d[d$seg == "web", ], weights = w)))[2])
  expect_equal(sw$pooled$slope, bw_all, tolerance = 1e-7)
  expect_equal(sw$table[segment == "web", slope_ratio], bw / bw_all, tolerance = 1e-7)
  # weighted rates and AUC; the error uses the unweighted class counts
  r <- d[d$seg == "web", ]
  expect_equal(sw$table[segment == "web", rate], stats::weighted.mean(r$y, r$w))
  idx <- data.table::frank(-r$score, ties.method = "dense")
  c1 <- as.numeric(tapply(r$w * r$y, factor(idx, seq_len(max(idx))), sum)); c1[is.na(c1)] <- 0
  c0 <- as.numeric(tapply(r$w * (1 - r$y), factor(idx, seq_len(max(idx))), sum)); c0[is.na(c0)] <- 0
  expect_equal(sw$table[segment == "web", auc], .auc_ks_counts(c1, c0)$auc)
  expect_equal(sw$table[segment == "web", auc_se], .study_delong_counts(c1, c0, sum(r$y), sum(1 - r$y)))
  # the Wald test of each slope against the slope of the other rows, from glm() on the rows
  wald <- function(dd, i, wt = NULL) {
    fit <- function(r) summary(suppressWarnings(stats::glm(y ~ score, stats::binomial(), r,
                                                         weights = if (is.null(wt)) NULL else r[[wt]],
                                                         control = stats::glm.control(epsilon = 1e-12))))$coefficients
    a <- fit(dd[i, ]); b <- fit(dd[!i, ])
    list(z = (a[2, 1] - b[2, 1]) / sqrt(a[2, 2]^2 + b[2, 2]^2), se = a[2, 2])
  }
  ps <- vapply(sg$table$segment, function(s) 2 * stats::pnorm(-abs(wald(d, d$seg == s)$z)), numeric(1))
  expect_equal(sg$table$p_slope, unname(ps), tolerance = 1e-6)
  expect_equal(sg$table$p_slope_adj, stats::p.adjust(sg$table$p_slope, "holm"))
  # the variance of the count fit is that of glm() on the rows
  cl <- .study_collapse(.study_hist(d$score[d$seg == "web"], d$y[d$seg == "web"]))
  sv <- .seg_slope(cl, mean(d$score), stats::sd(d$score))
  expect_equal(sqrt(sv$var), wald(d, d$seg == "web")$se, tolerance = 1e-6)
  # equal weights leave the test as it is: the variance runs on the Kish effective size
  d$c3 <- 3
  expect_equal(scr_segments(d, segment = "seg", weight = "c3")$table$p_slope, sg$table$p_slope, tolerance = 1e-6)
  # classes separated by the score, or a single score value: no slope
  sep <- data.frame(seg = "a", score = c(1, 2, 3, 4, 5, 6), y = c(1, 1, 1, 0, 0, 0))
  expect_true(is.na(scr_segments(sep, segment = "seg")$table$slope))
  expect_true(is.na(scr_segments(sep, segment = "seg")$table$p_slope))
  one <- data.frame(seg = "a", score = 5, y = c(1, 0, 0, 1))
  expect_true(is.na(scr_segments(one, segment = "seg")$table$slope))
})

test_that("the test of equal AUC is the inverse-variance chi-square", {
  d <- seg_df()
  sg <- scr_segments(d, segment = "seg")
  auc <- se <- numeric(3)
  sn <- c("app", "store", "web")
  for (k in 1:3) {
    r <- d[d$seg == sn[k], ]
    auc[k] <- .auc_ks(-r$score, r$y)$auc; se[k] <- ref_auc_se(r$score, r$y, higher_is_event = FALSE)
  }
  w <- 1 / se^2
  aw <- sum(w * auc) / sum(w)
  Q <- sum((auc - aw)^2 / se^2)
  te <- sg$test
  expect_equal(te$auc_w, aw); expect_equal(te$statistic, Q); expect_identical(te$df, 2L)
  expect_equal(te$p_value, stats::pchisq(Q, 2, lower.tail = FALSE)); expect_identical(te$segments, 3L)
  t <- sg$table
  expect_equal(t$auc_diff, auc - aw)
  # each segment against the mean that contains it: Var = se^2 - 1 / sum(w)
  pz <- 2 * stats::pnorm(-abs(auc - aw) / sqrt(se^2 - 1 / sum(w)))
  expect_equal(t$p_auc, pz); expect_equal(t$p_auc_adj, stats::p.adjust(pz, "holm"))
  # with two segments the two differences carry the same evidence as the chi-square
  two <- scr_segments(d[d$seg != "web", ], segment = "seg")
  expect_equal(two$table$p_auc, rep(two$test$p_value, 2))
  expect_equal(two$test$p_value, 2 * stats::pnorm(-abs(auc[1] - auc[2]) / sqrt(se[1]^2 + se[2]^2)))
  # one test told twice: no Holm adjustment, for the AUC and for the slope
  expect_identical(two$table$p_auc_adj, two$table$p_auc)
  expect_equal(two$table$p_slope[1], two$table$p_slope[2])
  expect_identical(two$table$p_slope_adj, two$table$p_slope)
  # a third segment left out of the AUC test (a single class) leaves two tests that are one
  z <- d; z$y[z$seg == "web"] <- 0L
  tz <- scr_segments(z, segment = "seg")$table
  expect_identical(tz$p_auc_adj[1:2], tz$p_auc[1:2]); expect_equal(tz$p_auc[1], tz$p_auc[2])
})

test_that("the actions on constructed segments", {
  set.seed(5)
  n <- 20000
  seg <- sample(c("base", "shifted", "weak", "tiny"), n, TRUE, c(0.696, 0.20, 0.10, 0.004))
  x <- stats::rnorm(n)
  lin <- -2 - x
  lin[seg == "shifted"] <- lin[seg == "shifted"] + 0.8                       # same ranking, higher level
  lin[seg == "weak"] <- -1.9 - 0.2 * x[seg == "weak"]                        # the score barely ranks
  d <- data.frame(seg = seg, score = round(600 + 50 * x), y = stats::rbinom(n, 1, stats::plogis(lin)))
  sg <- scr_segments(d, segment = "seg", min_events = 30)
  act <- stats::setNames(sg$table$action, sg$table$segment)
  expect_identical(act[["tiny"]], "too few events")
  expect_identical(act[["weak"]], "separate model")
  expect_identical(act[["shifted"]], "offset")
  expect_identical(act[["base"]], "shared")
  t <- sg$table
  # each action follows its documented rule
  expect_lt(t[segment == "tiny", events], 30)
  wk <- t[segment == "weak"]
  expect_true((abs(wk$auc_diff) > 0.03 && wk$p_auc_adj < 0.05) ||
                (abs(wk$slope_ratio - 1) > 0.25 && wk$p_slope_adj < 0.05))
  sh <- t[segment == "shifted"]
  expect_true(abs(sh$offset) > 0.25 && sh$oe_lo > 1)
  expect_lte(abs(sh$slope_ratio - 1), 0.25)
  # the tolerances decide: a tolerance nothing exceeds leaves only the event count
  loose <- scr_segments(d, segment = "seg", min_events = 30, auc_tol = 1, offset_tol = 100, slope_tol = 100)
  expect_identical(sort(unique(loose$table$action)), c("shared", "too few events"))
  # a tight offset tolerance flags the base, whose level sits below the pooled one
  tight <- scr_segments(d, segment = "seg", min_events = 30, offset_tol = 0.01, slope_tol = 100, auc_tol = 1)
  bs <- tight$table[segment == "base"]
  expect_identical(bs$action, if (bs$oe_hi < 1 || bs$oe_lo > 1) "offset" else "shared")
  # the slope alone can ask for a separate model
  sl <- scr_segments(d, segment = "seg", min_events = 30, auc_tol = 1, slope_tol = 0.25)
  expect_identical(sl$table[segment == "weak", action], "separate model")
  expect_gt(abs(sl$table[segment == "weak", slope_ratio] - 1), 0.25)
  expect_lt(sl$table[segment == "weak", p_slope_adj], 0.05)
  # a slope ratio away from 1 is not enough: with no tolerance at all, the test alone decides
  ns <- scr_segments(d, segment = "seg", min_events = 0, auc_tol = 1, slope_tol = 0)$table
  expect_identical(ns$action == "separate model", !is.na(ns$p_slope_adj) & ns$p_slope_adj < 0.05)
  expect_true(any(ns$action != "separate model"))
})

test_that("by-groups repeat the study within each group, on bands frozen on all rows", {
  d <- seg_df()
  sg <- scr_segments(d, segment = "seg", by = "per", n_bands = 6)
  expect_identical(sg$test$group, c("h1", "h2")); expect_equal(nrow(sg$table), 6L)
  expect_equal(sg$cuts, scr_segments(d, segment = "seg", n_bands = 6)$cuts)
  B <- length(sg$cuts) + 1L
  for (g in c("h1", "h2")) {
    dg <- d[d$per == g, ]
    one <- scr_segments(dg, segment = "seg", n_bands = 6)
    t <- sg$table[group == g]
    # what does not depend on the bands equals the study of the group alone
    for (cn in c("n", "events", "rate", "auc", "auc_se", "ks", "slope", "auc_diff", "p_auc", "p_slope")) {
      expect_equal(t[[cn]], one$table[[cn]], info = paste(g, cn))
    }
    expect_equal(sg$test[group == g, statistic], one$test$statistic)
    # the expected events use the band rates of the group
    ib <- factor(findInterval(dg$score, sg$cuts) + 1L, seq_len(B))
    Rb <- as.numeric(tapply(dg$y, ib, mean))
    expect_equal(t[segment == "web", expected], sum(as.numeric(table(ib[dg$seg == "web"])) * Rb))
    expect_equal(sum(t$expected), sum(dg$y))
  }
  expect_error(scr_segments(transform(d, per = ifelse(seq_along(per) == 1, NA, per)), segment = "seg", by = "per"),
               "missing values")
})

test_that("a scorecard scores the new rows and reads the segments", {
  sc <- sc_demo()
  sg <- scr_segments(sc, scr_demo, segment = "ds_channel")
  d <- data.frame(score = scr_apply(sc, scr_demo)$score, y = scr_demo$default, seg = scr_demo$ds_channel)
  ref <- scr_segments(d, segment = "seg")
  expect_equal(sg$table[, -"segment"], ref$table[, -"segment"])
  expect_identical(sg$table$segment, sort(unique(scr_demo$ds_channel)))
  expect_identical(sg$target, "default"); expect_identical(sg$direction, sc$direction)
  # a missing segment is a segment of its own
  so <- scr_segments(sc, scr_demo, segment = "ds_optin", by = "ref_date")
  expect_true("(missing)" %in% so$table$segment)
  expect_equal(so$table[segment == "(missing)", sum(n)], sum(is.na(scr_demo$ds_optin)))
  expect_identical(so$test$group, sort(unique(as.character(scr_demo$ref_date))))
  expect_error(scr_segments(sc, segment = "ds_channel"), "needs `newdata`")
  expect_error(scr_segments(sc, scr_demo, segment = "nope"), "not in `newdata`")
})

test_that("the bootstrap interval is reproducible and leaves the user's stream alone", {
  d <- seg_df(2000, 4)
  a <- scr_segments(d, segment = "seg", n_boot = 60, seed = 3)
  b <- scr_segments(d, segment = "seg", n_boot = 60, seed = 3)
  expect_identical(a$table$auc_lo, b$table$auc_lo); expect_identical(a$table$auc_hi, b$table$auc_hi)
  d0 <- scr_segments(d, segment = "seg")
  expect_false(isTRUE(all.equal(a$table$auc_lo, d0$table$auc_lo)))
  expect_true(all(a$table$auc_lo < a$table$auc & a$table$auc < a$table$auc_hi))
  # the point estimates and the DeLong error do not depend on the bootstrap
  expect_equal(a$table$auc, d0$table$auc); expect_equal(a$table$auc_se, d0$table$auc_se)
  set.seed(77); r1 <- stats::runif(3)
  set.seed(77); invisible(scr_segments(d, segment = "seg", n_boot = 30, seed = 2)); r2 <- stats::runif(3)
  expect_identical(r1, r2)
})

test_that("edge cases: a single class, one band, tied scores, zero weights, one segment", {
  d <- seg_df(1500, 6)
  # a segment with a single class: no AUC, no slope, too few events
  z <- d; z$y[z$seg == "web"] <- 0L
  t <- scr_segments(z, segment = "seg")$table
  web <- t[segment == "web"]
  expect_true(is.na(web$auc)); expect_true(is.na(web$auc_se)); expect_true(is.na(web$slope))
  expect_true(is.na(web$offset)); expect_equal(web$oe_ratio, 0); expect_identical(web$action, "too few events")
  expect_true(is.na(web$p_auc))
  # the test runs on the two segments that have an error
  expect_identical(scr_segments(z, segment = "seg")$test$segments, 2L)
  # one band: the expected rate of every segment is the pooled rate
  one <- scr_segments(d, segment = "seg", n_bands = 1)
  expect_equal(one$table$expected, one$table$n * mean(d$y))
  expect_true(all(is.na(one$table$psi)))
  # every score tied: AUC 0.5, no slope, no band
  tie <- d; tie$score <- 7
  tt <- scr_segments(tie, segment = "seg")
  expect_equal(tt$table$auc, rep(0.5, 3)); expect_true(all(is.na(tt$table$slope)))
  # nothing is ranked: no interval, and the segments are left out of the test
  expect_equal(tt$table$auc_se, rep(0, 3)); expect_true(all(is.na(tt$table$auc_lo))); expect_true(all(is.na(tt$table$auc_hi)))
  expect_identical(tt$test$segments, 0L)
  expect_equal(sum(grepl("0\\.5000 - ", utils::capture.output(print(tt)))), 3L)
  expect_true(all(is.na(scr_segments(tie, segment = "seg", n_boot = 20, seed = 1)$table$auc_lo)))
  expect_length(tt$cuts, 0L)
  # rows with a zero weight are not there
  w0 <- d; w0$w0 <- as.numeric(seq_len(nrow(d)) %% 3 > 0)
  expect_equal(scr_segments(w0, segment = "seg", weight = "w0")$table[, .(n, events, auc, auc_se, expected, slope)],
               scr_segments(w0[w0$w0 > 0, ], segment = "seg")$table[, .(n, events, auc, auc_se, expected, slope)])
  # one segment: nothing to compare with
  s1 <- scr_segments(d[d$seg == "app", ], segment = "seg")
  expect_true(is.na(s1$test$statistic)); expect_true(is.na(s1$test$df)); expect_true(is.na(s1$table$p_auc))
  expect_equal(s1$table$oe_ratio, 1); expect_equal(s1$table$slope_ratio, 1)
  # the segment is the whole pool: its PSI is 0 by construction and has no critical value
  expect_equal(s1$table$psi, 0); expect_true(is.na(s1$table$psi_critical)); expect_true(is.na(s1$table$p_slope))
  # a missing outcome counts in the volume only
  na <- d; na$y[1:100] <- NA
  tn <- scr_segments(na, segment = "seg")$table
  expect_equal(sum(tn$n), nrow(d)); expect_equal(sum(tn$events), sum(na$y, na.rm = TRUE))
  expect_error(scr_segments(d), "`segment` is needed")
  expect_error(scr_segments(d, segment = "nope"), "not in `x`")
  expect_error(scr_segments(d, segment = "seg", auc_tol = -1), "^scr_segments\\(\\): `auc_tol` is outside")
  expect_error(scr_segments(d, segment = "seg", slope_tol = "a"), "^scr_segments\\(\\): `slope_tol` must be")
  expect_error(scr_segments(d[0, ], segment = "seg"), "no row has a score")
  expect_error(scr_segments(d, segment = "seg", min_events = 1.5), "min_events")
  expect_error(scr_segments(d, segment = "seg", foo = 1), "unused argument")
  expect_error(scr_segments(transform(d, score = NA_real_), segment = "seg"), "no row has a score")
})

test_that("numeric groups and segments are listed in numeric order", {
  d <- seg_df(3000, 21)
  d$code <- c(app = 9, store = 10, web = 100)[d$seg]
  d$yr <- c(h1 = 2, h2 = 11)[d$per]
  sg <- scr_segments(d, segment = "code", by = "yr")
  expect_identical(sg$test$group, c("2", "11")); expect_identical(sg$pooled$group, c("2", "11"))
  expect_identical(sg$table$segment, rep(c("9", "10", "100"), 2))
  ref <- scr_segments(d, segment = "seg", by = "per")
  expect_equal(sg$table$auc, ref$table$auc); expect_equal(sg$table$expected, ref$table$expected)
  # a missing numeric segment comes last; text keeps the order of its labels
  d$code[1:40] <- NA
  expect_identical(scr_segments(d, segment = "code")$table$segment, c("9", "10", "100", "(missing)"))
  expect_identical(scr_segments(d, segment = "seg")$table$segment, c("app", "store", "web"))
})

test_that("a factor segment column is listed in the order of its levels", {
  d <- seg_df(3000, 21)
  ref <- scr_segments(d, segment = "seg")
  d$fac <- factor(d$seg, levels = c("web", "app", "store", "unused"))
  sg <- scr_segments(d, segment = "fac")
  # the levels in use, in their order; an unused level is not a segment
  expect_identical(sg$table$segment, c("web", "app", "store"))
  # the same rows as the text column, only reordered
  expect_equal(sg$table[match(ref$table$segment, segment), !"segment"], ref$table[, !"segment"])
  expect_equal(sg$test, ref$test); expect_equal(sg$pooled, ref$pooled)
  # within every group of `by`, and with the missing segment last
  d$fac[1:40] <- NA
  sb <- scr_segments(d, segment = "fac", by = "per")
  expect_identical(sb$table$segment, rep(c("web", "app", "store", "(missing)"), 2))
  # an ordered factor follows its levels too; a factor of numbers is not sorted as numbers
  d$ord <- factor(d$seg, levels = c("store", "web", "app"), ordered = TRUE)
  expect_identical(scr_segments(d, segment = "ord")$table$segment, c("store", "web", "app"))
  d$num <- factor(c(app = 9, store = 10, web = 100)[d$seg], levels = c(100, 9, 10))
  expect_identical(scr_segments(d, segment = "num")$table$segment, c("100", "9", "10"))
  # groups keep the order of their labels
  d$fper <- factor(d$per, levels = c("h2", "h1"))
  expect_identical(scr_segments(d, segment = "seg", by = "fper")$test$group, c("h1", "h2"))
})

test_that("under one model for all segments the rules keep their size", {
  skip_on_cran()
  # 200 samples of 4,000 rows from one model: a small segment (250 rows, about 38 events) and four large ones
  set.seed(33)
  res <- vapply(1:200, function(i) {
    seg <- c(rep("small", 250), sample(c("a", "b", "c", "d"), 3750, TRUE))
    x <- stats::rnorm(4000)
    t <- scr_segments(data.frame(seg = seg, score = round(600 + 50 * x), y = stats::rbinom(4000, 1, stats::plogis(-2 - x))),
                      segment = "seg")$table
    sm <- t[segment == "small"]
    c(ratio = abs(sm$slope_ratio - 1) > 0.25, flagged = sm$action == "separate model", wald = sm$p_slope < 0.05,
      psi_small = sm$psi > sm$psi_critical, psi_large = t[segment == "a", psi > psi_critical])
  }, numeric(5))
  rate <- rowMeans(res)
  # the slope ratio alone leaves its tolerance often on 38 events; the action needs the test as well.
  # 200 samples put a 5% rate within about 0.03 (two standard errors): bounds of 0.10 allow for that
  expect_gt(rate[["ratio"]], 0.15)
  expect_lte(rate[["flagged"]], 0.10)
  expect_lte(rate[["wald"]], 0.10); expect_gte(rate[["wald"]], 0.01)
  # the PSI critical value allows for the segment being part of the pool: its size is near 5%
  # whatever the share of the segment (the independent-samples value would almost never reject)
  expect_lte(rate[["psi_small"]], 0.10); expect_gte(rate[["psi_small"]], 0.01)
  expect_lte(rate[["psi_large"]], 0.10); expect_gte(rate[["psi_large"]], 0.01)
})

test_that("pooled score cells: exact counts, close discrimination and slope", {
  d <- seg_df(4000, 13)
  d$score <- d$score + stats::runif(nrow(d))
  ex <- scr_segments(d, segment = "seg", n_bands = 5)
  po <- scr_segments(d, segment = "seg", n_bands = 5, max_cells = 200)
  expect_true(po$quantized); expect_false(ex$quantized)
  expect_equal(po$table$n, ex$table$n); expect_equal(po$table$events, ex$table$events)
  expect_equal(po$table$auc, ex$table$auc, tolerance = 5e-3)
  expect_equal(po$table$slope_ratio, ex$table$slope_ratio, tolerance = 1e-2)
  # the bands hold whole cells: the standardization is exact for the pooled cuts
  ib <- factor(findInterval(d$score, po$cuts) + 1L, seq_len(length(po$cuts) + 1L))
  Rb <- as.numeric(tapply(d$y, ib, mean))
  expect_equal(po$table[segment == "web", expected], sum(as.numeric(table(ib[d$seg == "web"])) * Rb))
  # and within groups
  pb <- scr_segments(d, segment = "seg", by = "per", n_bands = 5, max_cells = 200)
  expect_equal(nrow(pb$table), 6L)
  expect_equal(pb$table[group == "h2" & segment == "app", n], sum(d$per == "h2" & d$seg == "app"))
  expect_equal(pb$pooled$n, as.numeric(table(d$per)))
})

test_that("print and export", {
  d <- seg_df(2000, 2)
  sg <- scr_segments(d, segment = "seg")
  out <- utils::capture.output(print(sg))
  expect_match(out[1], "^<scr_segments> target \"y\"")
  expect_true(any(grepl("equal AUC across 3 segments: chi-square", out)))
  # one line per segment, ending in its action
  expect_equal(sum(grepl("^  (app|store|web) .*(shared|offset|separate model|too few events)$", out)), 3L)
  skip_if_not_installed("openxlsx")
  old <- scr_verbose(FALSE); on.exit(scr_verbose(old), add = TRUE)
  dir <- file.path(tempdir(), "scr-segments-export")
  unlink(dir, recursive = TRUE)
  ex <- scr_export(sg, dir, stamp = FALSE)
  expect_identical(basename(ex$files$xlsx), "segments_y.xlsx")
  expect_identical(openxlsx::getSheetNames(ex$files$xlsx), c("Segments", "Test", "Pooled", "Cuts", "Settings"))
  expect_equal(nrow(openxlsx::read.xlsx(ex$files$xlsx, sheet = "Segments")), 3L)
})
