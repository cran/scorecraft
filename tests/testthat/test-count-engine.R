# Count engine: the band tables, the DeLong error and the bootstrap of
# scr_metrics() against the row-level implementations they replace or keep
# for scores with few ties. Those are kept below as test-only references.

# -- row-level references ------------------------------------------------------- #

# bands as a factor per row
ref_score_gains <- function(score, y, breaks, direction) {
  band <- cut(score, breaks = breaks, include.lowest = TRUE)
  d <- data.table::data.table(band = band, score = score, y = as.integer(y))[
    , .(n = .N, events = sum(y), min_score = min(score), mean_score = mean(score), max_score = max(score)),
    by = band]
  d <- if (identical(direction, "higher_is_safer")) d[order(band)] else d[order(-as.integer(band))]
  n_tot <- sum(d$n); e_tot <- sum(d$events); ne_tot <- n_tot - e_tot
  d[, `:=`(id = seq_len(.N), pct = n / n_tot, event_rate = events / n, non_events = n - events)]
  bw <- .band_woe(d$events, d$non_events)
  d[, `:=`(pct_event = bw$pct_event, pct_nonevent = bw$pct_nonevent, woe = bw$log_odds)]
  d[, `:=`(cum_pct = cumsum(pct), cum_event_pct = cumsum(events) / max(1, e_tot),
           cum_nonevent_pct = cumsum(non_events) / max(1, ne_tot))]
  safer <- identical(direction, "higher_is_safer")
  d[, `:=`(ks = abs(cum_event_pct - cum_nonevent_pct),
           lift = event_rate / (e_tot / n_tot),
           cum_lift = (cumsum(events) / cumsum(n)) / (e_tot / n_tot),
           odds = if (safer) (n - events + 0.5) / (events + 0.5) else (events + 0.5) / (n - events + 0.5))]
  d[, log_odds := log(odds)]
  data.table::setcolorder(d, c("id", "band", "n", "pct", "events", "non_events", "event_rate",
                               "pct_event", "pct_nonevent", "woe",
                               "min_score", "mean_score", "max_score", "cum_pct", "cum_event_pct",
                               "cum_nonevent_pct", "ks", "lift", "cum_lift", "odds", "log_odds"))
  d[, band := as.character(band)]
  d[]
}

ref_calibration <- function(s, breaks, al) {
  p <- .score_to_prob(al, s$score)
  y <- s$y
  band <- cut(s$score, breaks = breaks, include.lowest = TRUE)
  tb <- data.table::data.table(band = band, p = p, y = y)[
    , .(n = .N, expected = mean(p), observed = mean(y)), by = band][order(band)]
  tb[, band := as.character(band)]
  tb[, gap := observed - expected]
  ece <- sum(tb$n / sum(tb$n) * abs(tb$gap))
  lo <- suppressWarnings(stats::glm(y ~ stats::qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6)), family = stats::binomial()))
  cf <- stats::coef(lo)
  list(summary = data.table::data.table(sample = "holdout", n = length(y), brier = mean((p - y)^2),
                                        ece = ece, mce = max(abs(tb$gap)),
                                        intercept = unname(cf[1]), slope = unname(cf[2]),
                                        expected_rate = mean(p), observed_rate = mean(y)),
       table = tb[])
}

ref_strategy <- function(x, breaks = NULL, decisions = NULL, revenue_good = 1, loss_bad = 1,
                         sample = "holdout", rule = c("breakeven", "crossing")) {
  rule <- match.arg(rule)
  prop <- identical(x$config$objective, "propensity")
  breaks <- breaks %||% x$breaks
  s <- x$samples[[sample]]
  band <- cut(s$score, breaks = breaks, include.lowest = TRUE)
  d <- data.table::data.table(band = band, y = s$y, score = s$score)[
    , .(n = .N, events = sum(y), event_rate = mean(y), min_score = min(score), max_score = max(score)), by = band]
  d <- if (xor(identical(x$direction, "higher_is_safer"), prop)) d[order(-as.integer(band))] else d[order(band)]
  idx <- as.integer(d$band)
  d[, `:=`(id = seq_len(.N), pct = n / sum(n), band = as.character(band))]
  bw <- .band_woe(d$events, d$n - d$events)
  d[, `:=`(pct_event = bw$pct_event, pct_nonevent = bw$pct_nonevent, odds_event = bw$odds_event, log_odds = bw$log_odds)]
  p_bad <- if (prop) 1 - d$event_rate else d$event_rate
  be_bad <- revenue_good / (revenue_good + loss_bad)
  breakeven <- if (prop) loss_bad / (revenue_good + loss_bad) else be_bad
  d[, ep_per_account := (1 - p_bad) * revenue_good - p_bad * loss_bad]
  d[, band_profit := n * ep_per_account]
  cr <- .strategy_crossing(d$pct_event, d$pct_nonevent, d$log_odds, idx, d$band, d$min_score, d$max_score,
                           breaks, x$samples$train$score)
  lab <- if (prop) c("target", "review", "skip") else c("approve", "review", "decline")
  if (is.null(decisions)) {
    if (rule == "crossing") {
      if (is.na(cr$k)) stop("scr_strategy(): rule = \"crossing\" needs at least two bands and both classes in the sample.", call. = FALSE)
      d[, decision := data.table::fifelse(is.na(idx), NA_character_,
                       data.table::fifelse(seq_len(.N) <= cr$k, lab[1], lab[3]))]
    } else {
      d[, decision := data.table::fifelse(p_bad <= be_bad, lab[1],
                       data.table::fifelse(p_bad <= 1.25 * be_bad, lab[2], lab[3]))]
    }
  } else {
    d[, decision := as.character(decisions)]
  }
  d[, `:=`(cum_pct = cumsum(pct), cum_event_rate = cumsum(events) / cumsum(n), cum_profit = cumsum(band_profit))]
  data.table::setcolorder(d, c("id", "band", "min_score", "max_score", "n", "pct", "events", "event_rate",
                               "pct_event", "pct_nonevent", "odds_event", "log_odds",
                               "decision", "ep_per_account", "band_profit", "cum_pct", "cum_event_rate", "cum_profit"))
  structure(list(table = d[], breakeven = breakeven, revenue_good = revenue_good, loss_bad = loss_bad,
                 sample = sample, direction = x$direction, target = x$target,
                 objective = if (prop) "propensity" else "risk", rule = rule, crossing = cr$crossing),
            class = c("scr_strategy", "list"))
}

# coverage and sensitivity tables of scr_reject(), bands as factors
ref_reject <- function(x, population = NULL, accepted = NULL, multipliers = NULL, sample = "holdout") {
  multipliers <- multipliers %||% x$config$reject_multipliers
  s <- x$samples[[sample]]
  breaks <- x$breaks
  band_dev <- cut(s$score, breaks = breaks, include.lowest = TRUE)
  dev <- data.table::data.table(band = band_dev, y = s$y)[, .(n_dev = .N, events_dev = sum(y), rate_dev = mean(y)), by = band]
  pop_tb <- NULL
  if (!is.null(population)) {
    sp <- scr_apply(x, population)$score
    acc <- if (is.null(accepted)) rep(FALSE, length(sp)) else as.logical(accepted)
    band_pop <- cut(sp, breaks = breaks, include.lowest = TRUE)
    pop_tb <- data.table::data.table(band = band_pop, acc = acc)[, .(n_pop = .N, n_unknown = sum(!acc)), by = band]
  }
  lv <- levels(band_dev)
  cov <- data.table::data.table(band = factor(lv, levels = lv))
  cov <- merge(cov, dev, by = "band", all.x = TRUE)
  if (!is.null(pop_tb)) cov <- merge(cov, pop_tb, by = "band", all.x = TRUE) else cov[, `:=`(n_pop = NA_integer_, n_unknown = 0L)]
  for (cn in c("n_dev", "events_dev", "n_unknown")) cov[is.na(get(cn)), (cn) := 0L]
  cov[, coverage := if (all(is.na(n_pop))) NA_real_ else n_dev / pmax(1L, n_pop)]
  cov[, coverage_flag := data.table::fifelse(n_dev == 0L, "no_outcome",
                          data.table::fifelse(events_dev < 30L, "few_events", "ok"))]
  cov <- if (identical(x$direction, "higher_is_safer")) cov[order(-as.integer(band))] else cov[order(band)]
  cov[, band := as.character(band)]
  sens <- data.table::rbindlist(lapply(multipliers, function(m) {
    r <- data.table::copy(cov)
    r[, multiplier := m]
    r[, rate_unknown := pmin(1, rate_dev * m)]
    r[, events_implied := events_dev + n_unknown * data.table::fifelse(is.na(rate_unknown), 0, rate_unknown)]
    r[, rate_implied := events_implied / pmax(1L, n_dev + n_unknown)]
    tot <- data.table::data.table(band = "TOTAL", n_dev = sum(r$n_dev), events_dev = sum(r$events_dev),
                                  rate_dev = sum(r$events_dev) / max(1L, sum(r$n_dev)), n_pop = sum(r$n_pop),
                                  n_unknown = sum(r$n_unknown), coverage = NA_real_, coverage_flag = "",
                                  multiplier = m, rate_unknown = NA_real_, events_implied = sum(r$events_implied),
                                  rate_implied = sum(r$events_implied) / max(1L, sum(r$n_dev + r$n_unknown)))
    data.table::rbindlist(list(r, tot), use.names = TRUE, fill = TRUE)
  }))
  list(coverage = cov[, .(band, n_dev, events_dev, rate_dev, n_pop, n_unknown, coverage, coverage_flag)],
       sensitivity = sens[, .(multiplier, band, n_dev, events_dev, rate_dev, n_unknown, rate_unknown, events_implied, rate_implied)])
}

# numeric branch of scr_psi(), bands as factors
ref_psi <- function(base, compare, breaks = NULL, n_groups = 10L, alpha = 0.05, thresholds = c(0.10, 0.25)) {
  if (is.null(breaks)) {
    probs  <- seq(0, 1, length.out = n_groups + 1L)[-c(1L, n_groups + 1L)]
    breaks <- unique(c(-Inf, stats::quantile(base, probs = probs, na.rm = TRUE, names = FALSE), Inf))
  }
  gb <- cut(base, breaks = breaks, include.lowest = TRUE)
  gc <- cut(compare, breaks = breaks, include.lowest = TRUE)
  lv <- levels(gb)
  .psi_counts(tabulate(match(gb, lv), nbins = length(lv)), tabulate(match(gc, lv), nbins = length(lv)), lv, alpha, thresholds)
}

# the row-level DeLong error, ref_auc_se(), is in helper-scorecraft.R

# scr_metrics() as it was: the bootstrap on the rows for any score, with the
# seed of every resample drawn in the main process
ref_metrics <- function(score, y, higher_is_event = TRUE, ci = TRUE, n_boot = 200L,
                        level = 0.95, seed = NULL, nthread = 1L) {
  empty <- list(auc = NA_real_, ks = NA_real_, gini = NA_real_,
                auc_lo = NA_real_, auc_hi = NA_real_, ks_lo = NA_real_, ks_hi = NA_real_,
                gini_lo = NA_real_, gini_hi = NA_real_, n = 0L, events = 0L,
                n_boot = 0L, level = level)
  y <- .scr_y01(y, "scr_metrics")
  ok <- is.finite(score) & !is.na(y)
  if (!any(ok)) return(structure(empty, class = c("scr_metrics", "list")))
  score <- as.double(score[ok]); y <- y[ok]
  if (!isTRUE(higher_is_event)) score <- -score
  n1 <- sum(y == 1L); n0 <- sum(y == 0L)
  if (n1 == 0L || n0 == 0L) {
    empty$n <- length(y); empty$events <- n1
    return(structure(empty, class = c("scr_metrics", "list")))
  }
  idx <- data.table::frank(score, ties.method = "dense")
  K <- max(idx)
  idx1 <- idx[y == 1L]; idx0 <- idx[y == 0L]
  pt <- .auc_ks_counts(tabulate(idx1, K), tabulate(idx0, K))
  out <- empty
  out$auc <- pt$auc; out$ks <- pt$ks; out$gini <- pt$gini
  out$n <- length(y); out$events <- n1
  if (isTRUE(ci) && n_boot >= 2L) {
    .scr_local_seed(seed)
    seeds <- sample.int(.Machine$integer.max, n_boot)
    .scr_rng_guard()
    reps <- .scr_lapply(seeds, function(sd) {
      set.seed(sd)
      c1 <- tabulate(idx1[sample.int(n1, n1, replace = TRUE)], K)
      c0 <- tabulate(idx0[sample.int(n0, n0, replace = TRUE)], K)
      r <- .auc_ks_counts(c1, c0)
      c(r$auc, r$ks)
    }, nthread = nthread)
    b <- do.call(rbind, reps)
    a <- (1 - level) / 2
    q_auc <- stats::quantile(b[, 1], c(a, 1 - a), na.rm = TRUE, names = FALSE)
    q_ks  <- stats::quantile(b[, 2], c(a, 1 - a), na.rm = TRUE, names = FALSE)
    out$auc_lo <- q_auc[1]; out$auc_hi <- q_auc[2]
    out$ks_lo <- q_ks[1];   out$ks_hi <- q_ks[2]
    out$gini_lo <- 2 * q_auc[1] - 1; out$gini_hi <- 2 * q_auc[2] - 1
    out$n_boot <- as.integer(n_boot)
  }
  structure(out, class = c("scr_metrics", "list"))
}

# the bootstrap on the counts, as scr_metrics() runs it under many ties
ref_metrics_counts <- function(score, y, n_boot, level = 0.95, seed = NULL) {
  idx <- data.table::frank(score, ties.method = "dense"); K <- max(idx)
  .study_auc_boot(tabulate(idx[y == 1], K), tabulate(idx[y == 0], K), as.integer(n_boot), level, seed = seed, boot_cells = Inf)
}

# n rows over exactly K distinct scores, both classes present
ce_ties <- function(n, K, seed = 1) {
  set.seed(seed)
  s <- sample(c(seq_len(K), sample.int(K, n - K, replace = TRUE)))
  list(s = s / 10, y = stats::rbinom(n, 1, stats::plogis(-1 + 3 * s / K)))
}

# -- fixtures ---------------------------------------------------------------------- #

# the value of a call, or the message of its error
val <- function(expr) tryCatch(expr, error = function(e) conditionMessage(e))

ce_cards <- function() list(credit = sc_demo(), fraud = sc_fraud_demo(), propensity = sc_prop_demo())

# synthetic scores: heavy ties, continuous, a handful of values, missing and
# infinite scores, and scores sitting on the band edges
ce_scores <- function() {
  set.seed(4102)
  n <- 4000
  x <- stats::rnorm(n)
  y <- stats::rbinom(n, 1, stats::plogis(-2 - x))
  ties <- round(600 + 50 * x)
  miss <- ties; miss[sample.int(n, 80)] <- NA; miss[1:3] <- Inf; miss[4:6] <- -Inf
  edge <- c(500, 550, 600, 650, 700, 499.999, 700.001, 500 + 1e-9, 550 - 1e-9, rep(c(520, 575, 610, 690), 30))
  list(ties = list(s = ties, y = y), continuous = list(s = 600 + 50 * x + stats::rnorm(n), y = y),
       few = list(s = pmin(pmax(round(x), -2), 2), y = y), missing = list(s = miss, y = y),
       edge = list(s = edge, y = rep(c(0, 1, 0, 0, 1, 0, 1), length.out = length(edge))))
}

ce_breaks <- function(s) {
  list(quantile = .score_breaks(s, 10L), finite = c(500, 550, 600, 650, 700),
       unsorted = c(Inf, 620, -Inf, 580, 600), wide = c(-Inf, 1e-3, 599.5, 600, 1234567.891, Inf), one = c(-Inf, Inf))
}

# a scorecard carrying synthetic samples and breaks
ce_card <- function(base, s, y, breaks) {
  h <- seq_along(s) %% 2L == 0L
  base$samples$train <- data.table::data.table(score = s[!h], y = y[!h])
  base$samples$holdout <- data.table::data.table(score = s[h], y = y[h])
  base$breaks <- breaks
  base
}

# -- band assignment --------------------------------------------------------------- #

test_that(".score_band() gives the codes and the labels of cut(include.lowest = TRUE)", {
  set.seed(1)
  x <- c(-Inf, Inf, NA, NaN, 0, 1, 2, 3, 1 - 1e-12, 1 + 1e-12, 3 + 1e-9, -1e-300, -5, 10,
         stats::runif(2000, -1, 4), round(stats::runif(2000, -1, 4)))
  brs <- list(c(0, 1, 2, 3), c(-Inf, 1, 2, Inf), c(-Inf, Inf), c(3, 0, 2, 1), c(0, 1), c(-Inf, 0), c(0, Inf),
              c(-1, NA, 2, 5), c(-Inf, 1e-3, 599.5, 600, 1234567.891, Inf), c(0L, 2L, 4L))
  for (br in brs) {
    f <- cut(x, breaks = br, include.lowest = TRUE)
    b <- .score_band(x, br)
    expect_identical(b$idx, as.integer(f))
    expect_identical(b$labels, levels(f))
  }
  # the lowest edge is in the first band, the highest in the last, outside is NA
  b <- .score_band(c(0, 3, -1e-9, 3 + 1e-9, 1, 2), c(0, 1, 2, 3))
  expect_identical(b$idx, c(1L, 3L, NA, NA, 1L, 2L))
  expect_identical(b$labels, c("[0,1]", "(1,2]", "(2,3]"))
  # infinite edges hold infinite scores
  expect_identical(.score_band(c(-Inf, Inf), c(-Inf, 0, Inf))$idx, c(1L, 2L))
  # integer scores and no score at all
  expect_identical(.score_band(c(0L, 1L, 2L, 5L), c(0, 1, 3))$idx, c(1L, 1L, 2L, NA))
  expect_identical(.score_band(numeric(), c(0, 1, 3)), list(idx = integer(), labels = c("[0,1]", "(1,3]")))
  # a number of intervals takes its edges from the scores, as cut() does
  f <- cut(x[is.finite(x)], breaks = 4, include.lowest = TRUE)
  b <- .score_band(x[is.finite(x)], 4)
  expect_identical(b$idx, as.integer(f)); expect_identical(b$labels, levels(f))
  expect_error(.score_band(x, c(0, 1, 1)), "not unique")
})

test_that("gains and calibration equal the row-level band tables on the demo scorecards", {
  for (sc in ce_cards()) {
    for (smp in c("train", "holdout")) {
      s <- sc$samples[[smp]]
      expect_identical(.score_gains(s$score, s$y, sc$breaks, sc$direction), ref_score_gains(s$score, s$y, sc$breaks, sc$direction))
      # whole points: scores tied on the band edges
      bp <- .score_breaks(sc$samples$train$score_points, 10L)
      expect_identical(.score_gains(s$score_points, s$y, bp, sc$direction), ref_score_gains(s$score_points, s$y, bp, sc$direction))
      expect_identical(.calibration(s, sc$breaks, sc$alignment), ref_calibration(s, sc$breaks, sc$alignment))
    }
    # the tables stored in the scorecard are the ones of the references
    expect_identical(sc$gains[sample == "holdout", !"sample"],
                     ref_score_gains(sc$samples$holdout$score, sc$samples$holdout$y, sc$breaks, sc$direction))
    expect_identical(sc$calibration, ref_calibration(sc$samples$holdout, sc$breaks, sc$alignment))
  }
})

test_that("gains and calibration equal the row-level band tables on ties, edges and missing scores", {
  al <- sc_demo()$alignment
  sets <- ce_scores()
  for (nm in names(sets)) {
    s <- sets[[nm]]$s; y <- sets[[nm]]$y
    for (br in ce_breaks(s)) {
      for (dir in c("higher_is_safer", "higher_is_riskier")) {
        expect_identical(.score_gains(s, y, br, dir), ref_score_gains(s, y, br, dir))
      }
      if (nm != "missing") {
        d <- data.table::data.table(score = s, y = y)
        expect_identical(.calibration(d, br, al), ref_calibration(d, br, al))
      }
    }
  }
  # scores outside finite breaks form one band without a label
  g <- .score_gains(sets$edge$s, sets$edge$y, c(500, 550, 600, 650, 700), "higher_is_safer")
  expect_identical(g$band, c("[500,550]", "(550,600]", "(600,650]", "(650,700]", NA))
  expect_equal(g$n[5], 2L)
})

test_that("the strategy table equals the row-level one under both rules", {
  for (sc in ce_cards()) {
    ho <- sc$samples$holdout
    q <- stats::quantile(ho$score, c(0.1, 0.4, 0.7, 0.95), names = FALSE)
    for (smp in c("train", "holdout")) {
      expect_identical(scr_strategy(sc, sample = smp, revenue_good = 1080, loss_bad = 4500),
                       ref_strategy(sc, sample = smp, revenue_good = 1080, loss_bad = 4500))
      expect_identical(scr_strategy(sc, sample = smp, rule = "crossing"), ref_strategy(sc, sample = smp, rule = "crossing"))
    }
    # finite breaks (a row without band), unsorted breaks, a number of intervals, one band
    for (br in list(q, rev(c(-Inf, q, Inf)), 5, c(-Inf, Inf))) {
      expect_identical(scr_strategy(sc, breaks = br), ref_strategy(sc, breaks = br))
      expect_identical(val(scr_strategy(sc, breaks = br, rule = "crossing")), val(ref_strategy(sc, breaks = br, rule = "crossing")))
    }
    expect_error(scr_strategy(sc, breaks = c(-Inf, q[1], q[1], Inf)), "not unique")
  }
  sets <- ce_scores()
  for (sc in ce_cards()) for (nm in names(sets)) {
    s <- sets[[nm]]$s; y <- sets[[nm]]$y
    for (br in ce_breaks(s)) {
      z <- ce_card(sc, s, y, br)
      expect_identical(scr_strategy(z, revenue_good = 300, loss_bad = 1500), ref_strategy(z, revenue_good = 300, loss_bad = 1500))
      expect_identical(val(scr_strategy(z, rule = "crossing")), val(ref_strategy(z, rule = "crossing")))
    }
  }
})

test_that("the coverage and the sensitivity band of scr_reject() equal the row-level ones", {
  acc <- seq_len(nrow(scr_demo)) %in% res_demo()$split$holdout_idx
  for (sc in ce_cards()) {
    for (smp in c("train", "holdout")) {
      rj <- scr_reject(sc, sample = smp); ref <- ref_reject(sc, sample = smp)
      expect_identical(rj$coverage, ref$coverage); expect_identical(rj$sensitivity, ref$sensitivity)
    }
    rj <- scr_reject(sc, population = scr_demo, accepted = acc, multipliers = c(2, 4))
    ref <- ref_reject(sc, population = scr_demo, accepted = acc, multipliers = c(2, 4))
    expect_identical(rj$coverage, ref$coverage); expect_identical(rj$sensitivity, ref$sensitivity)
  }
  sets <- ce_scores()
  for (sc in ce_cards()) for (nm in names(sets)) {
    s <- sets[[nm]]$s; y <- sets[[nm]]$y
    for (br in ce_breaks(s)) {
      z <- ce_card(sc, s, y, br)
      rj <- scr_reject(z); ref <- ref_reject(z)
      expect_identical(rj$coverage, ref$coverage); expect_identical(rj$sensitivity, ref$sensitivity)
    }
  }
  # a band without a development row is kept, flagged as without outcome
  z <- ce_card(sc_demo(), sets$edge$s, sets$edge$y, c(-Inf, 400, 600, Inf))
  expect_identical(scr_reject(z)$coverage[band == "[-Inf,400]", coverage_flag], "no_outcome")
})

test_that("the numeric scr_psi() equals the row-level band counts", {
  for (sc in ce_cards()) {
    tr <- sc$samples$train; ho <- sc$samples$holdout
    expect_identical(scr_psi(tr$score, ho$score), ref_psi(tr$score, ho$score))
    expect_identical(scr_psi(tr$score, ho$score, breaks = sc$breaks), ref_psi(tr$score, ho$score, breaks = sc$breaks))
    expect_identical(scr_psi(tr$score_points, ho$score_points, n_groups = 20L), ref_psi(tr$score_points, ho$score_points, n_groups = 20L))
    expect_identical(sc$stability$score_psi, ref_psi(tr$score, ho$score, breaks = sc$breaks, alpha = sc$config$psi_alpha))
  }
  sets <- ce_scores()
  for (nm in names(sets)) {
    s <- sets[[nm]]$s
    h <- seq_along(s) %% 2L == 0L
    a <- s[!h]; b <- s[h] + 7
    expect_identical(scr_psi(a, b), ref_psi(a, b))
    expect_identical(scr_psi(a, b, n_groups = 4L, alpha = 0.01), ref_psi(a, b, n_groups = 4L, alpha = 0.01))
    for (br in ce_breaks(s)) expect_identical(scr_psi(a, b, breaks = br), ref_psi(a, b, breaks = br))
    # a number of intervals: each sample keeps its own edges, matched by label
    # (an error with infinite scores, the same as before)
    expect_identical(val(scr_psi(a, b, breaks = 4)), val(ref_psi(a, b, breaks = 4)))
  }
  # integer vectors, a sample without rows, a base without values
  expect_identical(scr_psi(1:50, 20:80), ref_psi(1:50, 20:80))
  expect_identical(scr_psi(numeric(), sets$ties$s), ref_psi(numeric(), sets$ties$s))
  expect_identical(scr_psi(rep(NA_real_, 5), 1:10), ref_psi(rep(NA_real_, 5), 1:10))
  expect_error(scr_psi(sets$ties$s, sets$ties$s, breaks = c(-Inf, 600, 600, Inf)), "not unique")
})

# -- DeLong ------------------------------------------------------------------------ #

test_that("the DeLong error from counts equals the mid-rank estimator on the rows", {
  for (sc in ce_cards()) {
    hie <- identical(sc$direction, "higher_is_riskier")
    for (s in sc$samples) {
      expect_equal(.pd_auc_se(s$score, s$y, higher_is_event = hie), ref_auc_se(s$score, s$y, higher_is_event = hie), tolerance = 1e-12)
      expect_equal(.pd_auc_se(s$score_points, s$y, higher_is_event = hie), ref_auc_se(s$score_points, s$y, higher_is_event = hie),
                   tolerance = 1e-12)
    }
  }
  sets <- ce_scores()
  for (nm in names(sets)) {
    s <- sets[[nm]]$s; y <- sets[[nm]]$y
    expect_equal(.pd_auc_se(s, y), ref_auc_se(s, y), tolerance = 1e-12)
    expect_equal(.pd_auc_se(s, y, higher_is_event = FALSE), ref_auc_se(s, y, higher_is_event = FALSE), tolerance = 1e-12)
  }
  # missing outcomes are dropped, a logical outcome is counted, every score tied gives 0
  y <- sets$ties$y; y[1:50] <- NA
  expect_equal(.pd_auc_se(sets$ties$s, y), ref_auc_se(sets$ties$s, y), tolerance = 1e-12)
  expect_equal(.pd_auc_se(sets$ties$s, sets$ties$y == 1), ref_auc_se(sets$ties$s, sets$ties$y), tolerance = 1e-12)
  expect_identical(.pd_auc_se(rep(1, 40), rep(0:1, 20)), 0)
  # fewer than two events or non-events, or no row: no estimate
  expect_true(is.na(.pd_auc_se(c(1, 2, 3), c(1, 0, 0))))
  expect_true(is.na(.pd_auc_se(c(1, 2, 3, NA), c(1, 1, 0, 0))))
  expect_true(is.na(.pd_auc_se(numeric(), integer())))
})

# -- bootstrap of scr_metrics() ---------------------------------------------------- #

bounds <- c("auc_lo", "auc_hi", "ks_lo", "ks_hi", "gini_lo", "gini_hi")

test_that("the point estimates of scr_metrics() are those of the previous version", {
  pt <- c("auc", "ks", "gini", "n", "events", "n_boot", "level")
  for (sc in ce_cards()) {
    hie <- identical(sc$direction, "higher_is_riskier")
    for (s in sc$samples) for (v in list(s$score, s$score_points)) {
      m <- scr_metrics(v, s$y, higher_is_event = hie, n_boot = 10, seed = 1)
      expect_identical(unclass(m)[pt], unclass(ref_metrics(v, s$y, higher_is_event = hie, n_boot = 10, seed = 1))[pt])
      expect_identical(scr_metrics(v, s$y, higher_is_event = hie, ci = FALSE), ref_metrics(v, s$y, higher_is_event = hie, ci = FALSE))
    }
  }
  sets <- ce_scores()
  for (nm in names(sets)) {
    s <- sets[[nm]]$s; y <- sets[[nm]]$y
    m <- scr_metrics(s, y, higher_is_event = FALSE, n_boot = 10, seed = 1)
    expect_identical(unclass(m)[pt], unclass(ref_metrics(s, y, higher_is_event = FALSE, n_boot = 10, seed = 1))[pt])
    expect_true(m$auc_lo <= m$auc && m$auc <= m$auc_hi && m$ks_lo <= m$ks_hi)
    expect_identical(m$gini_lo, 2 * m$auc_lo - 1); expect_identical(m$gini_hi, 2 * m$auc_hi - 1)
  }
  # degenerate inputs are untouched
  expect_identical(scr_metrics(1:5, rep(1, 5)), ref_metrics(1:5, rep(1, 5)))
  expect_identical(scr_metrics(numeric(), integer()), ref_metrics(numeric(), integer()))
})

test_that("the bootstrap runs on the counts up to one distinct score per ten rows, on the rows above", {
  for (n in c(200L, 205L)) {
    tenth <- n %/% 10L
    # K = floor(n / 10): the counts
    d <- ce_ties(n, tenth)
    expect_identical(data.table::uniqueN(d$s), tenth)
    m <- scr_metrics(d$s, d$y, n_boot = 30, seed = 6)
    expect_identical(unlist(unclass(m)[bounds]), unlist(ref_metrics_counts(d$s, d$y, 30, seed = 6)[bounds]))
    expect_false(identical(unclass(m)[bounds], unclass(ref_metrics(d$s, d$y, n_boot = 30, seed = 6))[bounds]))
    # one distinct score more: the rows, as before
    d <- ce_ties(n, tenth + 1L)
    m <- scr_metrics(d$s, d$y, n_boot = 30, seed = 6)
    expect_identical(m, ref_metrics(d$s, d$y, n_boot = 30, seed = 6))
    expect_false(identical(unlist(unclass(m)[bounds]), unlist(ref_metrics_counts(d$s, d$y, 30, seed = 6)[bounds])))
  }
  # between one distinct score per ten rows and one per two: the rows
  for (K in c(50L, 100L)) {
    d <- ce_ties(200L, K)
    expect_identical(scr_metrics(d$s, d$y, n_boot = 30, seed = 6), ref_metrics(d$s, d$y, n_boot = 30, seed = 6))
  }
  # the rows that count are those with a score and an outcome: 20 distinct
  # scores in 200 valid rows run on the counts, whatever the rows dropped
  d <- ce_ties(200L, 20L)
  s <- c(d$s, NA, Inf, 1001:1050 / 7); y <- c(d$y, 1, 0, rep(NA, 50))
  expect_identical(unlist(unclass(scr_metrics(s, y, n_boot = 30, seed = 6))[bounds]),
                   unlist(ref_metrics_counts(d$s, d$y, 30, seed = 6)[bounds]))
  # and 21 distinct scores in 200 valid rows run on the rows, although the
  # dropped rows would bring the share of distinct scores under one in ten
  d <- ce_ties(200L, 21L)
  s <- c(d$s, rep(NA, 100)); y <- c(d$y, rep(1, 100))
  expect_identical(unclass(scr_metrics(s, y, n_boot = 30, seed = 6))[bounds],
                   unclass(ref_metrics(d$s, d$y, n_boot = 30, seed = 6))[bounds])
  # the switch does not look at nthread
  d <- ce_ties(200L, 20L)
  expect_identical(scr_metrics(d$s, d$y, n_boot = 30, seed = 6, nthread = 2L), scr_metrics(d$s, d$y, n_boot = 30, seed = 6))
  # every score tied: one cell, on the counts
  tied <- scr_metrics(rep(1, 50), rep(0:1, 25), n_boot = 20, seed = 1)
  expect_identical(c(tied$auc_lo, tied$auc_hi, tied$ks_lo, tied$ks_hi), c(0.5, 0.5, 0, 0))
  # two rows, two scores: the rows
  expect_identical(scr_metrics(c(1, 2), c(0, 1), n_boot = 5, seed = 1), ref_metrics(c(1, 2), c(0, 1), n_boot = 5, seed = 1))
})

test_that("on scores with few ties the intervals are those of the previous version, bit for bit", {
  sets <- ce_scores()
  s <- sets$continuous$s; y <- sets$continuous$y
  expect_gt(data.table::uniqueN(s), length(s) / 10)
  for (sd in c(1, 7, 2026)) {
    expect_identical(scr_metrics(s, y, higher_is_event = FALSE, n_boot = 40, seed = sd),
                     ref_metrics(s, y, higher_is_event = FALSE, n_boot = 40, seed = sd))
  }
  expect_identical(scr_metrics(-s, y, n_boot = 25, seed = 3, level = 0.9), ref_metrics(-s, y, n_boot = 25, seed = 3, level = 0.9))
  # missing scores and outcomes, a logical outcome
  s2 <- s; s2[1:50] <- NA; y2 <- y == 1; y2[51:90] <- NA
  expect_identical(scr_metrics(s2, y2, higher_is_event = FALSE, n_boot = 25, seed = 3),
                   ref_metrics(s2, y2, higher_is_event = FALSE, n_boot = 25, seed = 3))
  # the exact score of the demo scorecards, where it has few ties
  for (sc in ce_cards()) {
    hie <- identical(sc$direction, "higher_is_riskier")
    for (smp in sc$samples) {
      if (data.table::uniqueN(smp$score) > nrow(smp) / 10) {
        expect_identical(scr_metrics(smp$score, smp$y, higher_is_event = hie, n_boot = 20, seed = 2),
                         ref_metrics(smp$score, smp$y, higher_is_event = hie, n_boot = 20, seed = 2))
      }
    }
  }
  # without a seed the same draws leave the user's stream at the same point
  set.seed(31); a <- scr_metrics(s, y, n_boot = 20); ra <- .Random.seed
  set.seed(31); b <- ref_metrics(s, y, n_boot = 20); rb <- .Random.seed
  expect_identical(a, b); expect_identical(ra, rb)
})

test_that("on the rows the bootstrap is seeded locally and does not depend on nthread", {
  set.seed(5)
  y <- stats::rbinom(600, 1, 0.3)
  s <- y + stats::rnorm(600)
  set.seed(42); before <- .Random.seed
  m1 <- scr_metrics(s, y, n_boot = 40, seed = 9)
  expect_identical(.Random.seed, before)
  expect_identical(m1, scr_metrics(s, y, n_boot = 40, seed = 9))
  expect_identical(m1, scr_metrics(s, y, n_boot = 40, seed = 9, nthread = 2L))
  expect_identical(m1, ref_metrics(s, y, n_boot = 40, seed = 9, nthread = 2L))
  expect_false(identical(m1$auc_lo, scr_metrics(s, y, n_boot = 40, seed = 10)$auc_lo))
  # the backend does not matter either
  withr::with_options(list(scorecraft.parallel = "serial"),
                      expect_identical(m1, scr_metrics(s, y, n_boot = 40, seed = 9, nthread = 2L)))
})

test_that("on the counts the interval bounds have the law of the row bootstrap", {
  set.seed(77)
  n <- 400
  s <- round(stats::rnorm(n) * 3)
  y <- stats::rbinom(n, 1, stats::plogis(-1 + 0.35 * s))
  expect_lte(data.table::uniqueN(s), n / 10)
  R <- 300; B <- 100
  bd <- c("auc_lo", "auc_hi", "ks_lo", "ks_hi")
  new <- vapply(seq_len(R), function(r) unlist(scr_metrics(s, y, n_boot = B, seed = r)[bd]), numeric(4))
  old <- vapply(seq_len(R), function(r) unlist(ref_metrics(s, y, n_boot = B, seed = 5000 + r)[bd]), numeric(4))
  for (k in bd) {
    sd_new <- stats::sd(new[k, ]); sd_old <- stats::sd(old[k, ])
    # two independent means of R bounds: four standard errors of their difference
    expect_lt(abs(mean(new[k, ]) - mean(old[k, ])), 4 * sqrt((sd_new^2 + sd_old^2) / R))
    # the sd of R draws has a relative error near 1 / sqrt(2 (R - 1)); the
    # bounds are not normal, so the ratio of two of them is given 0.25
    expect_lt(abs(sd_new / sd_old - 1), 0.25)
  }
  # the bounds bracket the point estimates
  m <- scr_metrics(s, y, ci = FALSE)
  expect_true(all(new["auc_lo", ] < m$auc & m$auc < new["auc_hi", ]))
})

test_that("on the counts the bootstrap is seeded locally and ignores nthread", {
  set.seed(5)
  y <- stats::rbinom(600, 1, 0.3)
  s <- round(2 * y + 3 * stats::rnorm(600))
  expect_lte(data.table::uniqueN(s), 60)
  set.seed(42); before <- .Random.seed
  m1 <- scr_metrics(s, y, n_boot = 40, seed = 9)
  expect_identical(.Random.seed, before)
  expect_identical(m1, scr_metrics(s, y, n_boot = 40, seed = 9))
  expect_identical(m1, scr_metrics(s, y, n_boot = 40, seed = 9, nthread = 2L))
  expect_identical(m1, scr_metrics(s, y, n_boot = 40, seed = 9, nthread = 8L))
  expect_false(identical(m1$auc_lo, scr_metrics(s, y, n_boot = 40, seed = 10)$auc_lo))
  # a call without the session's random stream leaves none behind
  rm(".Random.seed", envir = globalenv())
  invisible(scr_metrics(s, y, n_boot = 20, seed = 9))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  # without a seed the resamples come from the user's stream
  set.seed(42); a <- scr_metrics(s, y, n_boot = 40); after <- .Random.seed
  set.seed(42); b <- scr_metrics(s, y, n_boot = 40)
  expect_identical(a, b)
  expect_false(identical(after, before))
  # the level reaches the percentiles; fewer than two resamples give no interval
  w <- scr_metrics(s, y, n_boot = 200, seed = 9, level = 0.99); n9 <- scr_metrics(s, y, n_boot = 200, seed = 9, level = 0.8)
  expect_true(w$auc_lo < n9$auc_lo && n9$auc_hi < w$auc_hi)
  one <- scr_metrics(s, y, n_boot = 1, seed = 9)
  expect_true(is.na(one$auc_lo) && one$n_boot == 0L)
})

test_that("on the counts the bootstrap is exact above the pooling size of the studies", {
  set.seed(8)
  n <- 120000
  s <- sample.int(12000, n, replace = TRUE)
  y <- stats::rbinom(n, 1, stats::plogis(-1 + (s - 6000) / 3000))
  K <- data.table::uniqueN(s)
  expect_gt(K, 1e4); expect_lte(K, n / 10)
  m <- scr_metrics(s, y, n_boot = 30, seed = 3)
  ex <- ref_metrics_counts(s, y, 30, seed = 3)
  expect_identical(unlist(unclass(m)[c("auc", "ks", bounds)]), unlist(ex[c("auc", "ks", bounds)]))
  # the pooled resamples of the studies are another draw
  idx <- data.table::frank(s, ties.method = "dense")
  po <- .study_auc_boot(tabulate(idx[y == 1], K), tabulate(idx[y == 0], K), 30L, 0.95, seed = 3)
  expect_false(identical(m$auc_lo, po$auc_lo))
})

test_that("a chunk of one resample goes through the vector kernel with the same values", {
  set.seed(21)
  K <- 1000L
  c1 <- as.double(stats::rpois(K, 3)); c0 <- as.double(stats::rpois(K, 20))
  n1 <- sum(c1); n0 <- sum(c0)
  # 1001 resamples over 1000 cells run as chunks of 1000 and of 1
  full <- .study_auc_boot(c1, c0, 1001L, seed = 4, keep = TRUE)
  # the same draws by hand, the last one through the column kernel
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  set.seed(4)
  A1 <- stats::rmultinom(1000, n1, c1 / n1); A0 <- stats::rmultinom(1000, n0, c0 / n0)
  B1 <- stats::rmultinom(1, n1, c1 / n1);    B0 <- stats::rmultinom(1, n0, c0 / n0)
  if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv())
  first <- .study_auc_cols(A1, A0); last <- .study_auc_cols(B1, B0)
  expect_identical(full$boot_auc, c(first$auc, last$auc))
  expect_identical(full$boot_ks, c(first$ks, last$ks))
  # and the vector kernel agrees with the column kernel on that single draw
  v <- .auc_ks_counts(B1, B0)
  expect_identical(c(v$auc, v$ks), c(last$auc, last$ks))
})
