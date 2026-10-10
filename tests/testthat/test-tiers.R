# Tiers: the dynamic program against brute force, the fallback, anchors,
# rounding, stability and the production paths.

tier_df <- function(n = 6000, seed = 17, shift = 0) {
  set.seed(seed)
  x <- stats::rnorm(n)
  data.frame(score = round(600 + 40 * (x + shift)), y = stats::rbinom(n, 1, stats::plogis(-1.8 - 1.1 * x)),
             smp = sample(c("dev", "oot"), n, TRUE), stringsAsFactors = FALSE)
}

# Brute force over every contiguous segmentation, written independently of
# the C++ program (stats::fisher.test for the adjacency test)
tier_brute <- function(e, n, vol, er, nr, L, criterion, min_share, min_events, alpha) {
  M <- length(e)
  combos <- if (L == 1L) list(integer()) else utils::combn(M - 1L, L - 1L, simplify = FALSE)
  E <- sum(e); NE <- sum(n) - E
  vals <- vapply(combos, function(cb) {
    ends <- c(cb, M); starts <- c(1L, cb + 1L)
    seg <- function(v) mapply(function(a, b) sum(v[a:b]), starts, ends)
    se <- seg(e); sn <- seg(n); sv <- seg(vol); ser <- seg(er); snr <- seg(nr)
    if (any(sv < (min_share - 1e-12) * sum(vol)) || any(ser < min_events) || any(snr - ser < min_events)) return(-Inf)
    if (L > 1L) for (t in 2:L) {
      tab <- matrix(c(ser[t], snr[t] - ser[t], ser[t - 1], snr[t - 1] - ser[t - 1]), 2, byrow = TRUE)
      p <- if (ser[t] + ser[t - 1] == 0) 1 else stats::fisher.test(tab, alternative = "greater")$p.value
      if (!(p < alpha / (L - 1L))) return(-Inf)
    }
    if (criterion == 0L) {
      xl <- function(x, m) ifelse(x > 0, x * log(x / m), 0)
      sum(xl(se, sn) + xl(sn - se, sn))
    } else {
      pe <- se / E; pn <- (sn - se) / NE
      if (any(pe <= 0) || any(pn <= 0)) return(-Inf)
      sum((pe - pn) * log(pe / pn))
    }
  }, numeric(1))
  best <- max(vals)
  if (!is.finite(best)) return(list(feasible = FALSE))
  win <- which(abs(vals - best) <= 1e-9 * max(1, abs(best)))
  list(feasible = TRUE, objective = best, unique = length(win) == 1L, ends = c(combos[[win[1]]], M))
}

test_that("the tier program equals brute force on small cases, both criteria, with and without constraints", {
  set.seed(101)
  checked <- 0L
  for (M in c(5L, 8L, 12L)) for (rep in 1:3) {
    nr <- as.double(sample(30:200, M, TRUE))
    rate <- sort(stats::runif(M, 0.02, 0.6))
    er <- as.double(stats::rbinom(M, nr, rate))
    # weighted counts differ from the raw ones, as with case weights
    wt <- stats::runif(M, 0.5, 2)
    e <- er * wt; n <- nr * wt; vol <- n * stats::runif(M, 1, 1.1)
    for (L in 2:5) for (crit in 0:1) for (cons in c(FALSE, TRUE)) {
      if (L > M) next
      ms <- if (cons) 0.08 else 0; me <- if (cons) 8 else 0; al <- if (cons) 0.05 else 1e6
      got <- cpp_tier_dp(e, n, vol, er, nr, L, crit, ms, me, al)
      bf <- tier_brute(e, n, vol, er, nr, L, crit, ms, me, al)
      expect_identical(got$feasible, bf$feasible)
      if (bf$feasible) {
        expect_equal(got$objective, bf$objective, tolerance = 1e-10)
        if (bf$unique) expect_identical(as.integer(got$ends), as.integer(bf$ends))
        checked <- checked + 1L
      }
    }
  }
  expect_gt(checked, 50L)
  # fewer blocks than tiers, or an empty input: infeasible
  expect_false(cpp_tier_dp(c(1, 2), c(10, 10), c(10, 10), c(1, 2), c(10, 10), 3L, 0L, 0, 0, 1)$feasible)
  expect_error(cpp_tier_dp(c(1, 2), c(10), c(10, 10), c(1, 2), c(10, 10), 2L, 0L, 0, 0, 1), "differ in length")
})

test_that("pool adjacent violators matches an independent implementation", {
  pav_ref <- function(e, n) {
    groups <- as.list(seq_along(e))
    repeat {
      r <- vapply(groups, function(g) sum(e[g]) / sum(n[g]), 1)
      v <- which(diff(r) < 0)
      if (!length(v)) break
      i <- v[1]; groups[[i]] <- c(groups[[i]], groups[[i + 1L]]); groups[[i + 1L]] <- NULL
    }
    g <- integer(length(e)); for (k in seq_along(groups)) g[groups[[k]]] <- k
    g
  }
  set.seed(4)
  for (r in 1:20) {
    n <- as.double(sample(5:50, 15, TRUE)); e <- as.double(stats::rbinom(15, n, stats::runif(15)))
    expect_identical(.study_pav(e, n), pav_ref(e, n))
  }
  # a block without a known outcome joins its neighbor
  expect_identical(.study_pav(c(1, 0, 3), c(10, 0, 10)), c(1L, 1L, 2L))
})

test_that("scr_tiers optimal: monotone, distinct, labeled by event rate, event-richest first", {
  d <- tier_df()
  tr <- scr_tiers(d, sample = "smp", n_tiers = 5)
  expect_s3_class(tr, c("scr_study_tiers", "scr_study"))
  expect_identical(tr$n_tiers, 5L)
  expect_identical(tr$labels, c("very low", "low", "medium", "high", "very high"))
  ref <- tr$table[sample == "dev"]
  expect_identical(ref$tier, 5:1)
  expect_identical(ref$label, rev(tr$labels))
  expect_true(all(diff(ref$rate) < 0))
  expect_true(tr$summary[sample == "dev", monotone])
  expect_true(tr$summary[sample == "dev", all_distinct])
  # constraints hold on the reference
  expect_true(all(ref$pct >= 0.05))
  expect_true(all(ref$events >= 20 & ref$n - ref$events >= 20))
  # every cut is a pre-bin cut, hence never an observed reference score
  expect_false(any(tr$cuts %in% d$score[d$smp == "dev"]))
  # p_adjacent: one-sided Fisher exact test of a tier against the next lower one
  for (i in 1:4) {
    ft <- stats::fisher.test(matrix(c(ref$events[i], ref$n[i] - ref$events[i], ref$events[i + 1],
                                      ref$n[i + 1] - ref$events[i + 1]), 2, byrow = TRUE), alternative = "greater")$p.value
    expect_equal(ref$p_adjacent[i], ft, tolerance = 1e-10)
  }
  expect_true(is.na(ref$p_adjacent[5]))
  expect_identical(tr$measure, "risk")
  # three and seven tiers have their labels; IV gives a fit as well
  expect_identical(scr_tiers(d, sample = "smp", n_tiers = 3)$labels, c("low", "medium", "high"))
  t7 <- scr_tiers(d, sample = "smp", n_tiers = 7, criterion = "iv", min_pct = 0.02)
  expect_identical(t7$labels[c(1, 7)], c("extremely low", "extremely high"))
  expect_identical(t7$criterion, "iv")
  # even counts have no "medium"; 8 and 9 tiers are numbered
  expect_identical(scr_tiers(d, sample = "smp", n_tiers = 4)$labels, c("low", "medium low", "medium high", "high"))
  expect_identical(scr_tiers(d, sample = "smp", n_tiers = 2)$labels, c("low", "high"))
  expect_identical(.tier_labels(6), c("very low", "low", "medium low", "medium high", "high", "very high"))
  expect_identical(.tier_labels(8), paste0("T", 1:8))
  expect_identical(.tier_labels(9), paste0("T", 1:9))
  expect_identical(.tier_labels(1), "T1")
  t4 <- scr_tiers(d, sample = "smp", n_tiers = 4)
  expect_output(print(t4), "medium high")
  expect_identical(scr_tiers(d, sample = "smp", n_tiers = 3, labels = c("A", "B", "C"))$labels, c("A", "B", "C"))
  expect_error(scr_tiers(d, sample = "smp", n_tiers = 3, labels = c("A", "B")), "tiers were achieved")
  expect_output(print(tr), "scr_study_tiers")
})

test_that("an infeasible tier count falls back by two, is recorded and warns once", {
  set.seed(8)
  d <- data.frame(score = round(stats::rnorm(900) * 10))
  d$y <- stats::rbinom(900, 1, stats::plogis(-3 - 0.15 * d$score))
  expect_warning(tr <- scr_tiers(d, n_tiers = 5, min_events = 12), "5 tiers requested, 3 achieved")
  expect_identical(tr$n_tiers, 3L)
  expect_identical(tr$labels, c("low", "medium", "high"))
  led <- tr$ledger[step == "optimal"]
  expect_identical(led$n_tiers, c(5L, 3L))
  expect_identical(led$status, c("infeasible", "fitted"))
  # nothing feasible: one tier, recorded
  expect_warning(t1 <- scr_tiers(data.frame(score = 1:100, y = 0), n_tiers = 3, min_events = 0), "1 achieved")
  expect_identical(t1$n_tiers, 1L)
  expect_identical(t1$ledger$status[nrow(t1$ledger)], "single")
  expect_identical(t1$labels, "T1")
})

test_that("anchored cuts sit where the pooled block rate first reaches each anchor", {
  d <- tier_df(8000)
  anchors <- c(0.06, 0.25)
  ta <- scr_tiers(d, sample = "smp", method = "anchored", anchors = anchors, max_bins = 40)
  expect_identical(ta$n_tiers, 3L)
  # independent reconstruction: pre-bins, pooled until monotone, first block at or above the anchor
  h <- .study_cells(ta$hist, "dev")
  pre <- .tier_prebins(h, 40, "low")
  blk <- .tier_blocks(pre)
  rate <- blk$e / blk$n
  expect_true(all(diff(rate) >= 0))
  for (a in anchors) {
    j <- which(rate >= a)[1]
    expect_true(pre$bcut[blk$last[j - 1L]] %in% ta$cuts)
    expect_lt(rate[j - 1L], a)
  }
  # "overall" is the reference event rate
  ov <- scr_tiers(d, sample = "smp", method = "anchored", anchors = "overall", max_bins = 40)
  r0 <- mean(d$y[d$smp == "dev"])
  j <- which(rate >= r0)[1]
  expect_equal(ov$cuts, pre$bcut[blk$last[j - 1L]])
  # conservative: the lower Jeffreys bound reaches the anchor no earlier
  tc <- scr_tiers(d, sample = "smp", method = "anchored", anchors = anchors, max_bins = 40, conservative = TRUE)
  expect_true(all(tc$cuts <= ta$cuts))   # higher_is_safer: risk rises toward lower scores
  # an anchor above every rate gives no cut, recorded
  expect_warning(tn <- scr_tiers(d, sample = "smp", method = "anchored", anchors = c(0.1, 0.99)), "3 tiers requested, 2 achieved")
  expect_true(any(tn$ledger$status == "no cut"))
  expect_error(scr_tiers(d, method = "anchored"), "needs `anchors`")
  expect_error(scr_tiers(d, method = "anchored", anchors = 1.5), "anchors")
})

test_that("quantile tiers have equal shares; round_to re-evaluates every sample", {
  d <- tier_df()
  tq <- scr_tiers(d, sample = "smp", method = "quantile", n_tiers = 4)
  expect_equal(tq$table[sample == "dev", pct], rep(0.25, 4), tolerance = 0.03)
  expect_true(is.na(tq$criterion))
  tr <- scr_tiers(d, sample = "smp", n_tiers = 3, round_to = 25)
  expect_true(all(tr$cuts %% 25 == 0))
  expect_equal(tr$cuts, round(tr$cuts_raw / 25) * 25)
  for (nm in c("dev", "oot")) {
    s <- d$score[d$smp == nm]; y <- d$y[d$smp == nm]
    i <- findInterval(s, tr$cuts) + 1L
    t <- tr$table[sample == nm]
    tier <- tr$codes[i]
    expect_equal(t$n, as.numeric(vapply(t$tier, function(k) sum(tier == k), 1L)))
    expect_equal(t$events, as.numeric(vapply(t$tier, function(k) sum(y[tier == k]), 1L)))
  }
  expect_identical(tr$ledger$status[nrow(tr$ledger)], "rounded")
  # rounding that merges two cuts warns and relabels
  expect_warning(tm <- scr_tiers(d, sample = "smp", n_tiers = 5, round_to = 1000), "merged")
  expect_lt(tm$n_tiers, 5L)
})

test_that("rounded cuts on pooled scores: the table counts equal scr_apply()", {
  set.seed(19)
  n <- 4e5
  x <- stats::rnorm(n)
  d <- data.frame(score = 600 + 50 * x, y = stats::rbinom(n, 1, stats::plogis(-1.8 - 1.1 * x)),
                  smp = sample(c("dev", "oot"), n, TRUE))
  tr <- scr_tiers(d, sample = "smp", n_tiers = 5, max_cells = 1e4, round_to = 7)
  expect_true(tr$quantized)
  expect_true(any(tr$ledger$status == "recounted"))
  expect_false(any(tr$cuts %in% attr(.study_hist(d$score, d$y, max_cells = 1e4), "edges")))
  a <- scr_apply(tr, d)
  for (nm in c("dev", "oot")) {
    t <- tr$table[sample == nm]
    i <- a$smp == nm
    expect_equal(t$n, as.numeric(vapply(t$tier, function(k) sum(a$tier[i] == k), 1L)))
    expect_equal(t$events, as.numeric(vapply(t$tier, function(k) sum(d$y[i][a$tier[i] == k]), 1L)))
  }
  # the fit itself is on the pooled cells: raw cuts are cell edges
  expect_true(all(tr$cuts_raw %in% attr(.study_hist(d$score, d$y, max_cells = 1e4), "edges")))
})

test_that("tier stability is reproducible and leaves the user's stream alone", {
  d <- tier_df(3000)
  set.seed(1); r1 <- stats::runif(2)
  set.seed(1)
  s1 <- scr_tiers(d, sample = "smp", n_tiers = 3, max_bins = 30, n_boot = 15, seed = 9)
  r2 <- stats::runif(2)
  expect_identical(r1, r2)
  s2 <- scr_tiers(d, sample = "smp", n_tiers = 3, max_bins = 30, n_boot = 15, seed = 9)
  expect_identical(s1$stability, s2$stability)
  st <- s1$stability
  expect_identical(st$n_boot, 15L)
  expect_identical(nrow(st$cuts), 2L)
  expect_true(all(st$cuts$q25 <= st$cuts$median & st$cuts$median <= st$cuts$q75))
  expect_true(st$agreement > 0.5 && st$agreement <= 1)
  expect_null(scr_tiers(d, sample = "smp", n_tiers = 3, max_bins = 30)$stability)
})

test_that("tiers under propensity: tier 1 at the low scores, labels in the measure", {
  d <- tier_df()
  d$score <- -d$score
  tp <- scr_tiers(d, sample = "smp", objective = "propensity", n_tiers = 3)
  expect_identical(tp$measure, "propensity")
  expect_identical(tp$codes, 1:3)
  t <- tp$table[sample == "dev"]
  expect_true(all(diff(t$score_lo) < 0))
  expect_identical(t$tier, 3:1)
  expect_true(all(diff(t$rate) < 0))
})

test_that("scr_tiers on a scorecard reads the config keys", {
  sc <- sc_demo()
  sc2 <- sc; sc2$config$tier_min_pct <- 0.15; sc2$config$tier_max_bins <- 20L
  tr <- scr_tiers(sc2, n_tiers = 3)
  expect_equal(tr$min_pct, 0.15)
  expect_true(all(tr$table[sample == "train", pct] >= 0.15))
  expect_identical(unique(tr$table$sample), c("train", "holdout"))
  expect_error(scr_tiers(sc, n_tiers = 1), "n_tiers")
  expect_error(scr_tiers(sc, n_tiers = 10), "`n_tiers` must be at most 9")
  expect_error(scr_tiers(sc, max_bins = 1000), "`max_bins` must be at most 500")
  expect_error(scr_tiers(sc, foo = 1), "unused argument")
})

test_that("R and SQL tiers agree on DuckDB, NULL score included", {
  skip_if_not_installed("duckdb"); skip_if_not_installed("DBI")
  old <- scr_verbose(FALSE); on.exit(scr_verbose(old), add = TRUE)
  sc <- sc_demo()
  tr <- scr_tiers(sc, n_tiers = 5, round_to = 5)
  d <- data.frame(id = seq_len(nrow(scr_demo)), score = scr_apply(sc, scr_demo)$score)
  d$score[c(2, 50)] <- NA
  d$score[3] <- tr$cuts[2]
  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  DBI::dbWriteTable(con, "scores", d)
  sql <- scr_sql(tr, table = "scores", dialect = "duckdb")
  expect_true(any(grepl("IS NULL THEN NULL", sql)))
  got <- DBI::dbGetQuery(con, paste(sql, collapse = "\n"))
  exp <- scr_apply(tr, d)
  expect_identical(as.integer(got$tier), exp$tier)
  expect_identical(got$tier_label, exp$tier_label)
  expect_true(is.na(got$tier[2]))
  # the labels carry their order in both: credit with 5 tiers, tier 5 is "01.very high"
  expect_true(any(grepl("'01.very high'", sql, fixed = TRUE)))
  expect_true(any(grepl("ELSE '05.very low'", sql, fixed = TRUE)))
  expect_identical(got$tier_label[!is.na(got$tier)], tr$tier_labels[got$tier[!is.na(got$tier)]])
  expect_true(is.na(got$tier_label[2]))
  # numbered = FALSE: the plain labels, the same in R and SQL
  sql0 <- scr_sql(tr, table = "scores", dialect = "duckdb", numbered = FALSE)
  expect_false(any(grepl("'0[0-9]\\.", sql0)))
  expect_true(any(grepl("ELSE 'very low'", sql0, fixed = TRUE)))
  got0 <- DBI::dbGetQuery(con, paste(sql0, collapse = "\n"))
  exp0 <- scr_apply(tr, d, numbered = FALSE)
  expect_identical(got0$tier_label, exp0$tier_label)
  expect_identical(as.integer(got0$tier), exp$tier)
  expect_identical(exp0$tier_label, tr$code_labels[findInterval(d$score, tr$cuts) + 1L])
  # a score on a cut is on the upper side
  expect_identical(exp$tier[3], tr$codes[3])
  tmp <- tempfile(fileext = ".sql")
  expect_invisible(scr_sql(tr, file = tmp))
  expect_true(any(grepl("your_table", readLines(tmp))))
})

test_that("scr_export writes the tiers, the ledger and the stability", {
  skip_if_not_installed("openxlsx")
  old <- scr_verbose(FALSE); on.exit(scr_verbose(old), add = TRUE)
  tr <- scr_tiers(tier_df(2000), sample = "smp", n_tiers = 3, max_bins = 20, n_boot = 5, seed = 1)
  out <- file.path(tempdir(), "scr-tiers-export")
  unlink(out, recursive = TRUE)
  ex <- scr_export(tr, out, stamp = FALSE)
  expect_identical(openxlsx::getSheetNames(ex$files$xlsx),
                   c("Summary", "Tiers", "Cuts", "Settings", "Ledger", "Stability"))
})

test_that("production labels carry their order: 01 is the event-richest tier under every direction", {
  fx <- list(credit = sc_demo(), fraud = sc_fraud_demo(), propensity = sc_prop_demo())
  for (nm in names(fx)) {
    sc <- fx[[nm]]
    tr <- suppressWarnings(scr_tiers(sc, n_tiers = 5))
    L <- tr$n_tiers
    expect_gte(L, 2L)
    expect_identical(tr$tier_labels, sprintf("%02d.%s", rev(seq_len(L)), tr$labels))
    for (smp in c("train", "holdout")) {
      t <- tr$table[sample == smp]
      # the table lists the event-richest tier first: its rows read 01, 02, ..., and sort that way
      expect_identical(t$tier_label, sprintf("%02d.%s", seq_len(L), t$label))
      expect_identical(sort(t$tier_label), t$tier_label)
      expect_identical(t$tier, rev(seq_len(L)))   # the tier number still rises with the event rate
    }
    ref <- tr$table[sample == "train"]
    expect_identical(which.max(ref$rate), 1L)
    expect_match(ref$tier_label[1], "^01\\.")
    # scr_apply: numbered by default, plain on request, and joins to the table
    ho <- sc$samples$holdout
    a <- scr_apply(tr, ho)
    expect_identical(a$tier, tr$codes[findInterval(ho$score, tr$cuts) + 1L])
    expect_identical(a$tier_label, tr$tier_labels[a$tier])
    expect_identical(scr_apply(tr, ho, numbered = FALSE)$tier_label, tr$labels[a$tier])
    th <- tr$table[sample == "holdout"]
    j <- match(a$tier_label, th$tier_label)
    expect_false(anyNA(j))
    expect_identical(th$tier[j], a$tier)
    expect_equal(as.numeric(table(factor(a$tier_label, th$tier_label))), th$n)
    # the event-rich end of the score gets "01.": the low scores for credit, the high ones otherwise
    rich <- if (identical(tr$direction, "higher_is_safer")) min(ho$score) else max(ho$score)
    expect_match(scr_apply(tr, rich)$tier_label, "^01\\.")
    sql <- scr_sql(tr, dialect = "ansi")
    expect_true(any(grepl(sprintf("'%s'", tr$tier_labels[L]), sql, fixed = TRUE)))
  }
  # credit and fraud-like are mirror scales of the same model: the same labels, opposite score ends
  expect_identical(scr_tiers(fx$credit, n_tiers = 3)$tier_labels, c("03.low", "02.medium", "01.high"))
  expect_identical(scr_tiers(fx$fraud, n_tiers = 3)$tier_labels, c("03.low", "02.medium", "01.high"))
})

test_that("numbered labels for 2 to 9 tiers, user labels and older objects", {
  d <- tier_df()
  for (L in 2:9) {
    tq <- scr_tiers(d, sample = "smp", method = "quantile", n_tiers = L)
    expect_identical(tq$n_tiers, L)
    expect_identical(tq$tier_labels, sprintf("%02d.%s", L:1, .tier_labels(L)))
    t <- tq$table[sample == "dev"]
    expect_identical(t$tier_label[1], paste0("01.", .tier_labels(L)[L]))
    expect_identical(t$tier_label[L], sprintf("%02d.%s", L, .tier_labels(L)[1]))
    a <- scr_apply(tq, d)
    expect_setequal(unique(a$tier_label), t$tier_label)
    expect_identical(a$tier_label, tq$tier_labels[a$tier])
  }
  # the prefix is padded to at least two digits, more when the count needs it
  expect_identical(.tier_numbered(c("a", "b")), c("02.a", "01.b"))
  expect_identical(.tier_numbered(letters[1:12])[c(1, 12)], c("12.a", "01.l"))
  expect_identical(.tier_numbered(as.character(1:100))[c(1, 100)], c("100.1", "001.100"))
  expect_identical(.tier_numbered(character()), character())
  # user labels get the prefix too
  tu <- scr_tiers(d, sample = "smp", n_tiers = 3, labels = c("A", "B", "C"))
  expect_identical(tu$tier_labels, c("03.A", "02.B", "01.C"))
  expect_identical(tu$table[sample == "dev", tier_label], c("01.C", "02.B", "03.A"))
  expect_identical(sort(unique(scr_apply(tu, d)$tier_label)), c("01.C", "02.B", "03.A"))
  expect_true(any(grepl("ELSE '03.A'", scr_sql(tu), fixed = TRUE)))
  expect_setequal(unique(scr_apply(tu, d, numbered = FALSE)$tier_label), c("A", "B", "C"))
  expect_output(print(tu), "01.C")
  # a study fitted before the numbered labels existed still gets them
  old <- tu; old$tier_labels <- NULL; old$table <- old$table[, -"tier_label"]
  expect_identical(scr_apply(old, d)$tier_label, scr_apply(tu, d)$tier_label)
  expect_identical(scr_sql(old)[-3], scr_sql(tu)[-3])
  expect_output(print(old), "scr_study_tiers")
  # bands keep their interval labels, with or without the argument
  b <- scr_bands(d, sample = "smp", n_bands = 4, n_boot = 0)
  expect_identical(scr_apply(b, d)$tier_label, scr_apply(b, d, numbered = FALSE)$tier_label)
  expect_identical(scr_apply(b, d)$tier_label, b$code_labels[findInterval(d$score, b$cuts) + 1L])
  expect_identical(scr_sql(b)[-3], scr_sql(b, numbered = FALSE)[-3])
  expect_error(scr_apply(tu, d, numbered = NA), "`numbered`")
  expect_error(scr_sql(tu, numbered = "yes"), "`numbered`")
})
