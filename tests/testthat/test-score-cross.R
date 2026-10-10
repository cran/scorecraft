# Two scores on the same rows: the cross table against table(), the rank
# association against cor(), the overlap sets against set operations.

cross_df <- function(n = 3000, seed = 61) {
  set.seed(seed)
  z <- stats::rnorm(n)
  data.frame(a = round(600 + 40 * (-z + stats::rnorm(n, sd = 0.7))), b = round(450 + 30 * (0.5 * z + stats::rnorm(n))),
             y = stats::rbinom(n, 1, stats::plogis(-1.8 + z)), y2 = stats::rbinom(n, 1, stats::plogis(-0.8 + 0.4 * z)),
             w = stats::runif(n, 0.2, 3), stringsAsFactors = FALSE)
}

test_that("the cross table equals table() and tapply() over the bands", {
  d <- cross_df()
  cx <- scr_score_cross(d, "a", "b", y = "y", objective_b = "propensity", n_bands = 4)
  expect_s3_class(cx, "scr_score_cross")
  st <- cx$settings
  ia <- st$codes_a[findInterval(d$a, st$cuts_a) + 1L]
  ib <- st$codes_b[findInterval(d$b, st$cuts_b) + 1L]
  # band 1 is the event-richest: the low scores of A (higher_is_safer), the high scores of B
  expect_identical(st$codes_a, 1:4); expect_identical(st$codes_b, 4:1)
  t <- cx$table
  cell <- t[!is.na(band_a) & !is.na(band_b)]
  tab <- table(factor(ia, 1:4), factor(ib, 1:4))
  expect_equal(cell$n, as.numeric(tab[cbind(cell$band_a, cell$band_b)]))
  ev <- tapply(d$y, list(factor(ia, 1:4), factor(ib, 1:4)), sum)
  ev[is.na(ev)] <- 0
  expect_equal(cell$events, as.numeric(ev[cbind(cell$band_a, cell$band_b)]))
  expect_equal(cell$rate, ifelse(cell$n > 0, cell$events / cell$n, NA_real_))
  ok <- cell$n > 0
  expect_equal(cell$rate_lo[ok], ifelse(cell$events[ok] == 0, 0, stats::qbeta(0.025, cell$events[ok] + 0.5, cell$n[ok] - cell$events[ok] + 0.5)))
  expect_equal(cell$lift, cell$rate / mean(d$y))
  expect_equal(cell$pct, cell$n / nrow(d))
  # totals: the margins of table() and the grand total
  ra <- t[!is.na(band_a) & is.na(band_b)]; rb <- t[is.na(band_a) & !is.na(band_b)]; all_ <- t[is.na(band_a) & is.na(band_b)]
  expect_equal(ra$n, as.numeric(rowSums(tab))); expect_equal(rb$n, as.numeric(colSums(tab)))
  expect_identical(ra$label_b, rep("total", 4)); expect_identical(rb$label_a, rep("total", 4))
  expect_equal(all_$n, nrow(d)); expect_equal(all_$rate, mean(d$y)); expect_equal(all_$lift, 1)
  expect_identical(cell$label_a, .study_labels(st$cuts_a)[match(cell$band_a, st$codes_a)])
  # equal-share bands, tie-safe: the same cuts as scr_bands()
  expect_equal(st$cuts_a, scr_bands(d, score = "a", n_bands = 4, n_boot = 0)$cuts)
  expect_equal(st$cuts_b, scr_bands(d, score = "b", objective = "propensity", n_bands = 4, n_boot = 0)$cuts)
})

test_that("Kendall's tau-b and Spearman's rho equal cor(), ties included", {
  set.seed(5)
  for (r in 1:25) {
    n <- sample(c(2:12, 40, 150), 1)
    a <- sample(1:4, n, TRUE) + if (r %% 3 == 0) stats::rnorm(n) else 0
    b <- if (r %% 2 == 0) sample(1:3, n, TRUE) else round(a + stats::rnorm(n), 1)
    if (stats::sd(a) == 0 || stats::sd(b) == 0) next
    expect_equal(.scr_kendall_tau_b(a, b), stats::cor(a, b, method = "kendall"), tolerance = 1e-12)
  }
  expect_true(is.na(.scr_kendall_tau_b(c(1, 1, 1), c(1, 2, 3))))
  expect_true(is.na(.scr_kendall_tau_b(1, 2)))
  d <- cross_df(800)
  d$a <- round(d$a / 20) * 20   # heavy ties
  cx <- scr_score_cross(d, "a", "b", y = "y", n_bands = 3)
  as <- cx$association
  expect_equal(as$estimate[as$method == "kendall_tau_b"], stats::cor(d$a, d$b, method = "kendall"), tolerance = 1e-12)
  expect_equal(as$estimate[as$method == "spearman"], stats::cor(d$a, d$b, method = "spearman"), tolerance = 1e-12)
  expect_equal(as$n, c(800L, 800L))
  # oriented: A safer (event-rich at the low end), B riskier by default under risk with direction given
  cr <- scr_score_cross(d, "a", "b", direction_b = "higher_is_riskier")
  expect_equal(cr$association$oriented, -cr$association$estimate)
  cs <- scr_score_cross(d, "a", "b")
  expect_equal(cs$association$oriented, cs$association$estimate)
})

test_that("overlap and swap sets equal manual set operations at each depth", {
  d <- cross_df()
  cx <- scr_score_cross(d, "a", "b", y = "y", objective_b = "propensity", depths = c(0.2, 0.05, 0.1, 1))
  ov <- cx$overlap
  expect_equal(ov$depth, c(0.05, 0.1, 0.2, 1))
  # the share selected is the achievable share nearest the target (tie-safe), counted from the event-rich end
  ach <- function(s, high) {
    tb <- table(s); cnt <- as.numeric(tb)
    cumsum(if (high) rev(cnt) else cnt) / length(s)
  }
  sa <- ach(d$a, FALSE); sb <- ach(d$b, TRUE)
  for (k in 1:3) {
    t <- ov$depth[k]
    expect_equal(ov$share_a[k], sa[-length(sa)][which.min(abs(sa[-length(sa)] - t))])
    expect_equal(ov$share_b[k], sb[-length(sb)][which.min(abs(sb[-length(sb)] - t))])
  }
  r <- cx$overlap_rates
  for (k in seq_len(nrow(ov))) {
    A <- d$a < ov$cut_a[k]    # A is safer: its event-rich end is the low scores
    B <- d$b >= ov$cut_b[k]   # B is riskier: the high scores
    sets <- list(A = A, B = B, both = A & B, "A only" = A & !B, "B only" = B & !A)
    expect_equal(c(ov$n_a[k], ov$n_b[k], ov$n_both[k], ov$n_a_only[k], ov$n_b_only[k]),
                 vapply(sets, sum, 1, USE.NAMES = FALSE))
    expect_equal(ov$jaccard[k], sum(A & B) / sum(A | B))
    rk <- r[r$depth == ov$depth[k]]
    expect_identical(rk$set, names(sets))
    for (s in names(sets)) {
      i <- sets[[s]]
      expect_equal(rk$events[rk$set == s], sum(d$y[i]))
      expect_equal(rk$rate[rk$set == s], if (any(i)) mean(d$y[i]) else NA_real_)
      x <- sum(d$y[i]); n <- sum(i)
      if (n > 0 && x > 0 && x < n) expect_equal(rk$rate_hi[rk$set == s], stats::qbeta(0.975, x + 0.5, n - x + 0.5))
    }
  }
  # depth 1 selects every row with both scores
  expect_equal(ov$n_both[4], nrow(d)); expect_equal(ov$jaccard[4], 1)
  expect_equal(ov$cut_a[4], Inf); expect_equal(ov$cut_b[4], -Inf)
})

test_that("two outcomes give suffixed columns and one rate table per outcome", {
  d <- cross_df()
  cx <- scr_score_cross(d, "a", "b", y_a = "y", y_b = "y2", objective_b = "propensity", n_bands = 3)
  t <- cx$table
  expect_true(all(c("events_a", "rate_a", "rate_lo_a", "rate_hi_a", "lift_a", "events_b", "rate_b", "lift_b") %in% names(t)))
  expect_false("rate" %in% names(t))
  all_ <- t[is.na(band_a) & is.na(band_b)]
  expect_equal(all_$rate_a, mean(d$y)); expect_equal(all_$rate_b, mean(d$y2))
  st <- cx$settings
  ia <- st$codes_a[findInterval(d$a, st$cuts_a) + 1L]
  ra <- t[!is.na(band_a) & is.na(band_b)]
  expect_equal(ra$events_b, as.numeric(tapply(d$y2, factor(ia, sort(st$codes_a)), sum)))
  expect_setequal(unique(cx$overlap_rates$outcome), c("y", "y2"))
  # one outcome alone keeps its suffix; no outcome gives counts only
  expect_true("rate_b" %in% names(scr_score_cross(d, "a", "b", y_b = "y2")$table))
  c0 <- scr_score_cross(d, "a", "b")
  expect_null(c0$overlap_rates)
  expect_identical(names(c0$table), c("band_a", "label_a", "band_b", "label_b", "n", "pct"))
  expect_output(print(cx), "Event rate of \"y2\"")
  expect_output(print(cx), "swap-in")
})

test_that("cuts from a study object carry its numbers and labels; numeric cuts are used as they are", {
  d <- cross_df()
  tr <- scr_tiers(d, score = "a", n_tiers = 3)
  cx <- scr_score_cross(d, "a", "b", y = "y", cuts_a = tr, cuts_b = c(430, 470))
  st <- cx$settings
  expect_equal(st$cuts_a, tr$cuts); expect_identical(st$labels_a, tr$code_labels)
  expect_identical(st$bands_a, "tiers")
  expect_setequal(cx$table$label_a[!is.na(cx$table$band_a)], tr$labels)
  # the tier of every row as scr_apply() assigns it
  tier <- scr_apply(tr, d, score = "a")$tier
  ra <- cx$table[!is.na(band_a) & is.na(band_b)]
  expect_equal(ra$n, as.numeric(table(factor(tier, sort(tr$codes)))))
  expect_equal(st$cuts_b, c(430, 470))
  expect_equal(cx$table[is.na(band_a) & !is.na(band_b)]$n,
               as.numeric(table(factor(st$codes_b[findInterval(d$b, c(430, 470)) + 1L], 1:3))))
  expect_error(scr_score_cross(d, "a", "b", cuts_a = "x"), "cuts_a")
})

test_that("a study given as cuts sets the objective and the direction of its score", {
  d <- cross_df()
  d$fa <- 1000 - d$a   # the same score on a fraud-like scale: a higher score means more risk
  tf <- scr_tiers(d, score = "fa", direction = "higher_is_riskier", n_tiers = 3)
  cx <- scr_score_cross(d, "fa", "b", y = "y", cuts_a = tf, depths = c(0.05, 0.2))
  expect_identical(cx$settings$direction_a, "higher_is_riskier")
  expect_identical(cx$settings$objective_a, "risk")
  # the overlap selects the high scores of A, its event-rich end
  ov <- cx$overlap
  for (k in seq_len(nrow(ov))) {
    A <- d$fa >= ov$cut_a[k]
    expect_equal(ov$n_a[k], sum(A))
    expect_gt(mean(d$y[A]), mean(d$y))
    expect_equal(cx$overlap_rates[depth == ov$depth[k] & set == "A", rate], mean(d$y[A]))
  }
  # the mirror image of the same score read as higher_is_safer: the same rows selected
  cm <- scr_score_cross(d, "a", "b", y = "y", depths = c(0.05, 0.2))
  expect_equal(ov$n_a, cm$overlap$n_a); expect_equal(ov$n_both, cm$overlap$n_both)
  expect_equal(cx$association$oriented, cm$association$oriented)
  # an explicit argument that agrees is accepted; one that disagrees is an error
  expect_identical(scr_score_cross(d, "fa", "b", cuts_a = tf, direction_a = "higher_is_riskier",
                                   objective_a = "risk")$settings$direction_a, "higher_is_riskier")
  expect_error(scr_score_cross(d, "fa", "b", cuts_a = tf, direction_a = "higher_is_safer"),
               "`direction_a` is \"higher_is_safer\" but the study given as `cuts_a` was fitted with \"higher_is_riskier\"")
  expect_error(scr_score_cross(d, "fa", "b", cuts_a = tf, objective_a = "propensity"), "`objective_a` is \"propensity\"")
  # a propensity study as cuts of B: objective and direction of B come from it, not from A
  tp <- scr_tiers(d, score = "b", y = "y2", objective = "propensity", n_tiers = 3)
  cb <- scr_score_cross(d, "a", "b", y_a = "y", y_b = "y2", cuts_b = tp)
  expect_identical(cb$settings$objective_b, "propensity")
  expect_identical(cb$settings$direction_b, "higher_is_riskier")
  expect_equal(cb$overlap$n_b[1], sum(d$b >= cb$overlap$cut_b[1]))
  expect_error(scr_score_cross(d, "a", "b", cuts_b = tp, objective_b = "risk"), "`objective_b` is \"risk\"")
  expect_error(scr_score_cross(d, "a", "b", cuts_b = tp, direction_b = "higher_is_safer"), "`direction_b`")
})

test_that("weights act on the counts and the Kish size; the association stays unweighted", {
  d <- cross_df()
  d$w[1:40] <- 0
  cx <- scr_score_cross(d, "a", "b", y = "y", weight = "w", n_bands = 3, depths = 0.1)
  st <- cx$settings
  k <- d$w > 0
  ia <- st$codes_a[findInterval(d$a, st$cuts_a) + 1L]
  ra <- cx$table[!is.na(band_a) & is.na(band_b)]
  for (b in ra$band_a) {
    i <- k & ia == b
    ny <- sum(d$w[i]); e <- sum(d$w[i] * d$y[i]); neff <- ny^2 / sum(d$w[i]^2)
    j <- ra$band_a == b
    expect_equal(ra$n[j], ny); expect_equal(ra$events[j], e)
    expect_equal(ra$rate_lo[j], stats::qbeta(0.025, e / ny * neff + 0.5, neff - e / ny * neff + 0.5))
  }
  expect_equal(st$n_dropped, 40L)
  expect_equal(cx$association$estimate[2], stats::cor(d$a[k], d$b[k], method = "kendall"), tolerance = 1e-12)
})

test_that("edge cases: single class, all ties, missing scores, bad arguments", {
  d <- cross_df(500)
  one <- transform(d, y = 0)
  c1 <- scr_score_cross(one, "a", "b", y = "y", n_bands = 3)
  expect_true(all(c1$table$rate[c1$table$n > 0] == 0))
  expect_true(all(is.na(c1$table$lift)))
  expect_false(any(is.nan(c1$table$lift)))   # no event, no lift: missing, not NaN
  # all ties on A: one band, no depth boundary (nothing selected below depth 1), no association
  tie <- transform(d, a = 7)
  ct <- scr_score_cross(tie, "a", "b", y = "y", depths = c(0.1, 1))
  expect_identical(ct$settings$cuts_a, numeric())
  expect_true(is.na(ct$overlap$cut_a[1])); expect_equal(ct$overlap$n_a, c(0, 500))
  expect_equal(ct$overlap$n_both[1], 0); expect_equal(ct$overlap$jaccard[1], 0)
  expect_true(all(is.na(ct$association$estimate)))
  expect_true(is.na(ct$overlap_rates$rate[ct$overlap_rates$depth == 0.1 & ct$overlap_rates$set == "A"]))
  # rows with a missing or infinite score are left out
  dm <- d; dm$a[1:5] <- NA; dm$b[6] <- Inf
  cm <- scr_score_cross(dm, "a", "b", y = "y")
  expect_equal(cm$settings$n_dropped, 6L); expect_equal(cm$settings$n, 494)
  expect_output(print(cm), "6 left out")
  expect_error(scr_score_cross(d, "a", "zz"), "not in `x`")
  expect_error(scr_score_cross(d, "a", "b", y = "y", y_a = "y"), "not both")
  expect_error(scr_score_cross(d, "a", "b", depths = 0), "depths")
  expect_error(scr_score_cross(d, "a", "b", objective_a = "fraud"), "objective_a")
  expect_error(scr_score_cross(d, "a", "b", direction_b = "up"), "direction_b")
  expect_error(scr_score_cross(d, "a", "b", foo = 1), "unused argument")
  expect_error(scr_score_cross(transform(d, a = NA_real_), "a", "b"), "no row has both scores")
})

test_that("scr_export writes the cross workbook", {
  skip_if_not_installed("openxlsx")
  old <- scr_verbose(FALSE); on.exit(scr_verbose(old), add = TRUE)
  cx <- scr_score_cross(cross_df(600), "a", "b", y = "y", n_bands = 3)
  out <- file.path(tempdir(), "scr-cross-export")
  unlink(out, recursive = TRUE)
  ex <- scr_export(cx, out, stamp = FALSE)
  expect_identical(openxlsx::getSheetNames(ex$files$xlsx),
                   c("Cross_Table", "Overlap", "Overlap_Rates", "Association", "Settings"))
  expect_equal(nrow(openxlsx::read.xlsx(ex$files$xlsx, sheet = "Cross_Table")), nrow(cx$table))
})
