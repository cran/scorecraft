# Regression tests of the production triage: metrics, PSI, alignment,
# scorecard points, apply/SQL equivalence, cut-off, monitor, export, db and
# the coarse classing spec round trip.

# -- metrics ------------------------------------------------------------------ #

test_that("scr_metrics() refuses a factor or a non-0/1 outcome instead of miscounting it", {
  s <- c(0.1, 0.4, 0.35, 0.8)
  expect_error(scr_metrics(s, factor(c(0, 0, 1, 1)), ci = FALSE), "0/1 numeric or logical")
  expect_error(scr_metrics(s, c("a", "a", "b", "b"), ci = FALSE), "0/1 numeric or logical")
  expect_error(scr_metrics(s, c(0, 2, 1, 1), ci = FALSE), "must be 0/1")
  expect_equal(scr_metrics(s, c(FALSE, FALSE, TRUE, TRUE), ci = FALSE)$auc,
               scr_metrics(s, c(0, 0, 1, 1), ci = FALSE)$auc)
  expect_error(scr_align(s, c(0, 3, 1, 1)), "must be 0/1")
})

test_that("the rank-based AUC/KS equals the pairwise definition, ties counted as 1/2", {
  set.seed(11)
  y <- rbinom(300, 1, 0.3); s <- round(rnorm(300) + y, 1)
  d <- outer(s[y == 1], s[y == 0], "-")
  auc_pairs <- mean((d > 0) + 0.5 * (d == 0))
  m <- scr_metrics(s, y, ci = FALSE)
  expect_equal(m$auc, auc_pairs, tolerance = 1e-12)
  ks <- max(abs(ecdf(s[y == 1])(sort(unique(s))) - ecdf(s[y == 0])(sort(unique(s)))))
  expect_equal(m$ks, ks, tolerance = 1e-12)
  # the bootstrap is reproducible and brackets the point estimate
  b1 <- scr_metrics(s, y, n_boot = 40, seed = 3); b2 <- scr_metrics(s, y, n_boot = 40, seed = 3)
  expect_identical(b1, b2)
  expect_true(b1$auc_lo <= b1$auc && b1$auc <= b1$auc_hi)
})

test_that("scr_iv() does not tabulate up to a large integer code", {
  set.seed(2)
  y <- rbinom(500, 1, 0.3)
  g_int <- sample(c(20230101L, 20230201L, 20230301L), 500, replace = TRUE)
  expect_equal(scr_iv(g_int, y), scr_iv(as.character(g_int), y))
})

# -- PSI ---------------------------------------------------------------------- #

test_that("scr_psi(): a band empty in both samples changes neither the PSI nor the degrees of freedom", {
  set.seed(4)
  a <- sample(c("x", "y", "z"), 800, replace = TRUE, prob = c(0.5, 0.3, 0.2))
  b <- sample(c("x", "y", "z"), 700, replace = TRUE, prob = c(0.45, 0.35, 0.2))
  p0 <- scr_psi(a, b, levels = c("x", "y", "z"))
  p1 <- scr_psi(a, b, levels = c("x", "y", "z", "never"))
  expect_equal(p1$psi, p0$psi)
  expect_equal(p1$critical, p0$critical)
  expect_equal(p1$n_bins, 4L)
  expect_equal(p0$critical, (1 / 800 + 1 / 700) * stats::qchisq(0.95, 2))
  # the Yurdakul-Naranjo value quoted in the documentation
  base <- rep(1:10, each = 100); cmp <- rep(1:10, each = 100)
  expect_equal(round(scr_psi(base, cmp, breaks = c(-Inf, 1:9 + 0.5, Inf))$critical, 4),
               round(2 / 1000 * stats::qchisq(0.95, 9), 4))
  expect_equal(round(2 / 1000 * stats::qchisq(0.95, 9), 3), 0.034)
})

# -- alignment ---------------------------------------------------------------- #

test_that("scr_align() validates the weights", {
  set.seed(3)
  y <- rbinom(500, 1, 0.2); raw <- stats::qlogis(0.2) + y + rnorm(500)
  expect_error(scr_align(raw, y, weights = rep(1, 10)), "length of `raw`")
  expect_error(scr_align(raw, y, weights = c(-1, rep(1, 499))), "non-negative")
  expect_equal(scr_align(raw, y, weights = rep(1, 500))$b, scr_align(raw, y)$b)
})

# -- scorecard, apply and SQL ------------------------------------------------- #

test_that("a row in no fitted bin takes the points of WOE 0 in R and in SQL (distributed points)", {
  res <- res_demo()
  cats <- intersect(names(res$fit$results)[vapply(res$fit$results, function(r) identical(r$type, "categorical"), logical(1))],
                    scr_selected(res))
  skip_if(!length(cats), "no categorical variable in the demo shortlist")
  f <- cats[1]
  sc <- scr_scorecard(res, points_style = "distributed")
  expect_true(f %in% sc$features)
  new <- head(scr_demo, 20)
  new[[f]] <- as.character(new[[f]]); new[[f]][1:3] <- "LEVEL_NEVER_SEEN"
  ap <- scr_apply(sc, new, what = "all")
  unb <- .sc_unbinned_points(sc)
  k <- length(sc$features)
  base_raw <- sc$alignment$a + sc$alignment$b * unname(sc$coef["(Intercept)"])
  expect_equal(unb, round(base_raw / k))
  expect_true(all(ap[[paste0(f, "_points")]][1:3] == unb))
  expect_equal(ap[[paste0(f, "_woe")]][1:3], rep(0, 3))
  # whole points stay within rounding of the exact score on every row
  expect_true(all(abs(ap$score - ap$score_points) <= 0.5 * k + 1e-9))
  skip_if_not_installed("RSQLite"); skip_if_not_installed("DBI")
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con))
  d <- new; d$ref_date <- as.character(d$ref_date)
  DBI::dbWriteTable(con, "t", d)
  got <- DBI::dbGetQuery(con, paste(scr_sql(sc, table = "t", dialect = "sqlite"), collapse = "\n"))
  expect_equal(got$score, ap$score, tolerance = 1e-9)
  expect_equal(got$score_points, ap$score_points)
  expect_equal(got[[paste0(f, "_points")]], ap[[paste0(f, "_points")]])
})

test_that("the scorecard's own SQL quotes identifiers exactly as obwoe_sql() does", {
  ob_ident <- get0(".ob_sql_ident", envir = asNamespace("OptimalBinningWoE"), inherits = FALSE)
  ob_dial  <- get0(".ob_sql_dialect", envir = asNamespace("OptimalBinningWoE"), inherits = FALSE)
  skip_if(is.null(ob_ident) || is.null(ob_dial), "engine internals not available")
  nm <- c("income", "order", "User", "my var", "a\"b", "x.y", "_ok", "9lives", "limit_1", "Select")
  for (dl in c("ansi", "postgres", "duckdb", "sqlite", "databricks", "mysql", "sqlserver", "oracle", "snowflake")) {
    expect_identical(.sql_q(nm, dl), ob_ident(nm, ob_dial(dl), "auto"), info = dl)
  }
  cfg <- res_demo()$config; cfg$sql_keep_columns <- c("id", "order")
  sc <- sc_demo(); sc$config <- cfg
  sql <- paste(build_sql_score(sc), collapse = "\n")
  expect_match(sql, "\n    \"order\",", fixed = TRUE)
  expect_false(grepl("\n    order,", sql, fixed = TRUE))
})

test_that("a line break in a name cannot escape an SQL comment", {
  sc <- sc_demo(); sc$target <- "default\nDROP TABLE t"
  sql <- build_sql_score(sc)
  expect_false(any(grepl("^DROP TABLE", sql)))
})

test_that("scr_reasons() (vectorised) returns what the row-by-row ranking returned", {
  sc <- sc_demo()
  new <- head(scr_demo, 60)
  got <- scr_reasons(sc, new, k = 3)
  pts <- scr_apply(sc, new, what = "points")
  feats <- sc$features
  ref <- vapply(feats, function(f) { p <- sc$points[variable == f]; sum(p$points * p$pct_train) }, numeric(1))
  sgn <- if (identical(sc$direction, "higher_is_safer")) 1 else -1
  M <- as.matrix(pts[, paste0(feats, "_points"), with = FALSE])
  short <- sweep(-M, 2L, -ref) * sgn
  for (i in seq_len(nrow(M))) {
    o <- order(-short[i, ])[1:3]
    expect_identical(unname(unlist(got[i, c("reason_1", "reason_2", "reason_3")])), feats[o])
    expect_equal(unname(unlist(got[i, c("shortfall_1", "shortfall_2", "shortfall_3")])), unname(short[i, o]))
  }
})

test_that("the calibration table is in the numeric order of the bands", {
  sc <- sc_demo()
  lo <- as.numeric(sub("^[\\[(]([^,]+),.*$", "\\1", sc$calibration$table$band))
  expect_false(is.unsorted(lo))
})

# -- cut-off, strategy, reject ---------------------------------------------- #

test_that("the cut-off sweep (one sort per sample) equals a pass per cut", {
  sc <- sc_demo()
  ct <- scr_cutoff(sc, n_cuts = 12)
  for (nm in names(sc$samples)) {
    s <- sc$samples[[nm]]; tb <- ct$table[sample == nm]
    for (r in seq_len(nrow(tb))) {
      safe <- if (identical(sc$direction, "higher_is_safer")) s$score >= tb$cut[r] else s$score < tb$cut[r]
      expect_identical(tb$n_safe[r], sum(safe))
      expect_equal(tb$event_rate_safe[r], if (any(safe)) mean(s$y[safe]) else NA_real_)
      expect_equal(tb$events_avoided_pct[r], sum(s$y[!safe]) / max(1L, sum(s$y)))
    }
  }
})

test_that("scr_reject() validates its sample, accepted flags and multipliers", {
  sc <- sc_demo()
  expect_error(scr_reject(sc, sample = "nope"), "does not exist")
  acc <- rep(TRUE, 10); acc[2] <- NA
  expect_error(scr_reject(sc, population = head(scr_demo, 10), accepted = acc), "no NA")
  expect_error(scr_reject(sc, multipliers = c(2, -1)), "positive")
})

# -- monitor ------------------------------------------------------------------ #

test_that("scr_monitor() keeps undated rows as a period and matches scr_psi()", {
  sc <- sc_demo()
  d <- head(scr_demo, 600); d$ref_date[1:25] <- NA
  mo <- scr_monitor(sc, d, date_col = "ref_date", n_boot = 5)
  expect_equal(sum(mo$psi$n), nrow(d))
  expect_true(is.na(utils::tail(mo$psi$period, 1)))
  expect_equal(utils::tail(mo$psi$n, 1), 25L)
  sco <- scr_apply(sc, d)$score
  p1 <- mo$psi$period[1]
  i <- which(as.character(d$ref_date) == p1)
  r <- scr_psi(sc$samples$train$score, sco[i], breaks = sc$breaks, alpha = mo$thresholds$alpha,
               thresholds = mo$thresholds$psi)
  expect_equal(mo$psi$psi[1], r$psi)
  expect_equal(mo$psi$critical[1], r$critical)
  expect_error(scr_monitor(sc, d, alpha = 2), "alpha")
})

# -- export and db ------------------------------------------------------------ #

test_that("a target name never takes an export outside the given directory", {
  expect_identical(.file_tag("default"), "default")
  expect_identical(.file_tag("Bad.Flag"), "bad.flag")
  expect_false(grepl("/|\\.\\.", .file_tag("../../etc/passwd")))
  expect_identical(.file_tag(".."), "target")
})

test_that("scr_fetch() keeps a sampling fraction below 1e-6 (not printed as 0)", {
  skip_if_not_installed("RSQLite"); skip_if_not_installed("DBI")
  con <- scr_connect(driver = RSQLite::SQLite(), dbname = ":memory:")
  on.exit(DBI::dbDisconnect(con))
  DBI::dbWriteTable(con, "t", data.frame(a = 1:5))
  got <- scr_fetch(con, "t", sample_frac = 1e-7, sample_expr = "0.00000005", verbose = FALSE)
  expect_equal(nrow(got), 5L)
})

# -- classing spec round trip ------------------------------------------------- #

.spec_fixture <- function() {
  data.frame(
    target = "default",
    variable = c("v_num", "v_num", "v_num", "v_cat", "v_cat", "v_cat"),
    type = c("numeric", "numeric", "numeric", "categorical", "categorical", "categorical"),
    bin_id = c(1L, 2L, 3L, 1L, 2L, 3L),
    bin_label = c("a", "b", "c", "NA", "01%;%02", "-1%;%=x"),
    lower = c(NA, 10.5, 20, NA, NA, NA), upper = c(10.5, 20, NA, NA, NA, NA),
    categories = c(NA, NA, NA, "NA", "01%;%02", "-1%;%=x"),
    is_other = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
    reason = NA_character_, stringsAsFactors = FALSE)
}

test_that("a classing spec survives the CSV round trip: \"NA\" and \"01\" stay categories", {
  f <- tempfile(fileext = ".csv")
  utils::write.csv(.spec_fixture(), f, row.names = FALSE, na = "")
  sp <- scr_classing_read(f)
  expect_identical(sp$categories[sp$variable == "v_cat"], c("NA", "01%;%02", "-1%;%=x"))
  expect_equal(sp$upper[sp$variable == "v_num"], c(10.5, 20, NA))
})

test_that("a classing spec survives the xlsx round trip despite the formula guard", {
  skip_if_not_installed("openxlsx")
  f <- tempfile(fileext = ".xlsx")
  .scr_write_xlsx(list(Coarse_Classing = .spec_fixture()), f)
  sp <- scr_classing_read(f)
  expect_identical(sp$categories[sp$variable == "v_cat"], c("NA", "01%;%02", "-1%;%=x"))
})

test_that("scr_classing_read() refuses a bound on an open end", {
  d <- .spec_fixture(); d$upper[3] <- 99
  f <- tempfile(fileext = ".csv")
  utils::write.csv(d, f, row.names = FALSE, na = "")
  expect_error(scr_classing_read(f), "open ends")
})

test_that("the hardened xlsx writer accepts a legitimate all-NA row", {
  skip_if_not_installed("openxlsx")
  f <- tempfile(fileext = ".xlsx")
  d <- data.frame(a = c(1, NA, 3, NA), b = c("x", NA, "z", NA), stringsAsFactors = FALSE)
  expect_identical(.scr_write_xlsx(list(S = d, One = data.frame(z = NA)), f), f)
  expect_true(file.exists(f))
})

test_that("scr_metrics(seed = ) leaves the user's random stream untouched", {
  s <- stats::rnorm(200); y <- rep(0:1, 100)
  set.seed(9); before <- .Random.seed
  m <- scr_metrics(s, y, n_boot = 20, seed = 4)
  expect_identical(.Random.seed, before)
  expect_identical(m, scr_metrics(s, y, n_boot = 20, seed = 4))
})

test_that("the R pre-processing handles more features than the default column slots", {
  led <- data.table::data.table(kind = character(), output = character(), source = character(),
                                impute_value = numeric())
  nd <- as.data.frame(matrix(stats::rnorm(3 * 1100), 3))
  expect_identical(dim(.scr_preprocess(nd, led, names(nd), numeric())), c(3L, 1100L))
})
