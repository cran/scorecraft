# Regression tests for the core-stage audit (split, triage, bin, model,
# config, utils). One block per fixed defect.

test_that("a seed given to a stage is applied locally: the user's stream is restored", {
  set.seed(99); ref <- stats::runif(3)
  set.seed(99); a <- stats::runif(1)
  invisible(scr_split(scr_demo, "default", seed = 1, drop = "id"))
  expect_identical(c(a, stats::runif(2)), ref)

  set.seed(99); a <- stats::runif(1)
  i <- subsample_stratified(rep(0:1, 500), 100, seed = 3)
  expect_identical(c(a, stats::runif(2)), ref)
  expect_identical(i, subsample_stratified(rep(0:1, 500), 100, seed = 3))

  # no .Random.seed before the call -> none after it
  if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv())
  invisible(subsample_stratified(rep(0:1, 500), 100, seed = 3))
  expect_false(exists(".Random.seed", envir = globalenv()))
  set.seed(NULL)
})

test_that("a fractional target is refused instead of being truncated to 0/1", {
  d <- data.frame(y = c(0, 0.5, 1, 1), x = 1:4)
  expect_error(scr_split(d, "y"), "0/1")
  d$y <- c(0, 1, 1, 0.99)
  expect_error(scr_split(d, "y"), "0/1")
})

test_that("a text date column still gives the out-of-time split", {
  d <- data.table::as.data.table(scr_demo)
  ref <- scr_split(d, "default", date_col = "ref_date", drop = "id")
  d[, ref_date := as.character(ref_date)]
  tx <- scr_split(d, "default", date_col = "ref_date", drop = "id")
  expect_equal(tx$method, "out-of-time")
  expect_identical(tx$holdout_idx, ref$holdout_idx)
  d[, ref_date := format(as.Date(ref_date), "%Y%m")]
  expect_identical(scr_split(d, "default", date_col = "ref_date", drop = "id")$holdout_idx, ref$holdout_idx)
  d[, ref_date := format(as.Date(scr_demo$ref_date), "%Y-%m")]
  expect_identical(scr_split(d, "default", date_col = "ref_date", drop = "id")$holdout_idx, ref$holdout_idx)
})

test_that("rows with a missing date are reported and belong to neither side", {
  d <- data.table::as.data.table(scr_demo)
  d[1:10, ref_date := NA]
  old <- scr_verbose(TRUE); on.exit(scr_verbose(old))
  expect_message(sp <- scr_split(d, "default", date_col = "ref_date", drop = "id"), "missing 'ref_date'")
  expect_false(any(1:10 %in% c(sp$train_idx, sp$holdout_idx)))
  expect_length(intersect(sp$train_idx, sp$holdout_idx), 0L)
})

test_that("integer64 candidates are converted by value, not by bit pattern", {
  skip_if_not_installed("bit64")
  d <- head(scr_demo[, c("default", "vl_score_01")], 200)
  d$big <- bit64::as.integer64(seq_len(200))
  sp <- scr_split(d, "default")
  expect_true(is.double(sp$data$big) && !is.object(sp$data$big))
  expect_equal(sp$data$big, as.double(seq_len(200)))
})

test_that("sentinels are numeric-only: a text level equal to a sentinel is not special", {
  set.seed(5)
  y <- rep(0:1, 100)
  x <- ifelse(y == 1 & seq_along(y) %% 3 == 0, "-999", sample(c("a", "b"), 200, TRUE))
  x[seq(1, 200, by = 10)] <- NA
  set.seed(NULL)
  r <- .triage_one("x", x, y, FALSE, -999, scr_config(verbose = FALSE), length(y))
  expect_equal(r$row$woe_special, woe_subpop(is.na(x), y))
  expect_equal(r$row$n_special, 0L)
})

test_that("triage materialises more survivors than data.table over-allocates (1024)", {
  set.seed(11)
  n <- 120L; p <- 1100L
  y <- rep(0:1, length.out = n)
  d <- data.table::as.data.table(matrix(stats::rnorm(n * p), n, p))
  d[, target := y]
  set.seed(NULL)
  sp <- scr_split(d, "target", seed = 1)
  cfg <- scr_config(verbose = FALSE, nthread = 1, min_iv_quick = 0, check_duplicates = FALSE)
  tr <- scr_triage(sp, cfg)
  expect_gt(length(tr$keep), 1024L)
  expect_true(all(tr$keep %in% names(tr$clean)))
})

test_that("configuration keys of stages 0-7 are validated when the config is built", {
  expect_error(scr_config(cv_folds = 2), "cv_folds")
  expect_error(scr_config(corr_method = "pearsn"), "corr_method")
  expect_error(scr_config(verbose = "yes"), "verbose")
  expect_error(scr_config(special_values = NA), "special_values")
  expect_error(scr_config(rf_importance = "none"), "rf_importance")
  expect_error(scr_config(seed = 1.5), "seed")
  expect_error(scr_config(min_bins = 5, max_bins = 3), "max_bins")
  expect_error(scr_config(iv_min = 0.5, iv_max = 0.4), "iv_max")
  expect_error(scr_config(require_monotonic = "yes"), "require_monotonic")
  expect_error(scr_config(xgb_subsample = 0), "xgb_subsample")
  expect_error(scr_config(model_top_k = 0), "model_top_k")
  expect_error(scr_config(nthread = NA), "nthread")
  expect_s3_class(scr_config(special_values = numeric(), iv_max = Inf, model_max_rows = Inf), "scr_config")
})

test_that(".sql_str doubles the backslash only where it is an escape, and maps NA to NULL", {
  expect_identical(.sql_str("a\\b'c", dialect = "ansi"), "'a\\b''c'")
  expect_identical(.sql_str("a\\b'c", dialect = "databricks"), "'a\\\\b''c'")
  expect_identical(.sql_str("a\\b"), "'a\\\\b'")   # NULL dialect: historical behaviour
  expect_identical(.sql_str(c("x", NA)), c("'x'", "NULL"))
})

test_that("derived flags excluded by policy never prune a real column", {
  res <- res_demo()
  expect_false(any(res$derived_excluded %in% res$prune$dropped$feature))
  expect_false(any(res$prune$dropped$correlated_with %in% res$derived_excluded))
  expect_false(any(res$derived_excluded %in% res$prune$keep))
})

test_that("the xgboost adapter works on the installed API and leaves the RNG alone", {
  set.seed(21)
  n <- 600L
  x <- matrix(stats::rnorm(n * 3), n, 3, dimnames = list(NULL, c("a", "b", "c")))
  y <- as.integer(stats::runif(n) < stats::plogis(x[, 1]))
  cfg <- scr_config(verbose = FALSE, nthread = 1, xgb_rounds = 30L, xgb_early_stopping = 5L)
  set.seed(99); ref <- stats::runif(2); set.seed(99); a <- stats::runif(1)
  fx <- .fit_xgboost(x[1:400, ], y[1:400], x[401:600, ], y[401:600], cfg)
  expect_identical(c(a, stats::runif(1)), ref)
  expect_setequal(fx$importance$feature, c("a", "b", "c"))
  expect_length(fx$score, 200L)
  expect_true(all(fx$score > 0 & fx$score < 1))
  n_tr <- as.integer(sub(" trees", "", fx$note))
  expect_true(n_tr >= 1L && n_tr <= 30L)
  fx2 <- .fit_xgboost(x[1:400, ], y[1:400], x[401:600, ], y[401:600], cfg)
  expect_identical(fx$score, fx2$score)
  set.seed(NULL)
})

test_that("the lightgbm adapter maps min_child_weight to the hessian floor and leaves the RNG alone", {
  skip_if_not_installed("lightgbm")
  set.seed(22)
  n <- 600L
  x <- matrix(stats::rnorm(n * 3), n, 3, dimnames = list(NULL, c("a", "b", "c")))
  y <- as.integer(stats::runif(n) < stats::plogis(x[, 1]))
  cfg <- scr_config(verbose = FALSE, nthread = 1, xgb_rounds = 30L, xgb_early_stopping = 5L,
                    xgb_min_child_weight = 1)
  set.seed(99); ref <- stats::runif(2); set.seed(99); a <- stats::runif(1)
  fl <- .fit_lightgbm(x[1:400, ], y[1:400], x[401:600, ], y[401:600], cfg)
  expect_identical(c(a, stats::runif(1)), ref)
  expect_equal(fl$model$params$min_sum_hessian_in_leaf, 1)
  expect_null(fl$model$params$min_data_in_leaf)
  expect_length(fl$score, 200L)
  set.seed(NULL)
})
