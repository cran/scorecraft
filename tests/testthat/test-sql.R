test_that("the WOE SQL covers exactly the shortlist, in two blocks", {
  res <- res_demo()
  sql <- scr_sql(res)
  txt <- paste(sql, collapse = "\n")
  expect_match(txt, "WITH base_scr AS")
  expect_match(txt, "FROM your_table")
  for (f in scr_selected(res)) expect_match(txt, paste0(f, "_woe"), fixed = TRUE)
  for (f in setdiff(names(scr_demo), c(scr_selected(res), "default"))) {
    expect_false(grepl(paste0("AS ", f, "_woe"), txt, fixed = TRUE))
  }
  expect_match(paste(scr_sql(res, table = "prd.t", dialect = "databricks"), collapse = "\n"), "FROM prd.t")
  tmp <- tempfile(fileext = ".sql")
  expect_invisible(scr_sql(res, file = tmp))
  expect_true(file.exists(tmp))
})

test_that("R and SQL agree on WOE, bins, score and points (duckdb)", {
  skip_if_not_installed("duckdb"); skip_if_not_installed("DBI")
  res <- res_demo(); sc <- sc_demo()
  con <- DBI::dbConnect(duckdb::duckdb(), config = list(threads = "1"))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  DBI::dbWriteTable(con, "scr_demo", scr_demo)

  got <- DBI::dbGetQuery(con, paste(scr_sql(res, table = "scr_demo", dialect = "duckdb"), collapse = "\n"))
  exp <- scr_apply(res, scr_demo, what = "both")
  for (f in scr_selected(res)) {
    expect_equal(got[[paste0(f, "_woe")]], exp[[paste0(f, "_woe")]], tolerance = 1e-12)
    expect_identical(got[[paste0(f, "_bin")]], exp[[paste0(f, "_bin")]])
  }
  gs <- DBI::dbGetQuery(con, paste(scr_sql(sc, table = "scr_demo", dialect = "duckdb"), collapse = "\n"))
  es <- scr_apply(sc, scr_demo, what = "points")
  expect_equal(gs$score, es$score, tolerance = 1e-10)
  expect_equal(gs$score_points, es$score_points)
  for (f in sc$features) expect_equal(gs[[paste0(f, "_points")]], es[[paste0(f, "_points")]])
})

test_that("R and SQL agree on score and points (SQLite)", {
  skip_if_not_installed("RSQLite"); skip_if_not_installed("DBI")
  sc <- sc_demo()
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con))
  d <- scr_demo; d$ref_date <- as.character(d$ref_date)
  DBI::dbWriteTable(con, "scr_demo", d)
  gs <- DBI::dbGetQuery(con, paste(scr_sql(sc, table = "scr_demo", dialect = "sqlite"), collapse = "\n"))
  es <- scr_apply(sc, scr_demo, what = "points")
  expect_equal(gs$score, es$score, tolerance = 1e-10)
  expect_equal(gs$score_points, es$score_points)
})

test_that("keep columns and the WOE-only scorecard SQL work", {
  res <- res_demo()
  cfg <- res$config; cfg$sql_keep_columns <- "id"
  res2 <- res; res2$config <- cfg
  txt <- paste(scr_sql(res2, table = "t"), collapse = "\n")
  expect_match(txt, "\n    id,")
  sc <- sc_demo()
  w <- scr_sql(sc, what = "woe")
  expect_false(any(grepl("score_points", w)))
  expect_true(any(grepl("score_points", scr_sql(sc))))
})

test_that("the deployment SQL (what = \"all\") carries keys, bin labels, WOE, points and the score, equal to R", {
  skip_if_not_installed("DBI")
  sc <- sc_demo(); res <- res_demo()
  sql <- scr_sql(sc, table = "t", dialect = "duckdb", what = "all", keep_columns = c("id", "ref_date"))
  expect_true(any(grepl("bin label", sql)))
  feats <- sc$features
  for (f in feats) expect_true(any(grepl(sprintf("%s_bin,", f), sql, fixed = TRUE)))
  expect_error(scr_sql(sc, what = "all", keep_columns = c("id", NA)), "keep_columns")
  # the cached score SQL is untouched by keep_columns / what = "all"
  expect_identical(scr_sql(sc), sc$sql)
  exp_bin <- scr_apply(res, scr_demo, features = feats, what = "bin")
  exp_all <- scr_apply(sc, scr_demo, what = "all")
  run <- function(con) {
    DBI::dbWriteTable(con, "t", as.data.frame(scr_demo))
    got <- DBI::dbGetQuery(con, paste(sql, collapse = "\n"))
    expect_equal(got$id, scr_demo$id)
    for (f in feats) {
      expect_identical(got[[paste0(f, "_bin")]], exp_bin[[paste0(f, "_bin")]])
      expect_equal(got[[paste0(f, "_woe")]], exp_all[[paste0(f, "_woe")]], tolerance = 1e-12)
      expect_equal(got[[paste0(f, "_points")]], exp_all[[paste0(f, "_points")]])
    }
    expect_equal(got$score, exp_all$score, tolerance = 1e-9)
    expect_equal(got$score_points, exp_all$score_points)
  }
  if (requireNamespace("duckdb", quietly = TRUE)) {
    con <- DBI::dbConnect(duckdb::duckdb()); run(con); DBI::dbDisconnect(con, shutdown = TRUE)
  }
  if (requireNamespace("RSQLite", quietly = TRUE)) {
    sql <<- scr_sql(sc, table = "t", dialect = "sqlite", what = "all", keep_columns = c("id", "ref_date"))
    con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:"); run(con); DBI::dbDisconnect(con)
  }
})
