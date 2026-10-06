test_that("scr_export writes the selection deliverables and the four scorecard workbooks", {
  skip_if_not_installed("openxlsx")
  out <- file.path(tempdir(), "scr-export-test")
  unlink(out, recursive = TRUE)
  res <- scr_export(res_demo(), out, stamp = FALSE)
  expect_true(all(file.exists(unlist(res$files))))
  expect_setequal(names(res$files), c("xlsx", "sql", "md"))
  sheets <- openxlsx::getSheetNames(res$files$xlsx)
  expect_true(all(c("01_Funnel", "02_Gains", "07_Consensus", "08_Ledger") %in% sheets))
  fun <- openxlsx::read.xlsx(res$files$xlsx, sheet = "01_Funnel")
  expect_equal(nrow(fun), nrow(res$funnel))

  sc <- scr_export(sc_demo(), out, stamp = FALSE)
  expect_setequal(names(sc$files), c("scorecard", "validation", "strategy", "sql_score", "sql_woe"))
  expect_true(all(file.exists(unlist(sc$files))))
  expect_true(all(c("Score_Summary", "Final_Scorecard", "Model_Card", "Alignment") %in%
                    openxlsx::getSheetNames(sc$files$scorecard)))
  expect_true(all(c("Score_Gains_Frozen", "Discrimination_CI", "Stability_PSI_Timeline",
                    "Stability_CSI_Timeline", "Calibration", "Performance_By_Vintage",
                    "Rank_Order_Diagnostics") %in% openxlsx::getSheetNames(sc$files$validation)))
  expect_true(all(c("Population_Scope", "Cutoff_Sweep", "Strategy_Bands", "Reject_Sensitivity",
                    "Monitoring_Plan") %in% openxlsx::getSheetNames(sc$files$strategy)))
  ss <- openxlsx::read.xlsx(sc$files$scorecard, sheet = "Score_Summary")
  expect_true("odds_orientation" %in% ss$item)
  expect_true(any(grepl("score_points", readLines(sc$files$sql_score))))
})

test_that("the hardened writer sanitises formula injection and never fabricates a row", {
  skip_if_not_installed("openxlsx")
  f <- file.path(tempdir(), "scr-hardened.xlsx")
  sheets <- list(A = data.frame(x = c("=cmd()", "+1", "-1", "@x", "plain"), stringsAsFactors = FALSE),
                 B = data.frame())
  .scr_write_xlsx(sheets, f)
  a <- openxlsx::read.xlsx(f, sheet = "A")
  expect_true(all(startsWith(a$x[1:4], "'")))
  expect_equal(a$x[5], "plain")
  b <- openxlsx::read.xlsx(f, sheet = "B")
  expect_equal(b$availability, "not_available")
  expect_false(file.exists(file.path(dirname(f), paste0(".", basename(f), ".tmp"))))
})

test_that("stamp = TRUE writes into a timestamped subdirectory", {
  skip_if_not_installed("openxlsx")
  out <- file.path(tempdir(), "scr-stamp")
  res <- scr_export(res_demo(), out, stamp = TRUE)
  expect_match(dirname(res$files$xlsx), "[0-9]{8}_[0-9]{6}$")
})

test_that("the object's config$verbose governs the messages of scr_sql() and scr_export()", {
  old <- scr_verbose(TRUE); on.exit(scr_verbose(old), add = TRUE)
  res <- res_demo()
  res$config$verbose <- FALSE
  expect_no_message(scr_sql(res, file = withr::local_tempfile(fileext = ".sql")))
  expect_true(scr_verbose())
  res$config$verbose <- TRUE
  expect_message(scr_sql(res, file = withr::local_tempfile(fileext = ".sql")), "SQL written")
})
