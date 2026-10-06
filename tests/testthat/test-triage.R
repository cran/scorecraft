test_that("triage leaves no NA and decomposes sentinels into flags", {
  sp <- scr_split(scr_demo, "default", date_col = "ref_date", drop = c("id", "churn"))
  tr <- scr_triage(sp, cfg_test())
  expect_s3_class(tr, "scr_triage")
  expect_false(any(vapply(tr$clean, anyNA, logical(1))))
  expect_false(any(tr$clean$vl_hist_03 == -999))
  expect_true(any(grepl("__sp$", tr$derived)))
  expect_true(all(tr$derived %in% names(tr$clean)))
  expect_setequal(unique(tr$ledger$kind), c("num_impute", "num_flag", "cat_coalesce"))
  expect_output(print(tr), "survivors")
})

test_that("structural pathologies fail with the right reason", {
  sp <- scr_split(scr_demo, "default", date_col = "ref_date", drop = c("id", "churn"))
  p <- scr_triage(sp, cfg_test())$profile
  reason <- function(f) p[feature == f, triage_reason]
  expect_equal(reason("vl_constant"), "CONSTANT")
  expect_equal(reason("ds_constant"), "CONSTANT")
  expect_equal(reason("vl_near_const"), "NEAR_CONSTANT")
  expect_equal(reason("ds_high_card"), "HIGH_CARDINALITY")
  expect_match(reason("vl_duplicate"), "^DUPLICATE_OF:")
  expect_match(reason("vl_partial_03"), "TOO_MANY_MISSING")
  expect_false("vl_duplicate" %in% p[triage_status == "keep", feature])
})

test_that("the ledger is unconditional for survivors so hold-out and production agree", {
  sp <- scr_split(scr_demo, "default", date_col = "ref_date", drop = c("id", "churn"))
  tr <- scr_triage(sp, cfg_test())
  num_keep <- intersect(tr$keep, sp$cols$var_num)
  expect_true(all(num_keep %in% tr$ledger[kind == "num_impute", source]))
  cat_keep <- intersect(tr$keep, sp$cols$var_cat)
  expect_true(all(cat_keep %in% tr$ledger[kind == "cat_coalesce", source]))
})

test_that("triage validates its inputs", {
  expect_error(scr_triage(list(), cfg_test()), "scr_split")
  sp <- scr_split(head(scr_demo, 200), "default")
  expect_error(scr_triage(sp, list()), "scr_config")
})
