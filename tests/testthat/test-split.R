test_that("the out-of-time split uses whole periods and drops the date column", {
  sp <- scr_split(scr_demo, "default", date_col = "ref_date", drop = "id")
  expect_s3_class(sp, "scr_split")
  expect_equal(sp$method, "out-of-time")
  d_tr <- unique(scr_demo$ref_date[sp$train_idx]); d_ho <- unique(scr_demo$ref_date[sp$holdout_idx])
  expect_length(intersect(d_tr, d_ho), 0L)
  expect_gte(length(sp$holdout_idx) / nrow(scr_demo), 0.3)
  expect_true("ref_date" %in% sp$cols$dropped)
  expect_false("ref_date" %in% sp$cols$features)
  expect_output(print(sp), "out-of-time")
})

test_that("without a date column the split is stratified random and reproducible", {
  a <- scr_split(scr_demo, "default", seed = 1, drop = "id")
  b <- scr_split(scr_demo, "default", seed = 1, drop = "id")
  expect_equal(a$method, "stratified random")
  expect_identical(a$train_idx, b$train_idx)
  expect_equal(mean(scr_demo$default[a$train_idx]), mean(scr_demo$default[a$holdout_idx]), tolerance = 0.02)
  expect_message(scr_split(scr_demo, "default", date_col = "nope"), "does not exist")
})

test_that("event_level rewrites what is modelled and text targets are accepted", {
  inv <- scr_split(scr_demo, "default", event_level = 0)
  expect_true(inv$cols$event$inverted)
  expect_equal(inv$data$default, 1L - scr_demo$default)
  d <- head(scr_demo, 500)
  d$flag <- ifelse(d$default == 1, "BAD", "GOOD")
  tx <- scr_split(d, "flag", drop = c("default", "churn"))
  expect_equal(tx$cols$event$label, "GOOD")   # alphabetical default, reported
  expect_equal(scr_split(d, "flag", event_level = "BAD")$cols$event$label, "BAD")
  expect_error(scr_split(d, "flag", event_level = "UGLY"), "does not exist")
  d$lg <- as.logical(d$default)
  expect_equal(scr_split(d, "lg")$data$lg, d$default)
})

test_that("input errors are clear", {
  expect_error(scr_split(1:10, "y"), "data.frame")
  expect_error(scr_split(scr_demo, c("a", "b")), "single column")
  expect_error(scr_split(scr_demo, "missing_col"), "missing from the table")
  d <- data.frame(y = c(0, 1, 2), x = 1:3)
  expect_error(scr_split(d, "y"), "0/1")
  d$y <- c(0, 1, NA)
  expect_error(scr_split(d, "y"), "NA")
  expect_error(scr_split(scr_demo, "default", event_level = 3), "0 or 1")
})

test_that("copy = TRUE leaves the user's data untouched", {
  d <- head(scr_demo, 300)
  before <- vapply(d, function(x) class(x)[1], character(1))
  invisible(scr_split(d, "default"))
  expect_identical(vapply(d, function(x) class(x)[1], character(1)), before)
})
