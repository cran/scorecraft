test_that("presets change exactly four keys", {
  p <- scr_presets()
  expect_equal(p$preset, c("aggressive", "moderate", "lazy"))
  a <- scr_config("aggressive"); m <- scr_config("moderate"); l <- scr_config("lazy")
  diff_keys <- function(x, y) names(x)[!mapply(identical, x, y[names(x)])]
  expect_setequal(diff_keys(a, m), c("preset", "target_max", "min_votes", "corr_cutoff", "iv_min"))
  expect_setequal(diff_keys(l, m), c("preset", "target_max", "min_votes", "corr_cutoff"))
})

test_that("overrides are validated, NULL keeps the preset, duplicates fail", {
  expect_error(scr_config(iv_maximum = 1), "unknown key")
  expect_equal(scr_config("aggressive", iv_min = NULL)$iv_min, 0.03)
  expect_equal(scr_config(iv_min = 0.05)$iv_min, 0.05)
  expect_error(scr_config(iv_min = 0.1, iv_min = 0.2), "repeated")
  expect_error(scr_config(objective = "fraud"), "objective")
  expect_error(scr_config(direction = "up"), "direction")
  expect_error(scr_config(algorithm = "nope"), "algorithm")
  expect_error(scr_config(align_method = "magic"), "align_method")
  expect_error(scr_config(challenger = "ranger"), "challenger")
  expect_error(scr_config(pdo = -1), "pdo")
  expect_error(scr_config(holdout_ratio = 1.5), "holdout_ratio")
})

test_that("target bounds are enforced and direction derives from objective", {
  cfg <- scr_config(target_min = 30, target_max = 5)
  expect_equal(cfg$target_max, 30L)
  expect_equal(resolve_direction(scr_config()), "higher_is_safer")
  expect_equal(resolve_direction(scr_config(objective = "propensity")), "higher_is_riskier")
  expect_equal(resolve_direction(scr_config(objective = "propensity", direction = "higher_is_safer")), "higher_is_safer")
})

test_that("config keys dictionary matches the config object", {
  k <- scr_config_keys()
  cfg <- scr_config()
  expect_setequal(k$key, names(cfg))
  expect_true(all(scr_config_keys(stage = 5)$stage == 5))
  expect_output(print(cfg), "scr_config")
  expect_output(print(scr_config(objective = "propensity")), "higher_is_riskier")
})
