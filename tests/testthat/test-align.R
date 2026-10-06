test_that("a published calibration example is reproduced to the fifth decimal", {
  # ln(odds) = -1.1064394 + 1.29594211 * raw; PDO 15, 660 at 15:1
  al <- .scr_align_from(-1.1064394, 1.29594211, base_score = 660, base_odds = 15, pdo = 15,
                        direction = "higher_is_safer")
  expect_equal(al$factor, 21.64043, tolerance = 1e-6)
  expect_equal(al$offset, 601.39664, tolerance = 1e-6)
  expect_equal(al$a, 577.45282, tolerance = 1e-5)
  expect_equal(al$b, 28.04474, tolerance = 1e-5)
  expect_equal(al$odds_orientation, "safe:event")
})

test_that("direct alignment reproduces the PDO map and orientation", {
  al <- .scr_align_from(0, -1, 600, 50, 20, "higher_is_safer")
  expect_equal(al$factor, 20 / log(2))
  expect_equal(al$offset, 600 - 20 / log(2) * log(50))
  # a logit of ln(1/50) (odds safe:event = 50) must score exactly 600
  expect_equal(predict(al, log(1 / 50)), 600)
  # 20 points more halves the event odds
  expect_equal(predict(al, log(1 / 100)), 620)
  ar <- .scr_align_from(0, 1, 500, 1 / 9, 40, "higher_is_riskier")
  expect_equal(ar$odds_orientation, "event:safe")
  expect_equal(predict(ar, log(1 / 9)), 500)
  expect_equal(predict(ar, log(2 / 9)), 540)
})

test_that("regression alignment recovers a calibrated logit and warns when reversed", {
  set.seed(5)
  n <- 20000
  p <- stats::plogis(stats::rnorm(n, -1.5, 1.2))
  y <- stats::rbinom(n, 1, p)
  raw <- stats::qlogis(p)
  al <- scr_align(raw, y, direction = "higher_is_safer", n_bands = 10)
  expect_s3_class(al, "scr_align")
  expect_equal(al$calibration$method, "regression")
  expect_equal(al$calibration$slope, -1, tolerance = 0.08)
  expect_equal(al$calibration$intercept, 0, tolerance = 0.15)
  expect_gt(al$calibration$r2, 0.95)
  ar <- scr_align(raw, y, direction = "higher_is_riskier")
  expect_equal(ar$calibration$slope, 1, tolerance = 0.08)
  expect_warning(scr_align(-raw, y, direction = "higher_is_safer"), "opposite sign")
  expect_output(print(al), "calibration: ln\\(odds\\)")
})

test_that("predict(type = 'prob') inverts the scale and direct mode is honest about I and S", {
  set.seed(6)
  y <- stats::rbinom(3000, 1, 0.2)
  raw <- stats::qlogis(0.2) + y + stats::rnorm(3000)
  ad <- scr_align(raw, y, method = "direct")
  expect_equal(ad$calibration$intercept, 0)
  expect_equal(ad$calibration$slope, -1)
  expect_equal(predict(ad, raw, type = "prob"), stats::plogis(raw))
  ar <- scr_align(raw, y, method = "direct", direction = "higher_is_riskier")
  expect_equal(predict(ar, raw, type = "prob"), stats::plogis(raw))
  expect_output(print(ad), "direct")
})

test_that("alignment falls back to direct with a warning when bands are unusable", {
  y <- rep(0:1, 50)
  expect_warning(al <- scr_align(rep(0.3, 100), y), "not enough bands")
  expect_equal(al$calibration$method, "direct")
  expect_error(scr_align(1:3, 1:2), "same length")
  expect_error(scr_align(1:10, rep(0:1, 5), pdo = 0), "pdo")
})
