# Regression tests of the capital / ECL triage: hand-computed IRB values,
# numerical stability of the maturity adjustment, SME sales, maturity
# reporting, vectorised lookups and the ECL scenario guards.

k_ref <- function(pd, lgd, m, r) {
  b <- (0.11852 - 0.05478 * log(pd))^2
  lgd * (stats::pnorm((stats::qnorm(pd) + sqrt(r) * stats::qnorm(0.999)) / sqrt(1 - r)) - pd) *
    (1 + (m - 2.5) * b) / (1 - 1.5 * b)
}
r_corp <- function(pd) { w <- (1 - exp(-50 * pd)) / (1 - exp(-50)); 0.12 * w + 0.24 * (1 - w) }

test_that("scr_irb_rw() matches an independent implementation of CRE31", {
  p <- scr_irb_params("basel3_final")
  for (m in c(1, 2.5, 5)) {
    expect_equal(scr_irb_rw(0.02, 0.45, m = m, asset_class = "corporate", params = p)$k,
                 k_ref(0.02, 0.45, m, r_corp(0.02)), tolerance = 1e-12)
  }
  # SME firm-size adjustment at S = 20 (EUR m): R - 0.04 (1 - 15 / 45)
  expect_equal(scr_irb_rw(0.01, 0.45, m = 2.5, asset_class = "corporate_sme", sales = 20, params = p)$k,
               k_ref(0.01, 0.45, 2.5, r_corp(0.01) - 0.04 * (1 - 15 / 45)), tolerance = 1e-12)
  # 1 % / 45 % / 2.5y corporate: K = 0.0738534
  expect_equal(round(scr_irb_rw(0.01, 0.45, m = 2.5, asset_class = "corporate", params = p)$k, 7), 0.0738534)
})

test_that("the maturity adjustment stays finite and K increasing without a PD floor", {
  p <- scr_irb_params("basel3_final")
  pd <- c(0, 1e-8, 1e-7, 1e-6, 2.9e-6, 3e-6, 5e-6, 1e-5, 2e-5, 3e-5, 1e-4, 3e-4, 1e-3)
  for (m in c(1, 2.5, 5)) {
    r <- scr_irb_rw(pd, 0.45, m = m, asset_class = "sovereign", params = p)
    expect_true(all(is.finite(r$ma)) && all(r$ma > 0))
    expect_true(all(r$k >= 0))
    expect_true(all(diff(r$k) >= 0))
  }
  # unchanged at and above PD = 1e-5
  q <- c(1e-5, 3e-5, 1e-4)
  expect_equal(scr_irb_rw(q, 0.45, m = 5, asset_class = "sovereign", params = p)$k,
               k_ref(q, 0.45, 5, r_corp(q)), tolerance = 1e-12)
  # PD = 1 on a performing row gives K = 0, as a defaulted row with ELBE = LGD
  expect_equal(scr_irb_rw(1, 0.45, asset_class = "corporate", params = p)$k, 0)
})

test_that("missing SME sales take no firm-size adjustment", {
  p <- scr_irb_params("basel3_final")
  rc <- scr_irb_rw(0.01, 0.45, asset_class = "corporate", params = p)$r
  rs <- scr_irb_rw(0.01, 0.45, asset_class = "corporate_sme", sales = c(NA, 5, 50), params = p)$r
  expect_equal(rs[1], rc, tolerance = 1e-12)
  expect_equal(rs[2], rc - 0.04, tolerance = 1e-12)
  expect_equal(rs[3], rc, tolerance = 1e-12)
  expect_equal(scr_irb_rw(0.01, 0.45, asset_class = "corporate_sme", params = p)$r, rc, tolerance = 1e-12)
})

test_that("maturity is reported on wholesale rows only and F-IRB does not count clipping", {
  r <- scr_irb_rw(0.01, 0.2, m = 3, asset_class = c("retail_mortgage", "corporate"))
  expect_true(is.na(r$m[1])); expect_equal(r$m[2], 3)
  f <- scr_irb_rw(0.01, 0.45, m = c(0.2, 9), asset_class = "corporate", approach = "firb")
  expect_equal(f$m, c(2.5, 2.5))
  expect_equal(attr(f, "floors_hit")$n[3:4], c(0L, 0L))
  expect_error(scr_irb_rw(0.01, 0.45, m = -1, asset_class = "corporate"), "`m`")
})

test_that("the vectorised LGD-floor and LTV-band lookups keep their semantics", {
  p <- scr_irb_params("bcb")
  l <- scr_irb_rw(0.01, 0.01, asset_class = c("corporate", "corporate", "retail_other", "retail_mortgage", "qrre_revolver"),
                  collateral = c("receivables", "other_physical", "financial", "real_estate", "unsecured"), params = p)
  expect_equal(l$lgd_used, c(0.10, 0.15, 0.01, 0.05, 0.50))
  expect_equal(scr_sa_rw("retail_mortgage", ltv = c(-0.1, 0, 0.5, 0.500001, 0.6, 0.8, 0.9, 1, 1.01, Inf, NA)),
               c(0.20, 0.20, 0.20, 0.25, 0.25, 0.30, 0.40, 0.50, 0.70, 0.70, 0.70))
})

test_that("scr_pd_stress() rejects PDs outside [0, 1]", {
  expect_error(scr_pd_stress(1.2, 0.1, 0.9), "pd")
  expect_error(scr_pd_stress(-0.1, 0.1, 0.9), "pd")
  expect_true(is.na(scr_pd_stress(NA_real_, 0.1, 0.9)))
})

test_that("ECL scenario shocks cannot make the survival or the LGD negative", {
  cfg <- cfg_test(ecl_discount = "none")
  # hazard 0.9 after the shock plus prepayment 0.5: everything exits in month one
  e <- scr_ecl(0.3, 0.5, 100, prepay = 0.5, t_max = 6L, scenarios = list(bad = list(pd_mult = 3)), config = cfg)
  expect_equal(e$totals$ecl_12m, 0.9 * 0.5 * 100, tolerance = 1e-12)
  g <- scr_ecl(0.01, 0.1, 100, scenarios = list(good = list(lgd_add = -0.2)), config = cfg)
  expect_equal(g$totals$ecl, 0)
  expect_error(scr_ecl(0.01, 0.4, 100, scenarios = list(a = list(pd_mult = -1)), config = cfg), "pd_mult")
  expect_error(scr_ecl(0.01, 0.4, 100, scenarios = list(a = list(ead_mult = -1)), config = cfg), "ead_mult")
  expect_error(scr_ecl(0.01, 0.4, 100, scenarios = list(a = list(z = "x")), config = cfg), "z")
})

test_that("ECL input lengths and stages are validated", {
  cfg <- cfg_test(ecl_discount = "none")
  expect_error(scr_ecl(rep(0.01, 5), c(0.4, 0.5), 100, config = cfg), "length")
  expect_error(scr_ecl(0.01, 0.4, 100, stage = 2.5, config = cfg), "stage")
  expect_equal(scr_ecl(rep(0.01, 2), 0.4, 100, stage = c(1, 2), config = cfg)$stages$n, c(1L, 1L, 0L))
  # the vectorised cumulative PD is exact at the boundaries and for tiny hazards
  e <- scr_ecl(matrix(c(1e-12, 1, 0), 3, 12), 0.4, 100, config = cfg, keep_rows = TRUE)
  expect_equal(e$exposures$pd_12m, c(-expm1(12 * log1p(-1e-12)), 1, 0), tolerance = 1e-12)   # 1 - (1 - h)^12 loses 5 digits here
  expect_gt(e$exposures$pd_12m[1], 0)
})
