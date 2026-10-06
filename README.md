# scorecraft <img src="man/figures/logo.png" align="right" height="139" alt="scorecraft hex logo" />

<!-- badges: start -->
<!-- [![CRAN status](https://www.r-pkg.org/badges/version/scorecraft)](https://CRAN.R-project.org/package=scorecraft) -->
[![R-CMD-check](https://github.com/evandeilton/scorecraft/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/evandeilton/scorecraft/actions/workflows/R-CMD-check.yaml)
<!-- [![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/scorecraft)](https://CRAN.R-project.org/package=scorecraft) -->
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

A production-grade scorecard engine for binary targets (credit risk, fraud,
propensity), built on
[OptimalBinningWoE](https://CRAN.R-project.org/package=OptimalBinningWoE). It
selects variables through optimal binning and a multi-strategy consensus,
fits the points scorecard with an explicit and auditable **scale
alignment**, sweeps cut-offs with frozen cuts, performs honest reject
inference, monitors PSI/CSI with both the fixed and the sample-size-adjusted
threshold, and emits **production SQL** whose output is verified against R
by test. An IRB layer takes the scorecard to regulatory parameters: default
definition, PD calibration with rating grades and margins of conservatism,
workout LGD, EAD conversion factors, expected loss, risk weights, capital
and ECL, with the same ledgers, SQL and workbooks.

Every stage is an exported function; `scr_select()` and `scr_scorecard()`
are the shortcuts that chain them.

Documentation site: <https://evandeilton.github.io/scorecraft/>.

## Installation

```r
# from CRAN
install.packages("scorecraft")

# development version from GitHub
# install.packages("pak")
pak::pak("evandeilton/scorecraft")
```

`xgboost` and `OptimalBinningWoE` are required. `glmnet`, `ranger` and
`lightgbm` are optional consensus voters; `openxlsx` is needed by
`scr_export()`; `DBI` plus a driver by `scr_connect()`.

## The pipeline in one screen

```r
library(scorecraft)

cfg <- scr_config("moderate", objective = "risk", nthread = 4)

# Selection: split, triage, binning + screening, multi-strategy consensus
res <- scr_select(scr_demo, "default", config = cfg, drop = c("id", "churn"), date_col = "ref_date")
res
scr_selected(res)                 # the shortlist
scr_funnel(res)                   # every input column and the stage it died at
scr_leakage(res)                  # suspicious IV and degenerate bins

# Scorecard and alignment: points scorecard, aligned to 600 points at 50:1 with PDO 20
sc <- scr_scorecard(res, challenger = "xgboost")
sc
sc$alignment                      # ln(odds) = I + S * logit -> score = a + b * logit
scr_score_metrics(sc)             # AUC/KS/Gini with bootstrap CI, per sample
scr_score_gains(sc, "holdout")    # frozen bands, from the risky to the safe side

# Cut-off, strategy with marginal expected profit, honest reject inference
scr_cutoff(sc)
scr_strategy(sc, revenue_good = 1080, loss_bad = 4500)
scr_reject(sc)

# Coarse classing lab: manual bins and manual variable choice, with a reason
lab <- scr_coarse_classing(res)
p   <- scr_classing_propose(lab, "ds_region",
                            groups = list(edge = c("NORTH", "SOUTH"), core = c("EAST", "WEST", "CENTRE")))
lab <- scr_classing_accept(lab, p, reason = "edge/core is what pricing uses")
lab <- scr_classing_choose(lab, drop = "vl_score_10", reason = "not available at decision time")
res2 <- scr_classing_apply(lab)      # a new scr_result; scr_scorecard(res2) refits on it
scr_decisions(res2)                  # the append-only decision ledger

# Production: R and SQL give the same numbers
scr_apply(sc, newdata)
cat(scr_sql(sc, table = "prd.customers", dialect = "databricks"), sep = "\n")

# Monitoring, when the analyst decides to run it
scr_monitor(sc, newdata, date_col = "ref_date", target = "default")

# Deliverables: four hardened .xlsx workbooks, the SQL files, a Markdown summary
scr_export(sc, "output")
```

## From the scorecard to IRB risk parameters

Every regime-specific number is a table selected by a preset; the
functions read the tables, never the law.

```r
params <- scr_irb_params("bcb")          # or "basel3_final", "crr3"; editable tables

# Default flag from a monthly panel, then one-year default rates by cohort
d  <- scr_default(scr_demo_panel, id = "id", date = "ref_date", dpd = "dpd",
                  arrears = "arrears", exposure = "exposure", config = cfg)
dr <- scr_default_rate(d, by = "quarter")        # long-run average and its benchmark

# PD: calibrate the scorecard to the central tendency, cut grades, add MoC and the floor
cal <- scr_calibrate(sc, target = dr)             # a new alignment; the scorecard is untouched
gr  <- scr_grades(sc, calibration = cal, n_grades = 8)
gr  <- scr_moc(gr, category = "C", method = "ci_binomial")   # estimation error, computed
pd  <- scr_pd(gr, params = params, asset_class = "retail_other")
pnl <- merge(d$flags, scr_demo_panel[, c("id", "ref_date", "score")],
             by.x = c("id", "date"), by.y = c("id", "ref_date"))
scr_pd_validate(pd, pnl, id = "id", date = "date", default = "default", score = "score")

# LGD: workout cash flows -> realised LGD -> cure x severity -> pools -> downturn -> floor
wo  <- scr_workout(scr_demo_lgd, scr_demo_lgd_cashflows, rates = scr_demo_rates, config = cfg)
lgd <- scr_lgd(wo, drivers = c("product", "ltv", "months_on_book", "prior_dpd_max"), config = cfg)
lgd <- scr_lgd_downturn(lgd, periods = data.frame(start = as.Date("2022-01-01"), end = as.Date("2023-12-31")),
                        reason = "reference rate above 13% in 2022-2023")
lgd <- scr_lgd_floor(lgd, params = params, asset_class = "retail_other")

# EAD: facility snapshots -> realised conversion factors -> pools with the standardised floor
rds <- scr_ead_data(scr_demo_ead, facility_id = "facility_id", date_col = "ref_date",
                    limit = "limit", drawn = "drawn", defaulted = "defaulted",
                    drivers = c("product", "months_on_book"), config = cfg)
ead <- scr_ead(rds, drivers = c("product", "utilisation_ref", "months_on_book"), config = cfg)

# Expected loss, risk weights, capital, expected credit loss
cap <- scr_capital(scr_demo_portfolio, segment = "segment", asset_class = "asset_class",
                   provisions = "provision", params = params, config = cfg)
cap$totals; cap$segments
hz <- 1 - (1 - scr_demo_portfolio$pd)^(1 / 12)      # annual PD to a monthly hazard
scr_ecl(hz, lgd = scr_demo_portfolio$lgd, ead = scr_demo_portfolio$ead, eir = scr_demo_portfolio$eir)

# The same contracts as the scorecard: scoring, SQL, workbooks
scr_apply(pd, newdata); scr_apply(lgd, new_defaults); scr_apply(ead, new_facilities)   # each takes the columns its model was fitted on
scr_sql(pd, table = "prd.customers", dialect = "databricks")
scr_export(pd, "output"); scr_export(lgd, "output"); scr_export(ead, "output"); scr_export(cap, "output")
```

## Why another scorecard package

* **Alignment is a stage, not a footnote.** `scr_align()` takes the raw
  score of any engine to the declared scale by regressing empirical log-odds
  on score bands and composing with the PDO map, and records
  `odds_orientation`. Two scorecards aligned this way are comparable point
  for point.
* **The funnel is the deliverable.** No candidate disappears from the
  report: every input column carries the exact stage it failed at and why.
* **Production SQL is the core, not an export.** Two blocks (pre-processing
  CTE, WOE/BIN from the authoritative cut points), a third for the score,
  fourteen dialects, equivalence with `scr_apply()` verified by test in
  DuckDB and SQLite.
* **Both thresholds, always.** PSI/CSI report the market's 0.10/0.25 next to
  the sample-size-adjusted critical value; AUC/KS/Gini always come with a
  bootstrap interval.
* **Honest reject inference.** Population scope, band coverage and a
  sensitivity band instead of a single invented multiplier.
* **Manual binning with an audit trail.** `scr_coarse_classing()` opens a
  lab where the analyst proposes breaks or groupings, reads the comparison
  against the optimal bins on train and hold-out, accepts or discards with a
  mandatory reason, forces or drops variables, and commits to a new result
  that the scorecard, the R scoring and the SQL follow unchanged. The spec
  round-trips through CSV/xlsx for a business reviewer.
* **A challenger that never pretends.** A tree model aligned to the same
  scale for comparison, with `supports_scorecard = FALSE`: no points, no
  reason codes.

## Reading conventions

`objective = "risk"` means `target = 1` is the bad case and more points are
safer (`higher_is_safer`, odds `safe:event`). `objective = "propensity"`
means `target = 1` is the good case and more points are more likely
(`higher_is_riskier`, odds `event:safe`). `objective` never changes the
selection; `event_level` does.
