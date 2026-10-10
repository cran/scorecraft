# scorecraft 0.3.2

* New score studies, computed from one pass over the scored rows into a
  table of counts per score value:
  * `scr_bands()` cuts a score into percentile or tail bands frozen on a
    reference sample and reports, per band and sample, the event rate with
    a Jeffreys interval, lift, capture, KS, WOE, IV, PSI and a one-sided
    Fisher exact test of rank order, plus the AUC, Gini and KS with a
    bootstrap interval. Accepts a scorecard, a data.frame (with weights,
    a value column and sample labels) or pre-aggregated counts.
  * `scr_tiers()` fits 2 to 9 tiers (labeled from "low" to "high", with
    words for 2 to 7 tiers) by an exact dynamic program under share, event
    and distinctness constraints, by event-rate anchors or by equal shares;
    `round_to` gives policy-friendly cuts, and `n_boot` measures the
    stability of the cuts. An infeasible tier count falls back and is
    recorded in a ledger.
  * `scr_rag()` lights discrimination, calibration, stability and, for a
    scorecard, the variables red, amber or green, per period or segment
    with `by`; the thresholds are an editable table from `scr_rag_plan()`.
  * `scr_apply()`, `scr_sql()` and `scr_export()` assign the bands or tiers
    in R and in SQL and write the study to a workbook. Tier labels carry
    their order in production (`"01.very high"` for the tier with the
    highest event rate, down to `"05.very low"`), so they sort from the
    event-richest tier; the tiers table has the same value in `tier_label`,
    and `numbered = FALSE` gives the plain labels.
  * `scr_claims()` tests statements such as "rate >= 60% for scores of 625
    or more" on a band, a tier or a score range, with an exact one-sided
    binomial test, a one-sided Jeffreys bound and a Holm adjustment, and
    writes each verdict ("supported", "refuted" or "not proven") as one
    sentence. `type = "floor"` tests the weakest end of the group under a
    monotone fit of the rate instead of its average.
  * `scr_operating()` finds the cut that maximizes the value of targeting,
    alerting or approval under volume, share, budget, event-rate and daily
    capacity constraints, and reports the binding constraint and its
    shadow price.
  * `scr_score_cross()` crosses two scores on the same rows: the table of
    their bands with one or two outcomes, Spearman's rho and Kendall's
    tau-b, and the overlap of the rows each selects, with the swap-in and
    swap-out rates.
  * The three write a workbook with `scr_export()`; the vignette
    "Score studies" walks through them on credit, fraud and churn scores.
  * `scr_mix_shift()` splits the change of the event rate between two
    samples, or between each period and a base, into a mix effect and a
    rate effect per band; the effects add up to the change exactly.
  * `scr_segments()` reads one score on many segments: AUC with the DeLong
    standard error and a test of equal AUC, observed against expected
    events, the log-odds offset and the slope ratio, with a suggested
    action per segment.
  * `scr_maturity()` estimates the cumulative incidence of the event by
    band and horizon under censoring (Kaplan-Meier with Greenwood's
    variance), and the discrimination at each horizon.
  * `scr_uplift()` reads a score on a treated and a control group: uplift
    per band with the Newcombe interval, the Qini coefficient and the AUUC
    with bootstrap intervals, and a randomization check.
  * `scr_overlap()` compares rule flags with the alerts of a score: what
    each catches, the incremental recall of one over the other and the
    rules the score makes redundant.
  * `scr_detection()` measures, per alert threshold, the fraud episodes
    detected, the fraudulent transactions and the time before the first
    alert, and the loss prevented.
  * The six write a workbook with `scr_export()`.
* New configuration keys for the score studies (stage 13): `study_bands`,
  `study_level`, `tier_min_pct`, `tier_min_events` and `tier_max_bins`.
* `scr_pd_validate()` gives the light `"grey"` to a row without a testable
  result (a missing p-value) instead of `NA`, and the overall light is
  `"grey"`, not `"green"`, when no row has a testable result.
  `scr_ead_validate()` does the same and adds the overall `light`. In
  `scr_lgd_validate()`, a summary row over pools or drivers that are all
  grey is now grey, not green.
* `scr_monitor()` and the CSI timelines of `scr_export()` compute the CSI
  as `scr_psi()` does: a bin empty in both samples is left out of the
  smoothing and of the degrees of freedom of the adjusted threshold, and
  the CSI and its critical value are missing when the period has no rows
  in the bins or fewer than two bins are populated. Other results are
  unchanged.
* `scr_strategy()` freezes `crossing$cut` on the training scores, like the
  bands: midway between the training scores on either side of the band
  edge of the crossing. `score >= cut`, the convention of `scr_cutoff()`,
  then reproduces the split of the bands on train and on any score seen in
  training, so rows at a band edge seen in training no longer change side
  between the two functions. `breaks` given as a number of intervals now
  also gets a cut, taken on the evaluated sample.
* The rank-order diagnostics of `scr_scorecard()` test each band against
  the previous one with a one-sided Fisher exact test. The binomial test
  used before took the previous band's rate as known and flagged too many
  breaks when that band was small.
* `scr_scorecard()` says when ties in the training score give fewer score
  bands than `score_groups`.
* `scr_metrics()` draws its bootstrap on the counts per score value when
  the score has many ties (at most one distinct value for every ten rows:
  scorecard points, a grade scale, a WOE score on a large sample), so the
  cost of a resample follows the number of distinct scores, not of rows.
  The intervals of such scores keep their distribution but change for a
  given seed, wherever they are reported; scores with fewer ties keep the
  row resampling and their intervals. The band tables of the scorecard,
  `scr_strategy()`, `scr_reject()` and `scr_psi()` assign the bands as
  integer indices, without a factor per row. The DeLong standard error of
  `scr_pd_validate()` is computed from the counts and can differ in the
  last digits (about 1e-16). No other result changes.
* The cheat sheet gains a third page on the score studies.

# scorecraft 0.3.1

* `scr_strategy()` reports the event and non-event distributions of each
  band (`pct_event`, `pct_nonevent`, `odds_event` and `log_odds`, the band
  WOE) and the score where they cross (`crossing`); `rule = "crossing"`
  sets the decisions at that boundary.
* `scr_strategy()` follows the objective of the scorecard: under propensity
  the table starts at the most likely band, `revenue_good` is the revenue
  of an event and the decisions are `"target"`, `"review"` and `"skip"`.
  Credit and fraud results are unchanged.
* `scr_score_gains()` adds `pct_event`, `pct_nonevent` and `woe`. Under
  `higher_is_riskier`, `odds` and `log_odds` are now events per non-event,
  the orientation of the scale, so `log_odds` rises with the score in both
  directions.
* The points table and the `Variable_Gains_IV` sheet of `scr_export()` add
  `pct_event` and `pct_nonevent` per bin.
* Documentation, messages and comments use American English spelling.
  Column names, configuration keys and data values are unchanged (for
  example `utilisation`, `realised` and the `"grey"` traffic light).

# scorecraft 0.3.0

* `scr_lgd_downturn()` and `scr_ead_downturn()` estimate the observed
  downturn impact, and the LGD reference value, on the training rows only,
  like the long-run averages; the hold-out stays independent evidence.
* `scr_export()`, `scr_sql(file = )` and the classing lab functions follow
  the `verbose` key of the object's configuration, as the function that
  fitted it does.
* The classing lab raises `IV_RATIO_UNSTABLE`, an advisory warning, when the
  train IV of a proposal is below `iv_min` and the hold-out/train IV ratio
  carries little information.
* `scr_irb_params()` lists the regulatory texts behind the presets; users
  check the tables against the texts in force before any regulatory use.
* Three vignettes: Get started, coarse classing, and scaling, alignment and
  challengers. The PD, LGD/EAD and capital guides are articles on the
  package website.
* A two-page cheat sheet ships in `inst/cheatsheet` (PDF and its HTML
  source).

# scorecraft 0.2.0

The IRB layer: from the scorecard to regulatory risk parameters, with the
same contracts as the scorecard pipeline (one configuration, ledgers with
mandatory reasons, hold-out revalidation with frozen bins, hardened
workbooks, production SQL verified against DuckDB and SQLite). Regimes are
parameter tables selected by a preset, never prose.

* `scr_irb_params()` ships the numbers of three presets (`"bcb"`,
  `"basel3_final"`, `"crr3"`): PD floors, LGD input floors, foundation LGD,
  standardised CCFs, asset correlations, maturity rules, output floor and
  standardised risk weights; the tables are editable and edits are recorded.
* `scr_default()` builds the default flag from a monthly panel (days past
  due with absolute and relative materiality, unlikeliness to pay,
  probation, restructuring, obligor-level pulling effect);
  `scr_default_rate()` gives the default rates by cohort, grade, segment
  and exposure, with the long-run average and its benchmark.
* `scr_bin_continuous()` bins drivers against a bounded continuous target
  (LGD, CCF) and returns an object with the shape of the engine's, so
  `OptimalBinningWoE::obwoe_apply()` and `obwoe_sql()` reproduce the bin
  means in R and in every SQL dialect; hold-out revalidation with frozen
  cut points and PSI.
* `scr_config()` gains the keys of stages 8 to 12 (`default_*`, `pd_*`,
  `lgd_*`, `ccf_*`, `framework`, `capital_*`, `ecl_*`), all registered in
  `scr_config_keys()` and validated.
* `scr_demo_panel`, `scr_demo_lgd`, `scr_demo_lgd_cashflows`,
  `scr_demo_rates`, `scr_demo_ead` and `scr_demo_portfolio` are new
  demonstration data.
* PD: `scr_master_scale()`, `scr_calibrate()` (intercept shift, log-odds
  `(a, b)`, scaling, quasi-moment matching; a new alignment, the scorecard
  untouched), `scr_grades()` (geometric, quantile or supplied grades, merges
  below the minimum counts, monotone repair recorded), `scr_moc()`
  (estimation error computed; other categories with a mandatory reason),
  `scr_pd()` (floors from the preset), `scr_migration()`,
  `scr_pd_validate()` (Jeffreys, binomial, normal, Hosmer-Lemeshow,
  multi-period, AUC against the initial value, PSI, migration bandwidths,
  concentration; traffic lights), `scr_pd_pit_ttc()`, with `predict()`,
  `scr_apply()`, `scr_sql()` (grade and PD as a `CASE` on the score) and
  `scr_export()` methods.
* LGD: `scr_workout()` (discounted recoveries and costs, cures, merged
  re-defaults, extrapolated incomplete workouts, named funnel rules),
  `scr_lgd()` (cure stage on the binary engine, severity stage on the
  continuous binner with a fractional logit or a beta regression, hold-out
  revalidation, pools), `scr_lgd_downturn()`, `scr_lgd_floor()`,
  `scr_elbe()`, `scr_lgd_validate()`, with `scr_apply()`, `scr_sql()` (both
  stages, pool `CASE`, floored result) and `scr_export()` methods.
* EAD: `scr_ead_data()` (realised conversion factors under a fixed, cohort
  or variable horizon; conversion factor below and limit factor above a
  utilisation threshold; named funnel rules), `scr_ead()` (driver bins with
  admission rules, pools, estimation-error margin, standardised floor),
  `scr_ead_downturn()`, `scr_ead_validate()`, with `scr_apply()`,
  `scr_sql()` and `scr_export()` methods.
* Expected loss and capital: `scr_el()`, `scr_irb_rw()` (the risk-weight
  function with correlations, size adjustment, maturity, floors and the
  defaulted case), `scr_sa_rw()`, `scr_capital()` (reconciliation by
  segment, output floor, provisions shortfall and excess, floors impact,
  sensitivity grid, concentration), `scr_pd_stress()`, `scr_ecl()`
  (survival-weighted 12-month and lifetime expected credit loss with
  stages and scenarios), with `scr_sql()` (constants per pool, no normal
  quantile at run time) and `scr_export()` methods.
* `scr_irb_rw()` and `scr_capital()` read the supervisory LGD of the
  foundation approach from `params$lgd_firb` through a `claim` type;
  `scr_sa_rw()` and `scr_capital()` apply the non-granular retail weight
  with `granular = FALSE`; the Hosmer-Lemeshow light of
  `scr_pd_validate()` is green when every grade sits on the conservative
  side (the PD above the observed rate), since the statistic is two-sided.
* `betareg` (Suggests) powers the beta severity engine of `scr_lgd()`.
* `scr_sql()` on a scorecard gains `what = "all"` (bin label, WOE and
  points of every variable next to the exact score and the whole-points
  score) and `keep_columns` (key columns carried into the output), for a
  deployment that reports the band of each variable with the score.
* After the documentation audit: `scr_iv()` ignores `NA` for every group
  type; `scr_classing_read()` validates the separator and the spec carries
  it into `scr_classing_import()`; the `TOO_MANY_BINS` screening rule can
  fire (the screen reads `max_bins`); `scr_psi()` stores and prints its
  thresholds; `scr_default_rate()` reports one long-run mean and
  benchmarks an optional `lra_adjusted`; one `asset_class` configuration
  key replaces `pd_asset_class` and `capital_asset_class`;
  `scr_lgd_downturn()` always records a reason; the traffic-light
  convention is red at or below the first threshold in PD, LGD and EAD;
  `scr_apply()` on an `scr_ead` takes `what`; `scr_fetch()` gains
  `verbose` and `scr_run()` follows `config$verbose`; the `scr_demo`
  columns carry English names (`vl_partial_*`, `vl_noise_*`,
  `vl_constant`, `vl_near_const`, `vl_duplicate`, `vl_redundant`,
  `vl_late`, `ds_region`, `ds_band`, `ds_channel`, `ds_high_card`).

## Compiled kernels and big tables

* The package now compiles C++ code ('Rcpp' and 'RcppArmadillo'). Thin
  internal R wrappers validate the input and call the kernels, which take
  plain vectors and matrices (never a `data.table`) and read the columns
  without copying them. The kernels are tested against their reference R
  implementations.
* Redundancy pruning in `scr_bin()` (Pearson or Spearman): each WOE
  column is ranked once and the correlation matrix comes from one BLAS
  cross-product, instead of ranking both columns of every pair. The greedy
  sweep gives exactly the same result as
  `OptimalBinningWoE::obwoe_prune()`, and is about 50 times faster at 150
  columns.
* Somers' D of LGD and EAD models is counted exactly in
  `O(n log n)` (Knight's algorithm). It replaces an `O(n^2)` Kendall
  computation inside a 200-resample bootstrap. The EAD version is also
  exact now: it no longer groups the prediction into 60 quantile buckets.
* `scr_ecl()` streams the survival-weighted loss row by row and applies
  the scenario shocks on the fly. Memory is `O(n)` whatever the term (the
  matrix version built several `n x T` copies, about 2.9 GB each at
  `n = 1e6, T = 360`). Results are the same to 1e-15.
* `scr_metrics()` ranks the scores once; each bootstrap resample then
  re-tabulates counts, with no sort and no grouping by a double key. The
  cut-off sweep sorts once per sample. The monitor tabulates the base
  once for every period. Default rates by cohort use one rolling join.
  Every LGD cash-flow aggregation and EAD reference date is vectorised.
* Kernel threads follow `config$nthread`.

## Correctness and numerical stability

* Wide tables: `scr_triage()` and the R pre-processing of `scr_apply()`
  reserve column slots before adding columns, so they no longer fail past
  about 1024 columns.
* Random numbers: seeds are local to the call. `.Random.seed` is restored
  on exit, and a bootstrap advances the user's stream only by the
  replicate seeds it draws. Results for a given seed are
  unchanged.
* `scr_split()`:
  * the target is checked for 0/1 before integer coercion (0.5 was
    truncated to 0);
  * text dates (as DBI returns them) and `integer64` columns are read
    correctly;
  * rows with a missing date are reported.
* `scr_bin()`: under `allow_derived_final = FALSE`, derived flags leave
  before the redundancy pruning, so a flag can no longer remove a real
  column. The Rcpp subset-proxy warnings of
  `OptimalBinningWoE::obwoe_gains_score()` are muffled; its values are
  correct, and the fix belongs upstream.
* `scr_config()` validates every key of stages 0 to 7.
* LightGBM receives `min_sum_hessian_in_leaf`. The xgboost API is detected
  from `xgb.train()` (it works with xgboost 3).
* SQL string literals follow the dialect: backslashes are escaped only in
  MySQL, Spark/Hive/Databricks and BigQuery. A line break in a name can no
  longer escape an SQL comment. The scorecard's SQL quotes identifiers as
  `obwoe_sql()` does. A row that falls in no fitted bin takes the points
  of WOE 0, in R and in SQL.
* Metrics and monitoring:
  * `scr_metrics()` refuses a factor or a non-0/1 outcome, and counts are
    kept in double to avoid integer overflow;
  * `scr_psi()` leaves bands empty in both samples out of the index and
    out of the degrees of freedom;
  * `scr_monitor()` keeps undated rows as a period;
  * a classing spec survives the CSV and xlsx round trips;
  * an export cannot be written outside the given directory.
* PD:
  * `scr_default()` could assign one unit's flags to another under
    locales where the grouping order differed from the C-locale sort; it
    is fixed, and the state machine is now vectorised;
  * the `PSI_ACTION` gate of the LGD and EAD drivers compared against a
    flag `scr_psi()` never returns, so it never fired; it now does;
  * Hosmer-Lemeshow on fixed PDs uses `K` degrees of freedom;
  * the multi-period test is scaled by `sd(DR_t - PD_t)` (BCBS WP 14);
  * the AUC test uses the DeLong variance;
  * `scr_moc()` no longer edits the caller's ledger;
  * point-in-time PDs accept a grade at zero.
* LGD and EAD:
  * `scr_ead_data()` paired default dates with facilities after a sort,
    so unsorted input put defaults on the wrong facility; it is fixed;
  * `scr_ead()` no longer needs data.table 1.15;
  * `ccf_measure = "lf"` keeps its pools;
  * a cure recovers `ead + drawings - recovered`;
  * the `scr_elbe()` uplift starts at zero at the date of default;
  * `scr_lgd_downturn()`, `scr_lgd_floor()` and `scr_ead_validate()` no
    longer modify their input.
* Capital and ECL:
  * the maturity adjustment is held at `PD = 1e-5` below that point,
    where `1 - 1.5 b` approaches zero and the risk weight exploded or
    turned negative without a PD floor;
  * a missing SME sales figure gets no firm-size adjustment (CRE31.9);
  * the ECL exit probability is capped at one, the LGD floored at zero,
    and scenario shocks and input lengths are validated.

## Scorecard pipeline hardening

* `scr_monitoring_plan()` is the monitoring contract: created by
  `scr_scorecard()`, written to the `Monitoring_Plan` sheet, and read back by
  `scr_monitor(plan = )` (a table or the strategy workbook), which now takes
  its PSI/CSI thresholds, alpha and `min_events_per_period` from it.
* `scr_scorecard()` stores the hold-out bin index of every variable, so the
  `Stability_CSI_Timeline` sheet is a real timeline by vintage without a
  `scr_monitor()` object.
* `options(scorecraft.parallel = "fork" | "psock" | "serial")` selects
  the parallel backend; results are identical under the three.
* `scr_triage()` is parallel by column as well.
* Workers never fail silently on any backend: a worker error is re-thrown
  with the failing item, a worker killed by the system is reported as such
  (instead of a `NULL` that surfaces later as a subscript error), warnings
  raised in a worker are re-raised in the parent, PSOCK workers run with a
  single data.table thread, and a `data.table` returned by a worker is
  re-allocated so that `:=` works on it.
* `options(scorecraft.fork_mem_fraction = 0.75)` caps the fork workers by
  the memory available on Linux (`Inf` to disable), since forked workers
  duplicate the parent heap once the garbage collector runs.

# scorecraft 0.1.0

First release. A production-grade scorecard engine for binary targets, built
on 'OptimalBinningWoE': audit funnel, single configuration, named relaxation,
first-class scale alignment, cut-off strategy and hardened deliverables.

* Seven stages, each an exported function, chained by `scr_select()` and
  `scr_scorecard()`: `scr_split()`, `scr_triage()`, `scr_bin()`,
  `scr_model()`, `scr_align()`, `scr_cutoff()`/`scr_strategy()`/`scr_reject()`.
* `scr_align()` is a first-class stage: banded log-odds regression on the raw
  score composed with the PDO map, `odds_orientation` recorded, applied to
  any engine. It runs automatically inside `scr_scorecard()`.
* `scr_sql()` emits production SQL in fourteen dialects: a pre-processing
  CTE, the WOE/BIN transformation from the authoritative cut points and, for
  a scorecard, the exact score plus whole points from the bin index. R-SQL
  equivalence is verified by test against DuckDB and SQLite.
* `scr_bin()` is parallelised by column (`nthread`), with the fits merged;
  the result is identical to the serial one.
* `scr_metrics()` always reports a bootstrap confidence interval for
  AUC/KS/Gini; `scr_psi()` reports the fixed threshold next to the
  sample-size-adjusted critical value of Yurdakul and Naranjo (2020).
* `scr_scorecard(challenger = )` fits a tree challenger (`xgboost` or
  `lightgbm`) on the same WOE columns, aligned to the same scale, with
  `supports_scorecard = FALSE`.
* `scr_reject()` implements honest reject inference: population scope, band
  coverage and a 2x/4x/8x sensitivity band, never parcelling by default.
* `scr_monitor()` recomputes PSI/CSI (with the signed points shift) and the
  performance by vintage on new data; it never schedules itself.
* `scr_export()` writes four hardened `.xlsx` workbooks (selection,
  scorecard, validation, strategy), the SQL files and a Markdown summary.
* `scr_coarse_classing()` opens a manual binning lab: `scr_classing_view()`,
  `scr_classing_propose()` (breaks, groups, merge, split, missing_to,
  other_to, reset), `scr_classing_accept()`/`scr_classing_discard()` with a
  mandatory reason, `scr_classing_choose()` (keep/drop/force),
  `scr_classing_spec()`/`scr_classing_read()`/`scr_classing_import()` for a
  CSV/xlsx round trip, `scr_classing_apply()` to commit into a new
  `scr_result`, and `scr_decisions()` for the append-only ledger. Manual
  bins share the engine's contract, so `scr_scorecard()`, `scr_apply()` and
  `scr_sql()` follow them unchanged; the funnel gains `provenance`.
* `scr_connect()` accepts any DBI driver next to an ODBC DSN, and
  `scr_fetch()` samples server-side with a dialect-aware random expression.
