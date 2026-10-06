## ----include = FALSE----------------------------------------------------------
has_glmnet <- requireNamespace("glmnet", quietly = TRUE)
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", eval = has_glmnet)
old_dt <- data.table::setDTthreads(2)

## ----setup--------------------------------------------------------------------
library(scorecraft)
cfg <- scr_config(verbose = FALSE, nthread = 1, use_glmnet = TRUE, use_ranger = FALSE,
                  use_lightgbm = FALSE, xgb_rounds = 60, n_boot = 20)
res <- scr_select(scr_demo, "default", config = cfg,
                  drop = c("id", "churn"), date_col = "ref_date")
scr_selected(res)

## ----open---------------------------------------------------------------------
lab <- scr_coarse_classing(res, author = "analyst")
lab

## ----overview-----------------------------------------------------------------
ov <- scr_classing_view(lab)

## ----view-num-----------------------------------------------------------------
scr_classing_view(lab, "vl_score_01")

## ----view-cat-----------------------------------------------------------------
scr_classing_view(lab, "ds_region")

## ----propose-breaks-----------------------------------------------------------
p_breaks <- scr_classing_propose(lab, "vl_score_01", breaks = c(40, 55, 70))
p_breaks

## ----propose-merge-split------------------------------------------------------
p_merge <- scr_classing_propose(lab, "vl_score_01", merge = c(1, 2))
p_merge$entry$cutpoints
p_split <- scr_classing_propose(lab, "vl_score_01",
                                split = c(1, res$fit$results$vl_score_01$cutpoints[1] - 5))
p_split

## ----propose-groups-----------------------------------------------------------
p_groups <- scr_classing_propose(lab, "ds_region",
                                 groups = list(edge = c("NORTH", "SOUTH"),
                                               core = c("EAST", "WEST", "CENTRE")))
p_groups

## ----propose-other------------------------------------------------------------
p_other <- scr_classing_propose(lab, "ds_region",
                                groups = list(edge = c("NORTH", "SOUTH"), rest = "EAST"),
                                other_to = "rest")
p_other$entry$bin
p_other$entry$manual$is_other

## ----propose-missing----------------------------------------------------------
p_missing <- scr_classing_propose(lab, "ds_optin", missing_to = 1)
p_missing$entry$bin
p_missing$verdict
p_missing$warnings

## ----accept-------------------------------------------------------------------
lab <- scr_classing_accept(lab, p_breaks, reason = "policy bands 40/55/70 used by underwriting")

## ----accept-groups------------------------------------------------------------
lab <- scr_classing_accept(lab, p_groups, reason = "edge/core is what pricing uses")

## ----discard------------------------------------------------------------------
lab <- scr_classing_discard(lab, p_missing, reason = "folding MISSING into YES erases the signal")

## ----blocked------------------------------------------------------------------
p_blocked <- scr_classing_propose(lab, "vl_score_04", breaks = c(-5000, 50))
p_blocked$verdict
p_blocked$blocking

## ----blocked-accept, error = TRUE---------------------------------------------
try({
scr_classing_accept(lab, p_blocked, reason = "we need this band for the policy")
})

## ----override-----------------------------------------------------------------
lab_override <- scr_classing_accept(lab, p_blocked, reason = "deliberate policy floor at -5000",
                                    override = TRUE)
scr_decisions(lab_override)[variable == "vl_score_04", .(seq, action, proposal_id, verdict, warnings, reason)]

## ----choose-------------------------------------------------------------------
scr_funnel(res, cols = "all")[feature %in% c("vl_score_10", "vl_score_03"),
                              .(feature, exit_stage, screen_reason)]
lab <- scr_classing_choose(lab, drop = "vl_score_10", force = "vl_score_03",
                           reason = c(vl_score_10 = "not available at decision time",
                                      vl_score_03 = "policy: bureau band must be scored"))

## ----summary------------------------------------------------------------------
lab

## ----ledger-------------------------------------------------------------------
scr_decisions(lab)[, .(seq, variable, action, proposal_id, verdict, reason)]

## ----spec, message = FALSE----------------------------------------------------
spec <- scr_classing_spec(lab)
spec
spec_file <- file.path(tempdir(), "classing_default.csv")
scr_classing_spec(lab, file = spec_file)

## ----edit---------------------------------------------------------------------
sheet <- read.csv(spec_file, stringsAsFactors = FALSE)
i1 <- sheet$variable == "vl_score_01" & sheet$bin_id == 1
sheet$upper[i1]  <- 42
sheet$reason[sheet$variable == "vl_score_01"] <- "reviewer: first cut moved to 42 to match the bureau band"
write.csv(sheet, spec_file, row.names = FALSE, na = "")

## ----read-refused, error = TRUE-----------------------------------------------
try({
scr_classing_read(spec_file)
})

## ----read-ok------------------------------------------------------------------
sheet$lower[sheet$variable == "vl_score_01" & sheet$bin_id == 2] <- 42
write.csv(sheet, spec_file, row.names = FALSE, na = "")
spec_back <- scr_classing_read(spec_file)

## ----import-------------------------------------------------------------------
imported <- scr_classing_import(lab, spec_back)
names(imported)
imported$vl_score_01$imported_reason
imported$vl_score_01

## ----import-accept------------------------------------------------------------
lab <- scr_classing_accept(lab, imported$vl_score_01, reason = imported$vl_score_01$imported_reason)
scr_decisions(lab)[variable == "vl_score_01", .(seq, action, proposal_id, instruction, verdict)]

## ----apply--------------------------------------------------------------------
res2 <- scr_classing_apply(lab)
res2

## ----selected-----------------------------------------------------------------
scr_selected(res2)
scr_selected(res2, "consensus")
setdiff(scr_selected(res2), scr_selected(res2, "consensus"))
setdiff(scr_selected(res2, "consensus"), scr_selected(res2))

## ----funnel-------------------------------------------------------------------
touched <- c("vl_score_01", "ds_region", "vl_score_03", "vl_score_10")
scr_funnel(res2, cols = "all")[feature %in% touched,
                               .(feature, exit_stage, provenance, manual_reason)]

## ----frozen-------------------------------------------------------------------
res2$fit_auto$results$vl_score_01$cutpoints
res2$fit$results$vl_score_01$cutpoints
res2$fit$summary[res2$fit$summary$feature %in% c("vl_score_01", "ds_region"),
                 c("feature", "algorithm", "n_bins", "total_iv")]

## ----scorecard----------------------------------------------------------------
sc <- scr_scorecard(res2)
sc
sc$points[variable == "vl_score_01", .(variable, bin, woe, points)]
sc$points[variable == "ds_region", .(variable, bin, woe, points)]

## ----price--------------------------------------------------------------------
sc_auto <- scr_scorecard(res)
rbind(scr_score_metrics(sc_auto)[, .(card = "optimal", sample, auc, auc_lo, auc_hi, ks)],
      scr_score_metrics(sc)[, .(card = "manual", sample, auc, auc_lo, auc_hi, ks)])

## ----model-card---------------------------------------------------------------
str(sc$model_card[c("binning_algorithm", "shortlist_source", "n_manual_bins",
                    "manual_bins", "forced_in", "manual_dropped", "n_decisions")])
nrow(scr_decisions(sc))

## ----apply-new----------------------------------------------------------------
new <- head(scr_demo, 5)
new[, c("vl_score_01", "ds_region")]
scr_apply(sc, new, what = "points")[, .(score, score_points, vl_score_01_points, ds_region_points)]

## ----sql----------------------------------------------------------------------
sql <- scr_sql(sc, table = "prd.customers", dialect = "databricks")
sql_lines <- unlist(strsplit(sql, "\n", fixed = TRUE))
cat(grep("^-- Provenance", sql_lines, value = TRUE), sep = "\n")
cat(grep("WHEN (vl_score_01|ds_region) ", sql_lines, value = TRUE), sep = "\n")
cat(tail(sql, 8), sep = "\n")

## ----include = FALSE, eval = TRUE---------------------------------------------
data.table::setDTthreads(old_dt)

