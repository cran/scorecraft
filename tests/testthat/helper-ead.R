# Fixtures of the EAD/CCF module: the reference data set and the pools of the
# demo snapshots are built once per session and reused; a hand-made snapshot
# table exercises the identities and the funnel rules.

ead_cfg <- function(...) do.call(cfg_test, utils::modifyList(list(n_boot = 20L, ccf_min_defaults = 20L), list(...)))

ead_demo <- function() {
  if (is.null(.fx$ead_data)) {
    .fx$ead_data <- scr_ead_data(scr_demo_ead, facility_id = "facility_id", obligor_id = "obligor_id",
                                 date_col = "ref_date", limit = "limit", drawn = "drawn", defaulted = "defaulted",
                                 drivers = c("product", "months_on_book", "dpd"), config = ead_cfg())
  }
  .fx$ead_data
}

ead_model <- function() {
  if (is.null(.fx$ead_model)) {
    .fx$ead_model <- scr_ead(ead_demo(), drivers = c("utilisation_ref", "product", "months_on_book"), config = ead_cfg())
  }
  .fx$ead_model
}

# Hand-made snapshots: one facility per rule, 13 monthly snapshots, default in
# the last month (2024-01) unless stated. Amounts are chosen so that every
# realised value is a short decimal.
ead_hand <- function() {
  dates <- seq(as.Date("2023-01-01"), by = "month", length.out = 13L)
  mk <- function(id, obl, limit, drawn, from = 1L, def_from = 13L) {
    k <- length(dates) - from + 1L
    data.frame(facility_id = id, obligor_id = obl, ref_date = dates[from:13],
               limit = rep_len(limit, k), drawn = rep_len(drawn, k), product = "card",
               defaulted = as.integer(seq_len(k) + from - 1L >= def_from), stringsAsFactors = FALSE)
  }
  rbind(
    mk("A", "oA", 1000, c(rep(400, 12), 700)),                  # CCF 0.5, LF 0.7
    mk("B", "oB", 1000, c(rep(1000, 12), 1100)),                # zero undrawn -> LF 1.1
    mk("C", "oC", 1000, c(rep(800, 12), 500)),                  # CCF -1.5 -> floored to 0
    mk("D", "oD", 1000, c(rep(1200, 12), 1300)),                # over limit -> LF 1.3
    mk("E", "oE", 1000, c(200, 300, 800), from = 11L),          # fast default: 2 months of history
    mk("F", "oF", 1000, 900, from = 13L),                       # default on the first snapshot
    mk("G", "oG", 0, rep(0, 13)),                               # no limit -> not in scope
    mk("H", "oH", 1000, c(rep(960, 12), 1000)),                 # u = 0.96 >= u* -> LF 1.0
    mk("I", "oA", 2000, c(rep(500, 12), 500), def_from = 99L),  # sibling of A: never flags itself
    mk("J", "oJ", 1000, rep(300, 13), def_from = 99L),          # never defaults
    mk("L", "oL", 1000, c(rep(400, 12), 1100))                  # CCF 7/6 > 1
  )
}
