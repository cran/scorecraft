# Fixtures of the LGD module, fitted once per session (the `.fx` pattern of
# helper-scorecraft.R). The hand-made workout below has numbers that can be
# checked with a pocket calculator.

lgd_cfg <- function(...) {
  base <- list(verbose = FALSE, nthread = 1L, n_boot = 10L)
  do.call(scr_config, utils::modifyList(base, list(...)))
}

wo_demo <- function() {
  if (is.null(.fx$wo)) {
    .fx$wo <- scr_workout(scr_demo_lgd, scr_demo_lgd_cashflows, rates = scr_demo_rates, config = lgd_cfg())
  }
  .fx$wo
}

lgd_drivers <- function() c("product", "ltv", "prior_dpd_max", "months_on_book", "region")

lgd_demo <- function() {
  if (is.null(.fx$lgd)) .fx$lgd <- scr_lgd(wo_demo(), drivers = lgd_drivers(), config = lgd_cfg())
  .fx$lgd
}

# a downturn-and-floored copy of the demo model
lgd_final <- function() {
  if (is.null(.fx$lgd_final)) {
    m <- scr_lgd_downturn(lgd_demo(), periods = data.frame(start = as.Date("2022-01-01"), end = as.Date("2023-12-31")),
                          reason = "reference rate above 13% in 2022-2023")
    .fx$lgd_final <- scr_lgd_floor(m, asset_class = "retail_other", secured_share = 0.4)
  }
  .fx$lgd_final
}

# hand-made workout: three events, flat 12% rate + 0 add-on
hand_defaults <- function() {
  data.frame(
    default_id = c("A", "B", "C"), facility_id = c("f1", "f2", "f3"),
    default_date = as.Date(c("2024-01-15", "2024-01-15", "2024-01-15")),
    ead = c(1000, 2000, 1000), product = c("p", "p", "q"), status = c("closed", "closed", "cured"),
    close_date = as.Date(c("2025-01-15", "2025-01-15", "2024-04-15")), stringsAsFactors = FALSE)
}
hand_cashflows <- function() {
  data.frame(
    default_id = c("A", "A", "A", "B", "C"),
    date = as.Date(c("2024-07-15", "2025-01-15", "2024-02-15", "2024-07-15", "2024-02-15")),
    amount = c(600, 300, 50, 500, 100), type = c("recovery", "recovery", "direct_cost", "recovery", "recovery"),
    stringsAsFactors = FALSE)
}

# can a PSOCK worker see the internals of the namespace under development?
psock_has_dev_ns <- function() {
  # a PSOCK worker resolves the package namespace by name, so it sees the
  # INSTALLED scorecraft: the probe must look inside that namespace
  ok <- tryCatch(withr::with_options(list(scorecraft.parallel = "psock"),
                   .scr_lapply(1:2, function(i) {
                     ns <- tryCatch(asNamespace("scorecraft"), error = function(e) NULL)
                     !is.null(ns) && exists(".cbin_num", envir = ns, inherits = FALSE)
                   }, nthread = 2L)),
                 error = function(e) list(FALSE))
  all(unlist(ok))
}
