# Fixtures of the PD module: the calibration, the grades and the default-rate
# series of the demo panel are built once per session and reused.

pd_cal <- function() {
  if (is.null(.fx$pd_cal)) .fx$pd_cal <- scr_calibrate(sc_demo(), target = 0.06)
  .fx$pd_cal
}

pd_grades <- function() {
  if (is.null(.fx$pd_gr)) {
    .fx$pd_gr <- scr_grades(sc_demo(), pd_cal(), n_grades = 7L, min_obligors = 30L, min_defaults = 10L)
  }
  .fx$pd_gr
}

# the demo panel flagged with the default engine, plus its behavioural score
pd_panel <- function() {
  if (is.null(.fx$pd_panel)) {
    d <- scr_default(scr_demo_panel, "id", "ref_date", dpd = "dpd", arrears = "arrears", exposure = "exposure",
                     config = cfg_test())
    p <- merge(d$flags, data.table::as.data.table(scr_demo_panel)[, list(id, date = ref_date, score)], by = c("id", "date"))
    .fx$pd_panel <- p
  }
  data.table::copy(.fx$pd_panel)
}

# default-rate series keyed by the final grades of pd_grades()
pd_series <- function() {
  if (is.null(.fx$pd_dr)) {
    p <- pd_panel()
    p[, grade := predict(pd_grades(), score = score)]
    .fx$pd_dr <- scr_default_rate(p, id = "id", date = "date", default = "default", grade = "grade", by = "quarter",
                                  config = cfg_test())
  }
  .fx$pd_dr
}

pd_model <- function() {
  if (is.null(.fx$pd_model)) {
    gr <- scr_moc(pd_grades(), "C", method = "ci_timeseries", dr = pd_series())
    gr <- scr_moc(gr, "A", value = 0.001, reason = "unlikeliness-to-pay trigger not available before 2024")
    .fx$pd_model <- scr_pd(gr)
  }
  .fx$pd_model
}
