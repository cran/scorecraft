# Fixtures of the capital module: the portfolio run once per session.

cap_demo <- function() {
  if (is.null(.fx$cap)) {
    .fx$cap <- scr_capital(scr_demo_portfolio, segment = "segment", asset_class = "asset_class", m = "m",
                           defaulted = "defaulted", elbe = "elbe", provisions = "provision", ltv = "ltv",
                           rating = "rating", sales = "sales", transactor = "transactor", grade = "grade",
                           id = "id", config = cfg_test(), keep_rows = TRUE)
  }
  .fx$cap
}

cap_args <- function(...) {
  utils::modifyList(list(segment = "segment", asset_class = "asset_class", m = "m", defaulted = "defaulted",
                         elbe = "elbe", provisions = "provision", ltv = "ltv", rating = "rating", sales = "sales",
                         transactor = "transactor", grade = "grade", id = "id"), list(...))
}
