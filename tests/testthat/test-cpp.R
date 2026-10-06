# Compiled kernels against their reference R implementations

test_that("the correlation kernel matches stats::cor(), constant columns give NA", {
  set.seed(11)
  n <- 400
  x <- data.frame(a = round(rnorm(n), 1), b = sample(1:4, n, TRUE), c = rnorm(n))
  x$d <- x$a + rnorm(n, 0, 0.4)
  for (m in c("pearson", "spearman")) {
    expect_equal(.scr_cor_matrix(x, m), cor(x, method = m), tolerance = 1e-12)
  }
  x$k <- 1
  cm <- .scr_cor_matrix(x, "spearman")
  expect_true(all(is.na(cm["k", ])) && all(is.na(cm[, "k"])))
  expect_error(cpp_cor_matrix(list(c(1, NA, 3), c(1, 2, 3)), TRUE, 1L), "non-finite")
  expect_error(cpp_cor_matrix(list(c(1, 2, 3), c(1, 2)), TRUE, 1L), "length")
})

test_that("the greedy pruning is identical to obwoe_prune()", {
  set.seed(12)
  n <- 1500; p <- 25
  z <- matrix(rnorm(n * 5), n, 5)
  x <- as.data.frame(vapply(seq_len(p), function(j) round(z[, (j %% 5) + 1] + rnorm(n, 0, runif(1, 0.1, 1.2)), 1), numeric(n)))
  names(x) <- paste0("v", seq_len(p))
  rk <- sample(names(x))[-(1:3)]   # three variables outside the ranking
  for (m in c("pearson", "spearman")) for (cut in c(0.4, 0.7)) {
    ref <- OptimalBinningWoE::obwoe_prune(x, ranking = rk, cutoff = cut, method = m)
    got <- .scr_prune_matrix(.scr_cor_matrix(x, m), rk, cut)
    expect_identical(got$keep, ref$keep)
    expect_identical(got$dropped$variable, ref$dropped$variable)
    expect_identical(got$dropped$correlated_with, ref$dropped$correlated_with)
    expect_equal(got$dropped$correlation, ref$dropped$correlation, tolerance = 1e-12)
  }
})

test_that("Somers' D from the concordance kernel matches the O(n^2) definition", {
  brute <- function(p, r) {
    s <- outer(p, p, function(a, b) sign(a - b)) * outer(r, r, function(a, b) sign(a - b))
    tr <- outer(r, r, `==`)
    (sum(s[upper.tri(s)] > 0) - sum(s[upper.tri(s)] < 0)) / sum(!tr[upper.tri(tr)])
  }
  set.seed(13)
  for (k in 1:6) {
    p <- round(runif(150), sample(1:2, 1)); r <- round(runif(150) * p, 1)
    expect_equal(.scr_somers(p, r), brute(p, r), tolerance = 1e-14)
  }
  cc <- cpp_concordance(c(1, 2, 2, 3), c(1, 1, 2, 3))
  expect_equal(unname(cc), c(4, 6, 1, 1))              # C - D, pairs, ties on p, ties on r
  expect_true(is.na(.scr_somers(rep(1, 5), 1:5)))      # constant prediction: NA by default
  expect_equal(.scr_somers(rep(1, 5), 1:5, const_p = 0), 0)
  expect_true(is.na(.scr_somers(1:5, rep(2, 5))))      # no pair untied on the outcome
  expect_equal(.scr_somers(c(1, NA, 2, Inf, 3), c(1, 5, 2, 1, 3)), 1)   # non-finite pairs dropped
})

test_that("the ECL kernel matches the matrix formulation", {
  set.seed(14)
  n <- 30; H <- 24
  h <- matrix(runif(n * H, 0, 0.02), n, H); L <- runif(n, 0.2, 0.6); E <- runif(n, 100, 1000)
  P <- matrix(0.01, n, H); r <- runif(n, 0, 0.1); st3 <- rep(c(FALSE, FALSE, TRUE), 10)
  ref <- function(h, L, E, P, r, hz, z = NULL, mult = NULL, add = NULL, emult = NULL, rho = 0.15) {
    Lm <- matrix(L, n, H); Em <- matrix(E, n, H)
    if (!is.null(z)) h <- stats::pnorm((stats::qnorm(h) - sqrt(rho) * z) / sqrt(1 - rho))
    if (!is.null(mult)) h <- pmin(h * mult, 1)
    if (!is.null(add)) Lm <- pmax(Lm + add, 0)
    if (!is.null(emult)) Em <- Em * emult
    DF <- outer(1 + r, -(seq_len(H)) / 12, `^`)
    S <- rep(1, n); life <- numeric(n); m12 <- numeric(n)
    for (j in seq_len(H)) {
      life <- life + S * h[, j] * Lm[, j] * Em[, j] * DF[, j]
      if (j == hz) m12 <- life
      S <- S * (1 - pmin(h[, j] + P[, j], 1))
    }
    life[st3] <- m12[st3] <- Lm[st3, 1] * Em[st3, 1]
    list(ecl_12m = m12, ecl_life = life)
  }
  a <- cpp_ecl_paths(h, L, E, P, r, n, H, 12L, TRUE, st3, NULL, NULL, NULL, NULL, 0.15, 1L)
  b <- ref(h, L, E, P, r, 12L)
  expect_equal(a$ecl_12m, b$ecl_12m, tolerance = 1e-13)
  expect_equal(a$ecl_life, b$ecl_life, tolerance = 1e-13)
  expect_equal(a$pd_life, 1 - apply(1 - h, 1, prod), tolerance = 1e-13)
  a <- cpp_ecl_paths(h, L, E, P, r, n, H, 12L, TRUE, st3, -1.5, 3, -0.1, 1.2, 0.15, 2L)
  b <- ref(h, L, E, P, r, 12L, z = -1.5, mult = 3, add = -0.1, emult = 1.2)
  expect_equal(a$ecl_life, b$ecl_life, tolerance = 1e-12)
  expect_error(cpp_ecl_paths(h, L, E, NULL, r, n, H, 30L, TRUE, st3, NULL, NULL, NULL, NULL, 0.15, 1L), "dimensions")
  expect_error(cpp_ecl_paths(h[, 1:3], L, E, NULL, r, n, H, 1L, TRUE, st3, NULL, NULL, NULL, NULL, 0.15, 1L), "must be")
})

test_that("kernel thread counts honour the check limit", {
  withr::local_envvar(c("_R_CHECK_LIMIT_CORES_" = "TRUE"))
  expect_equal(.scr_threads(8L), 2L)
  withr::local_envvar(c("_R_CHECK_LIMIT_CORES_" = ""))
  expect_equal(.scr_threads(8L), 8L)
  expect_equal(.scr_threads(NA), 1L)
  expect_equal(.scr_threads(NULL), 1L)
})
