library(testthat)

# Pseudo-dyad generators. Convention source: rMEA shuffle() (Kleinbub &
# Ramseyer, 2020), start-aligned crop to the shorter series.

# --- Fixture helpers ------------------------------------------------------

# Build a dyad list whose every series is pairwise distinct (fresh rnorm
# draws), so each surrogate column has exactly one possible source.
# `forms` picks the element shape per dyad:
#   "df"      data frame, x then y
#   "named"   list(x = , y = )
#   "rev"     list(y = , x = ) -- names, not position, decide the roles
#   "unnamed" list(<x>, <y>)
make_dyads <- function(lengths, forms = rep("named", length(lengths))) {
  Map(function(n, form) {
    x <- stats::rnorm(n)
    y <- stats::rnorm(n)
    switch(form,
      df = data.frame(x = x, y = y),
      named = list(x = x, y = y),
      rev = list(y = y, x = x),
      unnamed = list(x, y)
    )
  }, lengths, forms)
}

# Every candidate source series as (dyad, role, values).
source_pool <- function(dyad_list) {
  out <- list()
  for (i in seq_along(dyad_list)) {
    d <- dyad_list[[i]]
    xy <- if (is.data.frame(d)) {
      list(x = d[[1]], y = d[[2]])
    } else if (!is.null(names(d)) && all(c("x", "y") %in% names(d))) {
      list(x = d[["x"]], y = d[["y"]])
    } else {
      list(x = d[[1]], y = d[[2]])
    }
    out[[length(out) + 1]] <- list(dyad = i, role = "x", values = xy$x)
    out[[length(out) + 1]] <- list(dyad = i, role = "y", values = xy$y)
  }
  out
}

# Identify the source of one column: every pool series whose first
# length(col) samples equal the column exactly.
identify_source <- function(col, pool) {
  n <- length(col)
  hits <- Filter(function(s) {
    length(s$values) >= n && identical(as.double(s$values[seq_len(n)]), col)
  }, pool)
  vapply(hits, function(s) paste0(s$dyad, s$role), character(1))
}

# Identify every column; each must have exactly one source.
column_sources <- function(mat, dyad_list) {
  pool <- source_pool(dyad_list)
  vapply(seq_len(ncol(mat)), function(j) {
    src <- identify_source(mat[, j], pool)
    expect_length(src, 1)
    src[1]
  }, character(1))
}

target_y <- function(dyad_list, dyad) source_pool(dyad_list)[[2 * dyad]]$values

# =========================================================================
# --- generate_surrogate_pseudo(): source identity ------------------------
# =========================================================================

test_that("keep_roles = TRUE: each column is another dyad's y, start-aligned", {
  set.seed(101)
  # N = 5 (odd); unequal lengths; every element form; dyad 3 is shorter
  # than the target and must be excluded.
  dl <- make_dyads(
    lengths = c(50, 60, 40, 70, 55),
    forms = c("df", "named", "unnamed", "rev", "df")
  )
  expect_warning(
    expect_message(
      mat <- generate_surrogate_pseudo(dl, dyad = 1),
      "Cropped 3 partner series"
    ),
    "Excluded 1 partner series"
  )
  expect_true(is.matrix(mat))
  expect_true(is.double(mat))
  expect_equal(nrow(mat), length(target_y(dl, 1)))
  expect_equal(ncol(mat), 3) # eligible y partners: dyads 2, 4, 5

  src <- column_sources(mat, dl)
  expect_setequal(src, c("2y", "4y", "5y"))
  expect_false(anyDuplicated(src) > 0)
})

test_that("keep_roles = FALSE: partners can be another dyad's x or y", {
  set.seed(102)
  dl <- make_dyads(
    lengths = c(50, 60, 40, 70, 55),
    forms = c("df", "named", "unnamed", "rev", "df")
  )
  expect_warning(
    expect_message(
      mat <- generate_surrogate_pseudo(dl, dyad = 1, keep_roles = FALSE),
      "Cropped 6 partner series"
    ),
    "Excluded 2 partner series"
  )
  src <- column_sources(mat, dl)
  expect_setequal(src, c("2x", "2y", "4x", "4y", "5x", "5y"))
  expect_false(anyDuplicated(src) > 0)
  # The target's own series never appear.
  expect_false(any(c("1x", "1y") %in% src))
})

test_that("the target's y is read by name, not position, for named lists", {
  set.seed(103)
  dl <- make_dyads(lengths = c(30, 30, 30), forms = c("rev", "rev", "named"))
  mat <- generate_surrogate_pseudo(dl, dyad = 1)
  expect_equal(nrow(mat), 30)
  expect_setequal(column_sources(mat, dl), c("2y", "3y"))
})

test_that("N = 2 with equal lengths gives one column, no warning or message", {
  set.seed(104)
  dl <- make_dyads(lengths = c(25, 25), forms = c("unnamed", "df"))
  expect_no_warning(expect_no_message(
    mat <- generate_surrogate_pseudo(dl, dyad = 2)
  ))
  expect_equal(dim(mat), c(25L, 1L))
  expect_equal(column_sources(mat, dl), "1y")
})

test_that("integer n_surrogates samples distinct eligible partners", {
  set.seed(105)
  dl <- make_dyads(lengths = c(40, 40, 45, 50, 40, 60, 41))
  mat <- suppressMessages(suppressWarnings(
    generate_surrogate_pseudo(dl, dyad = 3, n_surrogates = 2, keep_roles = FALSE)
  ))
  expect_equal(dim(mat), c(45L, 2L))
  src <- column_sources(mat, dl)
  expect_false(anyDuplicated(src) > 0)
  # Eligible partners: series of length >= 45 from dyads other than 3.
  expect_true(all(src %in% c("4x", "4y", "6x", "6y")))
})

test_that("n_surrogates = NULL uses every eligible partner", {
  set.seed(106)
  dl <- make_dyads(lengths = c(30, 31, 32, 33, 29))
  mat <- suppressMessages(suppressWarnings(
    generate_surrogate_pseudo(dl, dyad = 1, keep_roles = FALSE)
  ))
  # Eligible: x and y of dyads 2, 3, 4 (dyad 5 is shorter).
  expect_equal(ncol(mat), 6)
})

test_that("the sources attribute records each column's dyad and role", {
  set.seed(107)
  dl <- make_dyads(lengths = c(20, 20, 20, 20))
  mat <- generate_surrogate_pseudo(dl, dyad = 2, keep_roles = FALSE)
  srcs <- attr(mat, "sources")
  expect_s3_class(srcs, "data.frame")
  expect_equal(nrow(srcs), ncol(mat))
  expect_equal(paste0(srcs$dyad, srcs$role), column_sources(mat, dl))
})

# =========================================================================
# --- generate_surrogate_pseudo(): conditions ------------------------------
# =========================================================================

test_that("no warning fires when no partner is excluded", {
  set.seed(108)
  dl <- make_dyads(lengths = c(30, 40, 50))
  expect_no_warning(suppressMessages(generate_surrogate_pseudo(dl, dyad = 1)))
})

test_that("generate_surrogate_pseudo() aborts on invalid input", {
  set.seed(109)
  dl <- make_dyads(lengths = c(30, 30, 30))

  expect_error(
    generate_surrogate_pseudo(dl[1], dyad = 1),
    "must contain at least two dyads"
  )
  expect_error(
    generate_surrogate_pseudo("not a list", dyad = 1),
    "must contain at least two dyads"
  )
  for (bad in list(0, 4, c(1, 2), 1.5, "1", NA)) {
    expect_error(
      generate_surrogate_pseudo(dl, dyad = bad),
      "must be a single index into"
    )
  }
  for (bad in list(0, -1, 1.5, c(1, 2), "2", NA)) {
    expect_error(
      generate_surrogate_pseudo(dl, dyad = 1, n_surrogates = bad),
      "must be NULL or a single positive integer"
    )
  }
  expect_error(
    generate_surrogate_pseudo(dl, dyad = 1, n_surrogates = 3),
    "exceeds the 2 eligible partner series"
  )
  expect_error(
    generate_surrogate_pseudo(dl, dyad = 1, keep_roles = NA),
    "must be TRUE or FALSE"
  )
  bad_type <- list(list(x = 1:5, y = letters[1:5]), list(x = 1:5, y = 1:5))
  expect_error(
    generate_surrogate_pseudo(bad_type, dyad = 2),
    "must be numeric"
  )
})

test_that("no eligible partner aborts, and that check runs first", {
  set.seed(110)
  dl <- make_dyads(lengths = c(50, 20, 30))
  expect_error(
    suppressWarnings(generate_surrogate_pseudo(dl, dyad = 1)),
    "No partner series is at least as long"
  )
  # n_surrogates also exceeds the (zero) eligible count; the no-eligible
  # condition must be the one reported.
  expect_error(
    suppressWarnings(generate_surrogate_pseudo(dl, dyad = 1, n_surrogates = 5)),
    "No partner series is at least as long"
  )
})

# =========================================================================
# --- generate_surrogate_pseudo(): RNG ------------------------------------
# =========================================================================

test_that("generate_surrogate_pseudo() follows set.seed and never reseeds", {
  set.seed(111)
  dl <- make_dyads(lengths = rep(30, 10))

  set.seed(7)
  a <- generate_surrogate_pseudo(dl, dyad = 1, n_surrogates = 3)
  set.seed(7)
  b <- generate_surrogate_pseudo(dl, dyad = 1, n_surrogates = 3)
  expect_identical(a, b)

  # Consecutive calls without reseeding draw fresh partners.
  set.seed(8)
  c1 <- generate_surrogate_pseudo(dl, dyad = 1, n_surrogates = 3)
  c2 <- generate_surrogate_pseudo(dl, dyad = 1, n_surrogates = 3)
  expect_false(identical(c1, c2))
})

# =========================================================================
# --- generate_pseudo_dyads(): pairing set --------------------------------
# =========================================================================

# Independently enumerated expected pairings, as keys. keep_roles = TRUE:
# ordered "<i>x|<j>y" for i != j. FALSE: unordered pairs of the 2N series
# from different dyads, key = the two labels sorted.
expected_keys <- function(n, keep_roles) {
  if (keep_roles) {
    g <- expand.grid(i = seq_len(n), j = seq_len(n))
    g <- g[g$i != g$j, ]
    return(paste0(g$i, "x|", g$j, "y"))
  }
  labels <- paste0(rep(seq_len(n), each = 2), c("x", "y"))
  dyad_of <- rep(seq_len(n), each = 2)
  keys <- character(0)
  for (a in seq_along(labels)) {
    for (b in seq_along(labels)) {
      if (a < b && dyad_of[a] != dyad_of[b]) {
        keys <- c(keys, paste(sort(c(labels[a], labels[b])), collapse = "|"))
      }
    }
  }
  keys
}

recorded_keys <- function(pd, keep_roles) {
  s <- attr(pd, "sources")
  xl <- paste0(s$x_dyad, s$x_role)
  yl <- paste0(s$y_dyad, s$y_role)
  if (keep_roles) {
    paste0(xl, "|", yl)
  } else {
    mapply(function(a, b) paste(sort(c(a, b)), collapse = "|"), xl, yl,
      USE.NAMES = FALSE
    )
  }
}

# Each element's samples must equal the first min(length) samples of the
# two sources the attribute records.
expect_samples_match_sources <- function(pd, dyad_list) {
  pool <- source_pool(dyad_list)
  get_src <- function(d, r) pool[[2 * d - (r == "x")]]$values
  s <- attr(pd, "sources")
  for (k in seq_along(pd)) {
    sx <- get_src(s$x_dyad[k], s$x_role[k])
    sy <- get_src(s$y_dyad[k], s$y_role[k])
    m <- min(length(sx), length(sy))
    expect_identical(names(pd[[k]]), c("x", "y"))
    expect_identical(pd[[k]]$x, as.double(sx[seq_len(m)]))
    expect_identical(pd[[k]]$y, as.double(sy[seq_len(m)]))
  }
}

test_that("NULL n_pairs returns every pairing, for N = 2, 3, 4", {
  forms <- c("df", "named", "unnamed", "rev")
  for (n in 2:4) {
    for (kr in c(TRUE, FALSE)) {
      set.seed(200 + n)
      dl <- make_dyads(lengths = 20 + 3 * seq_len(n), forms = forms[seq_len(n)])
      pd <- suppressMessages(generate_pseudo_dyads(dl, keep_roles = kr))
      expect_length(pd, if (kr) n * (n - 1) else 2 * n * (n - 1))

      rec <- recorded_keys(pd, kr)
      expect_false(anyDuplicated(rec) > 0)
      expect_setequal(rec, expected_keys(n, kr))

      s <- attr(pd, "sources")
      expect_s3_class(s, "data.frame")
      expect_named(s, c("x_dyad", "x_role", "y_dyad", "y_role"))
      expect_true(all(s$x_dyad != s$y_dyad))
      if (kr) {
        expect_true(all(s$x_role == "x" & s$y_role == "y"))
      }
      expect_samples_match_sources(pd, dl)
    }
  }
})

test_that("pairs are cropped start-aligned and the crop is announced", {
  set.seed(210)
  dl <- make_dyads(lengths = c(30, 40, 30))
  # keep_roles = TRUE pairs with unequal lengths: (1,2), (2,1), (2,3), (3,2).
  expect_message(
    pd <- generate_pseudo_dyads(dl),
    "Cropped 4 pseudo-dyads"
  )
  expect_true(all(lengths(lapply(pd, `[[`, "x")) == 30))

  set.seed(211)
  same <- make_dyads(lengths = c(25, 25, 25))
  expect_no_message(generate_pseudo_dyads(same))
})

test_that("integer n_pairs samples distinct pairings", {
  set.seed(212)
  dl <- make_dyads(lengths = c(30, 30, 30, 30))
  for (kr in c(TRUE, FALSE)) {
    pd <- generate_pseudo_dyads(dl, n_pairs = 5, keep_roles = kr)
    expect_length(pd, 5)
    rec <- recorded_keys(pd, kr)
    expect_false(anyDuplicated(rec) > 0)
    expect_true(all(rec %in% expected_keys(4, kr)))
    expect_equal(nrow(attr(pd, "sources")), 5)
    expect_samples_match_sources(pd, dl)
  }
})

test_that("generate_pseudo_dyads() aborts on invalid input", {
  set.seed(213)
  dl <- make_dyads(lengths = c(30, 30, 30, 30))
  expect_error(generate_pseudo_dyads(dl[1]), "must contain at least two dyads")
  for (bad in list(0, -1, 1.5, c(1, 2), "2", NA)) {
    expect_error(
      generate_pseudo_dyads(dl, n_pairs = bad),
      "must be NULL or a single positive integer"
    )
  }
  expect_error(
    generate_pseudo_dyads(dl, n_pairs = 13),
    "exceeds the 12 possible pseudo-dyads"
  )
  expect_error(
    generate_pseudo_dyads(dl, n_pairs = 25, keep_roles = FALSE),
    "exceeds the 24 possible pseudo-dyads"
  )
  expect_error(
    generate_pseudo_dyads(dl, keep_roles = "yes"),
    "must be TRUE or FALSE"
  )
})

test_that("generate_pseudo_dyads() follows set.seed and never reseeds", {
  set.seed(214)
  dl <- make_dyads(lengths = rep(30, 6))

  set.seed(9)
  a <- generate_pseudo_dyads(dl, n_pairs = 4)
  set.seed(9)
  b <- generate_pseudo_dyads(dl, n_pairs = 4)
  expect_identical(a, b)

  set.seed(10)
  c1 <- generate_pseudo_dyads(dl, n_pairs = 4)
  c2 <- generate_pseudo_dyads(dl, n_pairs = 4)
  expect_false(identical(c1, c2))
})

# =========================================================================
# --- Integration with the surrogate wrappers -----------------------------
# =========================================================================

test_that("a pseudo-dyad matrix works in all four surrogate wrappers", {
  # Equal-length, NA-free dyads from sim_dyad's three axes (first 400
  # samples, 5 s at 80 Hz, to keep WDTW fast).
  idx <- 1:400
  dl <- list(
    list(x = sim_dyad$x_A[idx], y = sim_dyad$x_B[idx]),
    list(x = sim_dyad$y_A[idx], y = sim_dyad$y_B[idx]),
    list(x = sim_dyad$z_A[idx], y = sim_dyad$z_B[idx])
  )
  y_pseudo <- generate_surrogate_pseudo(dl, dyad = 3)
  expect_equal(dim(y_pseudo), c(400L, 2L))
  x <- dl[[3]]$x
  y <- dl[[3]]$y

  res_wcc <- wcc_surrogate(x, y, y_pseudo,
    window_size = 80, lag_max = 8, window_increment = 20
  )
  res_wdtw <- wdtw_surrogate(x, y, y_pseudo,
    window_size = 80, lag_max = 8, window_increment = 20
  )
  res_wgr <- wgranger_surrogate(x, y, y_pseudo,
    window_size = 80, window_increment = 20
  )
  res_wph <- wphase_surrogate(x, y, y_pseudo,
    window_size = 80, lag_max = 8, window_increment = 20
  )

  in_unit <- function(p) length(p) == 1 && !is.na(p) && p >= 0 && p <= 1
  for (res in list(res_wcc, res_wdtw, res_wph)) {
    expect_true(in_unit(res$p_value))
    expect_equal(res$n_surrogates, ncol(y_pseudo))
  }
  expect_true(in_unit(res_wgr$p_value_xy))
  expect_true(in_unit(res_wgr$p_value_yx))
  expect_equal(res_wgr$n_surrogates, ncol(y_pseudo))
})

test_that("pseudo-dyad list elements work as estimator input", {
  set.seed(215)
  dl <- make_dyads(lengths = c(120, 150, 130))
  pd <- suppressMessages(generate_pseudo_dyads(dl))
  res <- wcc(pd[[1]]$x, pd[[1]]$y, window_size = 30, lag_max = 5)
  expect_s3_class(res, "wcc_res")
})
