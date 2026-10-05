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
