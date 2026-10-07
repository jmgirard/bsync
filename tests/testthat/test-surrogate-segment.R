library(testthat)

# Segment-shuffling surrogates. Convention source: SUSY 0.1.0 susy()
# (Tschacher & Meier, 2020); source note cairn/references/tschacher2020.md.
# Order rule: uniform over all orders except the original (D-004; review
# cairn/reviews/archive/RR01-segment-surrogate-size.md, Q1).
#
# Oracle records (DESIGN.md section 13):
#   segment-closed-form (closed-form): segment_reference() below, a plain
#     reimplementation with explicit loops. It lists the orders other than
#     the original by stepping through every permutation in lexicographic
#     order (next_permutation()), not by the generator's matrix build, and
#     it uses the generator's random-number order: up to 8 segments, one
#     sample.int(k! - 1, n_surrogates) draw from that list; above 8, one
#     sample.int(k) per try until n_surrogates distinct non-identity orders
#     are accepted. Asserted in the "matches the plain reference" test.
#   segment-invariants (invariant): every column is a reordering of the
#     segments of y other than the original, the tail stays, NA values move
#     with their segments, and the segment cut equals the eps loop of SUSY's
#     susy() (tschacher2020). Asserted in the AC1 tests below.
#   segment-uniformity (simulation-coverage): chi-square tests that single
#     orders, pairs of orders, and the first-placed segment are equally
#     frequent. Asserted in the AC3 tests below.
#   segment-calibration (simulation-coverage): test-surrogate-calibration.R.

# --- Plain reference -------------------------------------------------------

# Lexicographic successor of a permutation, or NULL after the last one.
next_permutation <- function(p) {
  k <- length(p)
  i <- k - 1
  while (i >= 1 && p[i] >= p[i + 1]) {
    i <- i - 1
  }
  if (i < 1) {
    return(NULL)
  }
  j <- k
  while (p[j] <= p[i]) {
    j <- j - 1
  }
  tmp <- p[i]
  p[i] <- p[j]
  p[j] <- tmp
  p[(i + 1):k] <- rev(p[(i + 1):k])
  p
}

is_identity <- function(p) {
  for (i in seq_along(p)) {
    if (p[i] != i) {
      return(FALSE)
    }
  }
  TRUE
}

# Every order of 1..k except the original, one per row, in lexicographic
# order.
all_orders_plain <- function(k) {
  out <- list()
  p <- seq_len(k)
  while (!is.null(p)) {
    if (!is_identity(p)) {
      out[[length(out) + 1]] <- p
    }
    p <- next_permutation(p)
  }
  matrix(unlist(out), ncol = k, byrow = TRUE)
}

segment_reference <- function(y, segment_size, n_surrogates) {
  n <- length(y)
  k <- floor(n / segment_size)
  if (k <= 8) {
    orders <- all_orders_plain(k)
    if (n_surrogates < nrow(orders)) {
      orders <- orders[sample.int(nrow(orders), n_surrogates), , drop = FALSE]
    }
  } else {
    orders <- matrix(0L, nrow = 0, ncol = k)
    while (nrow(orders) < n_surrogates) {
      p <- sample.int(k)
      seen <- FALSE
      for (r in seq_len(nrow(orders))) {
        if (all(orders[r, ] == p)) {
          seen <- TRUE
        }
      }
      if (!is_identity(p) && !seen) {
        orders <- rbind(orders, p)
      }
    }
  }
  out <- matrix(0, nrow = n, ncol = nrow(orders))
  for (j in seq_len(nrow(orders))) {
    for (i in seq_len(k)) {
      for (t in seq_len(segment_size)) {
        out[(i - 1) * segment_size + t, j] <-
          y[(orders[j, i] - 1) * segment_size + t]
      }
    }
    if (k * segment_size < n) {
      for (r in (k * segment_size + 1):n) {
        out[r, j] <- y[r]
      }
    }
  }
  out
}

# For one column of the generator's output on y = 1:n, return the source
# segment of each position (NA if the block is not a whole segment of y).
source_segments <- function(col, segment_size, k) {
  vapply(seq_len(k), function(i) {
    block <- col[((i - 1) * segment_size + 1):(i * segment_size)]
    q <- (block[1] - 1) / segment_size + 1
    whole <- as.double(((q - 1) * segment_size + 1):(q * segment_size))
    if (q == round(q) && q >= 1 && q <= k && identical(block, whole)) {
      q
    } else {
      NA_real_
    }
  }, numeric(1))
}

# The segment order of every column of a surrogate matrix built from 1:n.
column_orders <- function(surr, segment_size) {
  k <- nrow(surr) %/% segment_size
  t(apply(surr, 2, source_segments, segment_size = segment_size, k = k))
}

expect_segment_shuffle <- function(surr, n, segment_size) {
  k <- floor(n / segment_size)
  expect_true(is.matrix(surr))
  expect_true(is.numeric(surr))
  expect_identical(nrow(surr), as.integer(n))
  for (j in seq_len(ncol(surr))) {
    q <- source_segments(surr[, j], segment_size, k)
    expect_false(anyNA(q))
    expect_identical(sort(q), as.double(seq_len(k)))
    expect_false(all(q == seq_len(k)))
    if (k * segment_size < n) {
      tail_rows <- (k * segment_size + 1):n
      expect_identical(surr[tail_rows, j], as.double(tail_rows))
    }
  }
  expect_identical(nrow(unique(t(surr))), ncol(surr))
}

# --- AC1: segment structure ------------------------------------------------

test_that("segment surrogates reorder whole segments, never the original", {
  # 60 = 6 * 10 (no tail); 64 = 6 * 10 + 4 (tail of 4); 97 = 12 * 8 + 1
  # (12 segments, the draw-and-reject branch).
  for (case in list(c(60, 10), c(64, 10), c(97, 8))) {
    n <- case[1]
    s <- case[2]
    set.seed(n)
    surr <- suppressMessages(
      generate_surrogate_segment(as.double(seq_len(n)), s, n_surrogates = 9)
    )
    expect_identical(ncol(surr), 9L)
    expect_segment_shuffle(surr, n, s)
  }
})

test_that("a segment can stay in place", {
  # Uniform non-identity orders leave about one of k segments in place on
  # average, so some of 50 columns of a 6-segment series keep a segment.
  set.seed(16)
  surr <- generate_surrogate_segment(as.double(1:60), 10, n_surrogates = 50)
  orders <- column_orders(surr, 10)
  expect_true(any(orders == matrix(1:6, 50, 6, byrow = TRUE)))
})

test_that("the segment cut matches the eps loop of SUSY's susy()", {
  # SUSY: numberEpochen = round(size/range - 0.499999) and the eps loop puts
  # rows (i0 - 1) * range + 1 to i0 * range in segment i0. For size = 23,
  # range = 5: round(4.6 - 0.499999) = round(4.100001) = 4 segments,
  # rows 1-5, 6-10, 11-15, 16-20; rows 21-23 are the tail (dropped by SUSY,
  # kept in place here).
  boundaries <- list(1:5, 6:10, 11:15, 16:20)
  y <- as.double(1:23)
  set.seed(7)
  surr <- suppressMessages(generate_surrogate_segment(y, 5, n_surrogates = 9))
  for (j in seq_len(ncol(surr))) {
    blocks <- lapply(boundaries, function(rows) surr[rows, j])
    sources <- vapply(blocks, function(b) {
      which(vapply(boundaries, function(rows) identical(b, y[rows]), TRUE))
    }, integer(1))
    expect_setequal(sources, 1:4)
    expect_false(all(sources == 1:4))
    expect_identical(surr[21:23, j], c(21, 22, 23))
  }
})

test_that("NA values move with their segments", {
  n <- 64
  y_id <- as.double(seq_len(n))
  y_na <- stats::rnorm(n)
  y_na[c(3, 4, 17, 40, 62)] <- NA
  set.seed(11)
  surr_id <- suppressMessages(generate_surrogate_segment(y_id, 10, 6))
  set.seed(11)
  surr_na <- suppressMessages(generate_surrogate_segment(y_na, 10, 6))
  expect_segment_shuffle(surr_id, n, 10)
  for (j in seq_len(ncol(surr_id))) {
    expect_identical(surr_na[, j], y_na[surr_id[, j]])
  }
  expect_identical(colSums(is.na(surr_na)), rep(5, 6))
})

test_that("segment surrogates accept integer input and return doubles", {
  set.seed(14)
  surr <- generate_surrogate_segment(1:40, 10, n_surrogates = 2)
  expect_type(surr, "double")
})

# --- AC2: closed-form oracle -----------------------------------------------

test_that("segment surrogates match the plain reference", {
  cases <- list(
    list(seed = 201, n = 60, s = 10, m = 20), # 6 segments, no tail
    list(seed = 202, n = 67, s = 8, m = 30), # 8 segments, tail of 3
    list(seed = 203, n = 101, s = 10, m = 25) # 10 segments, tail of 1
  )
  for (case in cases) {
    set.seed(case$seed)
    y <- stats::rnorm(case$n)
    set.seed(case$seed + 1000)
    want <- segment_reference(y, case$s, case$m)
    set.seed(case$seed + 1000)
    got <- suppressMessages(generate_surrogate_segment(y, case$s, case$m))
    expect_equal(got, want)
  }
})

test_that("the plain reference lists k! - 1 orders", {
  counts <- vapply(2:5, function(k) nrow(all_orders_plain(k)), integer(1))
  expect_identical(counts, c(1L, 5L, 23L, 119L))
})

# --- AC3: uniform orders -----------------------------------------------------

order_key <- function(o) paste(o, collapse = "")

test_that("with 4 segments, all 23 non-identity orders are returned", {
  surr <- suppressMessages(
    generate_surrogate_segment(as.double(1:40), 10, n_surrogates = 30)
  )
  got <- sort(apply(column_orders(surr, 10), 1, order_key))
  want <- sort(apply(all_orders_plain(4), 1, order_key))
  expect_identical(got, want)
  expect_length(want, 23)
})

test_that("single orders and pairs of orders are equally frequent", {
  y <- as.double(1:40)
  levels_1 <- apply(all_orders_plain(4), 1, order_key)

  set.seed(301)
  one <- vapply(seq_len(2300), function(i) {
    order_key(column_orders(generate_surrogate_segment(y, 10, 1), 10)[1, ])
  }, character(1))
  counts_1 <- table(factor(one, levels = levels_1))
  expect_identical(sum(counts_1), 2300L)
  expect_gte(stats::chisq.test(counts_1)$p.value, 0.001)

  pair_levels <- utils::combn(sort(levels_1), 2, paste, collapse = "|")
  set.seed(302)
  two <- vapply(seq_len(2530), function(i) {
    keys <- apply(column_orders(generate_surrogate_segment(y, 10, 2), 10), 1, order_key)
    paste(sort(keys), collapse = "|")
  }, character(1))
  counts_2 <- table(factor(two, levels = pair_levels))
  expect_length(counts_2, 253)
  expect_identical(sum(counts_2), 2530L)
  expect_gte(stats::chisq.test(counts_2)$p.value, 0.001)
})

test_that("with 9 segments, the first-placed segment is equally frequent", {
  # Excluding the original order changes P(segment 1 first) from 1/9 to
  # (8! - 1) / (9! - 1), a relative difference of 2.5e-5.
  y <- as.double(1:18)
  set.seed(303)
  first <- vapply(seq_len(2000), function(i) {
    column_orders(generate_surrogate_segment(y, 2, 1), 2)[1, 1]
  }, numeric(1))
  counts <- table(factor(first, levels = 1:9))
  expect_identical(sum(counts), 2000L)
  expect_gte(stats::chisq.test(counts)$p.value, 0.001)
})

# --- AC5: messages and errors ----------------------------------------------

test_that("the tail message names the tail length", {
  expect_message(
    generate_surrogate_segment(as.double(1:64), 10, n_surrogates = 3),
    "last 4 samples",
    class = "bsync_segment_tail"
  )
  expect_no_message(
    generate_surrogate_segment(as.double(1:60), 10, n_surrogates = 3),
    class = "bsync_segment_tail"
  )
})

test_that("3 segments return all 5 orders and warn; 4 segments do not", {
  set.seed(15)
  expect_warning(
    expect_message(
      surr <- generate_surrogate_segment(as.double(1:30), 10, n_surrogates = 19),
      "all 5 possible",
      class = "bsync_segment_all_orders"
    ),
    "no test .* can reach p <= .05",
    class = "bsync_segment_few_orders"
  )
  expect_identical(ncol(surr), 5L)
  expect_segment_shuffle(surr, 30, 10)

  set.seed(15)
  expect_no_warning(
    surr <- generate_surrogate_segment(as.double(1:40), 10, n_surrogates = 19)
  )
  expect_identical(ncol(surr), 19L)
  expect_segment_shuffle(surr, 40, 10)
})

test_that("generate_surrogate_segment() aborts on invalid input", {
  y <- as.double(1:40)
  expect_error(
    generate_surrogate_segment(letters, 2),
    "`y` must be a numeric vector"
  )
  expect_error(
    generate_surrogate_segment(matrix(y, 20), 2),
    "`y` must be a numeric vector"
  )
  for (bad in list(0, -3, 2.5, NA, c(5, 10), "5")) {
    expect_error(
      generate_surrogate_segment(y, 10, n_surrogates = bad),
      "`n_surrogates` must be a single positive integer"
    )
  }
  for (bad in list(1, 0, 2.5, NA, c(5, 10), "5")) {
    expect_error(
      generate_surrogate_segment(y, bad),
      "`segment_size` must be a single whole number of at least 2"
    )
  }
  expect_error(
    generate_surrogate_segment(y, 21),
    "`segment_size` \\(21\\) must be at most `floor\\(length\\(y\\) / 2\\)` \\(20\\)"
  )
  expect_error(
    generate_surrogate_segment(as.double(1:3), 2),
    "must be at most `floor\\(length\\(y\\) / 2\\)` \\(1\\)"
  )
  # 9 segments have 9! - 1 = 362879 orders other than the original.
  expect_error(
    generate_surrogate_segment(as.double(1:18), 2, n_surrogates = 362879),
    "`n_surrogates` \\(362879\\) must be below the 362879 possible orders"
  )
})
