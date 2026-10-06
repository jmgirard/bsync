# Extracted from test-surrogate-pseudo.R:346

# prequel ----------------------------------------------------------------------
library(testthat)
make_dyads <- function(lengths, forms = rep("named", length(lengths))) {
  Map(
    function(n, form) {
      x <- stats::rnorm(n)
      y <- stats::rnorm(n)
      switch(
        form,
        df = data.frame(x = x, y = y),
        named = list(x = x, y = y),
        rev = list(y = y, x = x),
        unnamed = list(x, y)
      )
    },
    lengths,
    forms
  )
}
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
identify_source <- function(col, pool) {
  n <- length(col)
  hits <- Filter(
    function(s) {
      length(s$values) >= n && identical(as.double(s$values[seq_len(n)]), col)
    },
    pool
  )
  vapply(hits, function(s) paste0(s$dyad, s$role), character(1))
}
column_sources <- function(mat, dyad_list) {
  pool <- source_pool(dyad_list)
  vapply(
    seq_len(ncol(mat)),
    function(j) {
      src <- identify_source(mat[, j], pool)
      expect_length(src, 1)
      src[1]
    },
    character(1)
  )
}
target_y <- function(dyad_list, dyad) source_pool(dyad_list)[[2 * dyad]]$values
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
    mapply(
      function(a, b) paste(sort(c(a, b)), collapse = "|"),
      xl,
      yl,
      USE.NAMES = FALSE
    )
  }
}
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

# test -------------------------------------------------------------------------
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
      expect_named(
        s,
        c("x_dyad", "x_role", "x_name", "y_dyad", "y_role", "y_name")
      )
      expect_true(all(s$x_dyad != s$y_dyad))
      if (kr) {
        expect_true(all(s$x_role == "x" & s$y_role == "y"))
      }
      expect_samples_match_sources(pd, dl)
    }
  }
