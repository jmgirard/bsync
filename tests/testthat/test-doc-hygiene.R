# User-facing documentation must not reference internal milestone numbers (M<n>).
# The milestone log lives in MILESTONES.md; user-facing docs should describe
# features in plain terms instead. Internal `#` code comments in R/ may reference
# milestones for development provenance, so R/ is deliberately NOT scanned here.

test_that("user-facing docs contain no milestone-number references", {
  skip_on_cran()

  root <- testthat::test_path("..", "..")
  files <- c(
    file.path(root, "NEWS.md"),
    file.path(root, "README.md"),
    list.files(file.path(root, "vignettes"), pattern = "\\.Rmd$", full.names = TRUE),
    list.files(file.path(root, "man"), pattern = "\\.Rd$", full.names = TRUE)
  )
  files <- files[file.exists(files)]
  skip_if(length(files) == 0L, "no user-facing docs found in this build")

  pattern <- "\\bM[0-9]+\\b"
  hits <- character()
  for (f in files) {
    lines <- readLines(f, warn = FALSE)
    idx <- grep(pattern, lines)
    if (length(idx) > 0L) {
      hits <- c(hits, sprintf("%s:%d: %s", basename(f), idx, trimws(lines[idx])))
    }
  }

  if (length(hits) > 0L) {
    fail(paste0(
      "Milestone-number references found in user-facing docs. Describe the ",
      "feature in plain terms; the milestone log lives in MILESTONES.md.\n",
      paste(hits, collapse = "\n")
    ))
  } else {
    succeed()
  }
})
