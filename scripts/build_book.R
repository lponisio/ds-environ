# Run from the project root (RStudio console or Rscript):
#   source("scripts/build_book.R")
#   Rscript scripts/build_book.R
if (!file.exists("index.Rmd") || !file.exists("_bookdown.yml")) {
  stop("Open the .Rproj from the book's root directory before compiling.", call. = FALSE)
}
source("scripts/check_inputs.R", local = TRUE)
required <- c("bookdown", "rmarkdown", "knitr", "tidyverse", "igraph",
              "bipartite", "networkD3", "vegan", "terra", "sf", "ggspatial",
              "RcppRoll", "RColorBrewer", "maps")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) {
  stop("Install the required packages first with source('scripts/install_packages.R'). ",
       "Missing: ", paste(missing, collapse = ", "), call. = FALSE)
}
message("Rendering a fresh HTML book (no previous compiled files required)...")
bookdown::render_book(input = "index.Rmd", output_format = "bookdown::gitbook",
                      clean = TRUE)
# Explicitly copy static teaching resources referenced by the child Rmd files.
# Markdown images are normally copied by bookdown, but this also ensures
# that the embedded PDF is present when serving the finished docs/ website.
unit1_resources <- c(
  "readings/reading-Alston2020.pdf",
  "fig/unit1-data-structures.png",
  "fig/unit1-logistic-trajectories.png",
  "fig/unit1-reproduction-and-growth.png"
)
for (resource in unit1_resources) {
  if (!file.exists(resource)) {
    warning("Unit 1 resource not found: ", resource)
    next
  }
  target <- file.path("docs", resource)
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(resource, target, overwrite = TRUE)) {
    warning("Could not copy unit 1 resource: ", resource)
  }
}
message("Book built: open docs/index.html in your browser.")
