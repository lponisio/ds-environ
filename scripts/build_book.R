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
message("Book built: open docs/index.html in your browser.")
