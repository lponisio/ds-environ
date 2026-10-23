# Run once from the project root in the RStudio console:
#   source("scripts/install_packages.R")
if (!file.exists("index.Rmd") || !file.exists("_bookdown.yml")) {
  stop("Open the book project in RStudio (or set working directory to its root).",
       call. = FALSE)
}
# HTML only: LaTeX/TinyTeX is NOT required.
book_packages <- c(
  "bookdown", "rmarkdown", "knitr", "tidyverse", "igraph", "bipartite",
  "networkD3", "vegan", "terra", "sf", "ggspatial", "RcppRoll",
  "RColorBrewer", "maps"
)
missing <- book_packages[!vapply(book_packages, requireNamespace,
                                   logical(1), quietly = TRUE)]
if (length(missing)) {
  message("Installing missing R packages: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
} else {
  message("All required book packages are already installed.")
}
still_missing <- book_packages[!vapply(book_packages, requireNamespace,
                                         logical(1), quietly = TRUE)]
if (length(still_missing)) {
  stop("Some packages did not install: ", paste(still_missing, collapse = ", "),
       ". For sf/terra, see platform-specific system requirements.", call. = FALSE)
}
message("Package setup complete. Run source('scripts/build_book.R').")
