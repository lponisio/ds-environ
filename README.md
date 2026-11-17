# Data Science for the Environment — reproducible bookdown sources

This is the editable, **student-facing source** for the University of Oregon
Data Science for the Environment course. It includes the original chapter
`.Rmd` documents, embedded demos/labs/readings, figures, bibliography, and
necessary bookdown configuration. 

## First setup (once per computer)

1. Install a current **R** and **RStudio**. RStudio includes Pandoc.
2. Clone or pull the course repository. 
3. In RStudio, open `data-science-for-environment.Rproj` (or set the R working
   directory to the folder with `index.Rmd`).
4. In the RStudio **Console**, run:

   ```r
   source("scripts/install_packages.R")
   source("scripts/build_book.R")
   ```

5. Open `docs/index.html` after the build completes. No PDF/TinyTeX is needed.

The first installation can take time, especially spatial packages `sf` and
`terra`. On some operating systems they may require platform tools/system
libraries; installing precompiled binaries when available is easiest.

## Every subsequent update

In RStudio, pull changes using the Git pane (or run `git pull` in a terminal),
then re-run:

```r
source("scripts/build_book.R")
```

The HTML pages are regenerated locally. 

## Source files that must remain in Git

- `index.Rmd` and `01-` through `10-*.Rmd`: chapters
- `_bookdown.yml`, `_output.yml`: bookdown configuration (**these are NOT outputs**)
- `_common.R`, `src/setup.R`, `scripts/`: R infrastructure
- `demos/`, `labs/`, `readings/`: embedded teaching content and support files
- `fig/`, `refs/`: teaching visuals and bibliography
- `.gitignore`, `.Rproj`, `README.md`, `LICENSE`


## File paths and rendering

Build from the **book root** with the provided script. Chapter paths are
relative to the book root; code in embedded labs and demos can use paths
relative to their own directories. Do not flatten these directories or move
`.Rmd` files around: that will break relative paths and images. E.g. the
embedded OBA network demo uses `../data/OBA_2018-2024.csv` relative to the
`demos/` folder. Missing input paths are reported by `check_inputs.R`.

The existing `new_session: true` bookdown option intentionally renders each
chapter in a separate R session; teaching examples can reuse object names and
chunk labels without colliding across chapters.

**Note:** Some pedagogical code in the original course intentionally produces
errors (particularly the debugging material). The book still displays those
examples. The preflight script catches missing source/data and packages before
rendering so missing inputs do not get silently buried in chapter output.
