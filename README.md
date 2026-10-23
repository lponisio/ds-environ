# Data Science for the Environment — reproducible bookdown sources

This is the editable, **student-facing source** for the University of Oregon
Data Science for the Environment course. It includes the original chapter
`.Rmd` documents, embedded demos/labs/readings, figures, bibliography, and
necessary bookdown configuration. It intentionally **does not include**
previous compiled pages, caches, the instructor's `.git` history, or the
large top-level `data/` folder.

## First setup (once per computer)

1. Install a current **R** and **RStudio**. RStudio includes Pandoc.
2. Clone or pull the course repository. The instructor will restore the
   **top-level `data/`** folder alongside `index.Rmd` before distributing
   the complete course (or supply it separately if it is too large for Git).
   Keep the subdirectory structure and original capitalizations.
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

The HTML pages are regenerated locally. You do **not** need to commit or pull
`docs/`, `_book/`, `_bookdown_files/`, `*.rds` build state, or HTML. `data/`
must also be available at the root: it can be tracked in Git, or distributed
separately if it is too large for a normal repository.
When a new R package is added to a chapter, run
`source("scripts/install_packages.R")` again; update its package vector as well.

## Source files that must remain in Git

- `index.Rmd` and `01-` through `10-*.Rmd`: chapters
- `_bookdown.yml`, `_output.yml`: bookdown configuration (**these are NOT outputs**)
- `_common.R`, `src/setup.R`, `scripts/`: R infrastructure
- `demos/`, `labs/`, `readings/`: embedded teaching content and support files
- `fig/`, `refs/`: teaching visuals and bibliography
- `.gitignore`, `.Rproj`, `README.md`, `LICENSE`

The package preserves the complete `labs/data/` and `demos/data/` directories
because they belong to the embedded lessons. They are **not** the separately
supplied top-level `data/` directory.

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

## Instructor: remove previously tracked build artifacts

`.gitignore` only prevents *new* untracked files being added. If you already
committed generated files to Git, remove them **from the Git index only** (the
files can stay on your computer):

```bash
git rm -r --cached --ignore-unmatch docs _book _bookdown_files
git rm --cached --ignore-unmatch data-science-for-environment.rds
# If the large data folder is distributed separately and should remain untracked:
# git rm -r --cached --ignore-unmatch data
git add .gitignore README.md *.Rmd _bookdown.yml _output.yml _common.R   scripts src demos labs readings fig refs data-science-for-environment.Rproj LICENSE
git commit -m "Make bookdown sources independently buildable"
git push
```

If you are adding the omitted data folder to Git, explicitly `git add data/`
(after checking GitHub file-size limits) before committing.

Keep any current `docs/` publishing setup in mind: if GitHub Pages is configured
to publish from tracked `docs/`, ignoring `docs/` will stop updates to that
published site. Use a GitHub Actions build/deployment workflow if you also want
an automatically published website without storing build artifacts in source.
