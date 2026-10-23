# Source/build audit — October 2026

## Confirmed source issues addressed

- `index.Rmd` called ten `install.packages()` operations during each render
  using obsolete HTTP CRAN URLs. Removed. Dependency installation is now a
  separate opt-in step in `scripts/install_packages.R` using HTTPS.
- A nested programming demo (`demos/01-demo-fundamentals.Rmd`) attempted to
  install `lme4` while compiling. The installation demonstration is now
  displayed without being evaluated; students can run it manually.
- The nested first programming demo invoked `fahr_to_celsius()` after defining
  `convert_fahr_to_celsius()`. Calls have been corrected.
- That demo also called `ggplot` without attaching the package and handed it a
  matrix. Corrected to an explicit `ggplot2` call on a data frame.
- `04-simulation-hypothesis-testing.Rmd` linked to image filenames with
  incorrect capitalization (`Fdiv_high`, `Fdiv_low`, `Fdis_low`). Corrected
  to actual files `FDiv_high`, `FDiv_low`, `FDis_low`.
- The optional blank population-growth demo linked to a nonexistent
  `demos/figures/r0.png`. It now links to the image in `fig/r0.png` from the
  `demos/` directory.
- The chapter 03 CO2 setup used the omitted top-level `data/` even though
  the same climate data is bundled in `labs/data/`. It now uses the bundled
  CO2 file.
- A lab setup guarded `library(bipartite)` behind a check for the unrelated
  `networks3D` package. Removed the mistaken check.
- `src/setup.R` would automatically fetch large NEON data archives and
  Natural Earth files as a side effect of knitting chapters 07–08.
  It now checks required local inputs and never downloads on compilation.
- The original `.gitignore` did not exclude compiled outputs. It now ignores
  regenerable HTML, intermediate build directories, and local RStudio state
  while keeping bookdown YAML source files tracked.
- `README.md` incorrectly instructed students to open nonexistent `index.R`
  rather than the bookdown project/`index.Rmd`. Corrected.

## Input/layout checks performed here

- Parsed `_bookdown.yml` successfully and confirmed all **11** included
  chapter source paths exist.
- Confirmed every listed `child=` path in the main chapters resolves to an
  included `demos/`, `labs/`, or `readings/` document.
- Confirmed all non-placeholder Markdown image references resolve in the
  provided source package after the capitalization fixes (allowing project-
  root and child-local paths).
- Confirmed the top-level `data/` folder is deliberately absent and that
  `labs/data/` and `demos/data/` are present in the source package.
- Confirmed all package-install calls that would execute during the main
  book render have been removed.
- Kept `_bookdown.yml`, `_output.yml`, `_common.R`, bibliography, scripts,
  figures, local teaching datasets, and standalone teaching documents.

## What was NOT tested

**A real R/bookdown render was not possible in the preparation environment
because no R runtime is installed and package installation was unavailable.**
The top-level `data/` folder was not present in the upload, so its exact
files/contents cannot be verified. As a result, this package has passed
static source/reference checks, *not an end-to-end successful knit*.

Once the instructor supplies the original top-level `data/` folder, run:

```r
source("scripts/install_packages.R")
source("scripts/build_book.R")
```

If `check_inputs.R` identifies any additional missing paths, check the
folder structure/names in the original data archive. Some pedagogical code
in the original lessons intentionally throws errors (notably debugging).

## Git publishing caveat

The new `.gitignore` ignores the `docs/` build output; it does not untrack
files already committed. Use `git rm --cached` as shown in `README.md`.
If the existing public book is served from GitHub Pages `docs/`, switch to a
separate build/deployment action before removing tracked `docs/` output.
