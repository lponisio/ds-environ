# Check the minimal data/source inputs before knitting.
# Called by scripts/build_book.R, and safe to source independently.
if (!file.exists("index.Rmd") || !file.exists("_bookdown.yml")) {
  stop("Working directory must be the book root, containing index.Rmd and _bookdown.yml.",
       call. = FALSE)
}
needed_files <- c(
  "_output.yml", "_common.R", "src/setup.R",
  "refs/ecol_primer.bib",
  "fig/funcDiv/FDiv_high.png", "fig/funcDiv/FDiv_low.png",
  "fig/funcDiv/Fdis_high.png", "fig/funcDiv/FDis_low.png",
  "demos/data/co2_mm_mlo.txt",
  "labs/data/co2_mm_mlo.txt",
  "labs/data/specimens-complete.csv",
  "labs/data/egypt_data.txt", "labs/data/egypt_mass.txt",
  "data/OBA_2018-2024.csv",
  "data/NEON-DS-Site-Layout-Files/HARV/HARVtower_UTM18N.shp",
  "data/NEON-DS-Site-Layout-Files/HARV/HARV_roads.shp",
  "data/NEON-DS-Site-Layout-Files/HARV/HarClip_UTMZ18.shp",
  "data/NEON-DS-Site-Layout-Files/HARV/PlotLocations_HARV.shp",
  "data/NEON-DS-Site-Layout-Files/HARV/HARV_PlotLocations.csv",
  "data/NEON-DS-Site-Layout-Files/HARV/HARV_2NewPhenPlots.csv",
  "data/NEON-DS-Site-Layout-Files/US-Boundary-Layers/US-Boundary-Dissolved-States.shp",
  "data/NEON-DS-Airborne-Remote-Sensing/HARV/DSM/HARV_dsmCrop.tif",
  "data/NEON-DS-Airborne-Remote-Sensing/HARV/CHM/HARV_chmCrop.tif"
)
missing <- needed_files[!file.exists(needed_files)]
if (length(missing)) {
  stop("Missing input files (most should arrive in the instructor's separately supplied ",
       "top-level data/ folder):\n  - ", paste(missing, collapse = "\n  - "),
       "\nDo not copy docs/, _book/ or _bookdown_files/; those are generated.",
       call. = FALSE)
}
message("Source and required data files are present.")
