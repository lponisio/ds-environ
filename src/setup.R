# Preflight for the geospatial lessons. Data is provided separately.
# This file is sourced by chapters 07 and 08 during rendering.
# Do not install packages or download datasets here.
needed <- c(
  "data/NEON-DS-Site-Layout-Files/HARV/HarClip_UTMZ18.shp",
  "data/NEON-DS-Site-Layout-Files/HARV/HARVtower_UTM18N.shp",
  "data/NEON-DS-Airborne-Remote-Sensing/HARV/DSM/HARV_dsmCrop.tif",
  "data/NEON-DS-Airborne-Remote-Sensing/HARV/CHM/HARV_chmCrop.tif"
)
missing <- needed[!file.exists(needed)]
if (length(missing)) {
  stop("Geospatial lesson data are missing. Restore the original top-level data/ ",
       "folder (including its NEON subdirectories) before building. Missing:\n - ",
       paste(missing, collapse = "\n - "), call. = FALSE)
}
rm(needed, missing)
