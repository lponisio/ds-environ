




# Geospatial data: data types (raster and vector)

###    Lecture summary

####    Learning Objectives

- Describe the fundamental attributes of a raster dataset.
- Explore raster attributes and metadata using R.
- Import rasters into R using the `terra` package.
- Plot a raster file in R using the `ggplot2` package.
- Describe the difference between single- and multi-band rasters.
- Build customized plots for a single band raster using the `ggplot2` package.
- Layer a raster dataset on top of a hillshade to create an elegant basemap.
- Reproject a raster

####    Introduction to geospatial data

**Summary**

- Describe the fundamental attributes of a raster dataset.
- Explore raster attributes and metadata using R.
- Import rasters into R using the `terra` package.
- Plot a raster file in R using the `ggplot2` package.
- Describe the difference between single- and multi-band rasters.

We will introduce the fundamental principles, packages and
metadata/raster attributes that are needed to work with raster data in R. We will discuss some of the core metadata elements that we need to understand to work with rasters in R, including CRS and resolution. We will also explore missing and bad data values as stored in a raster and how R handles these elements.

**About Raster Data**

Raster data is any pixelated (or gridded) data where each pixel is associated with a specific geographical location. The value of a pixel can be continuous (e.g. elevation) or categorical (e.g. land use). If this sounds familiar, it is because this data structure is very common: it's how we represent any digital image. A geospatial raster is only different
from a digital photo in that it is accompanied by spatial information
that connects the data to a particular location. This includes the
raster's extent and cell size, the number of rows and columns, and
its coordinate reference system (or CRS).

![](fig/dc-spatial-raster/raster_concept.png){alt='Raster Concept'}

**Continuous rasters**
Some examples of continuous rasters include:

1. Precipitation maps.
2. Maps of tree height derived from LiDAR data.
3. Elevation values for a region.

**Categorial rasters**
Some rasters contain categorical data where each pixel represents a discrete class such as a landcover type (e.g., "forest" or "grassland") rather than a continuous value such as elevation or temperature. Some examples of classified maps include:

1. Landcover / land-use maps.
2. Tree height maps classified as short, medium, and tall trees.
3. Elevation maps classified as low, medium, and high elevation.

![](fig/USA_landcover_classification.png){alt='USA landcover classification'}

#### Reading in raster data and understanding its attributes 

**The Data**
The data  was collected by [National Ecological Observatory Network (NEON)](https://www.neonscience.org/). 
The data uses data from Harvard Forest (HARV) - Massachusetts, USA - [fieldsite description](https://www.neonscience.org/field-sites/harv)

**View Raster File Attributes**

We will be working with a series of GeoTIFF files. The
GeoTIFF format contains a set of embedded tags with metadata about the raster data. We can use the function `describe()` to get information about our raster data before we read that data into R. It is ideal to do this before importing your data.

We will be using the digital elevation model (DSM) data from the forest


``` r
describe("data/NEON-DS-Airborne-Remote-Sensing/HARV/DSM/HARV_dsmCrop.tif")
```

```
##  [1] "Driver: GTiff/GeoTIFF"                                                                                                                                                                                                                                                          
##  [2] "Files: data/NEON-DS-Airborne-Remote-Sensing/HARV/DSM/HARV_dsmCrop.tif"                                                                                                                                                                                                          
##  [3] "Size is 1697, 1367"                                                                                                                                                                                                                                                             
##  [4] "Coordinate System is:"                                                                                                                                                                                                                                                          
##  [5] "PROJCRS[\"WGS 84 / UTM zone 18N\","                                                                                                                                                                                                                                             
##  [6] "    BASEGEOGCRS[\"WGS 84\","                                                                                                                                                                                                                                                    
##  [7] "        DATUM[\"World Geodetic System 1984\","                                                                                                                                                                                                                                  
##  [8] "            ELLIPSOID[\"WGS 84\",6378137,298.257223563,"                                                                                                                                                                                                                        
##  [9] "                LENGTHUNIT[\"metre\",1]]],"                                                                                                                                                                                                                                     
## [10] "        PRIMEM[\"Greenwich\",0,"                                                                                                                                                                                                                                                
## [11] "            ANGLEUNIT[\"degree\",0.0174532925199433]],"                                                                                                                                                                                                                         
## [12] "        ID[\"EPSG\",4326]],"                                                                                                                                                                                                                                                    
## [13] "    CONVERSION[\"UTM zone 18N\","                                                                                                                                                                                                                                               
## [14] "        METHOD[\"Transverse Mercator\","                                                                                                                                                                                                                                        
## [15] "            ID[\"EPSG\",9807]],"                                                                                                                                                                                                                                                
## [16] "        PARAMETER[\"Latitude of natural origin\",0,"                                                                                                                                                                                                                            
## [17] "            ANGLEUNIT[\"degree\",0.0174532925199433],"                                                                                                                                                                                                                          
## [18] "            ID[\"EPSG\",8801]],"                                                                                                                                                                                                                                                
## [19] "        PARAMETER[\"Longitude of natural origin\",-75,"                                                                                                                                                                                                                         
## [20] "            ANGLEUNIT[\"degree\",0.0174532925199433],"                                                                                                                                                                                                                          
## [21] "            ID[\"EPSG\",8802]],"                                                                                                                                                                                                                                                
## [22] "        PARAMETER[\"Scale factor at natural origin\",0.9996,"                                                                                                                                                                                                                   
## [23] "            SCALEUNIT[\"unity\",1],"                                                                                                                                                                                                                                            
## [24] "            ID[\"EPSG\",8805]],"                                                                                                                                                                                                                                                
## [25] "        PARAMETER[\"False easting\",500000,"                                                                                                                                                                                                                                    
## [26] "            LENGTHUNIT[\"metre\",1],"                                                                                                                                                                                                                                           
## [27] "            ID[\"EPSG\",8806]],"                                                                                                                                                                                                                                                
## [28] "        PARAMETER[\"False northing\",0,"                                                                                                                                                                                                                                        
## [29] "            LENGTHUNIT[\"metre\",1],"                                                                                                                                                                                                                                           
## [30] "            ID[\"EPSG\",8807]]],"                                                                                                                                                                                                                                               
## [31] "    CS[Cartesian,2],"                                                                                                                                                                                                                                                           
## [32] "        AXIS[\"(E)\",east,"                                                                                                                                                                                                                                                     
## [33] "            ORDER[1],"                                                                                                                                                                                                                                                          
## [34] "            LENGTHUNIT[\"metre\",1]],"                                                                                                                                                                                                                                          
## [35] "        AXIS[\"(N)\",north,"                                                                                                                                                                                                                                                    
## [36] "            ORDER[2],"                                                                                                                                                                                                                                                          
## [37] "            LENGTHUNIT[\"metre\",1]],"                                                                                                                                                                                                                                          
## [38] "    USAGE["                                                                                                                                                                                                                                                                     
## [39] "        SCOPE[\"Navigation and medium accuracy spatial referencing.\"],"                                                                                                                                                                                                        
## [40] "        AREA[\"Between 78°W and 72°W, northern hemisphere between equator and 84°N, onshore and offshore. Bahamas. Canada - Nunavut; Ontario; Quebec. Colombia. Cuba. Ecuador. Greenland. Haiti. Jamaica. Panama. Turks and Caicos Islands. United States (USA). Venezuela.\"],"
## [41] "        BBOX[0,-78,84,-72]],"                                                                                                                                                                                                                                                   
## [42] "    ID[\"EPSG\",32618]]"                                                                                                                                                                                                                                                        
## [43] "Data axis to CRS axis mapping: 1,2"                                                                                                                                                                                                                                             
## [44] "Origin = (731453.000000000000000,4713838.000000000000000)"                                                                                                                                                                                                                      
## [45] "Pixel Size = (1.000000000000000,-1.000000000000000)"                                                                                                                                                                                                                            
## [46] "Metadata:"                                                                                                                                                                                                                                                                      
## [47] "  AREA_OR_POINT=Area"                                                                                                                                                                                                                                                           
## [48] "Image Structure Metadata:"                                                                                                                                                                                                                                                      
## [49] "  COMPRESSION=LZW"                                                                                                                                                                                                                                                              
## [50] "  INTERLEAVE=BAND"                                                                                                                                                                                                                                                              
## [51] "Corner Coordinates:"                                                                                                                                                                                                                                                            
## [52] "Upper Left  (  731453.000, 4713838.000) ( 72d10'52.71\"W, 42d32'32.18\"N)"                                                                                                                                                                                                      
## [53] "Lower Left  (  731453.000, 4712471.000) ( 72d10'54.71\"W, 42d31'47.92\"N)"                                                                                                                                                                                                      
## [54] "Upper Right (  733150.000, 4713838.000) ( 72d 9'38.40\"W, 42d32'30.35\"N)"                                                                                                                                                                                                      
## [55] "Lower Right (  733150.000, 4712471.000) ( 72d 9'40.41\"W, 42d31'46.08\"N)"                                                                                                                                                                                                      
## [56] "Center      (  732301.500, 4713154.500) ( 72d10'16.56\"W, 42d32' 9.13\"N)"                                                                                                                                                                                                      
## [57] "Band 1 Block=1697x1 Type=Float64, ColorInterp=Gray"                                                                                                                                                                                                                             
## [58] "  Min=305.070 Max=416.070 "                                                                                                                                                                                                                                                     
## [59] "  Minimum=305.070, Maximum=416.070, Mean=359.853, StdDev=17.832"                                                                                                                                                                                                                
## [60] "  NoData Value=-9999"                                                                                                                                                                                                                                                           
## [61] "  Metadata:"                                                                                                                                                                                                                                                                    
## [62] "    STATISTICS_MAXIMUM=416.06997680664"                                                                                                                                                                                                                                         
## [63] "    STATISTICS_MEAN=359.85311802914"                                                                                                                                                                                                                                            
## [64] "    STATISTICS_MINIMUM=305.07000732422"                                                                                                                                                                                                                                         
## [65] "    STATISTICS_STDDEV=17.83169335933"
```

If you wish to store this information in R, you can do the following:


``` r
HARV_dsmCrop_info <- capture.output(
  describe("data/NEON-DS-Airborne-Remote-Sensing/HARV/DSM/HARV_dsmCrop.tif")
)
```

Each line of text that was printed to the console is now stored as an element of
the character vector `HARV_dsmCrop_info`. We will be exploring this data throughout this
demo. By the end of this demo, you will be able to explain and understand the output above.

**BP programming tip - Object names**

To improve code readability, file and object names should be used that make it clear what is in the file. The data for this demo were collected from Harvard Forest so we'll use a naming convention of `datatype_HARV`.

**Open a Raster in R**
Now that we've previewed the metadata for our GeoTIFF, let's import this
raster dataset into R and explore its metadata more closely. We can use the `rast()` function to open a raster in R.


``` r
DSM_HARV <-
  rast("data/NEON-DS-Airborne-Remote-Sensing/HARV/DSM/HARV_dsmCrop.tif")

DSM_HARV
```

```
## class       : SpatRaster 
## size        : 1367, 1697, 1  (nrow, ncol, nlyr)
## resolution  : 1, 1  (x, y)
## extent      : 731453, 733150, 4712471, 4713838  (xmin, xmax, ymin, ymax)
## coord. ref. : WGS 84 / UTM zone 18N (EPSG:32618) 
## source      : HARV_dsmCrop.tif 
## name        : HARV_dsmCrop 
## min value   :       305.07 
## max value   :       416.07
```

The information above includes a report of min and max values, but no other data range statistics. Similar to other R data structures like vectors and data frame columns, descriptive statistics for raster data can be retrieved like


``` r
summary(DSM_HARV)
```

```
## Warning: [summary] used a sample
```

```
##   HARV_dsmCrop  
##  Min.   :305.6  
##  1st Qu.:345.6  
##  Median :359.6  
##  Mean   :359.8  
##  3rd Qu.:374.3  
##  Max.   :414.7
```

but note the warning - unless you force R to calculate these statistics using every cell in the raster, it will take a random sample of 100,000 cells and calculate from that instead. To force calculation all the values, you can use the function `values`:


``` r
summary(values(DSM_HARV))
```

```
##   HARV_dsmCrop  
##  Min.   :305.1  
##  1st Qu.:345.6  
##  Median :359.7  
##  Mean   :359.9  
##  3rd Qu.:374.3  
##  Max.   :416.1
```

To visualise this data in R using `ggplot2`, we need to convert it to a
dataframe.  The `terra` package has an built-in function for conversion to a plot-able dataframe. We need to tell R that the data is geo-referenced using the argument xy=TRUE.


``` r
DSM_HARV_df <- as.data.frame(DSM_HARV, xy = TRUE)
```

Now when we view the structure of our data, we will see a standard
dataframe format.


``` r
str(DSM_HARV_df)
```

```
## 'data.frame':	2319799 obs. of  3 variables:
##  $ x           : num  731454 731454 731456 731456 731458 ...
##  $ y           : num  4713838 4713838 4713838 4713838 4713838 ...
##  $ HARV_dsmCrop: num  409 408 407 407 409 ...
```

We can use `ggplot()` to plot this data. 
- We will set the color scale to  `scale_fill_viridis_c` which is a color-blindness friendly color scale. 


``` r
ggplot() +
    geom_raster(data = DSM_HARV_df , aes(x = x, y = y, fill = HARV_dsmCrop)) +
    scale_fill_viridis_c() 
```

<div class="figure">
<img src="08-geospatial-I_files/figure-html/ggplot-raster-1.png" alt="Raster plot with ggplot2 using the viridis color scale" width="672" />
<p class="caption">(\#fig:ggplot-raster)Raster plot with ggplot2 using the viridis color scale</p>
</div>

**Plotting Tip**

More information about the Viridis palette used above at
[R Viridis package documentation](https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html).

**What do we see?** 
This map shows the digital surface model, so the elevation of the surface of our study site in Harvard Forest. From the legend, we can see that the maximum elevation is ~400, but we can't tell whether this is 400 feet or 400 meters because the legend doesn't show us the units. We
can look at the metadata of our object to see what the units are. Much of the metadata that we're interested in is part of the Coordinate Reference System (CRS). 

Now we will see how features of the CRS appear in our data file and what
meanings they have.

#### Raster Coordinate Reference System (CRS) 

We can view the CRS string associated with our R object using the`crs()`
function.


``` r
crs(DSM_HARV, proj = TRUE)
```

```
## [1] "+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs"
```


What units are our data in?
`+units=m` tells us that our data is in meters.

**Understanding CRS in Proj4 Format**

The CRS for our data is given to us by R in `proj4` format. Let's break down the pieces of `proj4` string. The string contains all of the individual CRS elements that R or another GIS might need. Each element is specified with a `+` sign, similar to how a `.csv` file is delimited or broken up by a `,`. After each `+` we see the CRS element being defined. For example projection (`proj=`) and datum (`datum=`).

**UTM Proj4 String**

A projection string (like the one of `DSM_HARV`) specifies the UTM projection  as follows:

`+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs'

- **proj=utm:** the projection is UTM, UTM has several zones.
- **zone=18:** the zone is 18
- **datum=WGS84:** the datum is WGS84 (the datum refers to the  0,0 reference for
  the coordinate system used in the projection)
- **units=m:** the units for the coordinates are in meters

Sometimes you will also see
+ellps=WGS84 +towgs84=0,0,0`

- **ellps=WGS84:** the ellipsoid (how the earth's  roundness is calculated) for
  the data is WGS84 (not in our crs above, but it can be.)
- **+towgs84=0,0,0**  parameter shifts (translation + rotation + scaling) realted to the datum

Note that the zone is unique to the UTM projection. Not all CRSs will have a zone. Image source: Chrismurf at English Wikipedia, via [Wikimedia Commons](https://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system#/media/File:Utm-zones-USA.svg) (CC-BY).


![The UTM zones across the continental United States. From: https://upload.wikimedia.org/wikipedia/commons/8/8d/Utm-zones-USA.svg](fig/Utm-zones-USA.svg){alt='UTM zones in the USA.'}

#### Extent

The spatial extent is the geographic area that the raster data covers.
The spatial extent of an R spatial object represents the geographic edge or
location that is the furthest north, south, east and west. In other words, extent
represents the overall geographic coverage of the spatial object.

![](fig/dc-spatial-raster/spatial_extent.png){alt='Spatial extent image'}

(Image Source: National Ecological Observatory Network (NEON))
{: .text-center}


``` r
ext(DSM_HARV)
```

```
## SpatExtent : 731453, 733150, 4712471, 4713838 (xmin, xmax, ymin, ymax)
```

#### Resolution

A resolution of a raster represents the area on the ground that each
pixel of the raster covers. The image below illustrates the effect
of changes in resolution.

![](fig/dc-spatial-raster/raster_resolution.png){alt='Resolution image'}

(Source: National Ecological Observatory Network (NEON))
{: .text-center}


``` r
res(DSM_HARV)
```

```
## [1] 1 1
```

#### Raster Bands

The Digital Surface Model object (`DSM_HARV`) that we've been working with is a single band raster. This means that there is only one dataset stored in the raster: surface elevation in meters for one time period.

![](fig/dc-spatial-raster/single_multi_raster.png){alt='Multi-band raster image'}

A raster dataset can contain one or more bands. We can use the `rast()`
function to import one single band from a single or multi-band raster. We can view the number of bands in a raster using the `nly()` function.


``` r
nlyr(DSM_HARV)
```

```
## [1] 1
```

However, raster data can also be multi-band, meaning that one raster file
contains data for more than one variable or time period for each cell. By
default the `raster()` function only imports the first band in a raster
regardless of whether it has one or more bands. 

**Calculate Raster Min and Max Values**

It is useful to know the minimum or maximum values of a raster dataset. In this case, given we are working with elevation data, these values represent the min/max elevation range at our site.

Raster statistics are often calculated and embedded in a GeoTIFF for us. We can view these values:

We can see that the elevation at our site ranges from 305.0700073m to
416.0699768m.



``` r
minmax(DSM_HARV)
```

```
##     HARV_dsmCrop
## min       305.07
## max       416.07
```

``` r
min(values(DSM_HARV))
```

```
## [1] 305.07
```

``` r
max(values(DSM_HARV))
```

```
## [1] 416.07
```
#### "Classifying" a raster

**Binning Raster Data and Plotting Binned data**

Above, we viewed our data using a continuous color ramp. For
clarity and visibility of the plot, we may prefer to view the data in discrete bins or colored according to ranges of values. This is comparable to a "classified" map. To do this, we need to tell `ggplot` how many groups to break our data  into, and where those breaks should be. To make these decisions, it is useful  to first explore the distribution of the data using a bar plot. To begin with,  we will use `dplyr`'s `mutate()` function combined with `cut()` to split the  data into 3 bins.

We want to set cutoff values for these groups to have the elevation ranges of 301–350 m, 351–400 m, and 401–450 m. To implement this we will give `mutate()` a numeric vector of break points instead of the number of breaks we want.


``` r
custom_bins <- c(300, 350, 400, 450)

DSM_HARV_df <- DSM_HARV_df %>%
  mutate(fct_elevation_2 = cut(HARV_dsmCrop, breaks = custom_bins))

unique(DSM_HARV_df$fct_elevation_2)
```

```
## [1] (400,450] (350,400] (300,350]
## Levels: (300,350] (350,400] (400,450]
```

We can use those groups to plot our raster data, with each group being a 
different color:


``` r
ggplot() +
  geom_raster(data = DSM_HARV_df , aes(x = x, y = y, fill = fct_elevation_2)) 
```

<img src="08-geospatial-I_files/figure-html/raster-with-breaks-1.png" width="672" />

The plot above uses the default colors inside `ggplot` for raster objects.
We can specify our own colors to make the plot look a little nicer.
R has a built in set of colors for plotting terrain, which are built in
to the `terrain.colors()` function.
Since we have three bins, we want to create a 3-color palette:


``` r
terrain.colors(3)
```

```
## [1] "#00A600" "#ECB176" "#F2F2F2"
```

The `terrain.colors()` function returns *hex colors* -
each of these character strings represents a color.
To use these in our map, we pass them across using the
`scale_fill_manual()` function.


``` r
ggplot() +
 geom_raster(data = DSM_HARV_df , aes(x = x, y = y,
                                      fill = fct_elevation_2)) + 
    scale_fill_manual(values = terrain.colors(3)) 
```

<img src="08-geospatial-I_files/figure-html/ggplot-breaks-customcolors-1.png" width="672" />

#### Layering Rasters

We can layer a raster on top of a hillshade raster for the same area, and use a transparency factor to create a 3-dimensional shaded effect. 

The **hillshade layer** maps the terrain using light and shadow to create a  3D-looking image, based on a **hypothetical illumination of the ground level**.A hillshade is a raster that maps the shadows and texture that you would see from above when viewing terrain. We will add a custom color, making the plot grey.

First we need to read in our DSM hillshade data and view the structure:


``` r
DSM_hill_HARV <-
  rast("data/NEON-DS-Airborne-Remote-Sensing/HARV/DSM/HARV_DSMhill.tif")

DSM_hill_HARV
```

```
## class       : SpatRaster 
## size        : 1367, 1697, 1  (nrow, ncol, nlyr)
## resolution  : 1, 1  (x, y)
## extent      : 731453, 733150, 4712471, 4713838  (xmin, xmax, ymin, ymax)
## coord. ref. : WGS 84 / UTM zone 18N (EPSG:32618) 
## source      : HARV_DSMhill.tif 
## name        : HARV_DSMhill 
## min value   :   -0.7136298 
## max value   :    0.9999997
```

Next we convert it to a dataframe, so that we can plot it using `ggplot2`:


``` r
DSM_hill_HARV_df <- as.data.frame(DSM_hill_HARV, xy = TRUE) 

str(DSM_hill_HARV_df)
```

```
## 'data.frame':	2313675 obs. of  3 variables:
##  $ x           : num  731454 731456 731456 731458 731458 ...
##  $ y           : num  4713836 4713836 4713836 4713836 4713836 ...
##  $ HARV_DSMhill: num  -0.15567 0.00743 0.86989 0.9791 0.96283 ...
```

Now we can plot the hillshade data:


``` r
ggplot() +
  geom_raster(data = DSM_hill_HARV_df,
              aes(x = x, y = y, alpha = HARV_DSMhill)) + 
  scale_alpha(range =  c(0.15, 0.65), guide = "none") 
```

<img src="08-geospatial-I_files/figure-html/raster-hillshade-1.png" width="672" />

We can layer another raster on top of our hillshade by adding another call to
the `geom_raster()` function. Let's overlay `DSM_HARV` on top of the `hill_HARV`.


``` r
ggplot() +
  geom_raster(data = DSM_HARV_df , 
              aes(x = x, y = y, 
                  fill = HARV_dsmCrop)) + 
  geom_raster(data = DSM_hill_HARV_df, 
              aes(x = x, y = y, 
                  alpha = HARV_DSMhill)) +  
  scale_fill_viridis_c() +  
  scale_alpha(range = c(0.15, 0.65), guide = "none") +  
  ggtitle("Elevation with hillshade") 
```

<img src="08-geospatial-I_files/figure-html/overlay-hillshade-1.png" width="672" />

# 6. Reprojection

Sometimes we encounter raster datasets that do not "line up" when plotted or analyzed. Rasters that don't line up are most often in different Coordinate Reference Systems (CRS).  We will work through how to deal with rasters in  different, known CRSs. It will walk though reprojecting rasters in R using 
the `project()` function in the `terra` package.

#### Raster Projection in R

We will read in the Harvard Forest Digital Terrain
Model data. This differs from the surface model data we've been working with so far in that the digital surface model (DSM) includes the tops of trees, while the digital terrain model (DTM) shows the ground level.

Here, we will create a map of the Harvard ForestDigital Terrain Model (`DTM_HARV`)  layered on top of the hillshade 
(`DTM_hill_HARV`).

First, we need to import the DTM and DTM hillshade data.


``` r
DTM_HARV <- 
    rast("data/NEON-DS-Airborne-Remote-Sensing/HARV/DTM/HARV_dtmCrop.tif")

DTM_hill_HARV <- 
    rast("data/NEON-DS-Airborne-Remote-Sensing/HARV/DTM/HARV_DTMhill_WGS84.tif")
```

Next, we will convert each of these datasets to a dataframe for
plotting with `ggplot`.


``` r
DTM_HARV_df <- as.data.frame(DTM_HARV, xy = TRUE)

DTM_hill_HARV_df <- as.data.frame(DTM_hill_HARV, xy = TRUE)
```

Now we can create a map of the DTM layered over the hillshade.


``` r
ggplot() +
     geom_raster(data = DTM_HARV_df , 
                 aes(x = x, y = y, 
                  fill = HARV_dtmCrop)) + 
     geom_raster(data = DTM_hill_HARV_df, 
                 aes(x = x, y = y, 
                   alpha = HARV_DTMhill_WGS84)) +
     scale_fill_gradientn(name = "Elevation", colors = terrain.colors(10)) 
```

<img src="08-geospatial-I_files/figure-html/unnamed-chunk-12-1.png" width="672" />

Our results are curious - neither the Digital Terrain Model (`DTM_HARV_df`)
nor the DTM Hillshade (`DTM_hill_HARV_df`) plotted.
Let's try to plot the DTM on its own to make sure there are data there.


``` r
ggplot() +
geom_raster(data = DTM_HARV_df,
    aes(x = x, y = y,
    fill = HARV_dtmCrop)) +
scale_fill_gradientn(name = "Elevation", colors = terrain.colors(10))
```

<img src="08-geospatial-I_files/figure-html/plot-DTM-1.png" width="672" />

Our DTM seems to contain data and plots just fine.

Next we plot the DTM Hillshade on its own to see whether everything is OK.


``` r
ggplot() +
geom_raster(data = DTM_hill_HARV_df,
    aes(x = x, y = y,
    alpha = HARV_DTMhill_WGS84))
```

<img src="08-geospatial-I_files/figure-html/plot-DTM-hill-1.png" width="672" />

If we look at the axes, we can see that the projections of the two rasters are  different. When this is the case, `ggplot` won't render the image. It won't even throw an  error message to tell you something has gone wrong. We can look at Coordinate  Reference Systems (CRSs) of the DTM and the hillshade data to see how they 
differ.


``` r
# view crs for DTM
crs(DTM_HARV, proj = TRUE)
```

```
## [1] "+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs"
```

``` r
# view crs for hillshade
crs(DTM_hill_HARV, proj = TRUE)
```

```
## [1] "+proj=longlat +datum=WGS84 +no_defs"
```

`DTM_HARV` is in the UTM projection, with units of meters.
`DTM_hill_HARV` is in `Geographic WGS84` - which is represented by latitude and longitude values.

Because the two rasters are in different CRSs, they don't line up when plotted in R. We need to reproject (or change the projection of) `DTM_hill_HARV` into the UTM CRS. Alternatively, we could reproject `DTM_HARV` into WGS84.

**Reproject Rasters**

We can use the `project()` function to reproject a raster into a new CRS.
Keep in mind that reprojection only works when you first have a defined CRS for the raster object that you want to reproject. It cannot be used if no CRS is defined. Lucky for us, the `DTM_hill_HARV` has a defined CRS.

When we reproject a raster, we move it from one "grid" to another. Thus, we are  modifying the data! Keep this in mind as we work with raster data.

To use the `project()` function, we need to define two things:

1. the object we want to reproject and
2. the CRS that we want to reproject it to.

The syntax is `project(RasterObject, crs)`

We want the CRS of our hillshade to match the `DTM_HARV` raster. We can thus assign the CRS of our `DTM_HARV` to our hillshade within the `project()` function as follows: `crs(DTM_HARV)`.
Note that we are using the `project()` function on the raster object,
not the `data.frame()` we use for plotting with `ggplot`.

First we will reproject our `DTM_hill_HARV` raster data to match the `DTM_HARV`  raster CRS:


``` r
DTM_hill_UTMZ18N_HARV <- project(DTM_hill_HARV,
                                 crs(DTM_HARV))
```

Now we can compare the CRS of our original DTM hillshade and our new DTM 
hillshade, to see how they are different.


``` r
crs(DTM_hill_UTMZ18N_HARV, proj = TRUE)
```

```
## [1] "+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs"
```

``` r
crs(DTM_HARV, proj = TRUE)
```

```
## [1] "+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs"
```

We can also compare the extent of the two objects.


``` r
ext(DTM_hill_UTMZ18N_HARV)
```

```
## SpatExtent : 731402.31567604, 733200.221994349, 4712407.19751409, 4713901.78222079 (xmin, xmax, ymin, ymax)
```

``` r
ext(DTM_hill_HARV)
```

```
## SpatExtent : -72.1819236223343, -72.1606102223342, 42.5294079700285, 42.5423355900285 (xmin, xmax, ymin, ymax)
```

Notice in the output above that the `crs()` of `DTM_hill_UTMZ18N_HARV` is now UTM. However, the extent values of `DTM_hillUTMZ18N_HARV` are different from `DTM_hill_HARV`.

## In class challange 1: Extent Change with CRS Change

Why do you think the two extents differ?

**Answers**

The extent for DTM\_hill\_UTMZ18N\_HARV is in UTMs so the extent is in meters.  The extent for DTM\_hill\_HARV is in lat/long so the extent is expressed in  decimal degrees.

**Deal with Raster Resolution**

Let's next have a look at the resolution of our reprojected hillshade versus our original data.


``` r
res(DTM_hill_UTMZ18N_HARV)
```

```
## [1] 1.001061 1.001061
```

``` r
res(DTM_HARV)
```

```
## [1] 1 1
```

These two resolutions are different, but they're representing the same data. We  can tell R to force our newly reprojected raster to be 1m x 1m resolution by  adding a line of code `res=1` within the `project()` function. In the  example below, we ensure a resolution match by using `res(DTM_HARV)` as a variable.


``` r
  DTM_hill_UTMZ18N_HARV <- project(DTM_hill_HARV, 
                                   crs(DTM_HARV), 
                                   res = res(DTM_HARV)) 
```

Now both our resolutions and our CRSs match, so we can plot these two data sets together. Let's double-check our resolution to be sure:


``` r
res(DTM_hill_UTMZ18N_HARV)
```

```
## [1] 1 1
```

``` r
res(DTM_HARV)
```

```
## [1] 1 1
```

For plotting with `ggplot()`, we will need to create a dataframe from our newly reprojected raster.


``` r
DTM_hill_HARV_2_df <- as.data.frame(DTM_hill_UTMZ18N_HARV, xy = TRUE)
```

We can now create a plot of this data.


``` r
ggplot() +
     geom_raster(data = DTM_HARV_df , 
                 aes(x = x, y = y, 
                  fill = HARV_dtmCrop)) + 
     geom_raster(data = DTM_hill_HARV_2_df, 
                 aes(x = x, y = y, 
                   alpha = HARV_DTMhill_WGS84)) +
     scale_fill_gradientn(name = "Elevation", colors = terrain.colors(10))
```

<img src="08-geospatial-I_files/figure-html/plot-projected-raster-1.png" width="672" />

We have now successfully draped the Digital Terrain Model on top of our
hillshade to produce a nice looking, textured map!



## In class challenge 2: Reproject, then Plot a Digital Terrain Model

Create a map of the
[San Joaquin Experimental Range](https://www.neonscience.org/field-sites/field-sites-map/SJER)
field site using the `SJER_DSMhill_WGS84.tif` and `SJER_dsmCrop.tif` files.

Reproject the data as necessary to make things line up!


**Answers**


``` r
# import DSM
DSM_SJER <- 
    rast("data/NEON-DS-Airborne-Remote-Sensing/SJER/DSM/SJER_dsmCrop.tif")
# import DSM hillshade
DSM_hill_SJER_WGS <-
    rast("data/NEON-DS-Airborne-Remote-Sensing/SJER/DSM/SJER_DSMhill_WGS84.tif")

# reproject raster
DSM_hill_UTMZ18N_SJER <- project(DSM_hill_SJER_WGS,
                                 crs(DSM_SJER),
                                 res = 1)

# convert to data.frames
DSM_SJER_df <- as.data.frame(DSM_SJER, xy = TRUE)

DSM_hill_SJER_df <- as.data.frame(DSM_hill_UTMZ18N_SJER, xy = TRUE)

ggplot() +
     geom_raster(data = DSM_hill_SJER_df, 
                 aes(x = x, y = y, 
                   alpha = SJER_DSMhill_WGS84)
                 ) +
     geom_raster(data = DSM_SJER_df, 
             aes(x = x, y = y, 
                  fill = SJER_dsmCrop,
                  alpha=0.8)
             ) + 
     scale_fill_gradientn(name = "Elevation", colors = terrain.colors(10))
```

<img src="08-geospatial-I_files/figure-html/challenge-code-reprojection-1.png" width="672" />


#### Create A Histogram of Raster Values

We can explore the distribution of values contained within our raster using the `geom_histogram()` function which produces a histogram. Histograms are often useful in identifying outliers and bad data values in our raster data.


``` r
ggplot() +
    geom_histogram(data = DSM_HARV_df, aes(HARV_dsmCrop), bins = 40)
```

<img src="08-geospatial-I_files/figure-html/view-raster-histogram2-1.png" width="672" />



