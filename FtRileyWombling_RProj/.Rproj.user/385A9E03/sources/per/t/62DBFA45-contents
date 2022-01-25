##############################################################################
########### Tracking spatial regimes in animal communities: ##################
############ implications for resilience-based management ####################
############# Roberts et al. Ecological Indicators, 2022 #####################
##############################################################################

###### Description of data:
# 'GWR_Bird_outputs_AllYears.csv'
# GWRr2   = wombling R^2 values
# GWRbx   = wombling beta coefficient for x coordinates
# GWRby   = wombling beta coefficient for y coordinates 
# GWR_xSE = standard error for wombling beta coefficient for x coordinates 
# GWR_ySE = standard error for wombling beta coefficient for y coordinates 
# DMCE    = UTM zone 14 easting coordinate
# DMCN    = UTM zone 14 northing coordinate
# Year    = year of point-count survey
# PlotID  = identifier for the point-count survey transect
# long    = longitude
# lat     = latitude
#
  # NOTE: see Diniz-Filho et al., 2016 in Genetica and 'womblingFunctions.R'
  # script for details about computing wombling values
  
#=============================================================================
## Preparation

# # Clean environment?
# rm(list=ls())

require("librarian")

# List of packages necessary to run this script:
shelf(stringr, dplyr, ggplot2, car, ggsn, vegan, rgeos,
      raster, data.table, here, scales, rgdal, sp, 
      gratia, cowplot, mgcv, 
      lib = tempdir())

# Set paths
FigurePath <- "SpCov_Figures"
ResultsPath <- "SpCov_Results"
DataPath <- "SpCov_Data"

# Source some functions
source(here("womblingFunctions.R"))

# Set random numbers
set.seed(41001)

#=============================================================================
## Load data

# Load wombling outputs + point-count survey locations + dates
outputs <- fread(here("GWR_Bird_outputs_AllYears.csv"))

# Load grass:woody spatial regime boundary (spatial covariance) data
projRasts_139 <- read_project(here("SpCov_RAPv2/R2_NLCD2016_FtRiley_30m_PFGC_TREE_139W"))

# Load Ft Riley boundary and maneuver areas
riley_boundary <- readOGR(here(DataPath), "FtRiley_Boundary")
riley <- readOGR(dsn = here("FtRiley_Boundaries"), 
                 layer = "FtRiley_ManeuverAreas")
riley <- spTransform(riley, CRS("+proj=utm +zone=14 +datum=WGS84 +units=m"))

#=============================================================================
## Predict wombling values across a continuous surface

# Run hierarchical generalized additive model (HGAM)
fit1 <- gam(GWRr2 ~ te(Year, DMCE, DMCN),
            method = "REML", 
            data = outputs)
# summary(fit1)
# gratia::appraise(fit1)

# Convert the Ft Riley SpatialPolygons to a raster and remove the Impact Area
rileyRast <- mask(crop(rasterize(riley, 
                                 projRasts_139[[1]]),
                       subset(riley, Maneuver != "IA")),
                  subset(riley, Maneuver != "IA"))

# Create a data.frame for the newdata argument in predict()
nd <- data.table(as.data.frame(rileyRast, xy = TRUE))
nd <- nd[!is.na(layer_Maneuver)]
setnames(nd, c("x", "y"), c("DMCE", "DMCN"))
pred_nd <- data.table(Year = rep(c(1991, 2005, 2017), 
                                 each = nrow(nd)),
                      rbind(nd, nd, nd))

# Predict HGAM across all pixels in the Ft Riley raster
predFit <- predict(fit1, 
                   newdata = pred_nd,
                   type = "response")
predFit_dt <- data.table(fit = predFit,
                         pred_nd)

#=============================================================================
## Fire data

fileFun <- function(theDir,
                    searchKey = "burn_bndy") {
  ## Look for files (directories included for now)
  allFiles <- list.files(theDir, no.. = TRUE, full.names = TRUE)
  ## Look for directory names
  allDirs <- list.dirs(theDir, full.names = FALSE, recursive = FALSE)
  ## If there are any directories,
  if(length(allDirs)) {
    ## then call this function again
    moreFiles <- lapply(file.path(theDir, allDirs), fileFun)
    ## Set names for the new list
    names(moreFiles) <- allDirs
    ## Determine files found, excluding directory names
    outFiles <- allFiles[!allFiles %in% allDirs]
    ## Combine appropriate results for current list
    if(length(outFiles)) {
      allFiles <- c(outFiles, moreFiles)
    } else {
      allFiles <- moreFiles
    }
  }
  return(allFiles)
}
# debug(fileFun)

## Try with your directory
fireBoundaries <- unlist(fileFun(here("FtRiley_Fire_Data")))

# Load all fire boundary polygons
firePolys <- lapply(unique(str_remove(fireBoundaries[str_detect(fireBoundaries, "burn_bndy")], "\\.([a-z])+")), 
                    function(X) readOGR(str_sub(X, 1, 92), str_sub(X, 94, str_length(X))))
firePolys <- do.call(rbind, firePolys)
firePolys@data$ID <- as.numeric(as.factor((firePolys@data$Ig_Date)))
firePolys@data$Year <- as.numeric(str_sub(firePolys@data$Ig_Date, 1, 4))
firePolys <- spTransform(firePolys, CRS(proj4string(projRasts_139[[1]])))

# Rasterize fire polygons and convert to a data.table
fire_dt <- rbindlist(lapply(c(1991, 2005, 2017),
                            function(X) {
                              fireRast <- mask(rasterize(subset(firePolys, Year < X), 
                                                         rileyRast,
                                                         field = "Year",
                                                         fun= "count"),
                                               subset(riley, Maneuver != "IA"))
                              fire_dt <- na.omit(as.data.table(as.data.frame(fireRast, xy = TRUE)))
                              # fire_dt <- fire_dt[!is.na(layer_Maneuver)]
                              fire_dt$Year <- X
                              return(fire_dt)
                            }))
setnames(fire_dt, "layer", "Number of Fires")

#=============================================================================
## Conduct hypothesis tests

# Select focal years from GAM-predicted wombling values
womble <- predFit_dt[Year %in% c(1991, 2005, 2017)]

# select the same focal years from spatial covariance rasters and convert to
# data.table
spcov <- rbindlist(lapply(list("R2_NLCD2016_FtRiley_30m_PFGC_TREE_139W_1991.tif",
                               "R2_NLCD2016_FtRiley_30m_PFGC_TREE_139W_2005.tif",
                               "R2_NLCD2016_FtRiley_30m_PFGC_TREE_139W_2017.tif"),
                          function(X) {
                            temp <- data.table(as.data.frame(projRasts_139[[X]], xy = TRUE))
                            setnames(temp, names(temp), c("DMCE", "DMCN", "spcov"))
                            temp$Year <- as.numeric(str_sub(X, -8, -5))
                            return(na.omit(temp))
                          }))

# Join wombling, spatial covariance, and fire data
tst_cor <- left_join(na.omit(womble[ , c("fit", "DMCE", "DMCN", "Year")]),
                     na.omit(spcov))
setnames(fire_dt, c("x", "y", "Number of Fires"), c("DMCE", "DMCN", "n_fires"))

# Remove NAs, constrain GAM-predicted wombling values to be between 0-1, replace
# NAs in the fire data with 0, and make year a factor.
tst_cor <- na.omit(left_join(tst_cor, fire_dt))
tst_cor[fit > 1]$fit <- 1
tst_cor[fit < 0]$fit <- 0
tst_cor[is.na(n_fires)]$n_fires <- 0
tst_cor[ , "Year_factor" := as.factor(Year)]

# Randomly select 10% of pixels in each year
tst_cor_sample <- 
  tst_cor %>% 
  group_by(Year) %>% 
  sample_n(nrow(tst_cor)/3 * 0.1)

# Run HGAM to determine how number of fires and spatial covariance affect
# predicted wombling values
fit_womble <- gam(fit ~ s(spcov) + s(n_fires, k = 5) + s(Year_factor, bs = "re"), 
                  data = tst_cor_sample, 
                  method = "REML")
# summary(fit_womble)
# appraise(fit_womble)
# draw(fit_womble)

#=============================================================================
## Make plots

source(here("FtRileyWombling_MakePlots.R"))
