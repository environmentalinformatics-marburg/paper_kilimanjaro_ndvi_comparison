### environmental stuff

## functions (located in repository, not on external hdd)
source("R/aggregateGimms.R")

## working directory
library(Orcs)
setwdOS(path_lin = "/media/fdetsch/XChange/", path_win = "D:/", 
        path_ext = "kilimanjaro/ndvi_comparison")

## install old 'MODIS' version (bug in MODIS_0.10-33::whittaker.raster has not 
## been fixed yet)
install.packages("http://download.r-forge.r-project.org/src/contrib/MODIS_0.10-18.tar.gz", 
                 repos = NULL)

## packages
lib <- c("Rsenal", "doParallel", "MODIS", "remote")
jnk <- sapply(lib, function(x) library(x, character.only = TRUE, quietly = TRUE))

## 'MODIS' global settings
MODISoptions(localArcPath = paste0(getwd(), "/data/MODIS_ARC/"), 
             outDirPath = paste0(getwd(), "/data/MODIS_ARC/PROCESSED/"), 
             MODISserverOrder = c("LAADS","LPDAAC"), quiet = TRUE)

## parallelization
cl <- makeCluster(4)
registerDoParallel(cl)

## geographic extent
rst_kili <- kiliAerial(rasterize = TRUE, minNumTiles = 20)

## download available GIMMS data
fls_gimms <- downloadGimms(dsn = "data/GIMMS")

# Rearrange GIMMS files according to timestamp
fls_gimms <- rearrangeGimms(dsn = "data/GIMMS", 
                            pattern = "^geo",
                            rename_yearmon = TRUE,
                            full.names = TRUE)

# Create .hdr companion file
fls_hdr <- createHdr(file = "data/gimms3g.hdr")


## Data processing
rst_gimms <- 
  foreach(i = fls_gimms, .packages = c("Rsenal", "zoo"), .combine = "stack", 
          .export = ls(envir = globalenv())) %dopar% {
            
    # rasterize
    fls_rst <- paste0("data/rst/GIMMS3g/rst/", basename(i))
    rst <- rasterizeGimms(file = i, 
                          headerfile = fls_hdr, 
                          file_out = fls_rst, 
                          format = "GTiff", overwrite = TRUE)
    
    # crop
    fls_crp <- paste0("data/rst/GIMMS3g/crp/CRP_", basename(fls_rst), ".tif")
    rst_crp <- crop(rst, extent(c(37, 37.72, -3.4, -2.84)), snap = "out", 
                    filename = fls_crp, format = "GTiff", overwrite = TRUE)
    
    # project
    fls_prj <- paste0("data/rst/GIMMS3g/prj/PRJ_", basename(fls_crp))
    rst_crp_prj <- projectRaster(rst_crp, crs = "+init=epsg:21037", 
                                 filename = fls_prj, 
                                 format = "GTiff", overwrite = TRUE)
    # return processed raster
    return(rst_crp_prj)
}


## remove white margins

# trim
fls_prj <- list.files("data/rst/GIMMS3g/prj", pattern = "^PRJ_.*.tif$", 
                      full.names = TRUE)
rst_prj <- stack(fls_prj)

spy_prj <- rasterToPolygons(rst_prj[[1]])

dir_trm <- "data/rst/GIMMS3g/trm/"
fls_trm <- paste0(dir_trm, "TRM_", basename(fls_prj))

rst_trm <- foreach(i = (1:nlayers(rst_prj))[c(283, 284, 291)], .packages = c("raster", "rgdal"), 
                   .combine = "stack") %dopar% {
                     rst <- crop(rst_prj[[i]], spy_prj)
                     rst <- writeRaster(rst, filename = fls_trm[i], format = "GTiff", 
                                        overwrite = TRUE)
                   }

# reimport files
fls_trm <- list.files(dir_trm, pattern = "^TRM_.*.tif$", full.names = TRUE)
rst_trm <- stack(fls_trm)

# ## test
# st <- grep("2000", fls_trm)[1]
# nd <- grep("2010", fls_trm)[length(grep("2012", fls_trm))]
# rst_trm <- rst_trm[[st:nd]]
# fls_prj <- fls_prj[st:nd]

## gap-filling

# setup `orgTime` object -> replace %Y%m with %Y%m%d (compatible to `as.Date` in
# `orgTime`)
org_gimms <- basename(fls_prj)
for (i in 1:length(org_gimms)) {
  dt_yrmn <- substr(org_gimms[i], 12, 17)
  dt_yrmndy <- paste0(dt_yrmn, ifelse(substr(org_gimms[i], 18, 20) == "15a", "01", "15"))
  org_gimms[i] <- gsub(dt_yrmn, dt_yrmndy, org_gimms[i])
}

org_gimms <- 
  orgTime(org_gimms, pillow = 0, pos1 = 12, pos2 = 19, format = "%Y%m%d")

# whittaker.raster
rst_wht <- 
  whittaker.raster(rst_trm, timeInfo = org_gimms, lambda = 6000, nIter = 3, 
                   removeOutlier = TRUE, outlierThreshold = 0.2, groupYears = FALSE,
                   outDirPath = "data/rst/GIMMS3g/whittaker", overwrite = TRUE)

# store
rst_wht <- stack("data/rst/GIMMS3g/whittaker/NDVI_YearlyLambda6000_fullPeriod.tif")

dir_wht <- "data/rst/GIMMS3g/whittaker/"
fls_wht <- paste0(dir_wht, "WHT_", basename(fls_trm))
ls_rst_wht <- foreach(i = 1:nlayers(rst_wht), j = fls_wht, 
                      .packages = c("raster", "rgdal")) %dopar% {
  writeRaster(rst_wht[[i]], filename = j, format = "GTiff", overwrite = TRUE)
}
rst_wht <- stack(ls_rst_wht)

## Monthly aggregation

# available files
fls_wht <- list.files("data/rst/GIMMS3g/whittaker", 
                      pattern = "^WHT_TRM.*.tif$", full.names = TRUE)

# aggregate
dir_mvc <- "data/rst/GIMMS3g/whittaker_mvc/"
fls_mvc <- paste0(dir_mvc, "MVC_", unique(substr(basename(fls_wht), 1, 25)))
rst_agg <- aggregateGimms(files = fls_wht, files_out = fls_mvc, 
                          start = 20, stop = 25, nodes = 3L, 
                          format = "GTiff", overwrite = TRUE)

## deseason

# import data
fls_agg <- list.files("data/rst/GIMMS3g/whittaker_mvc/", 
                      pattern = "^MVC_.*.tif$", full.names = TRUE)
fls_agg <- fls_agg[7:length(fls_agg)] # remove jul-dec 1981
rst_agg <- stack(fls_agg)

# deseason
rst_dsn <- deseason(rst_agg, use.cpp = TRUE)

# store
dir_dsn <- "data/rst/GIMMS3g/whittaker_dsn/"
fls_dsn <- paste0(dir_dsn, "DSN_", basename(fls_agg))

ls_rst_dsn <- foreach(i = 1:nlayers(rst_dsn), j = fls_dsn, 
                   .packages = c("raster", "rgdal")) %dopar% {
  writeRaster(rst_dsn[[i]], filename = j, format = "GTiff", overwrite = TRUE)
}
rst_dsn <- stack(ls_rst_dsn)

## deregister parallel backend
stopImplicitCluster()

## re-install most recent 'MODIS' version
detach("package:MODIS")
install.packages("http://download.r-forge.r-project.org/src/contrib/MODIS_0.10-33.tar.gz", 
                 repos = NULL)
