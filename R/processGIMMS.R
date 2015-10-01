### environmental stuff

# working directory
library(Orcs)
setwdOS(path_lin = "/media/fdetsch/XChange/", path_win = "D:/", 
        path_ext = "kilimanjaro/ndvi_comparison")

## packages
lib <- c("Rsenal", "doParallel", "MODIS", "remote")
jnk <- sapply(lib, function(x) library(x, character.only = TRUE, quietly = TRUE))

## functions
source("src/aggregateGimms.R")

## parallelization
cl <- makeCluster(3)
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

fls_trm <- paste0("data/rst/GIMMS3g/trm/TRM_", basename(fls_prj))

rst_trm <- foreach(i = 1:nlayers(rst_prj), .packages = c("raster", "rgdal"), 
                   .combine = "stack") %dopar% {
  rst <- crop(rst_prj[[i]], spy_prj)
  rst <- writeRaster(rst, filename = fls_trm[i], format = "GTiff", 
                     overwrite = TRUE)
}

## Whittaker smoothing

# Setup `orgTime` object -> replace %Y%m with %Y%m%d (compatible to `as.Date` in
# `orgTime`)
org_gimms <- basename(fls_gimms)
for (i in 1:length(org_gimms)) {
  dt_yrmn <- substr(org_gimms[i], 12, 17)
  dt_yrmndy <- paste0(dt_yrmn, ifelse(i %% 2 == 1, "01", "15"))
  org_gimms[i] <- gsub(dt_yrmn, dt_yrmndy, org_gimms[i])
}

# org_gimms <- sapply(org_gimms[seq(1, length(org_gimms), 2)], function(i) {
#   gsub(substr(i, 12, 17), paste0(substr(i, 12, 17), "01"), i)
# })

org_gimms <- 
  orgTime(org_gimms, pillow = 0, pos1 = 12, pos2 = 19, format = "%Y%m%d")

# Perform Whittaker smoothing
rst_gimms_wht <- 
  whittaker.raster(fls_gimms, timeInfo = org_gimms, lambda = 6000, nIter = 3, 
                   removeOutlier = TRUE, threshold = 0.2, groupYears = FALSE,
                   outDirPath = "data/rst/GIMMS3g/whittaker", overwrite = TRUE)

## Monthly aggregation

# available files
fls_gimms_wht <- list.files("data/rst/GIMMS3g/whittaker", 
                            pattern = "^MCD.*.tif$", full.names = TRUE)

# aggregate
rst_gimms_agg <- aggregateGimms(files = fls_gimms_wht, 
                                start = 5, stop = 11)

# Output storage
fls_ext <- strftime(as.Date(substr(basename(fls_gimms_wht), 5, 11), format = "%Y%j"), 
                    format = "%Y%m")
fls_out <- paste0("data/rst/GIMMS3g/whittaker_mvc/MVC_MCD.", fls_ext, ".yL6000.ndvi")

ls_gimms_agg <- lapply(seq(nlayers(rst_gimms_agg)), function(i) {
  writeRaster(rst_gimms_agg[[i]], filename = fls_out[i], 
              format = "GTiff", overwrite = TRUE)
})

# fls_gimms_agg <- list.files("data/rst/", pattern = "aggmax.tif$", 
#                             full.names = TRUE)
# ls_gimms_agg <- lapply(fls_gimms_agg, raster)
rst_gimms_agg <- stack(ls_gimms_agg)


## Deseasoning

fls_gimms_aggmax <- list.files("data/rst/whittaker", pattern = "_aggmax.tif$", 
                            full.names = TRUE)
fls_gimms_aggmax <- 
  fls_gimms_aggmax[grep("1982", fls_gimms_aggmax)[1]:length(fls_gimms_aggmax)]
rst_gimms_aggmax <- stack(fls_gimms_aggmax)

rst_gimms_aggmax_dsn <- deseason(rst_gimms_aggmax)

# Outnames
dir_out <- unique(dirname(fls_gimms_aggmax))
fls_out <- paste0(dir_out, "/", names(rst_gimms_aggmax), "_dsn")

ls_gimms_aggmax_dsn <- lapply(1:nlayers(rst_gimms_aggmax_dsn), function(i) {
  writeRaster(rst_gimms_aggmax_dsn[[i]], filename = fls_out[i], 
              format = "GTiff", overwrite = TRUE) 
})

rst_gimms_aggmax_dsn <- stack(ls_gimms_aggmax_dsn)

# Deregister parallel backend
stopCluster(cl)