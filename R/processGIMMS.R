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


## Whittaker smoothing

# Avl files
fls_gimms <- list.files("data/rst/", pattern = "_crp_utm.tif$", 
                        full.names = TRUE)

# Setup `orgTime` object -> replace %Y%m with %Y%m%d (compatible to `as.Date` in
# `orgTime`)
org_gimms <- basename(fls_gimms)
org_gimms <- sapply(org_gimms, function(i) {
  gsub(substr(i, 4, 9), paste0(substr(i, 4, 9), "01"), i)
})
org_gimms <- 
  orgTime(org_gimms, pillow = 0, pos1 = 4, pos2 = 11, format = "%Y%m%d")

# Perform Whittaker smoothing
rst_gimms_wht <- 
  whittaker.raster(fls_gimms, timeInfo = org_gimms, lambda = 6000, nIter = 3, 
                   removeOutlier = TRUE, threshold = 0.2, groupYears = FALSE,
                   outDirPath = "data/rst/whittaker", overwrite = TRUE)

# Save Whittaker files separately
rst_gimms_wht <- 
  foreach(i = rst_gimms_wht[[1]], j = fls_gimms, .combine = "stack") %do% {
            outdir <- paste0(dirname(j), "/whittaker")
            outfile <- paste0(substr(basename(j), 1, nchar(basename(j))-4), "_wht")
            file_out <- paste(outdir, outfile, sep = "/")
            writeRaster(i, filename = file_out, format = "GTiff", overwrite = TRUE)
          }


## Monthly aggregation

# Avl files
fls_gimms_wht <- list.files("data/rst/whittaker", pattern = "_crp_utm_wht.tif$", 
                            full.names = TRUE)

# Outnames
fls_out <- paste0(substr(basename(fls_gimms_wht), 1, 9),
                  substr(basename(fls_gimms_wht), 17, (nchar(basename(fls_gimms_wht))-4)), 
                  "_aggmax")
fls_out <- paste(unique(dirname(fls_gimms_wht)), unique(fls_out), sep = "/")

# Aggregation, `fun = max`
rst_gimms_agg <- aggregateGimms(files = fls_gimms_wht, 
                                start = 4, stop = 9)

# Output storage
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