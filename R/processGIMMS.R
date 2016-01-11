### environmental stuff

# working directory
library(Orcs)
setwdOS(path_lin = "/media/fdetsch/XChange/", path_win = "D:/",
        path_ext = "kilimanjaro/ndvi_comparison")

## packages
lib <- c("Rsenal", "gimms", "doParallel", "MODIS", "remote")
jnk <- sapply(lib, function(x) library(x, character.only = TRUE, quietly = TRUE))

## functions
source("src/aggregateGimms.R")

## parallelization
cl <- makeCluster(3)
registerDoParallel(cl)

## geographic extent
rst_kili <- kiliAerial(minNumTiles = 20)

## download available GIMMS data
fls_gimms <- downloadGimms(dsn = "data/GIMMS")


## Data processing

## ndvi

# base directory
dir_base <- "data/rst/GIMMS3g/"

lst_gimms <- foreach(type = c("ndvi", "flag"), flag = c(FALSE, TRUE)) %do% {

  # rasterize ndvi
  dir_rst <- paste0(dir_base, "rst/", type, "/")
  fls_rst <- paste0(dir_rst, basename(fls_gimms), "_", type, ".tif")
  rst <- rasterizeGimms(fls_gimms, cores = 3L, flag = flag,
                        filename = fls_rst, format = "GTiff", overwrite = TRUE)

  # crop
  dir_crp <- paste0(dir_base, "crp/", type, "/")
  rst_crp <- crop(rst, extent(c(37, 37.72, -3.4, -2.84)), snap = "out",
                  filename = paste0(dir_crp, "CRP"), format = "GTiff",
                  overwrite = TRUE, bylayer = TRUE, suffix = names(rst))

  # project
  dir_prj <- paste0(dir_base, "prj/", type, "/")
  rst_prj <- projectRaster(rst_crp, crs = "+init=epsg:21037", format = "GTiff",
                           filename = paste0(dir_prj, "PRJ"), overwrite = TRUE,
                           bylayer = TRUE, suffix = names(rst_crp))

  # remove white margins
  spy_prj <- rasterToPolygons(rst_prj[[1]])

  dir_trm <- paste0(dir_base, "trm/", type, "/")
  rst_trm <- crop(rst_prj, spy_prj, filename = paste0(dir_trm, "TRM"),
                  format = "GTiff", overwrite = TRUE, bylayer = TRUE,
                  suffix = names(rst_prj))

  return(rst_trm)
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
