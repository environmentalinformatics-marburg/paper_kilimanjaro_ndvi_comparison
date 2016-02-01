### environmental stuff

# working directory
path_pkg <- getwd()

library(Orcs)
setwdOS(path_lin = "/media/fdetsch/XChange/", path_win = "D:/",
        path_ext = "kilimanjaro/ndvi_comparison")

## install old 'MODIS' version (bug in MODIS_0.10-33::whittaker.raster has not 
## been fixed yet)
install.packages(paste0(path_pkg, "inst/ex/MODIS_0.10-18.tar.gz"), repos = NULL)

## packages
lib <- c("Rsenal", "gimms", "doParallel", "MODIS", "remote")
jnk <- sapply(lib, function(x) library(x, character.only = TRUE, quietly = TRUE))

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
  rst <- rasterizeGimms(fls_gimms, cores = 2L, flag = flag,
                        filename = fls_rst, format = "GTiff", overwrite = TRUE)
  
  # crop
  dir_crp <- paste0(dir_base, "crp/", type, "/")
  rst_crp <- crop(rst, extent(c(37, 37.72, -3.4, -2.84)), snap = "out",
                  filename = paste0(dir_crp, "CRP"), format = "GTiff",
                  overwrite = TRUE, bylayer = TRUE, suffix = names(rst))
  
  # project
  dir_prj <- paste0(dir_base, "prj/", type, "/")
  rst_prj <- projectRaster(rst_crp, crs = "+init=epsg:21037", format = "GTiff",
                           method = ifelse(type == "ndvi", "bilinear", "ngb"), 
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

## re-import files
lst_gimms <- foreach(type = c("ndvi", "flag")) %do% {
  dir_trm <- paste0(dir_base, "trm/", type, "/")
  fls_trm <- rearrangeFiles(dsn = dir_trm, pattern = "^TRM_.*.tif$", 
                            pos = c(4, 6, 11) + 12, full.names = TRUE)
  rst_trm <- stack(fls_trm)
  return(rst_trm)
}


### flag correction ------------------------------------------------------------

# ## flag value counts
# mat <- as.matrix(lst_gimms[[2]])
# table(mat)
# 
# ## pixels with more than 0.25
# mat_log <- mat > 1
# (rowSums(mat_log) / ncol(mat_log)) > .25

## overlay ndvi with flags and keep 'good' values only (flags 1 and 2)
dir_qc <- paste0(dir_base, "qc/")
rst_qc <- overlay(lst_gimms[[1]], lst_gimms[[2]], fun = function(x, y) {
  x[y[] > 2] <- NA
  return(x)
}, filename = paste0(dir_qc, "QC"), format = "GTiff", bylayer = TRUE, 
suffix = names(lst_gimms[[1]]), overwrite = TRUE)

## re-import files
fls_qc <- rearrangeFiles(dsn = dir_qc, pattern = "^QC_.*.tif$", 
                         pos = c(4, 6, 11) + 15, full.names = TRUE)
rst_qc <- stack(fls_qc)


### whittaker smoothing --------------------------------------------------------

systime_locale <- Sys.getlocale(category = "LC_TIME")
if (Sys.info()[["sysname"]] == "Windows") {
  invisible(Sys.setlocale(category = "LC_TIME", locale = "C"))
} else {
  invisible(Sys.setlocale(category = "LC_TIME", locale = "en_US.UTF-8"))
}

# Setup `orgTime` object -> replace %Y%m with %Y%m%d (compatible to `as.Date` in
# `orgTime`)
org_gimms <- basename(fls_gimms)
for (i in 1:length(org_gimms)) {
  dt_yrmn <- substr(org_gimms[i], 4, 8)
  dt_yrmndy <- paste0(dt_yrmn, ifelse(i %% 2 == 1, "01", "15"))
  org_gimms[i] <- gsub(dt_yrmn, dt_yrmndy, org_gimms[i])
}

# org_gimms <- sapply(org_gimms[seq(1, length(org_gimms), 2)], function(i) {
#   gsub(substr(i, 12, 17), paste0(substr(i, 12, 17), "01"), i)
# })

org_gimms <-
  orgTime(org_gimms, pillow = 0, pos1 = 4, pos2 = 10, format = "%y%b%d")

# Perform Whittaker smoothing
dir_wht <- paste0(dir_base, "whittaker2")
rst_wht <-
  whittaker.raster(rst_qc, timeInfo = org_gimms, lambda = 6000, nIter = 3, 
                   removeOutlier = TRUE, threshold = 0.2, groupYears = FALSE,
                   outDirPath = dir_wht, overwrite = TRUE)

# store
fls_wht <- paste0(dir_wht, "/NDVI_YearlyLambda6000_fullPeriod.tif")
rst_wht <- stack(fls_wht)

fls_wht <- paste0(dir_wht, "/WHT_", basename(fls_qc))
ls_rst_wht <- foreach(i = unstack(rst_wht), j = fls_wht, 
                      .packages = c("raster", "rgdal")) %dopar% {
                        writeRaster(i, filename = j, 
                                    format = "GTiff", overwrite = TRUE)
                      }
rst_wht <- stack(ls_rst_wht)


### monthly aggregation --------------------------------------------------------

# available files
fls_wht <- rearrangeFiles(dsn = "data/rst/GIMMS3g/whittaker", 
                          pos = c(4, 6, 11) + 19, 
                          pattern = "^WHT.*.tif$", full.names = TRUE)

dir_mvc <- paste0(dir_base, "whittaker_mvc/")
fls_mvc <- paste0(dir_mvc, "MVC_", unique(substr(basename(fls_wht), 1, 27)))
rst_mvc <- monthlyComposite(fls_wht, pos1 = 23L, pos2 = 27L, cores = 3L,
                            filename = fls_mvc, format = "GTiff", 
                            overwrite = TRUE)

## re-import files
rst_mvc <- stack(paste0(fls_mvc, ".tif"))


### deseasoning ----------------------------------------------------------------

## create temporal subset (jan 82 to dec 12)
st <- grep("82jan", fls_mvc)
nd <- grep("12dec", fls_mvc)
fls_mvc <- fls_mvc[st:nd]
rst_mvc <- rst_mvc[[st:nd]]

rst_dsn <- deseason(rst_mvc, use.cpp = TRUE)

## write to file
dir_dsn <- paste0(dir_base, "whittaker_dsn/")
fls_dsn <- paste0(dir_dsn, "DSN_", basename(fls_mvc))

ls_rst_dsn <- foreach(i = unstack(rst_dsn), j = fls_dsn, 
                      .packages = c("raster", "rgdal")) %dopar% {
                        writeRaster(i, filename = j, format = "GTiff", 
                                    overwrite = TRUE)
                      }

rst_dsn <- stack(ls_rst_dsn)

## deregister parallel backend
stopCluster(cl)

## re-install most recent 'MODIS' version
detach("package:MODIS")
install.packages("http://download.r-forge.r-project.org/src/contrib/MODIS_0.10-33.tar.gz", 
                 repos = NULL)

# revoke locale time adjustment
Sys.setlocale(category = "LC_TIME", locale = systime_locale)
