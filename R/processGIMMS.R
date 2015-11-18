### environmental stuff

## working directory
# library(Orcs)
# setwdOS(path_lin = "/media/fdetsch/XChange/", path_win = "D:/", 
#         path_ext = "kilimanjaro/ndvi_comparison")

library(Orcs)
setwdOS(path_lin = "/media/fdetsch/modis_data/", path_win = "D:/", 
        path_ext = "kilimanjaro/ndvi_comparison")

## install old 'MODIS' version (bug in MODIS_0.10-33::whittaker.raster has not 
## been fixed yet)
install.packages("http://download.r-forge.r-project.org/src/contrib/MODIS_0.10-18.tar.gz", 
                 repos = NULL)

## packages
lib <- c("Rsenal", "doParallel", "MODIS", "remote", "gimms")
jnk <- sapply(lib, function(x) library(x, character.only = TRUE, quietly = TRUE))

## 'MODIS' global settings
MODISoptions(localArcPath = paste0(getwd(), "/data/MODIS_ARC/"), 
             outDirPath = paste0(getwd(), "/data/MODIS_ARC/PROCESSED/"), 
             MODISserverOrder = c("LAADS","LPDAAC"), quiet = TRUE)

## parallelization
cl <- makeCluster(3)
registerDoParallel(cl)

## geographic extent
rst_kili <- kiliAerial(minNumTiles = 20)

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
lst_gimms <- foreach(i = fls_gimms, .packages = c("Rsenal", "zoo")) %dopar% {
            
  # rasterize
  fls_rst <- paste0("data/rst/GIMMS3g/rst/", basename(i), ".tif")
  
  if (file.exists(fls_rst)) {
    rst <- raster(fls_rst)
  } else {
    rst <- rasterizeGimms(file = i, 
                          headerfile = fls_hdr, 
                          file_out = fls_rst, 
                          format = "GTiff", overwrite = TRUE)
  }
  
  # crop
  fls_crp <- paste0("data/rst/GIMMS3g/crp/CRP_", basename(fls_rst), ".tif")
  rst_crp <- crop(rst, extent(c(37, 37.72, -3.4, -2.84)), snap = "out", 
                  filename = fls_crp, format = "GTiff", overwrite = TRUE)
  
  # project
  fls_prj <- paste0("data/rst/GIMMS3g/prj/PRJ_", basename(fls_crp))
  rst_crp_prj <- projectRaster(rst_crp, crs = "+init=epsg:21037", 
                               filename = fls_prj, 
                               format = "GTiff", overwrite = TRUE)
  
  # remove margins
  rst_crp_prj <- trim(rst_crp_prj, filename = fls_prj, 
                      format = "GTiff", overwrite = TRUE)
  
  # return processed raster
  return(rst_crp_prj)
}

rst_gimms <- stack(lst_gimms)

# reimport files
fls_prj <- rearrangeFiles(dsn = "data/rst/GIMMS3g/prj", pos = c(4, 6, 11) + 8, 
                          pattern = "^PRJ.*.tif$", full.names = TRUE)
rst_prj <- stack(fls_prj)

## gap-filling

systime_locale <- Sys.getlocale(category = "LC_TIME")
if (Sys.info()[["sysname"]] == "Windows") {
  invisible(Sys.setlocale(category = "LC_TIME", locale = "C"))
} else {
  invisible(Sys.setlocale(category = "LC_TIME", locale = "en_US.UTF-8"))
}

# setup `orgTime` object -> replace %Y%m with %Y%m%d (compatible to `as.Date` in
# `orgTime`)
org_gimms <- basename(fls_prj)
for (i in 1:length(org_gimms)) {
  dt_yrmn <- substr(org_gimms[i], 12, 16)
  dt_yrmndy <- paste0(dt_yrmn, ifelse(substr(org_gimms[i], 17, 19) == "15a", "01", "15"))
  org_gimms[i] <- gsub(dt_yrmn, dt_yrmndy, org_gimms[i])
}

org_gimms <- 
  orgTime(org_gimms, pillow = 0, pos1 = 12, pos2 = 19, format = "%y%b%d")

# whittaker.raster
rst_wht <- 
  whittaker.raster(rst_prj, timeInfo = org_gimms, lambda = 6000, nIter = 3, 
                   removeOutlier = TRUE, threshold = 0.2, groupYears = FALSE,
                   outDirPath = "data/rst/GIMMS3g/whittaker", overwrite = TRUE)

# store
rst_wht <- stack("data/rst/GIMMS3g/whittaker/NDVI_YearlyLambda6000_fullPeriod.tif")

dir_wht <- "data/rst/GIMMS3g/whittaker/"
fls_wht <- paste0(dir_wht, "WHT_", basename(fls_prj))
ls_rst_wht <- foreach(i = 1:nlayers(rst_wht), j = fls_wht, 
                      .packages = c("raster", "rgdal")) %dopar% {
  writeRaster(rst_wht[[i]], filename = j, format = "GTiff", overwrite = TRUE)
}
rst_wht <- stack(ls_rst_wht)

## Monthly aggregation

# available files
fls_wht <- rearrangeFiles(dsn = "data/rst/GIMMS3g/whittaker", 
                          pos = c(4, 6, 11) + 12, 
                          pattern = "^WHT.*.tif$", full.names = TRUE)

rst_mvc <- 
  monthlyComposite(fls_wht, pos1 = 16L, pos2 = 20L, 
                   filename = "data/rst/GIMMS3g/whittaker_mvc/MVC", 
                   format = "GTiff", overwrite = TRUE, bylayer = TRUE, 
                   suffix = unique(substr(basename(fls_wht), 1, 20)))

## deseason

# import data
fls_mvc <- list.files("data/rst/GIMMS3g/whittaker_mvc/", 
                      pattern = "^MVC_.*.tif$", full.names = TRUE)

dts_mvc <- as.Date(paste0(substr(basename(fls_mvc), 20, 24), "01"), 
                   format = "%y%b%d")
fls_mvc <- fls_mvc[order(dts_mvc)]
fls_mvc <- fls_mvc[-(1:6)]
rst_mvc <- stack(fls_mvc)

# deseason
rst_dsn <- deseason(rst_mvc, use.cpp = TRUE)

# store
drs_dsn <- "data/rst/GIMMS3g/whittaker_dsn/"
fls_dsn <- paste0(drs_dsn, "DSN_", basename(fls_mvc))

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

# revoke locale time adjustment
Sys.setlocale(category = "LC_TIME", locale = systime_locale)
