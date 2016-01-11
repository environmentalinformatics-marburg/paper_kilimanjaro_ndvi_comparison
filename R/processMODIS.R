### environmental stuff

# workspace clearance
rm(list = ls(all = TRUE))

# working directory
# library(Orcs)
# setwdOS(path_lin = "/media/fdetsch/XChange/", path_win = "D:/", 
#         path_ext = "kilimanjaro/ndvi_comparison")

library(Orcs)
setwdOS(path_lin = "/media/fdetsch/XChange/", path_win = "D:/", 
        path_ext = "kilimanjaro/ndvi_comparison")

# packages
# install.packages("MODIS", repos="http://R-Forge.R-project.org")
lib <- c("raster", "rgdal", "MODIS", "foreach", "Kendall", "RColorBrewer", 
         "reshape2", "ggplot2", "zoo", "GSODTools", "remote", "Rsenal")
loadPkgs(lib)

## functions
source("/media/permanent/repositories/paper_kilimanjaro_ndvi_comparison/R/qcMCD13.R")
source("src/aggregateNDVICells.R")


### Data import

# ## re-organize files
# orgStruc(move = TRUE)

## 'MODIS' global settings
MODISoptions(localArcPath = paste0(getwd(), "/data/MODIS_ARC/"), 
             outDirPath = paste0(getwd(), "/data/MODIS_ARC/PROCESSED/"), 
             MODISserverOrder = c("LAADS","LPDAAC"), quiet = TRUE)

## Extract .hdf container files for further processing
for (i in c("MOD13Q1", "MYD13Q1")) {
  
  runGdal(i, tileH = 21, tileV = 9, job = paste0(i, ".006"), 
          collection = "006",
          SDSstring = "000000000010", outProj = "EPSG:21037")
}

## geographic extent
rst_kili <- kiliAerial(minNumTiles = 20L)

## plots
shp_plots <- suppressWarnings(
  readOGR(dsn = "/media/permanent/kilimanjaro/coordinates/coords/", 
          layer = "PlotPoles_ARC1960_mod_20140807_final", 
          p4s = "+init=epsg:21037"))
shp_plots_amp <- subset(shp_plots, PoleType == "AMP")

## DEM
rst_dem <- raster("/media/permanent/kilimanjaro/coordinates/coords/DEM_ARC1960_30m_Hemp.tif")

## NDVI data
lst_ndvi_qc <- lapply(c("MOD13Q1", "MYD13Q1"), function(product) {
  qcMCD13(product, ref_ext = rst_kili, 
          inpath = paste0(options()$MODIS_outDirPath, product, ".005"), 
          dsn = paste0("data/rst/", product, ".005/"), 
          cores = 3L, apply_crop = FALSE, apply_qc = FALSE, 
          apply_tso = FALSE, apply_adj = FALSE)
})  

### Gap filling

lst_ndvi <- foreach(product = list("MOD13Q1", "MYD13Q1")) %do% {
  
  ## initial files (for date information)
  ndvi.fls.init <- list.files(paste0("data/MODIS_ARC/PROCESSED/", product, ".006"),
                              pattern = paste(product, "NDVI.tif$", sep = ".*"), 
                              full.names = TRUE, recursive = TRUE)
  
  ## application of whittaker smoothing algorithm
  drs_wht <- paste0("data/rst/", product, ".006/whittaker")
  #             rst.wht <- whittaker.raster(vi = layers, removeOutlier = TRUE, 
  #                                         threshold = 2000,
  #                                         timeInfo = orgTime(ndvi.fls.init, pillow = 0), 
  #                                         lambda = 6000, nIter = 3, groupYears = FALSE, 
  #                                         outDirPath = drs_wht, 
  #                                         overwrite = TRUE, format = "raster")
  
  
  # Application of scale factor and removal of inconsistent values
  fls_wht <- list.files(drs_wht, pattern = "^MCD.*yL6000.ndvi.tif$", 
                        full.names = TRUE)
  rst_wht <- raster::stack(fls_wht)
  dir_scl <- paste0("data/rst/", product, ".006/whittaker_scl")
  fls_scl <- paste0(dir_scl, "/SCL_", names(rst_wht))
  
  #           lst_scl <- foreach(i = unstack(rst_wht), j = as.list(fls_scl), 
  #                              .packages = c("raster", "rgdal"), 
  #                              .export = ls(envir = globalenv())) %dopar% {  
  #                                
  #                                # scale factor                   
  #                                rst <- i
  #                                rst <- rst / 10000
  #                                
  #                                # rejection of inconsistent values
  #                                id <- which(rst[] < -1 | rst[] > 1)
  #                                
  #                                if (length(id) > 0) {
  #                                  rst[id] <- NA
  #                                }
  #                                
  #                                # store
  #                                rst <- writeRaster(rst, filename = j, format = "GTiff", overwrite = TRUE)
  #                                
  #                                return(rst)
  #                              }
  #           
  #           rst_scl <- stack(lst_scl)
  
  fls_scl <- substr(fls_scl, 1, nchar(fls_scl) - 5)
  fls_scl <- paste0(fls_scl, ".tif")
  rst_scl <- raster::stack(fls_scl)
  
  ## monthly aggregation (shifted scaling and mvc creation; remember to 
  ## adjust this when running the script the next time)
  
  dates <- MODIS::orgTime(ndvi.fls.init)$inputLayerDates
  dates_agg <- dates + 8
  yearmon_agg <- zoo::as.yearmon(dates_agg)
  indices_agg <- as.numeric(as.factor(yearmon_agg))
  
  outdir <- paste0("data/rst/", product, ".006/whittaker")
  
  # composite day of the year
  fls_doy <- list.files(paste0("data/rst/", product, ".006/crp"), 
                        full.names = TRUE,
                        pattern = paste0("^CRP_", toupper(product), ".*composite"))
  rst_doy <- raster::stack(fls_doy)
  
  # output file
  dates_doy <- as.numeric(substr(basename(fls_doy), 14, 20))
  months_doy <- unique(strftime(as.Date(as.character(dates_doy), 
                                        format = "%Y%j"), "%Y%m"))
  
  file_out <- paste0("data/rst/", product, ".006/whittaker_mvc/MVC")
  
  #           id_start <- grep("2013", fls_doy)[1]
  #           
  #           id_end <- grep("2014", fls_doy)
  #           id_end <- id_end[length(id_end)]
  #           
  #           rst_scl <- rst_scl[[id_start:id_end]]
  #           rst_doy <- rst_doy[[id_start:id_end]]
  #           dates_doy <- dates_doy[id_start:id_end]
  #           months_doy <- months_doy[127:(127+23)]
  
  # aggregate to maximum value composites (mvc)
  cat("Initializing creation of maximum value composites for", product, "...\n")
  rst_wht_mvc <- aggregateNDVICells(rst = rst_scl, 
                                    rst_doy = rst_doy, 
                                    dates = dates_doy, 
                                    cores = 3, 
                                    save_output = TRUE, filename = file_out, 
                                    bylayer = TRUE, suffix = months_doy, 
                                    format = "GTiff", overwrite = TRUE)
  
  # reimport mvc layers
  fls_wht_mvc <- list.files(paste0("data/rst/", product, ".006/whittaker_mvc/"), 
                            pattern = "^MVC.*.tif$", full.names = TRUE)
  
  
  ## deseasoning
  fls_scl <- list.files(dir_scl, pattern = "^SCL_MVC", full.names = TRUE)
  
  # start and end month
  st <- grep("200301", fls_scl)
  nd <- grep("201412", fls_scl)
  
  # deseason
  rst_scl <- raster::stack(fls_scl[st:nd])
  rst_dsn <- remote::deseason(rst_scl)
  
  # store
  dir_dsn <- paste0("data/rst/", h, "/whittaker_dsn")
  fls_dsn <- paste0(dir_dsn, "/DSN_", names(rst_scl))
  lst_dsn <- foreach(i = raster::unstack(rst_dsn), j = as.list(fls_dsn)) %do% {
                       raster::writeRaster(i, filename = j, format = "GTiff", overwrite = TRUE)
  }
  
  rst_dsn <- raster::stack(lst_dsn)
}

# Store percentage information about significant NDVI pixels
write.csv(do.call("rbind", stats), "out/mk_na_stats.csv", row.names = FALSE)

# Remove white margins from output images
system("cd out/; for file in *.png; do convert -trim $file $file; done")

# Deregister parallel backend
stopCluster(cl)


## Visualization of Whittaker gap-filling performance

# Import quality-controlled raster data of both Terra and Aqua
rst.orig <- lapply(c("MOD13Q1", "MYD13Q1"), function(h) {
  
  fls <- list.files("data/processed/", full.names = TRUE, 
                    pattern = paste("^CRP_", h, "NDVI.tif$", sep = ".*"))
  
  #   st <- ifelse(h == "MOD13Q1", "2001", "2003")
  st <- "2003"
  nd <- "2013"
  
  fls <- fls[grep(st, fls)[1]:grep(nd, fls)[length(grep(nd, fls))]]
  rst <- stack(fls)
  
  return(rst)
})

# Import quality-controlled raster data of both Terra and Aqua
rst.qa <- lapply(c("MOD13Q1", "MYD13Q1"), function(h) {
  
  fls <- list.files("data/processed/", full.names = TRUE, 
                    pattern = paste("^BF_", h, sep = ".*"))
  
  #   st <- ifelse(h == "MOD13Q1", "2001", "2003")
  st <- "2003"
  nd <- "2013"
  
  fls <- fls[grep(st, fls)[1]:grep(nd, fls)[length(grep(nd, fls))]]
  rst <- stack(fls)
  
  return(rst)
})

# Import Whittaker-filled raster data and corresponding dates
rst.wht <- lapply(c("MOD13Q1", "MYD13Q1"), function(h) {
  fls <- list.files(paste0("data/processed/whittaker_", tolower(h)), 
                    pattern = "NDVI_Year.*_year.*.tif$", full.names = TRUE)
  
  #   st <- ifelse(h == "MOD13Q1", "2001", "2003")
  st <- "2003"
  nd <- "2013"
  
  fls <- fls[grep(st, fls)[1]:grep(nd, fls)[length(grep(nd, fls))]]
  rst <- stack(fls)
  
  return(rst)
})

dat.wht <- lapply(c("MOD13Q1", "MYD13Q1"), function(h) {
  fls.crp <- list.files("data/processed/",
                        pattern = paste("^CRP", h, "NDVI.tif$", sep = ".*"))
  
  #   st <- ifelse(h == "MOD13Q1", "2001", "2003")
  st <- "2003"
  nd <- "2013"
  
  fls.crp <- fls.crp[grep(st, fls.crp)[1]:
                       grep(nd, fls.crp)[length(grep(nd, fls.crp))]]
  
  return(substr(fls.crp, 14, 20))
})

# Extract cell numbers for each KiLi plot, and extract corresponding time
# series from Terra and Aqua RasterStack objects
cell.numbers <- data.frame(PlotID = plt.wgs.apoles@data$PlotID, 
                           Cell = cellFromXY(rst.wht[[1]], plt.wgs.apoles))

# Original, quality-checked and gap-filled data as matrices
mat.orig.terra <- as.matrix(rst.orig[[1]])
mat.qa.terra <- as.matrix(rst.qa[[1]])
mat.wht.terra <- as.matrix(rst.wht[[1]])

mat.orig.aqua <- as.matrix(rst.orig[[2]])
mat.qa.aqua <- as.matrix(rst.qa[[2]])
mat.wht.aqua <- as.matrix(rst.wht[[2]])

index <- grep("cof2", cell.numbers[, 1])

ts.orig.terra <- mat.orig.terra[cell.numbers[index, 2], ] / 10000
ts.gappy.terra <- mat.qa.terra[cell.numbers[index, 2], ] / 10000
ts.filled.terra <- mat.wht.terra[cell.numbers[index, 2], ] / 10000

ts.orig.aqua <- mat.orig.aqua[cell.numbers[index, 2], ] / 10000
ts.gappy.aqua <- mat.qa.aqua[cell.numbers[index, 2], ] / 10000
ts.filled.aqua <- mat.wht.aqua[cell.numbers[index, 2], ] / 10000

# Terra data
dat.terra <- data.frame(date = as.Date(dat.wht[[1]], format = "%Y%j"),
                        "Original" = ts.orig.terra,
                        "Filtered" = ts.gappy.terra, 
                        "Imputed" = ts.filled.terra)

dat.terra <- melt(dat.terra, id.vars = 1)
dat.terra$sensor <- "Terra"

# Aqua data
dat.aqua <- data.frame(date = as.Date(dat.wht[[2]], format = "%Y%j"),
                       "Original" = ts.orig.aqua,
                       "Filtered" = ts.gappy.aqua, 
                       "Imputed" = ts.filled.aqua)

dat.aqua <- melt(dat.aqua, id.vars = 1)
dat.aqua$sensor <- "Aqua"

# Plot Terra and Aqua data separately, 
# quality-controlled data = grey, gap-filled data = black
# ggplot(aes(x = date, y = value, group = variable, colour = variable), 
#        data = dat.terra) + 
#   geom_line(size = .6) + 
#   facet_wrap(~ sensor) +
#   scale_color_manual("", values = c("Original" = "grey75", "Imputed" = "black")) +
#   labs(x = "Time [16-day intervals]", y = "NDVI") +
#   theme_bw()

png("out/whittaker_gapfill_performance_cof2.png", width = 25, height = 15, 
    units = "cm", res = 300, pointsize = 16)
ggplot(aes(x = date, y = value, group = variable, colour = variable), 
       data = dat.aqua) + 
  geom_line(size = 1) + 
  scale_color_manual("", values = c("Original" = "grey75", 
                                    "Filtered" = "darkolivegreen2", 
                                    "Imputed" = "darkred")) +
  labs(x = "Time [16-day intervals]", y = "NDVI") +
  theme_bw()
dev.off()

# Plot Terra and Aqua data in combination
dat <- rbind(dat.terra, dat.aqua)
dat$sensor <- factor(dat$sensor, levels = c("Terra", "Aqua"))

# ggplot(aes(x = date, y = value, group = variable, colour = variable), 
#        data = dat) + 
#   geom_line(size = .6) + 
# #   stat_smooth(method = "lm", se = FALSE, size = 2) + 
#   facet_wrap(~ sensor, ncol = 1) + 
#   scale_color_manual("", values = c("Original" = "grey75", "Imputed" = "black")) +
#   labs(x = "Time [16-day intervals]", y = "NDVI") +
#   theme_bw()

# Compare Terra and Aqua (e.g., COF3)
cof3.terra <- data.frame(date = as.Date(dat.wht[[1]], format = "%Y%j"), 
                         ndvi = mat.wht.terra[cell.numbers[index, 2], ] / 10000, 
                         sensor = "Terra")
cof3.aqua <- data.frame(date = as.Date(dat.wht[[2]], format = "%Y%j"), 
                        ndvi = mat.wht.aqua[cell.numbers[index, 2], ] / 10000, 
                        sensor = "Aqua")
cof3 <- rbind(cof3.terra, cof3.aqua)
cof3$plot <- cell.numbers[index, 1]

ggplot(aes(x = date, y = ndvi, group = sensor, colour = sensor), data = cof3) + 
  geom_line(size = .6) + 
  stat_smooth(method = "lm", se = FALSE, size = 2) + 
  facet_wrap(~ plot) + 
  scale_color_manual("", values = c("Terra" = "brown", "Aqua" = "blue")) +
  labs(x = "Time [16-day intervals]", y = "NDVI") +
  theme_bw()


