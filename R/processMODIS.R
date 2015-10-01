### environmental stuff

# workspace clearance
rm(list = ls(all = TRUE))

# working directory
library(Orcs)
setwdOS(path_lin = "/media/fdetsch/XChange/", path_win = "D:/", 
        path_ext = "kilimanjaro/ndvi_comparison")

# packages
# install.packages("MODIS", repos="http://R-Forge.R-project.org")
lib <- c("raster", "rgdal", "MODIS", "doParallel", "Kendall", "RColorBrewer", 
         "reshape2", "ggplot2", "zoo", "GSODTools", "remote", "Rsenal")
sapply(lib, function(...) require(..., character.only = TRUE))

## functions
source("src/aggregateNDVICells.R")

## parallelization
cl <- makeCluster(3)
registerDoParallel(cl)


### Data import

# ## re-organize files
# orgStruc(move = TRUE)

## 'MODIS' global settings
MODISoptions(localArcPath = paste0(getwd(), "/data/MODIS_ARC/"), 
             outDirPath = paste0(getwd(), "/data/MODIS_ARC/PROCESSED/"), 
             MODISserverOrder = c("LAADS","LPDAAC"), quiet = TRUE)

## Extract .hdf container files for further processing
for (i in c("MOD13Q1", "MYD13Q1")) {

  #   ## rename collection 6 files (currently not supported by 'MODIS' package);
  #   ## remember to change folders (e.g. 'MOD13Q1.006') manually beforehand
  #   ch_fls_old <- list.files(getOption("MODIS_localArcPath"), pattern = i, 
  #                            recursive = TRUE, full.names = TRUE)
  #   ch_fls_new <- paste(dirname(ch_fls_old), 
  #                       gsub(".006.", ".005.", basename(ch_fls_old)), 
  #                       sep = "/")
  #   file.rename(ch_fls_old, ch_fls_new)
  
  runGdal(i, tileH = 21, tileV = 9, job = paste0(i, ".006"), begin = "2013-05-09", end = "2013-09-14", collection = "005",
          SDSstring = "100000000011", outProj = "EPSG:21037")
}

# 2012-04-12

## geographic extent
rst_kili <- kiliAerial(rasterize = TRUE, minNumTiles = 20)
# fls_gimms <- list.files("../gimms3g/gimms3g/data/rst/", pattern = "_crp_utm.tif$", 
#                         full.names = TRUE)
# rst_gimms <- raster(fls_gimms[[1]])
# kili <- rasterToPolygons(rst_gimms)

## plots
shp_plots <- suppressWarnings(
  readOGR(dsn = "data/shp/", 
          layer = "PlotPoles_ARC1960_mod_20140807_final", 
          p4s = "+init=epsg:21037"))
shp_plots_amp <- subset(shp_plots, PoleType == "AMP")

## DEM
rst_dem <- raster("data/dem/DEM_ARC1960_30m_Hemp.tif")

## NDVI data
stats <- lapply(c("MOD13Q1", "MYD13Q1"), function(h) {
  
  pttrn <- paste(h, c("NDVI.tif$", "pixel_reliability.tif$", 
                      "composite_day_of_the_year.tif$"), sep = ".*")
  
  ## crop
  ndvi.rst <- foreach(i = pttrn, .packages = c("raster", "rgdal")) %dopar% {                                      
    
    # available files
    fls <- list.files(paste0("data/MODIS_ARC/PROCESSED/", h, ".006"), 
                      pattern = i, full.names = TRUE)  
    
    # stack, crop and store
    rst <- stack(fls)
    rst.crp <- crop(rst, extent(rst_kili))
    rst.crp <- writeRaster(rst.crp, format = "GTiff", overwrite = TRUE,
                           filename = paste0("data/rst/", h, "/crp/CRP"), 
                           bylayer = TRUE, suffix = names(rst))
    return(rst.crp)
  }
  
  #   ndvi.rst <- ndvi.rst[1:2]
  #   
  #   ndvi.rst <- lapply(pttrn, function(i) {
  #     fls <- list.files("data/processed/", full.names = TRUE, 
  #                       pattern = paste("^CRP_", i, sep = ".*"))
  #     stack(fls)
  #   })
  
  ## quality assurance
  
  # overlay ndvi and qa layers
  ndvi.rst.qa <- overlay(ndvi.rst[[1]], ndvi.rst[[2]], fun = function(x, y) {
    x[!y[] %in% c(0:2)] <- NA
    return(x)
  })
  
  # store
  ndvi.rst.qa <- writeRaster(ndvi.rst.qa, format = "GTiff", overwrite = TRUE,
                             filename = paste0("data/rst/", h, "/qa/QA"), 
                             bylayer = TRUE, suffix = names(ndvi.rst[[1]]))
  
  #     ndvi.fls.qa <- list.files("data/processed/", full.names = TRUE, 
  #                               pattern = paste("^QA_", pttrn[1], sep = ".*"))
  #     ndvi.rst.qa <- stack(ndvi.fls.qa)
  
  ## additional outlier check
  
  # calc, tsOutliers
  ndvi.rst.qa.sd <- calc(ndvi.rst.qa, fun = function(x) {
    id <- tsOutliers(x, lower_quantile = .4, upper_quantile = .9, index = TRUE)
    x[id] <- NA
    return(x)
  })
  
  # store
  ndvi.rst.qa.sd <- writeRaster(ndvi.rst.qa.sd, format = "GTiff", overwrite = TRUE,
                                filename = paste0("data/rst/", h, "/sd/SD"), 
                                bylayer = TRUE, suffix = names(ndvi.rst.qa))
  
  #   ndvi.fls.qa.sd <- list.files("data/processed/", full.names = TRUE, 
  #                                pattern = paste("^SD_", pttrn[1], sep = ".*"))
  #   ndvi.rst.qa.sd <- stack(ndvi.fls.qa.sd)
  
  ## reject neighboring pixels
  
  # adjacent
  ndvi.rst.qa.sd.fc <- foreach(i = unstack(ndvi.rst.qa.sd), .combine = "stack",
                               .packages = c("raster", "rgdal")) %dopar% {
    cells <- which(is.na(i[]))
    id <- adjacent(i, cells = cells, directions = 8, pairs = FALSE)
    i[id] <- NA
    return(i)
  }

  # store
  ndvi.rst.qa.sd.fc <- writeRaster(ndvi.rst.qa.sd.fc, format = "GTiff",
                                   filename = paste0("data/rst/", h, "/adj/ADJ"),  
                                   bylayer = TRUE, suffix = names(ndvi.rst.qa.sd), 
                                   overwrite = TRUE)

  #   ndvi.fls.qa.sd.fc <- list.files("data/processed/", full.names = TRUE, 
  #                                   pattern = paste("^BF_SD_QA_CRP", pttrn[1], sep = ".*"))
  #   ndvi.rst.qa.sd.fc <- stack(ndvi.fls.qa.sd.fc)
  # 
  # dates <- orgTime(ndvi.fls.init)$inputLayerDates
  # dates_agg <- dates + 8
  # yearmon_agg <- as.yearmon(dates_agg)
  # indices_agg <- as.numeric(as.factor(yearmon_agg))
  # 
  # outdir <- paste0("data/processed")
  # rst_qa_sd_fc_aggmax <- 
  #   stackApply(ndvi.rst.qa.sd.fc, indices = indices_agg, fun = max, bylayer = TRUE,
  #              filename = paste0(outdir, "/AGGMAX_BF_SD_QA_", h), format = "GTiff",  
  #              suffix = strftime(unique(yearmon_agg), format = "%Y%m"), 
  #              overwrite = TRUE)
  
  
  ### Gap filling
  
  ## initial files (for date information)
  ndvi.fls.init <- list.files(paste0("data/MODIS_ARC/PROCESSED/", h, ".006"),
                              pattern = paste(h, "NDVI.tif$", sep = ".*"), 
                              full.names = TRUE, recursive = TRUE)
  
  # org_agg <- sapply(strsplit(names(rst_qa_sd_fc_aggmax), "_"), "[[", 6)
  # org_agg <- paste0(org_agg, "01")
  # org_agg <- orgTime(org_agg, nDays = "1 month", pillow = 0, 
  #                    pos1 = 1, pos2 = 8, format = "%Y%m%d")
  
  ## application of whittaker smoothing algorithm
  rst.wht <- whittaker.raster(vi = ndvi.rst.qa.sd.fc, removeOutlier = TRUE, 
                              threshold = 2000,
                              timeInfo = orgTime(ndvi.fls.init, pillow = 0), 
                              lambda = 6000, nIter = 3, groupYears = FALSE, 
                              outDirPath = paste0("data/rst/", h, "/whittaker"), 
                              overwrite = TRUE, format = "raster")
  
  ## monthly aggregation
  
  #   # import files
  #   fls.wht <- list.files(paste0("data/processed/whittaker_agg1km_", tolower(h)), 
  #                         pattern = "^WHT.*.tif$", full.names = TRUE)
  #   rst.wht <- stack(fls.wht)
  
  dates_agg <- dates + 8
  yearmon_agg <- as.yearmon(dates_agg)
  indices_agg <- as.numeric(as.factor(yearmon_agg))
  
  outdir <- paste0("data/rst/", h, "/whittaker")
  # rst_wht_aggmax <- 
  #   stackApply(rst.wht, indices = indices_agg, fun = max, bylayer = TRUE,
  #              filename = paste0(outdir, "/AGGMAX_WHT"), format = "GTiff",  
  #              suffix = strftime(unique(yearmon_agg), format = "%Y%m"), 
  #              overwrite = TRUE)
  # 
  # rst_wht_aggmin <- 
  #   stackApply(rst.wht, indices = indices_agg, fun = min, bylayer = TRUE,
  #              filename = paste0(outdir, "/AGGMIN_WHT"), format = "GTiff",  
  #              suffix = strftime(unique(yearmon_agg), format = "%Y%m"), 
  #              overwrite = TRUE)
  
  fls_doy <- list.files(paste0("data/rst/", h, "/crp"), full.names = TRUE,
                        pattern = paste0("^CRP_", toupper(h), ".*composite"))
  rst_doy <- ndvi.rst[[3]]
  
  dates_doy <- as.numeric(substr(basename(fls_doy), 14, 20))
  months_doy <- unique(strftime(as.Date(as.character(dates_doy), format = "%Y%j"), "%Y%m"))
  
  file_out <- paste0("data/rst/", h, "/whittaker_mvc/MVC")
  rst_wht_mvc <- aggregateNDVICells(rst = stack(rst.wht), 
                                    rst_doy = rst_doy, 
                                    dates = dates_doy, 
                                    n_cores = 3, 
                                    save_output = TRUE, filename = file_out, 
                                    bylayer = TRUE, suffix = months_doy, 
                                    format = "GTiff", overwrite = TRUE)
  
  
  # Application of scale factor and removal of inconsistent values
  # fls_wht_aggmax <- list.files(outdir, pattern = "^AGGMAX_WHT", full.names = TRUE)
  # rst_wht_aggmax <- stack(fls_wht_aggmax)
  dir_scl <- paste0("data/rst/", h, "/whittaker_scl")
  fls_scl <- paste0(dir_scl, "/SCL_", names(rst_wht_mvc))
  
  rst_scl <- foreach(i = unstack(rst_wht_mvc), j = as.list(fls_scl), 
                     .packages = c("raster", "rgdal"), .combine = "stack") %dopar% {  
                       
    # scale factor                   
    rst <- i
    rst <- rst / 10000
    
    # rejection of inconsistent values
    id <- which(rst[] < -1 | rst[] > 1)
    
    if (length(id) > 0) {
      rst[id] <- NA
    }
    
    # store
    rst <- writeRaster(rst, filename = j, format = "GTiff", overwrite = TRUE)
    
    return(rst)
  }
  
  ## deseasoning
  fls_scl <- list.files(dir_scl, pattern = "^SCL_MVC", full.names = TRUE)
  
  # start and end month
  st <- grep("200301", fls_scl)
  nd <- grep("201212", fls_scl)
  
  # deseason
  rst_scl <- stack(fls_scl[st:nd])
  rst_dsn <- deseason(rst_scl)
  
  # store
  dir_dsn <- paste0("data/rst/", h, "/whittaker_dsn")
  fls_dsn <- paste0(dir_dsn, "/DSN_", names(rst_scl))
  rst_dsn <- foreach(i = unstack(rst_dsn), j = as.list(fls_dsn), 
                     .packages = c("raster", "rgdal"), .combine = "stack") %dopar% {
    writeRaster(i, filename = j, format = "GTiff", overwrite = TRUE)
  }
  
  # # 2011-2013
  # fls_scl <- list.files(outdir, pattern = "^SCL", full.names = TRUE)
  # 
  # # st <- grep("201101", fls_scl)
  # # nd <- grep("201312", fls_scl)
  # # fls_scl <- fls_scl[st:nd]
  # 
  # rst_scl <- stack(fls_scl)
  # 
  # mat_plt_val <- extract(rst_scl, plt_apoles)
  # df_plt_val <- data.frame(PlotID = plt_apoles@data[, 1], mat_plt_val)
  # names(df_plt_val)[2:ncol(df_plt_val)] <- 
  #   sapply(strsplit(names(df_plt_val)[2:ncol(df_plt_val)], "_"), "[[", 4)
  # write.csv(df_plt_val, "out/csv/ndvi_aggmin_200207_201409.csv", row.names = FALSE)
  # 
  # months <- substr(sapply(strsplit(fls_scl, "_"), "[[", 5), 5, 6)
  # indices <- as.numeric(as.factor(months))
  # 
  # rst_scl_mmonth <- stackApply(rst_scl, indices = indices, fun = mean)
  # 
  # plt_scl_mmonth <- data.frame(PlotId = plt_apoles$PlotID, 
  #                              extract(rst_scl_mmonth, plt_apoles))
  # 
  # id <- grep("mai2", plt_scl_mmonth$PlotId)
  # plot(unlist(plt_scl_mmonth[id, 2:ncol(plt_scl_mmonth)]), type = "l", 
  #      xlab = "Month", ylab = "NDVI")
  # 
  #   st <- ifelse(h == "MOD13Q1", "2001", "2003")
  #   st <- "2003"
  #   nd <- "2013"
  # 
  #   fls.wht <- fls.wht[grep(st, fls.wht):grep(nd, fls.wht)]
  #   rst.wht <- stack(fls.wht)
  #   
  #   rst.mk <- overlay(rst.wht, fun = function(x) MannKendall(as.numeric(x))$tau, 
  #                     filename = paste("out/MK", toupper(h), unique(substr(names(rst.wht), 1, 21)), st, nd, sep = "_"), 
  #                     format = "GTiff", overwrite = TRUE)
  #   rst.mk.p <- overlay(rst.wht, fun = function(x) MannKendall(as.numeric(x))$sl, 
  #                       filename = paste("out/MK_p", toupper(h), unique(substr(names(rst.wht), 1, 21)), st, nd, sep = "_"), 
  #                       format = "GTiff", overwrite = TRUE)
  # 
  #   fls.mk <- list.files("out", pattern = paste0("MK_", toupper(h), ".*.tif$"), 
  #                        full.names = TRUE)
  #   rst.mk <- raster(fls.mk)
  #   
  #   png(paste0(substr(fls.mk, 1, nchar(fls.mk)-4), ".png"), units = "mm", 
  #       width = 300, res = 300, pointsize = 20)
  #   print(spplot(rst.mk, scales = list(draw = TRUE), xlab = "x", ylab = "y", 
  #          col.regions = colorRampPalette(brewer.pal(11, "BrBG")), 
  #          sp.layout = list("sp.lines", rasterToContour(dem)), 
  #          par.settings = list(fontsize = list(text = 15)), at = seq(-.9, .9, .1)))
  #   dev.off()
  #   
  #   stats <- lapply(c(.01, .001), function(i) {
  #     rst.mk <- overlay(rst.wht, fun = function(x) {
  #       mk <- MannKendall(as.numeric(x))
  #       if (mk$sl >= i) return(NA) else return(mk$tau)
  #     }, filename = paste("out/MK", i, toupper(h), 
  #                         unique(substr(names(rst.wht), 1, 21)), sep = "_"), 
  #     format = "GTiff", overwrite = TRUE)
  #     
  #     fls.mk <- list.files("out", pattern = paste(i, toupper(h), ".tif$", sep = ".*"), 
  #                          full.names = TRUE)
  #     rst.mk <- raster(fls.mk)
  #     
  #     png(paste0(substr(fls.mk, 1, nchar(fls.mk)-4), ".png"), units = "mm", 
  #         width = 300, res = 300, pointsize = 20)
  #     print(spplot(rst.mk, scales = list(draw = TRUE), xlab = "x", ylab = "y", 
  #                  col.regions = colorRampPalette(brewer.pal(11, "BrBG")), 
  #                  sp.layout = list("sp.lines", rasterToContour(dem)), 
  #                  par.settings = list(fontsize = list(text = 15)), 
  #                  at = seq(-.9, .9, .1)))
  #     dev.off()
  #     
  #     val <- round(sum(!is.na(rst.mk[]))/ncell(rst.mk), digits = 3)
  #     
  #     val.pos <- round(sum(rst.mk[] > 0, na.rm = TRUE) / sum(!is.na(rst.mk[])), 3)
  #     val.neg <- round(sum(rst.mk[] < 0, na.rm = TRUE) / sum(!is.na(rst.mk[])), 3)
  #     
  #     return(data.frame(sensor = h, p = as.character(i), nona = val, 
  #                       nona_pos = val.pos, nona_neg = val.neg))
  #   })
  #   
  #   return(do.call("rbind", stats))
})

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

