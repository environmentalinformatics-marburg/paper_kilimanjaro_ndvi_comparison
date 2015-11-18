## clear workspace
rm(list = ls(all = TRUE))

## packages
library(raster)
library(bfast)
library(doParallel)
library(plyr)

## functions
source("R/visDEM.R")

## parallelization
cl <- makeCluster(3)
registerDoParallel(cl)


### data processing

## folders
ch_dir_extdata <- "/media/fdetsch/XChange/kilimanjaro/ndvi_comparison/data/rst/"
ch_dir_outdata <- "/media/fdetsch/XChange/kilimanjaro/ndvi_comparison/out/"

## start and end year
st_year <- 2003
nd_year <- 2012

## reimport data
fls_bfast <- list.files("data/", pattern = "_bfast.rds$", full.names = TRUE)
fls_bfast <- fls_bfast[c(1, 2, 4, 3, 5)]
lst_bfast <- vector("list", length(fls_bfast))
for (i in 1:length(fls_bfast)) {
  
  # import .rds file
  dat_tmp <- readRDS(fls_bfast[[i]])
  
  # identify and remove missing data (due to missing values)
  log_isdf <- sapply(dat_tmp, function(j) {
    is.data.frame(j)
  })
  if (any(!log_isdf))
    dat_tmp <- dat_tmp[log_isdf]
  
  # write to list
  lst_bfast[[i]] <- do.call("rbind.fill", dat_tmp)
  rm(dat_tmp)
}


### bfast value insertion

## raster templates
rst_tmp_gimms <- raster(paste0(ch_dir_extdata, "GIMMS3g/whittaker_dsn/DSN_MVC_WHT_TRM_PRJ_CRP_geo198201.tif"))
rst_tmp_gimms[] <- NA
rst_bfast_gimms <- rst_bfast_gimms_maxtime <- rst_bfast_gimms_maxmagn <- rst_tmp_gimms

rst_tmp_modis <- raster(paste0(ch_dir_extdata, "MYD13Q1.006/whittaker_dsn/DSN_SCL_MVC_200301.tif"))
rst_tmp_modis[] <- NA
rst_bfast_modis <- rst_bfast_modis_maxtime <- rst_bfast_modis_maxmagn <- rst_tmp_modis

## reference extent
fls_ndvi <- paste0(ch_dir_extdata, "MOD13Q1.006/whittaker_dsn/DSN_SCL_MVC_200301.tif")

rst_ndvi <- raster(fls_ndvi[1])
rst_ndvi <- projectRaster(rst_ndvi, crs = "+init=epsg:4326")
rst_ndvi <- trim(rst_ndvi)

num_xmin <- xmin(rst_ndvi)
num_xmax <- xmax(rst_ndvi)
num_ymin <- ymin(rst_ndvi)
num_ymax <- ymax(rst_ndvi)

## digital elevation model (dem)
ch_fls_dem <- paste0(ch_dir_extdata, "../dem/DEM_ARC1960_30m_Hemp.tif")
rst_dem <- raster(ch_fls_dem)
rst_dem <- aggregate(rst_dem, fact = 10)
rst_dem <- projectRaster(rst_dem, crs = "+init=epsg:4326")
p_dem <- visDEM(rst_dem, labcex = .6, cex = 1.6, col = "grey70")

## insert values
lst_p <- foreach(i = lst_bfast, j = 1:length(lst_bfast)) %do% {
  
  # templates
  if (j == 1) {
    rst_bp <- rst_maxmagn <- rst_maxtime <- rst_bfast_gimms
  } else {
    rst_bp <- rst_maxmagn <- rst_maxtime <- rst_bfast_modis
  }
  
  # insert values
  rst_bp[i$cell] <- i$bp_vt_ct
  rst_bp <- projectRaster(rst_bp, crs = "+init=epsg:4326", method = "ngb")
  
  rst_maxmagn[i$cell] <- i$bp_vt_maxmagn
  rst_maxmagn <- projectRaster(rst_maxmagn, crs = "+init=epsg:4326", 
                               method = "ngb")
  
  rst_maxtime[i$cell] <- i$bp_vt_maxtime
  rst_maxtime <- projectRaster(rst_maxtime, crs = "+init=epsg:4326", 
                               method = "ngb")
  
  cat("Total counts:", maxValue(rst_bp), 
      "\tMagnitudes:", minValue(rst_maxmagn), ",", maxValue(rst_maxmagn), 
      "\tTimes:", minValue(rst_maxtime), ",", maxValue(rst_maxtime), "\n")
  
  # # reject pixels with low magnitude
  # int_lowmagn <- which(rst_maxmagn[] < .1 & rst_maxmagn[] > -.1)
  # rst_bp[int_lowmagn] <- rst_maxmagn[int_lowmagn] <- rst_maxtime[int_lowmagn] <- NA
  
  # number of breakpoints
  col_counts <- brewer.pal(5, "YlGnBu")
  p_counts <- spplot(rst_bp, col.regions = col_counts, 
                     at = seq(.5, 5.5, 1), scales = list(draw = TRUE), 
                     colorkey = FALSE, xlim = c(num_xmin, num_xmax), 
                     ylim = c(num_ymin, num_ymax)) + 
    as.layer(p_dem)
  
  # maximum magnitude
  col_maxmagn <- colorRampPalette(brewer.pal(9, "BrBG"))
  # brks <- quantile(rst_breaks_maxmagn, seq(0, 1, length.out = 256))
  p_maxmagn <- spplot(rst_maxmagn, col.regions = col_maxmagn(1000),
                      at = seq(-.6, .6, .01), scales = list(draw = TRUE), 
                      xlim = c(num_xmin, num_xmax), ylim = c(num_ymin, num_ymax)) + 
    as.layer(p_dem)
  
  # timing of maximum magnitude
  col_maxtime <- envinmrPalette(1000)[100:1000]
  p_maxtime <- spplot(rst_maxtime, col.regions = col_maxtime, 
                      at = seq(2003, 2012, 1), scales = list(draw = TRUE), 
                      colorkey = FALSE, xlim = c(num_xmin, num_xmax), 
                      ylim = c(num_ymin, num_ymax)) + 
    as.layer(p_dem)
  
  # combine and return figures
  lst_p_tmp <- list(p_counts, p_maxtime, p_maxmagn)
  latticeCombineGrid(lst_p_tmp, layout = c(3, 1))
}

png(paste0(ch_dir_outdata, "figure_03.png"), width = 20, height = 25, 
    units = "cm", res = 500)
plot.new()

vp0 <- viewport(x = 0, y = 0, width = 1, height = .9, 
                just = c("left", "bottom"), name = "vp_count")
pushViewport(vp0)
print(latticeCombineGrid(lst_p, layout = c(3, 5)), newpage = FALSE)

# key left
downViewport(trellis.vpname("figure"))
vp1 <- viewport(x = 0, y = 1.04, width = 1/3, height = .1, 
                just = c("left", "bottom"), name = "vp_key1")
pushViewport(vp1)
draw.colorkey(key = list(col = col_counts, width = 1, height = .5,
                         at = seq(.5, 5.5, 1), 
                         space = "bottom"), draw = TRUE)
grid.text("Trend breakpoints", x = 0.5, y = 1.12, just = c("centre", "top"), 
          gp = gpar(font = 2, cex = .85))

# key right
upViewport()
vp2 <- viewport(x = 2/3, y = 1.04, width = 1/3, height = .1, 
                just = c("left", "bottom"), name = "vp_key2")
pushViewport(vp2)
draw.colorkey(key = list(col = col_maxmagn(1000), width = 1, height = .5,
                         at = seq(-.6, .6, .01), 
                         space = "bottom"), draw = TRUE)
grid.text("Magnitude of biggest change", x = 0.5, y = 1.12, 
          just = c("centre", "top"), gp = gpar(font = 2, cex = .85))

# key center
upViewport()
vp3 <- viewport(x = 1/3, y = 1.04, width = 1/3, height = .1, 
                just = c("left", "bottom"), name = "vp_key3")
pushViewport(vp3)
draw.colorkey(key = list(col = col_maxtime, width = 1, height = .5,
                         at = 2003:2012, space = "bottom"), draw = TRUE)
grid.text("Timing of biggest change", x = 0.5, y = 1.12, 
          just = c("centre", "top"), gp = gpar(font = 2, cex = .85))

dev.off()

