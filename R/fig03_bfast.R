## clear workspace
rm(list = ls(all = TRUE))

## packages
library(raster)
library(bfast)
library(doParallel)
library(plyr)

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

## insert values
rst_bfast_gimms[lst_bfast[[1]]$cell] <- lst_bfast[[1]]$bp_vt_ct

lst_p <- foreach(i = lst_bfast, j = 1:length(lst_bfast)) %do% {
  
  # templates
  if (j == 1) {
    rst_bp <- rst_maxmagn <- rst_maxtime <- rst_bfast_gimms
  } else {
    rst_bp <- rst_maxmagn <- rst_maxtime <- rst_bfast_modis
  }
  
  # insert values
  rst_bp[i$cell] <- i$bp_vt_ct
  rst_maxmagn[i$cell] <- i$bp_vt_maxmagn
  rst_maxtime[i$cell] <- i$bp_vt_maxtime
  
  # reject pixels with low magnitude
  int_lowmagn <- which(rst_maxmagn[] < .1 & rst_maxmagn[] > -.1)
  rst_bp[int_lowmagn] <- rst_maxmagn[int_lowmagn] <- rst_maxtime[int_lowmagn] <- NA
  
  # number of breakpoints
  cols <- colorRampPalette(brewer.pal(5, "RdPu"))
  p_bp <- spplot(rst_bp, scales = list(draw = TRUE), at = 0.5:5.5, 
                 col.regions = cols(10), colorkey = list(space = "top"))

  # maximum magnitude
  cols <- colorRampPalette(rev(brewer.pal(6, "RdBu")))
  p_maxmagn <- spplot(rst_maxmagn, scales = list(draw = TRUE), 
                      at = seq(-.6, .6, .1), col.regions = cols(12), 
                      colorkey = list(space = "top"))
  
  # timing of maximum magnitude
  col_month_max <- data.frame(cell = 1:10, 
                              h = -360/12 + 1:10 * 360/12, 
                              l = 80,
                              c = 65)
  
  p_maxtime <- spplot(rst_maxtime, scales = list(draw = TRUE), at = 2003:2012, 
                      col.regions = hcl(h = col_month_max$h, c = 70, l = 65), 
                      colorkey = list(space = "top"))
  
  # return artworks
  return(list(p_bp, p_maxmagn, p_maxtime))
}

p_bp <- latticeCombineGrid(lapply(lst_p, "[[", 1), layout = c(1, 5))
p_maxmagn <- latticeCombineGrid(lapply(lst_p, "[[", 2), layout = c(1, 5))
p_maxtime <- latticeCombineGrid(lapply(lst_p, "[[", 3), layout = c(1, 5))
