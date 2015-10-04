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

## start processing single products...
ch_dir_ndvi <- list.dirs(ch_dir_extdata, full.names = TRUE, recursive = FALSE)
ls_bfast <- lapply(ch_dir_ndvi, function(h) {

  # adjust search pattern (gimms required no scaling)
  pttrn <- if (basename(h) == "GIMMS3g") "^MVC_.*.tif$" else "^SCL_.*.tif$"

  # search for available files
  fls_scl <- list.files(h, recursive = TRUE,
                        pattern = pttrn, full.names = TRUE)

  # subset available files by start and end year
  st <- grep(st_year, fls_scl)[1]
  nd <- grep(nd_year, fls_scl)[length(grep(nd_year, fls_scl))]
  fls_scl <- fls_scl[st:nd]

  # import data
  rst_ndvi <- stack(fls_scl)
  mat_ndvi <- as.matrix(rst_ndvi)

  # breakpoint detection
  df_bp <- foreach(i = 1:nrow(mat_ndvi), .packages = "bfast", 
                   .combine = "rbind.fill", .export = ls(envir = globalenv())) %dopar% {
                     
    # bfast                 
    num_val <- mat_ndvi[i, ]
    ts_val <- ts(num_val, start = c(st_year, 1), end = c(nd_year, 12), frequency = 12)
    bf_val <- bfast(ts_val, season = "harmonic", max.iter = 10)
    
    # trend (vt) breakpoints
    num_bp_vt <- bf_val$output[[length(bf_val$output)]]$Vt.bp
    num_bp_vt_ct <- ifelse(!bf_val$nobp$Vt, 
                           length(num_bp_vt), 
                           NA)
    
    if (!is.na(num_bp_vt_ct)) {
      
      # times of trend breakpoints
      num_bp_vt_dt <- time(ts_val)[num_bp_vt]
      df_bp_vt_dt <- data.frame(matrix(num_bp_vt_dt, byrow = TRUE, nrow = 1))
      names(df_bp_vt_dt) <- paste0("bp_vt_", formatC(1:length(num_bp_vt_dt), width = 2, flag = "0"))
      df_bp_vt_dt <- data.frame("bp_vt_ct" = num_bp_vt_ct, df_bp_vt_dt)
      
      # time and magnitude of biggest change detected
      num_bp_vt_maxtime <- bf_val$Time
      num_bp_vt_maxtime <- time(ts_val)[num_bp_vt_maxtime]
      num_bp_vt_maxmagn <- bf_val$Magnitude
      df_bp_vt_max <- data.frame("bp_vt_maxtime" = num_bp_vt_maxtime, 
                                 "bp_vt_maxmagn" = num_bp_vt_maxmagn)
      
      df_bp_vt <- cbind(df_bp_vt_dt, df_bp_vt_max)
      
    } else {
      df_bp_vt <- data.frame(row.names = 1)
    } 
    
    # seasonal (wt) breakpoints
    num_bp_wt <- bf_val$output[[length(bf_val$output)]]$Wt.bp
    num_bp_wt_ct <- ifelse(!bf_val$nobp$Wt, 
                           length(num_bp_wt), 
                           NA)
    
    if (!is.na(num_bp_wt_ct)) {
      
      # times of seasonal breakpoints
      num_bp_wt_dt <- time(ts_val)[num_bp_wt]
      df_bp_wt_dt <- data.frame(matrix(num_bp_wt_dt, byrow = TRUE, nrow = 1))
      names(df_bp_wt_dt) <- paste0("bp_wt_", formatC(1:length(num_bp_wt_dt), width = 2, flag = "0"))
      df_bp_wt_dt <- data.frame("bp_wt_ct" = num_bp_wt_ct, df_bp_wt_dt)
      
      # time and magnitude of biggest change detected
      num_bp_wt_maxtime <- bf_val$Time
      num_bp_wt_maxtime <- time(ts_val)[num_bp_wt_maxtime]
      num_bp_wt_maxmagn <- bf_val$Magnitude
      df_bp_wt_max <- data.frame("bp_wt_maxtime" = num_bp_wt_maxtime, 
                                 "bp_wt_maxmagn" = num_bp_wt_maxmagn)
      
      df_bp_wt <- cbind(df_bp_wt_dt, df_bp_wt_max)
      
    } else {
      df_bp_wt <- data.frame(row.names = 1)
    }
    
    # return data
    data.frame(cell = i, df_bp_vt, df_bp_wt, 
               row.names = "")
  }
  
  # return concatenated data
  return(df_bp)
})