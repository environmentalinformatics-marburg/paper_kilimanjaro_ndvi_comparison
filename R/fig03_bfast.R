## packages
library(raster)
library(bfast)
library(doParallel)

## parallelization
cl <- makeCluster(3)
registerDoParallel(cl)

st_year <- "2003"
nd_year <- "2012"

fls_scl <- list.files(paste0(ch_dir_extdata, "MOD13Q1.005"), recursive = TRUE,
                      pattern = "^SCL_.*.tif", full.names = TRUE)

st <- grep(st_year, fls_scl)[1]
nd <- grep(nd_year, fls_scl)[length(grep(nd_year, fls_scl))]
fls_scl <- fls_scl[st:nd]
rst_ndvi <- stack(fls_scl)
mat_ndvi <- as.matrix(rst_ndvi)

?bfast

df_bp <- foreach(i = 1:nrow(mat_ndvi), .packages = "bfast", 
                 .combine = "rbind") %dopar% {
                   
  # bfast                 
  num_val <- mat_ndvi[i, ]
  ts_val <- ts(num_val, start = c(2003, 1), end = c(2012, 12), frequency = 12)
  bf_val <- bfast(ts_val, season = "harmonic", max.iter = 10)
  
  # trend (vt) breakpoints
  num_bp_vt <- ifelse(!bf_val$nobp$Vt, 
                      length(bf_val$output[[length(bf_val$output)]]$Vt.bp), 
                      NA)
  
  # seasonal (wt) breakpoints
  num_bp_wt <- ifelse(!bf_val$nobp$Wt, 
                      length(bf_val$output[[length(bf_val$output)]]$Wt.bp), 
                      NA)
  
  ## return data
  data.frame(cell = i, bp_vt = num_bp_vt, bp_wt = num_bp_wt, row.names = "")
}