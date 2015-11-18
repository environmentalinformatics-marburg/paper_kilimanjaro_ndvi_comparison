visIOA <- function(x, y, 
                   cores = 1,
                   fun = c("mean", "median", "max"),
                   sp.layout = NULL, scales = list(draw = TRUE), 
                   at = seq(.2, 1, .05),
                   col.regions = NULL,
                   alpha.regions = 1,
                   crs = NULL, 
                   xlab = "x", ylab = "y",
                   xlim = bbox(rst)[1, ], ylim = bbox(rst)[2, ],
                   filename = "", ...) {
  
  ## packages
  lib <- c("Rsenal", "RColorBrewer", "reshape2", "plyr")
  jnk <- sapply(lib, function(x) library(x, character.only = TRUE))
  
  ## if file already exists, skip processing
  if (file.exists(filename)) {
    cat("File found.\n")
    rst_ioa <- raster(filename)
  } else {
    cat("File not found.\n")
    
    ## resample data if x and y have different spatial resolution
    if (any(res(x) != res(y))) {
      
      # initialize parallel backend
      library(doParallel)
      cl <- makeCluster(cores)
      registerDoParallel(cl)
      
      # template used for value extraction
      spy_template <- rasterToPolygons(y[[1]])
      
      # function to be applied
      fun <- fun[1]
      
      # value extraction
      ls_vals <- 
        foreach(i = 1:nlayers(x), .packages = c("raster", "rgdal")) %dopar% {
          val <- extract(x[[i]], spy_template)
          val_stats <- lapply(val, function(j) {
            mn <- mean(j, na.rm = TRUE)
            qn <- quantile(j, probs = c(.1, .5, .9), na.rm = TRUE)
            mx <- max(j, na.rm = TRUE)
            matrix(c(mn, qn), ncol = 5, byrow = TRUE)
          })
          val_stats <- do.call("rbind", val_stats)
          return(val_stats)
        }
      
      # function selection
      int_id_col <- if (fun == "mean") 1L else if (fun == "median") 3L else 5L
      mat_x <- do.call("cbind", lapply(ls_vals, function(i) {
        i[, int_id_col]
      }))
      mat_x <- data.frame(cell = 1:nrow(mat_x), mat_x)
      names(mat_x)[2:ncol(mat_x)] <- formatC(1:(ncol(mat_x)-1), width = 3, flag = 0)
      mat_x <- melt(mat_x, id.vars = 1)
      
      # check out parallel backend
      stopCluster(cl)
    }
    
    # index of association (ioa)
    if (!exists("mat_x")) {
      mat_x <- getValues(x)
      mat_x <- data.frame(cell = 1:ncell(x), mat_x)
      names(mat_x)[2:ncol(mat_x)] <- formatC(1:(ncol(mat_x)-1), width = 3, flag = 0)
      mat_x <- melt(mat_x, id.vars = 1)
    }
    
    mat_y <- getValues(y)
    mat_y <- data.frame(cell = 1:ncell(y), mat_y)
    names(mat_y)[2:ncol(mat_y)] <- formatC(1:(ncol(mat_y)-1), width = 3, flag = 0)
    mat_y <- melt(mat_y, id.vars = 1)
    
    mat_mrg <- data.frame(mat_x, value_y = mat_y[, 3])
    
    num_ioa <- ddply(mat_mrg, .(mat_mrg$cell), 
                     summarise, ioa = ioa(value, value_y))
    
    # insert values
    rst_ioa <- y[[1]]
    rst_ioa[] <- num_ioa[, 2]
    
    # output storage (optional)
    if (filename != "")
      rst_ioa <- writeRaster(rst_ioa, filename = filename, ...)
    
  }
  
  ## reprojection (optional)
  if (!is.null(crs)) {
    rst_ioa <- projectRaster(rst_ioa, crs = crs, method = "ngb")
    rst_ioa <- trim(rst_ioa)
  }
  
  ## colors
  if (is.null(col.regions))
    col.regions <- colorRampPalette(brewer.pal(3, "Reds"))
  
  ## visualize
  spplot(rst_ioa, scales = scales, xlab = xlab, ylab = ylab, 
         col.regions = col.regions, 
         sp.layout = sp.layout, xlim = xlim, ylim = ylim,
         par.settings = list(fontsize = list(text = 15)), 
         at = at, alpha.regions = alpha.regions)
  
}