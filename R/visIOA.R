visIOA <- function(x, y, 
                   cores = 1, 
                   sp.layout = NULL,
                   at = seq(0, 1, .1),
                   col.regions = NULL,
                   alpha.regions = 1,
                   crs = NULL, 
                   xlab = "x", ylab = "y",
                   filename = "", ...) {
  
  ## packages
  lib <- c("Rsenal", "RColorBrewer", "reshape2", "plyr")
  jnk <- sapply(lib, function(x) library(x, character.only = TRUE))
  
  ## resample data if x and y have different spatial resolution
  if (any(res(x) != res(y))) {
    
    # initialize parallel backend
    library(doParallel)
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    
    # template used for value extraction
    spy_template <- rasterToPolygons(y[[1]])
    
    # value extraction
    ls_vals <- foreach(i = 1:nlayers(x), 
                       .packages = c("raster", "rgdal")) %dopar% {
      val <- extract(x[[i]], spy_template)
      val_stats <- lapply(val, function(j) {
        mn <- mean(j, na.rm = TRUE)
        qn <- quantile(j, probs = c(.1, .5, .9), na.rm = TRUE)
        matrix(c(mn, qn), ncol = 4, byrow = TRUE)
      })
      val_stats <- do.call("rbind", val_stats)
      return(val_stats)
    }
    
    # checkout parallel backend
    stopCluster(cl)
  }
  
  # Colors
  if (is.null(col.regions))
    col.regions <- colorRampPalette(brewer.pal(11, "BrBG"))
  
  if (file.exists(filename)) {
    cat("File found.\n")
    rst_ioa <- raster(filename)
  } else {
    cat("File not found.\n")
    
    # index of association (ioa)
    mat_x <- getValues(x)
    mat_x <- data.frame(cell = 1:ncell(x), mat_x)
    names(mat_x)[2:ncol(mat_x)] <- formatC(1:(ncol(mat_x)-1), width = 3, flag = 0)
    mat_x <- melt(mat_x, id.vars = 1)

    mat_y <- getValues(y)
    mat_y <- data.frame(cell = 1:ncell(y), mat_y)
    names(mat_y)[2:ncol(mat_y)] <- formatC(1:(ncol(mat_y)-1), width = 3, flag = 0)
    mat_y <- melt(mat_y, id.vars = 1)
    
    mat_mrg <- data.frame(mat_x, value_y = mat_y[, 3])
    
    num_ioa <- ddply(mat_mrg, .(mat_mrg$cell), 
                     summarise, ioa = round(ioa(value, value_y), 2))
    
    # insert values
    rst_ioa <- x[[1]]
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
  
  # Plotting
  p <- spplot(rst_ioa, scales = list(draw = TRUE), xlab = xlab, ylab = ylab, 
              col.regions = col.regions, 
              sp.layout = sp.layout, 
              par.settings = list(fontsize = list(text = 15)), 
              at = at, alpha.regions = alpha.regions)
  
  return(p)
}