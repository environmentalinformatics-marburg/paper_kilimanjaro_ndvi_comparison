visMannKendall <- function(rst, dem, max_ele,
                           p_value = NULL, 
                           sp.layout = NULL,
                           at = seq(-1, 1, .2),
                           col.regions = NULL,
                           alpha.regions = 1,
                           keycex = .8,
                           crs = NULL, 
                           xlab = "x", ylab = "y", scales = list(draw = TRUE),
                           xlim = bbox(rst)[1, ], ylim = bbox(rst)[2, ],
                           filename = "", rewrite = FALSE, cores = 1L, ...) {
  
  ## parallelize
  if (cores > 1L) {
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
  }
    
  # output files
  file_tau <- paste0(dirname(filename), Orcs::pureBasename(filename, TRUE), "_tau.tif")
  file_p <- paste0(dirname(filename), Orcs::pureBasename(filename, TRUE), "_p.tif")
  
  if (file.exists(file_tau) & file.exists(file_p) & !rewrite) {
    cat("Files found, working with existing files.\n")
    ndvi.mk <- raster::stack(file_tau, file_p)
  } else {
    cat("Files in, start processing.\n")
    
    # Mann-Kendall tau
    mat <- raster::as.matrix(rst)
    mat_mk <- foreach::foreach(i = 1:nrow(mat), .packages = "gimms", 
                               .combine = "rbind") %dopar% {
                                 gimms::significantTau(mat[i, ], df = TRUE, 
                                                       prewhitening = TRUE, 
                                                       conf.intervals = FALSE)
                               }
    
    if (nrow(mat_mk) != raster::ncell(rst))
      stop("Input raster and output matrix are of different size!\n")
    
    # raster templates
    rst_tau <- rst_p <- rst[[1]]
    rst_tau[] <- rst_p[] <- NA
    
    # insert values
    rst_tau[] <- mat_mk[, 1]
    rst_p[] <- mat_mk[, 2]
    
    # reject values above user-defined elevation level
    if (!missing(dem) & !missing(max_ele)) {
      rst_tau <- overlay(rst_tau, dem, fun = function(x, y) {
        x[y[] > max_ele] <- NA
        return(x)
      })
      rst_p <- overlay(rst_p, dem, fun = function(x, y) {
        x[y[] > max_ele] <- NA
        return(x)
      })    
    }

    rst_out <- raster::stack(rst_tau, rst_p)
    
    # write to file
    lst_out <- foreach::foreach(i = 1:2, j = c(file_tau, file_p)) %do% {
      raster::writeRaster(rst_out[[i]], filename = j, ...)
    }
    
    ndvi.mk <- raster::stack(lst_out)
  }

  
  # reject non-significant pixels (optional)    
  if (!is.null(p_value)) {
    file_sign <- paste0(dirname(filename), Orcs::pureBasename(filename, TRUE), 
                        "_tau", substr(p_value, 3, nchar(p_value)))
    ndvi.mk <- raster::overlay(ndvi.mk[[1]], ndvi.mk[[2]], 
                               fun = function(x, y) {
                                 x[y[] >= p_value] <- NA
                                 return(x)
                               }, filename = file_sign, 
                               format = "GTiff", overwrite = TRUE)
  }
  
  # reproject (optional)
  if (!is.null(crs)) {
    ndvi.mk <- raster::projectRaster(ndvi.mk, crs = crs, method = "ngb")
    if (!all(is.na(ndvi.mk[])))
      ndvi.mk <- raster::trim(ndvi.mk)
  }
  
  ## deregister parallel backend
  if (cores > 1L) {
    parallel::stopCluster(cl)
  }
  
  ## create figure
  
  # colors
  if (is.null(col.regions))
    col.regions <- colorRampPalette(RColorBrewer::brewer.pal(11, "BrBG"))
  
  # only missing values
  if (all(is.na(ndvi.mk[]))) {
    
    ndvi.mk[] <- 0
    sp::spplot(ndvi.mk, col.regions = "transparent", scales = scales, 
               colorkey = FALSE, xlab = xlab, ylab = ylab, 
               xlim = xlim, ylim = ylim) + 
      latticeExtra::layer(sp::sp.polygons(raster::rasterToPolygons(ndvi.mk), 
                                          lty = "dotted", col = "grey40"), 
                          data = list(ndvi.mk = ndvi.mk)) + 
      latticeExtra::layer(sp::sp.text(loc = c(37.04, -3.35), 
                                      txt = expression(bold("b) NDVI"["3g"])), 
                                      font = 2, cex = .6, adj = c(.1, 1)))
    
  # valid values  
  } else {
    
    sp::spplot(ndvi.mk, scales = scales, xlab = xlab, ylab = ylab, 
               col.regions = col.regions, 
               sp.layout = sp.layout, xlim = xlim, ylim = ylim,
               par.settings = list(fontsize = list(text = 15)), 
               at = at, alpha.regions = alpha.regions, 
               colorkey = list(labels = list(cex = keycex)))
  }
  
}