visMannKendall <- function(rst, 
                           p_value = NULL, 
                           sp.layout = NULL,
                           at = seq(-1, 1, .2),
                           col.regions = NULL,
                           alpha.regions = 1,
                           keycex = .8,
                           crs = NULL, 
                           xlab = "x", ylab = "y", scales = list(draw = TRUE),
                           xlim = bbox(rst)[1, ], ylim = bbox(rst)[2, ],
                           filename = "", rewrite = FALSE, ...) {
  
  lib <- c("raster", "Kendall", "rasterVis", "RColorBrewer")
  jnk <- sapply(lib, function(x) library(x, character.only = TRUE))
  
  # Colors
  if (is.null(col.regions))
    col.regions <- colorRampPalette(brewer.pal(11, "BrBG"))
  
  if (file.exists(filename) & !rewrite) {
    cat("Files found, working with existing files.\n")
    ndvi.mk <- raster(filename)
  } else {
    cat("Files in, start processing.\n")
    
    # Mann-Kendall tau
    if (is.null(p_value)) {
      ndvi.mk <- overlay(rst, fun = function(x) {
        Rsenal::significantTau(x, p = 1, prewhitening = TRUE)
      }, filename = filename, overwrite = TRUE, ...)
      
    } else {
      # Mann-Kendall tau of significant pixels only
      ndvi.mk <- overlay(rst, fun = function(x) {
        Rsenal::significantTau(x, p = p_value, prewhitening = TRUE)
      }, filename = filename, overwrite = TRUE, ...)
    }
  }
  
  ## status message
  cat("minValue:", minValue(ndvi.mk), "maxValue:", maxValue(ndvi.mk), "\n")
  
  if (!is.null(crs)) {
    ndvi.mk <- projectRaster(ndvi.mk, crs = crs, method = "ngb")
    ndvi.mk <- trim(ndvi.mk)
  }
  
  # Plotting
  spplot(ndvi.mk, scales = scales, xlab = xlab, ylab = ylab, 
         col.regions = col.regions, 
         sp.layout = sp.layout, xlim = xlim, ylim = ylim,
         par.settings = list(fontsize = list(text = 15)), 
         at = at, alpha.regions = alpha.regions, 
         colorkey = list(labels = list(cex = keycex)))
  
}