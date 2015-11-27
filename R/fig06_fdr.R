################################################################################
### environmental stuff

## clear workspace
rm(list = ls(all = TRUE))

## packages
lib <- c("Rsenal", "foreach", "RColorBrewer", "rasterVis")
Orcs::loadPkgs(lib)

## functions
source("R/panel.smoothconts.R")

## folders
ch_dir_outdata <- "/media/fdetsch/XChange/kilimanjaro/ndvi_comparison/out/"


################################################################################
### data processing

## loop over significance levels
lst_p_fdr <- lapply(c(.05, .001), function(p_value) {
  
  # status message
  cat("Processing significance level p =", p_value, "\n")
  
  # products
  products <- c("GIMMS3g", 
                "MOD13Q1.006", "MYD13Q1.006")
  
  if (p_value == .001)
    products <- products[-1]
  
  # import mann-kendall layers (p < 0.05)
  pattern <- paste0("0312_tau", substr(p_value, 3, nchar(p_value)), ".tif$")
  
  fls_mk <- list.files(ch_dir_outdata, pattern = pattern, 
                       full.names = TRUE)[c(1, 2, 4, 3, 5)]
  
  id_prod <- sapply(products, function(i) grep(i, fls_mk))
  fls_mk <- fls_mk[id_prod]
  rst_mk <- lapply(fls_mk, raster)
  
  # number of modis cells
  n <- ncell(rst_mk[[2]])
  
  # compute fdr for all values of prevalence...
  lst_fdr <- lapply(seq(.295, .005, -.001), function(prevalence) {
    # ...and sensitivity
    sapply(seq(.505, .995, .01), function(sensitivity) {
      statFDR(n, prevalence, sensitivity, p_value)
    })
  })    
  
  # insert values
  mat_fdr <- do.call("rbind", lst_fdr)
  rst_fdr <- raster(mat_fdr, 
                    xmn = .5, xmx = 1, 
                    ymn = 0, ymx = .3) 
  
  # create levelplot
  p_heatmap <- 
    levelplot(rst_fdr, margin = FALSE, col.regions = envinmrPalette(501),
              colorkey = list(space = "top", width = .6, height = .5), 
              main = list("FDR", hjust = -.5, cex = .9), at = seq(0, 1, .002),
              xlab = list("Power", cex = .8, vjust = -.5), 
              ylab = list("P", cex = .8),
              scales = list(cex = .7))
  
  # create contourplot
  rst_fdr_flipped <- flip(rst_fdr, "y")
  x <- coordinates(rst_fdr_flipped)[, 1]
  y <- coordinates(rst_fdr_flipped)[, 2]
  z <- rst_fdr_flipped[]
  
  p_contours <- if (p_value == .05) {
    levelplot(z ~ x * y, colorkey = FALSE, 
              panel = function(...) {
                panel.smoothconts(zlevs.conts = seq(.1, .9, .1), 
                                  # labels = c(0.2, "", 0.4, "", 0.6, "", 0.8, ""), 
                                  labels = c(seq(.1, .7, .1), rep("", 2)),
                                  col = "black", cex = 1.8, labcex = .8, method = "edge", ...)
              }, scales = list(cex = .7))
  } else {
    levelplot(z ~ x * y, colorkey = FALSE, 
              panel = function(...) {
                panel.smoothconts(zlevs.conts = seq(.01, .13, .02), 
                                  # labels = c(0.2, "", 0.4, "", 0.6, "", 0.8, ""), 
                                  labels = c(seq(.01, .09, .02), rep("", 2)),
                                  col = "black", cex = 1.8, labcex = .8, method = "edge", ...)
              }, scales = list(cex = .7))
  }
  
  # combine 
  p <- p_heatmap + 
    latticeExtra::as.layer(p_contours) + 
    latticeExtra::layer(sp.text(txt = ifelse(p_value == .05, 
                                             "a) p < 0.05", 
                                             "b) p < 0.001"), 
                                loc = c(.75, .275), font = 2, cex = .8), 
                        data = list(p_value = p_value))
  
  return(p)
})

## combine figures
p_fdr <- latticeCombineGrid(lst_p_fdr, layout = c(2, 1))
p_fdr <- update(p_fdr, scales = list(alternating = 1))

## in-text png version
png(paste0(ch_dir_outdata, "figure06.png"), width = 15, height = 10, 
    units = "cm", res = 500)
plot.new()    
print(p_fdr, newpage = FALSE)
dev.off()

## stand-alone tiff version
tiff(paste0(ch_dir_outdata, "figure06.tiff"), width = 15, height = 10, 
     units = "cm", res = 500, compression = "lzw")
plot.new()    
print(p_fdr, newpage = FALSE)
dev.off()
