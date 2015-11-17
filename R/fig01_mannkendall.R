### environmental stuff

## clear workspace
rm(list = ls(all = TRUE))

## packages
lib <- c("grid", "Rsenal", "foreach", "latticeExtra", "ggplot2", 
         "Orcs", "stargazer")
Orcs::loadPkgs(lib)

## functions
source("R/visKili.R")
source("R/visDEM.R")
source("R/visMannKendall.R")
source("R/visDensity.R")

## folders
ch_dir_extdata <- "/media/fdetsch/XChange/kilimanjaro/ndvi_comparison/data/rst/"
ch_dir_outdata <- "/media/fdetsch/XChange/kilimanjaro/ndvi_comparison/out/"

### data processing

## digital elevation model (dem)
ch_fls_dem <- paste0(ch_dir_extdata, "../dem/DEM_ARC1960_30m_Hemp.tif")
rst_dem <- raster(ch_fls_dem)
rst_dem <- aggregate(rst_dem, fact = 10)
rst_dem <- projectRaster(rst_dem, crs = "+init=epsg:4326")
p_dem <- visDEM(rst_dem, labcex = .6, cex = 1.6, col = "black")

## reference extent
fls_ndvi <- paste0(ch_dir_extdata, "MOD13Q1.006/whittaker_dsn/DSN_SCL_MVC_200301.tif")

rst_ndvi <- raster(fls_ndvi[1])
rst_ndvi <- projectRaster(rst_ndvi, crs = "+init=epsg:4326")
rst_ndvi <- trim(rst_ndvi)

num_xmin <- xmin(rst_ndvi)
num_xmax <- xmax(rst_ndvi)
num_ymin <- ymin(rst_ndvi)
num_ymax <- ymax(rst_ndvi)

## gimms grid
rst_gimms <- raster(paste0(ch_dir_outdata, "/GIMMS3g_mk_0312_tau.tif"))
rst_gimms <- projectRaster(rst_gimms, crs = "+init=epsg:4326")
spy_gimms <- rasterToPolygons(rst_gimms)

## study area
rst_kili <- kiliAerial(upperLeft = c(num_ymax, num_xmin), 
                       lowerRight = c(num_ymin, num_xmax),
                       minNumTiles = 12L, projection = "+init=epsg:4326")

# create figure
p_bing <- spplot(rst_kili[[1]], col.regions = NA, colorkey = FALSE, 
                 sp.layout = list(rgb2spLayout(rst_kili, quantiles = c(.005, .9775)), 
                                  list("sp.text", loc = c(37.04, -2.86), 
                                       txt = "a)", font = 2, cex = .6, 
                                       adj = c(.1, 1), col = "grey90"), 
                                  list("sp.text", loc = c(37.6, -3.4), 
                                       txt = "\uA9 OpenStreetMap contributors", 
                                       font = 2, cex = .4, col = "grey90")),
                 maxpixels = ncell(rst_kili), xlim = c(num_xmin, num_xmax), 
                 ylim = c(num_ymin, num_ymax), 
                 scales = list(draw = TRUE, cex = .5, 
                               y = list(at = seq(-2.9, -3.3, -.2))))

## topographic map
p_topo <- visKili(cex = .5, lwd = .05, ext = rst_ndvi)

## mann-kendall trend tests (2003-2012; p < 0.05)
st_year <- "2003"
nd_year <- "2012"

products <- list("GIMMS3g", 
                 "MOD13Q1.005", "MYD13Q1.005", 
                 "MOD13Q1.006", "MYD13Q1.006")

# ## breaks (works only when mann-kendall trend layers already exist)
# fls_mk <- list.files(ch_dir_outdata, pattern = "mk_0312_tau05", 
#                      full.names = TRUE)
# lst_mk <- lapply(fls_mk, raster)
# 
# sapply(lst_mk, function(i) {
#   cat("Minimum:", minValue(i), "\tMaximum:", maxValue(i), "\n")
#   return(invisible(NULL))
# })

## create and visualize mann-kendall trend layers
lst_p_mk <- lapply(c(.05, .001), function(p_value) {
  # status message
  cat("Processing p =", p_value, "\n")
  
  if (p_value == 0.05) {
    labels <- list(expression(bold("b) NDVI"["3g"])), 
                   expression(bold("c) NDVI"["Terra-C5"])), 
                   expression(bold("d) NDVI"["Aqua-C5"])), 
                   expression(bold("e) NDVI"["Terra-C6"])), 
                   expression(bold("f) NDVI"["Aqua-C6"])))
  } else {
    labels <- list(expression(bold("x) NDVI"["3g"])), 
                   expression(bold("a) NDVI"["Terra-C5"])), 
                   expression(bold("b) NDVI"["Aqua-C5"])), 
                   expression(bold("c) NDVI"["Terra-C6"])), 
                   expression(bold("d) NDVI"["Aqua-C6"])))
  }
  
  foreach(i = products, txt = labels, 
          .export = ls(envir = environment())) %do% {
    
    fls_ndvi <- list.files(paste0(ch_dir_extdata, i), recursive = TRUE,
                           pattern = "^DSN_.*.tif$", full.names = TRUE)
    
    st <- grep(st_year, fls_ndvi)[1]
    nd <- grep(nd_year, fls_ndvi)[length(grep(nd_year, fls_ndvi))]
    
    fls_ndvi <- fls_ndvi[st:nd]
    rst_ndvi <- stack(fls_ndvi)
    
    p <- visMannKendall(rst = rst_ndvi, dem = rst_dem, max_ele = 4000,
                        
                        xlab = "", ylab = "",
                        p_value = p_value, crs = "+init=epsg:4326",
                        filename = paste0(ch_dir_outdata, i, "_mk_0312.tif"), 
                        at = seq(-.55, .55, .01), cores = 3L, 
                        format = "GTiff", overwrite = TRUE, 
                        xlim = c(num_xmin, num_xmax), 
                        ylim = c(num_ymin, num_ymax), 
                        scales = list(draw = TRUE, cex = .5, 
                                      y = list(at = seq(-2.9, -3.3, -.2))), 
                        rewrite = FALSE)
   
    # add contour lines and text 
    p <- p + 
      latticeExtra::layer(sp.polygons(spy_gimms, lty = 3, col = "grey75"), 
                          data = list(i = i)) + 
      latticeExtra::layer(sp.text(loc = c(37.04, -2.86), txt = txt, font = 2, 
                                  cex = .6, adj = c(.1, 1), col = "black"), 
                          data = list(txt = txt)) 
    
    p <- envinmrRasterPlot(p, rot = 90, height = .5, width = .4, key.cex = .7)
    
    return(p)
  }
})
  
  
################################################################################
### visualization ##############################################################
################################################################################

## combination final figure, p < 0.05
p_mk_comb <- latticeCombineGrid(append(list(p_bing), lst_p_mk[[1]]), 
                                layout = c(2, 3))

p_mk_comb <- p_mk_comb + 
  latticeExtra::as.layer(p_dem)

## density plot
p_dens <- visDensity(p = .05, dsn = ch_dir_outdata, combined = FALSE)

# in-text png version
ch_fls_png <- paste0(ch_dir_outdata, "figure01.png")
png(ch_fls_png, width = 20.5, height = 25, units = "cm", res = 500)
plot.new()

print(p_mk_comb, newpage = FALSE)

# add key
vp_key <- viewport(x = .5, y = .905,
                   height = 0.1, width = .8,
                   just = c("center", "bottom"),
                   name = "key.vp")
pushViewport(vp_key)
draw.colorkey(key = list(col = colorRampPalette(brewer.pal(11, "BrBG")), 
                         width = .6, height = .5,
                         at = seq(-.55, .55, .01), 
                         space = "bottom"), draw = TRUE)
grid.text(expression("Kendall's " ~ tau), x = 0.5, y = .9, just = c("centre", "top"), 
          gp = gpar(font = 2, cex = .85))

# add density plot #1
upViewport(n = 0)
downViewport(trellis.vpname("figure"))
vp_dens1 <- viewport(x = .4925, y = .340, height = .13, width = .2, 
                     just = c("center", "bottom"), name = "dens1.vp")
pushViewport(vp_dens1)
print(p_dens[[1]], newpage = FALSE)
grid.rect(gp = gpar(col = "black", fill = "transparent", cex = 1.1))

# add density plot #2
upViewport()
vp_dens2 <- viewport(x = .4925, y = .004, height = .13, width = .2, 
                     just = c("center", "bottom"), name = "dens2.vp")
pushViewport(vp_dens2)
print(p_dens[[2]], newpage = FALSE)
grid.rect(gp = gpar(col = "black", fill = "transparent", cex = 1.1))

# add topographic map
upViewport(n = 0)
vp_rect <- viewport(x = .365, y = .6995, height = .315, width = .15, 
                    just = c("left", "bottom"))
pushViewport(vp_rect)
print(p_topo, newpage = FALSE)

# add equator label
downViewport(trellis.vpname("figure"))
grid.text(x = .05, y = .38, just = c("left", "bottom"), label = "Eq.", 
          gp = gpar(cex = .3))

dev.off()

# standalone tiff version
ch_fls_tif <- paste0(ch_dir_outdata, "figure01.tiff")
tiff(ch_fls_tif, width = 20.5, height = 25, units = "cm", res = 500, 
     compression = "lzw")
plot.new()

print(p_mk_comb, newpage = FALSE)

# add key
vp_key <- viewport(x = .5, y = .905,
                   height = 0.1, width = .8,
                   just = c("center", "bottom"),
                   name = "key.vp")
pushViewport(vp_key)
draw.colorkey(key = list(col = colorRampPalette(brewer.pal(11, "BrBG")), 
                         width = .6, height = .5,
                         at = seq(-.55, .55, .01), 
                         space = "bottom"), draw = TRUE)
grid.text(expression("Kendall's " ~ tau), x = 0.5, y = .9, just = c("centre", "top"), 
          gp = gpar(font = 2, cex = .85))

# add density plot #1
upViewport(n = 0)
downViewport(trellis.vpname("figure"))
vp_dens1 <- viewport(x = .4925, y = .340, height = .13, width = .2, 
                     just = c("center", "bottom"), name = "dens1.vp")
pushViewport(vp_dens1)
print(p_dens[[1]], newpage = FALSE)
grid.rect(gp = gpar(col = "black", fill = "transparent", cex = 1.1))

# add density plot #2
upViewport()
vp_dens2 <- viewport(x = .4925, y = .004, height = .13, width = .2, 
                     just = c("center", "bottom"), name = "dens2.vp")
pushViewport(vp_dens2)
print(p_dens[[2]], newpage = FALSE)
grid.rect(gp = gpar(col = "black", fill = "transparent", cex = 1.1))

# add topographic map
upViewport(n = 0)
vp_rect <- viewport(x = .365, y = .6995, height = .315, width = .15, 
                    just = c("left", "bottom"))
pushViewport(vp_rect)
print(p_topo, newpage = FALSE)

# add equator label
downViewport(trellis.vpname("figure"))
grid.text(x = .05, y = .38, just = c("left", "bottom"), label = "Eq.", 
          gp = gpar(cex = .3))

dev.off()

################################################################################
## combination final figure, p < 0.001
################################################################################
p_mk_comb <- latticeCombineGrid(lst_p_mk[[2]][2:5], 
                                layout = c(2, 2))

p_mk_comb <- p_mk_comb + 
  latticeExtra::as.layer(p_dem)

p_dens <- visDensity(p = 0.001, dsn = ch_dir_outdata, combined = FALSE)

# in-text png version
ch_fls_png <- paste0(ch_dir_outdata, "figure03.png")
png(ch_fls_png, width = 20.5, height = 18, units = "cm", res = 500)
plot.new()

print(p_mk_comb, newpage = FALSE)

# add key caption
grid.text(bquote(bold("Kendall's " ~ tau)), x = 0.5, y = 1, 
          just = c("centre", "top"), gp = gpar(font = 2, cex = .85))

# add density plot #1
upViewport(n = 0)
downViewport(trellis.vpname("figure"))
vp_dens1 <- viewport(x = .4925, y = .341 * 3/2, height = .195, width = .2, 
                     just = c("center", "bottom"), name = "dens1.vp")
pushViewport(vp_dens1)
print(p_dens[[1]], newpage = FALSE)
grid.rect(gp = gpar(col = "black", fill = "transparent", cex = 1.1))

# add density plot #2
upViewport()
vp_dens2 <- viewport(x = .4925, y = .004 * 3/2, height = .195, width = .2, 
                     just = c("center", "bottom"), name = "dens2.vp")
pushViewport(vp_dens2)
print(p_dens[[2]], newpage = FALSE)
grid.rect(gp = gpar(col = "black", fill = "transparent", cex = 1.1))

dev.off()

# standalone tiff version
ch_fls_tif <- paste0(ch_dir_outdata, "figure03.tiff")
tiff(ch_fls_tif, width = 20.5, height = 18, units = "cm", res = 500, 
     compression = "lzw")
plot.new()

print(p_mk_comb, newpage = FALSE)

# add key caption
grid.text(bquote(bold("Kendall's " ~ tau)), x = 0.5, y = 1, 
          just = c("centre", "top"), gp = gpar(font = 2, cex = .85))

# add density plot #1
upViewport(n = 0)
downViewport(trellis.vpname("figure"))
vp_dens1 <- viewport(x = .4925, y = .341 * 3/2, height = .195, width = .2, 
                     just = c("center", "bottom"), name = "dens1.vp")
pushViewport(vp_dens1)
print(p_dens[[1]], newpage = FALSE)
grid.rect(gp = gpar(col = "black", fill = "transparent", cex = 1.1))

# add density plot #2
upViewport()
vp_dens2 <- viewport(x = .4925, y = .004 * 3/2, height = .195, width = .2, 
                     just = c("center", "bottom"), name = "dens2.vp")
pushViewport(vp_dens2)
print(p_dens[[2]], newpage = FALSE)
grid.rect(gp = gpar(col = "black", fill = "transparent", cex = 1.1))

dev.off()

