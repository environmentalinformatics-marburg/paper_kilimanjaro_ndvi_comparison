### environmental stuff

## packages
lib <- c("grid", "Rsenal", "doParallel", "latticeExtra", "ggplot2", 
         "OpenStreetMap", "Orcs")
jnk <- sapply(lib, function(x) library(x, character.only = TRUE))

## functions
source("R/visKili.R")
source("R/visDEM.R")
source("R/visMannKendall.R")

## parallelization
cl <- makeCluster(3)
registerDoParallel(cl)

## folders
ch_dir_extdata <- "/media/fdetsch/XChange/kilimanjaro/ndvi_comparison/data/rst/"
ch_dir_outdata <- "/media/fdetsch/XChange/kilimanjaro/ndvi_comparison/out/"

### data processing

## digital elevation model (dem)
ch_fls_dem <- paste0(ch_dir_extdata, "../dem/DEM_ARC1960_30m_Hemp.tif")
rst_dem <- raster(ch_fls_dem)
rst_dem <- aggregate(rst_dem, fact = 10)
rst_dem <- projectRaster(rst_dem, crs = "+init=epsg:4326")
p_dem <- visDEM(rst_dem, labcex = .6, cex = 1.6)

## reference extent
fls_ndvi <- paste0(ch_dir_extdata, "MOD13Q1.006/whittaker_dsn/DSN_SCL_MVC_200301.tif")

rst_ndvi <- raster(fls_ndvi[1])
rst_ndvi <- projectRaster(rst_ndvi, crs = "+init=epsg:4326")
rst_ndvi <- trim(rst_ndvi)

num_xmin <- xmin(rst_ndvi)
num_xmax <- xmax(rst_ndvi)
num_ymin <- ymin(rst_ndvi)
num_ymax <- ymax(rst_ndvi)

## mann-kendall trend tests (2003-2012; p < 0.05)
st_year <- "2003"
nd_year <- "2012"

products <- c("GIMMS3g", 
              "MOD13Q1.005", "MYD13Q1.005", 
              "MOD13Q1.006", "MYD13Q1.006")

label <- c("b)", "c)", "d)", "e)", "f)")
label <- sapply(1:length(label), function(i) paste(label[i], products[i]))

p_mk <- foreach(i = products, txt = label, .packages = lib) %dopar% {
                  
  fls_ndvi <- list.files(paste0(ch_dir_extdata, i), recursive = TRUE,
                         pattern = "^DSN_.*.tif$", full.names = TRUE)
  
  st <- grep(st_year, fls_ndvi)[1]
  nd <- grep(nd_year, fls_ndvi)[length(grep(nd_year, fls_ndvi))]
  
  fls_ndvi <- fls_ndvi[st:nd]
  rst_ndvi <- stack(fls_ndvi)
  
  p <- visMannKendall(rst = rst_ndvi, xlab = "", ylab = "",
                      p_value = .05, crs = "+init=epsg:4326",
                      filename = paste0(ch_dir_outdata, i, "_mk05_0312.tif"), 
                      format = "GTiff", overwrite = TRUE, 
                      sp.layout = list("sp.text", loc = c(37.04, -3.35), 
                                       txt = txt, font = 2, cex = .6, adj = c(.1, 1)), 
                      xlim = c(num_xmin, num_xmax), 
                      ylim = c(num_ymin, num_ymax), 
                      scales = list(draw = TRUE, cex = .5, 
                                    y = list(at = seq(-2.9, -3.3, -.2))))
  
  p <- p + as.layer(p_dem)
  p <- envinmrRasterPlot(p, rot = 90, height = .5, width = .4, key.cex = .7)
  
  return(p)
}

## statistics
fls_mk05 <- list.files(ch_dir_outdata, pattern = "mk05_0312.tif$", 
                        full.names = TRUE)
rst_mk05 <- lapply(fls_mk05, raster)

df_mk_stats <- foreach(i = rst_mk05, .combine = "rbind") %do% mkStats(i)
df_mk_stats <- cbind(product = products, df_mk_stats)


# mean difference
meanDifference(rst_mk05[[2]], rst_mk05[[4]]) # MOD13Q1.005 vs. MYD13Q1.005
meanDifference(rst_mk05[[3]], rst_mk05[[5]]) # MOD13Q1.006 vs. MYD13Q1.006

## study area
osm_kili <- openproj(openmap(upperLeft = c(num_ymax, num_xmin), 
                             lowerRight = c(num_ymin, num_xmax), 
                             type = "bing", minNumTiles = 12L), 
                     projection = "+init=epsg:4326")

# visualize bing image
rst_kili <- raster(osm_kili)

p_bing <- spplot(rst_kili[[1]], col.regions = NA, colorkey = FALSE, 
                 sp.layout = list(rgb2spLayout(rst_kili), 
                                  list("sp.text", loc = c(37.04, -3.35), 
                                       txt = "a)", font = 2, cex = .6, adj = c(.1, 1))),
                 xlim = c(num_xmin, num_xmax), 
                 ylim = c(num_ymin, num_ymax), 
                 scales = list(draw = TRUE, cex = .5, 
                               y = list(at = seq(-2.9, -3.3, -.2))))

# visualize topographic map
p_topo <- visKili(cex = 1, lwd = 1)


################################################################################
### visualization ##############################################################
################################################################################

## combination final figure
p_mk_comb <- latticeCombineGrid(append(list(p_bing), p_mk), layout = c(2, 3))

# in-text png version
ch_fls_png <- paste0(ch_dir_outdata, "figure01.png")
png(ch_fls_png, width = 15.5, height = 19.5, units = "cm", res = 500)
plot.new()

print(p_mk_comb, newpage = FALSE)

# add key
vp_key <- viewport(x = .5, y = .885,
                height = 0.1, width = .8,
                just = c("center", "bottom"),
                name = "key.vp")
pushViewport(vp_key)
draw.colorkey(key = list(col = colorRampPalette(brewer.pal(11, "BrBG")), 
                         width = .6, height = .5,
                         at = seq(-1, 1, .2), 
                         space = "bottom"), draw = TRUE)
grid.text(expression("Kendall's " ~ tau), x = 0.5, y = 1, just = c("centre", "top"), 
          gp = gpar(font = 2, cex = .85))

# add topographic map
upViewport()
vp_rect <- viewport(x = .365, y = .6635, height = .315, width = .15, 
                    just = c("left", "bottom"))
pushViewport(vp_rect)
print(p_topo, newpage = FALSE)

dev.off()

# standalone tiff version
ch_fls_tif <- paste0(ch_dir_outdata, "figure01.tiff")
tiff(ch_fls_tif, width = 15.5, height = 19.5, units = "cm", res = 500, 
     compression = "lzw")
plot.new()

print(p_mk_comb, newpage = FALSE)

# add key
vp_key <- viewport(x = .5, y = .885,
                   height = 0.1, width = .8,
                   just = c("center", "bottom"),
                   name = "key.vp")
pushViewport(vp_key)
draw.colorkey(key = list(col = colorRampPalette(brewer.pal(11, "BrBG")), 
                         width = .6, height = .5,
                         at = seq(-1, 1, .2), 
                         space = "bottom"), draw = TRUE)
grid.text(expression("Kendall's " ~ tau), x = 0.5, y = 1, just = c("centre", "top"), 
          gp = gpar(font = 2, cex = .85))

# add topographic map
upViewport()
vp_rect <- viewport(x = .365, y = .6635, height = .315, width = .15, 
                    just = c("left", "bottom"))
pushViewport(vp_rect)
print(p_topo, newpage = FALSE)

dev.off()

## deregister parallel backend
stopCluster(cl)
