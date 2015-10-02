## packages
lib <- c("raster", "rgdal", "OpenStreetMap")
jnk <- sapply(lib, function(x) library(x, character.only = TRUE))

## functions
source("R/visKili.R")

## bing aerial image
osm_kili <- openproj(openmap(upperLeft = c(num_ymax, num_xmin), 
                             lowerRight = c(num_ymin, num_xmax), 
                             type = "bing", minNumTiles = 12L), 
                             projection = "+init=epsg:4326")

## visualize bing image
rst_kili <- raster(osm_kili)
p_bing <- spplot(rst_kili[[1]], col.regions = NA, colorkey = FALSE, 
                 sp.layout = rgb2spLayout(rst_kili), 
                 scales = list(draw = TRUE, cex = .6, 
                               y = list(at = seq(-2.9, -3.3, -.2))))

## visualize topographic map
p_topo <- visKili(cex = 1.5)

# ## manuscript version
# png("vis/fig01__map_w|o_grid.png", units = "cm", width = 9, 
#     height = 7.2, res = 500)
# 
# # satellite image
# grid.newpage()
# print(p_bing, newpage = FALSE)
# 
# # topo map
# vp_cont <- viewport(x = .735, y = .635, just = c("left", "bottom"), 
#                     width = .3, height = .4)
# pushViewport(vp_cont)
# print(p_topo, newpage = FALSE)
# 
# dev.off()
# 
# ## standalone version
# tiff("vis/figure_01.tiff", res = 500, width = 9, height = 7.2, 
#      units = "cm", compression = "lzw")
# 
# grid.newpage()
# print(p_bing, newpage = FALSE)
# 
# vp_cont <- viewport(x = .735, y = .635, just = c("left", "bottom"), 
#                     width = .3, height = .4)
# pushViewport(vp_cont)
# print(p_topo, newpage = FALSE)
# 
# dev.off()
