### environmental stuff

## packages
lib <- c("grid", "Rsenal", "foreach", "latticeExtra", "ggplot2")
jnk <- sapply(lib, function(x) library(x, character.only = TRUE))

## functions
source("R/visDEM.R")
source("R/visIOA.R")

## folders
ch_dir_extdata <- "/media/fdetsch/XChange/kilimanjaro/ndvi_comparison/data/"
ch_dir_outdata <- "/media/fdetsch/XChange/kilimanjaro/ndvi_comparison/out/"

### data processing

## digital elevation model (dem)
ch_fls_dem <- paste0(ch_dir_extdata, "dem/DEM_ARC1960_30m_Hemp.tif")
rst_dem <- raster(ch_fls_dem)
rst_dem <- aggregate(rst_dem, fact = 10)
rst_dem <- projectRaster(rst_dem, crs = "+init=epsg:4326")
p_dem <- visDEM(rst_dem, labcex = .8, cex = 1.6)

## reference extent
fls_ndvi <- list.files(paste0(ch_dir_extdata, "mod13q1"), 
                       pattern = ".tif$", full.names = TRUE)

rst_ndvi <- raster(fls_ndvi[1])
rst_ndvi <- projectRaster(rst_ndvi, crs = "+init=epsg:4326")
rst_ndvi <- trim(rst_ndvi)

num_xmin <- xmin(rst_ndvi)
num_xmax <- xmax(rst_ndvi)
num_ymin <- ymin(rst_ndvi)
num_ymax <- ymax(rst_ndvi)

## index of association (ioa; 2003-2013)
st_year <- "2003"
nd_year <- "2012"

ls_ndvi <- lapply(c("mod13q1", "myd13q1", "gimms"), function(i) {
                  
  fls_ndvi <- list.files(paste0(ch_dir_extdata, "ioa/", i), 
                         pattern = ".tif$", full.names = TRUE)
  
  st <- grep(st_year, fls_ndvi)[1]
  nd <- grep(nd_year, fls_ndvi)[length(grep(nd_year, fls_ndvi))]
  
  fls_ndvi <- fls_ndvi[st:nd]
  rst_ndvi <- stack(fls_ndvi)
  
  if (i == "gimms") {
    spy <- rasterToPolygons(rst_ndvi[[1]])
    rst_ndvi <- crop(rst_ndvi, spy)
  }
  
  return(rst_ndvi)
})

p_ioa <- foreach(i = list(ls_ndvi[[1]], ls_ndvi[[1]], ls_ndvi[[2]]), 
                 j = list(ls_ndvi[[2]], ls_ndvi[[3]], ls_ndvi[[3]]), 
                 k = list("a)", "b)", "c)"), fun = list("median", "median", "median"),
                 l = list("terra|aqua", "terra|gimms", "aqua|gimms")) %do% {
          
  p <- visIOA(i, j, cores = 3, crs = "+init=epsg:4326", xlab = "", ylab = "",
              sp.layout = list("sp.text", loc = c(37.04, -3.35), 
                               txt = k, font = 2, cex = .7), 
              fun = fun,
              xlim = c(num_xmin, num_xmax), 
              ylim = c(num_ymin, num_ymax), 
              filename = paste0(ch_dir_outdata, l, "_ioa_", fun, "_0312.tif"), 
              scales = list(draw = TRUE, cex = .6, 
                            y = list(at = seq(-2.9, -3.3, -.2))))

  p <- p + as.layer(p_dem)
  p <- envinmrRasterPlot(p, rot = 90, height = .5, width = .4, key.cex = .7)
  
  return(p)
}

## statistics
fls_ioa <- list.files(ch_dir_outdata, pattern = "ioa_median_0312.tif$", 
                        full.names = TRUE)
rst_ioa <- lapply(fls_ioa, raster)

# terra vs. aqua
val_ioa_modis <- rst_ioa[[2]][]
val_ioa_modis_mu <- mean(val_ioa_modis)
val_ioa_modis_quan <- quantile(val_ioa_modis, probs = seq(.25, .75, .25))

# terra vs. gimms
val_ioa_terra_gimms <- rst_ioa[[3]][]
val_ioa_terra_gimms_mu <- mean(val_ioa_terra_gimms)
val_ioa_terra_gimms_quan <- quantile(val_ioa_terra_gimms, 
                                     probs = seq(.25, .75, .25))

# aqua vs. gimms
val_ioa_aqua_gimms <- rst_ioa[[1]][]
val_ioa_aqua_gimms_mu <- mean(val_ioa_aqua_gimms)
val_ioa_aqua_gimms_quan <- quantile(val_ioa_aqua_gimms, 
                                     probs = seq(.25, .75, .25))

## trend values
fls_ioa <- list.files(ch_dir_outdata, pattern = paste0("ioa_median_0312.tif$"), 
                     full.names = TRUE)
val_ioa <- lapply(fls_ioa, function(i) {
  rst_ioa <- raster(i)
  values(rst_ioa)
})

## density plot
cols <- c("black", "black", "grey50")
ltys <- c("longdash", "solid", "solid")
names(ltys) <- names(cols) <- c("Aqua|GIMMS", "Terra|Aqua", "Terra|GIMMS")

p_dens <- ggplot() + 
  geom_line(aes(x = val_ioa[[2]], y = ..count../sum(..count..), 
                linetype = "Terra|Aqua", colour = "Terra|Aqua"), 
            stat = "density", lwd = .8) + 
  geom_line(aes(x = val_ioa[[3]], y = ..count../sum(..count..), 
                linetype = "Terra|GIMMS", colour = "Terra|GIMMS"), 
            stat = "density", lwd = .8) + 
  geom_line(aes(x = val_ioa[[1]], y = ..count../sum(..count..), 
                linetype = "Aqua|GIMMS", colour = "Aqua|GIMMS"), 
            stat = "density", lwd = .8) + 
  geom_text(aes(x = .275, y = .008), label = "d)", fontface = "bold", size = 4) + 
  scale_linetype_manual("", breaks = c("Terra|Aqua", "Terra|GIMMS", "Aqua|GIMMS"), values = ltys) + 
  scale_colour_manual("", breaks = c("Terra|Aqua", "Terra|GIMMS", "Aqua|GIMMS"), values = cols) +
  labs(x = "IOA", y = "Density") + 
  theme_bw() + 
  theme(text = element_text(size = 9.6), 
        panel.grid = element_blank(),
        legend.key.height = unit(.5, "cm"),
        legend.margin = unit(0, "mm"),
        legend.key.size = unit(1, "cm"), 
        legend.key = element_rect(colour = "transparent"),
        legend.background = element_rect(fill = "transparent"),
        legend.text = element_text(size = 8),
        legend.position = c(.435, .8), legend.justification = c("center", "center"), 
        plot.margin = unit(rep(0, 4), units = "mm"), 
        panel.border = element_rect(colour = "black"), 
        axis.title.y = element_text(vjust = 1), 
        axis.title.x = element_text(vjust = -.2))


## combination final figure
p_ioa_comb <- latticeCombineGrid(p_ioa, layout = c(2, 2))

# png version (deprecated)
ch_fls_png <- paste0(ch_dir_outdata, "figure02.png")
png(ch_fls_png, width = 15.5, height = 15, units = "cm", res = 300)
plot.new()

print(p_ioa_comb, newpage = FALSE)

vp_rect <- viewport(x = .4975, y = .1, height = .36, width = .1, 
                    just = c("left", "bottom"))
pushViewport(vp_rect)
grid.rect(gp = gpar(col = "white"))

upViewport()
vp_dens <- viewport(x = .52, y = 0.1, width = .43, height = .325, 
                    just = c("left", "bottom"))
pushViewport(vp_dens)
print(p_dens, newpage = FALSE)
dev.off()

# standalone tiff version
ch_fls_tif <- paste0(ch_dir_outdata, "figure02.tiff")
tiff(ch_fls_tif, width = 15.5, height = 15, units = "cm", res = 500, 
     compression = "lzw")
plot.new()

print(p_ioa_comb, newpage = FALSE)

vp_rect <- viewport(x = .4975, y = .1, height = .36, width = .1, 
                    just = c("left", "bottom"))
pushViewport(vp_rect)
grid.rect(gp = gpar(col = "white"))

upViewport()
vp_dens <- viewport(x = .52, y = 0.1, width = .43, height = .325, 
                    just = c("left", "bottom"))
pushViewport(vp_dens)
print(p_dens, newpage = FALSE)
dev.off()

