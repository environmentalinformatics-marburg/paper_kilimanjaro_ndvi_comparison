### environmental stuff

## packages
lib <- c("grid", "Rsenal", "foreach", "latticeExtra", "ggplot2")
jnk <- sapply(lib, function(x) library(x, character.only = TRUE))

## functions
source("R/visDEM.R")
source("R/visMannKendall.R")

## folders
ch_dir_extdata <- "/media/fdetsch/XChange/kilimanjaro/ndvi_comparison/data/"
ch_dir_outdata <- "/media/fdetsch/XChange/kilimanjaro/ndvi_comparison/out/"

### data processing

## digital elevation model (dem)
ch_fls_dem <- paste0(ch_dir_extdata, "dem/DEM_ARC1960_30m_Hemp.tif")
rst_dem <- raster(ch_fls_dem)
rst_dem <- projectRaster(rst_dem, crs = "+init=epsg:4326")
p_dem <- visDEM(rst_dem)

## mann-kendall trend tests (2003-2013; p < 0.001)
st_year <- "2003"
nd_year <- "2012"

p_mk <- foreach(i = c("mod13q1", "myd13q1", "gimms"), 
                txt = c("a)", "b)", "c)")) %do% {
                  
  fls_ndvi <- list.files(paste0(ch_dir_extdata, i), 
                         pattern = ".tif$", full.names = TRUE)
  
  st <- grep(st_year, fls_ndvi)[1]
  nd <- grep(nd_year, fls_ndvi)[length(grep(nd_year, fls_ndvi))]
  
  fls_ndvi <- fls_ndvi[st:nd]
  rst_ndvi <- stack(fls_ndvi)
  
  if (i == "gimms") {
    spy <- rasterToPolygons(rst_ndvi[[1]])
    rst_ndvi <- crop(rst_ndvi, spy)
  }

  p <- visMannKendall(rst = rst_ndvi, xlab = "", ylab = "",
                      p_value = .05, crs = "+init=epsg:4326",
                      filename = paste0(ch_dir_outdata, i, "_mk05_0312.tif"), 
                      format = "GTiff", overwrite = TRUE, 
                      sp.layout = list("sp.text", loc = c(37.04, -3.35), 
                                       txt = txt, font = 2, cex = 1.1))
  
  p <- p + as.layer(p_dem)
  p <- envinmrRasterPlot(p)
  
  return(p)
}

## statistics
fls_mk001 <- list.files(ch_dir_outdata, pattern = "mk05_0312.tif$", 
                        full.names = TRUE)
rst_mk001 <- lapply(fls_mk001, raster)

# terra
val_mk001_terra <- rst_mk001[[2]][]
val_mk001_terra_abs <- sum(!is.na(rst_mk001[[2]][]))
val_mk001_terra_rel <- val_mk001_terra_abs / ncell(rst_mk001[[2]])
val_mk001_terra_abs_pos <- sum(val_mk001_terra > 0, na.rm = TRUE)
val_mk001_terra_rel_pos <- val_mk001_terra_abs_pos / val_mk001_terra_abs
val_mk001_terra_abs_neg <- sum(val_mk001_terra < 0, na.rm = TRUE)
val_mk001_terra_rel_neg <- val_mk001_terra_abs_neg / val_mk001_terra_abs

# aqua
val_mk001_aqua <- rst_mk001[[3]][]
val_mk001_aqua_abs <- sum(!is.na(rst_mk001[[3]][]))
val_mk001_aqua_rel <- val_mk001_aqua_abs / ncell(rst_mk001[[3]])
val_mk001_aqua_abs_pos <- sum(val_mk001_aqua > 0, na.rm = TRUE)
val_mk001_aqua_rel_pos <- val_mk001_aqua_abs_pos / val_mk001_aqua_abs
val_mk001_aqua_abs_neg <- sum(val_mk001_aqua < 0, na.rm = TRUE)
val_mk001_aqua_rel_neg <- val_mk001_aqua_abs_neg / val_mk001_aqua_abs

# mean difference
val_mk001 <- cbind(val_mk001_terra, val_mk001_aqua)
val_mk001_cc <- val_mk001[complete.cases(val_mk001), ]
mean(val_mk001_cc[, 1] - val_mk001_cc[, 2])


# gimms
val_mk001_gimms <- rst_mk001[[1]][]
val_mk001_gimms_abs <- sum(!is.na(rst_mk001[[1]][]))
val_mk001_gimms_rel <- val_mk001_gimms_abs / ncell(rst_mk001[[1]])
val_mk001_gimms_abs_pos <- sum(val_mk001_gimms > 0, na.rm = TRUE)
val_mk001_gimms_rel_pos <- val_mk001_gimms_abs_pos / val_mk001_gimms_abs
val_mk001_gimms_abs_neg <- sum(val_mk001_gimms < 0, na.rm = TRUE)
val_mk001_gimms_rel_neg <- val_mk001_gimms_abs_neg / val_mk001_gimms_abs

## trend values
fls_mk <- list.files(ch_dir_outdata, pattern = paste0("mk001.*.tif$"), 
                     full.names = TRUE)
val_mk <- lapply(fls_mk, function(i) {
  rst_mk <- raster(i)
  values(rst_mk)
})

## density plot
cols <- c("black", "black", "grey50")
ltys <- c("longdash", "solid", "solid")
names(ltys) <- names(cols) <- c("GIMMS", "Terra", "Aqua")

p_dens <- ggplot() + 
  geom_line(aes(x = val_mk[[2]], y = ..count../sum(..count..), 
                linetype = "Terra", colour = "Terra"), 
            stat = "density", lwd = .8) + 
  geom_line(aes(x = val_mk[[3]], y = ..count../sum(..count..), 
                linetype = "Aqua", colour = "Aqua"), 
            stat = "density", lwd = .8) + 
  geom_line(aes(x = val_mk[[1]], y = ..count../sum(..count..), 
                linetype = "GIMMS", colour = "GIMMS"), 
            stat = "density", lwd = .8) + 
  geom_text(aes(x = -.675, y = .016), label = "d)", fontface = "bold", size = 6) + 
  scale_linetype_manual("", breaks = c("Terra", "Aqua", "GIMMS"), values = ltys) + 
  scale_colour_manual("", breaks = c("Terra", "Aqua", "GIMMS"), values = cols) +
  labs(x = expression(atop("", "Kendall's " * tau)), y = "Density\n") + 
  theme_bw() + 
  theme(text = element_text(size = 13), 
        panel.grid = element_blank(),
        legend.key.size = unit(1, "cm"), 
        legend.key = element_rect(colour = "transparent"),
        legend.text = element_text(size = 11),
        legend.position = c(.8, .65), legend.justification = c("center", "center"), 
        plot.margin = unit(rep(0, 4), units = "mm"), 
        panel.border = element_rect(colour = "black"))


## combination final figure
p_mk_comb <- latticeCombineGrid(p_mk, layout = c(2, 2))

ch_fls_fig <- paste0(ch_dir_outdata, "fig01__mannkendall05.png")
png(ch_fls_fig, width = 24, height = 26, units = "cm", 
    res = 300, pointsize = 15)
plot.new()

print(p_mk_comb, newpage = FALSE)

vp_dens <- viewport(x = .52, y = 0.15, width = .4575, height = .3, 
                    just = c("left", "bottom"))
pushViewport(vp_dens)
print(p_dens, newpage = FALSE)
dev.off()

