### environmental stuff

## packages
lib <- c("grid", "Rsenal", "doParallel", "latticeExtra", "ggplot2")
jnk <- sapply(lib, function(x) library(x, character.only = TRUE))

## functions
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
p_dem <- visDEM(rst_dem, labcex = .8, cex = 1.6)

## reference extent
fls_ndvi <- paste0(ch_dir_extdata, "MOD13Q1.006/whittaker_dsn/DSN_SCL_MVC_200301.tif")

rst_ndvi <- raster(fls_ndvi[1])
rst_ndvi <- projectRaster(rst_ndvi, crs = "+init=epsg:4326")
rst_ndvi <- trim(rst_ndvi)

num_xmin <- xmin(rst_ndvi)
num_xmax <- xmax(rst_ndvi)
num_ymin <- ymin(rst_ndvi)
num_ymax <- ymax(rst_ndvi)

## study area


## mann-kendall trend tests (2003-2012; p < 0.05)
st_year <- "2003"
nd_year <- "2012"

p_mk <- foreach(i = c("GIMMS3g", 
                      "MOD13Q1.005", "MYD13Q1.005", 
                      "MOD13Q1.006", "MYD13Q1.006"), 
                txt = c("b)", "c)", "d)", "e)", "f)"), .packages = lib) %dopar% {
                  
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
                                       txt = txt, font = 2, cex = .7), 
                      xlim = c(num_xmin, num_xmax), 
                      ylim = c(num_ymin, num_ymax), 
                      scales = list(draw = TRUE, cex = .6, 
                                    y = list(at = seq(-2.9, -3.3, -.2))), 
                      keycex = .4)
  
  # p <- p + as.layer(p_dem)
  p <- envinmrRasterPlot(p, rot = 90, height = .5, width = .4, key.cex = .7)
  
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
  geom_text(aes(x = -.675, y = .016), label = "d)", fontface = "bold", size = 4) + 
  scale_linetype_manual("", breaks = c("Terra", "Aqua", "GIMMS"), values = ltys) + 
  scale_colour_manual("", breaks = c("Terra", "Aqua", "GIMMS"), values = cols) +
  labs(x = expression("Kendall's " * tau), y = "Density") + 
  theme_bw() + 
  theme(text = element_text(size = 9.6), 
        panel.grid = element_blank(),
        legend.key.height = unit(.5, "cm"),
        legend.key.size = unit(1, "cm"), 
        legend.key = element_rect(colour = "transparent"),
        legend.text = element_text(size = 8),
        legend.position = c(.75, .7), legend.justification = c("center", "center"), 
        legend.background = element_rect(fill = "transparent"),
        plot.margin = unit(rep(0, 4), units = "mm"), 
        panel.border = element_rect(colour = "black"), 
        axis.title.y = element_text(vjust = 1), 
        axis.title.x = element_text(vjust = -.2))


## combination final figure
p_mk_comb <- latticeCombineGrid(p_mk, layout = c(2, 2))

## visualization

# in-text png version
ch_fls_png <- paste0(ch_dir_outdata, "figure01.png")
png(ch_fls_png, width = 15.5, height = 15, units = "cm", res = 500)
plot.new()

print(p_mk_comb, newpage = FALSE)

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
ch_fls_tif <- paste0(ch_dir_outdata, "figure01.tiff")
tiff(ch_fls_tif, width = 15.5, height = 15, units = "cm", res = 500, 
     compression = "lzw")
plot.new()

print(p_mk_comb, newpage = FALSE)

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

## deregister parallel backend
stopCluster(cl)