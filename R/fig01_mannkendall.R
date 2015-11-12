### environmental stuff

## clear workspace
rm(list = ls(all = TRUE))

## packages
lib <- c("grid", "Rsenal", "doParallel", "latticeExtra", "ggplot2", 
         "Orcs", "stargazer")
jnk <- sapply(lib, function(x) library(x, character.only = TRUE))

## functions
source("R/mkStats.R")
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
p_dem <- visDEM(rst_dem, labcex = .6, cex = 1.6, col = "grey70")

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
rst_kili <- kiliAerial(upperLeft = c(num_ymax, num_xmin), 
                       lowerRight = c(num_ymin, num_xmax),
                       minNumTiles = 12L, rasterize = TRUE, 
                       projection = "+init=epsg:4326")


# create figure
p_bing <- spplot(rst_kili[[1]], col.regions = NA, colorkey = FALSE, 
                 sp.layout = list(rgb2spLayout(rst_kili, alpha = .8), 
                                  list("sp.text", loc = c(37.02, -3.35), 
                                       txt = "a)", font = 2, cex = .6, adj = c(.1, 1))),
                 xlim = c(num_xmin, num_xmax), 
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

labels <- list(expression(bold("b) NDVI"["3g"])), 
               expression(bold("c) NDVI"["Terra-C5"])), 
               expression(bold("d) NDVI"["Aqua-C5"])), 
               expression(bold("e) NDVI"["Terra-C6"])), 
               expression(bold("f) NDVI"["Aqua-C6"])))

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
                        sp.layout = list("sp.text", loc = c(37.04, -3.35), 
                                         txt = txt, font = 2, cex = .6, adj = c(.1, 1)), 
                        xlim = c(num_xmin, num_xmax), 
                        ylim = c(num_ymin, num_ymax), 
                        scales = list(draw = TRUE, cex = .5, 
                                      y = list(at = seq(-2.9, -3.3, -.2))), 
                        rewrite = FALSE)
    
    p <- p + as.layer(p_dem)
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

# add topographic map
upViewport()
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

# add topographic map
upViewport()
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

# in-text png version
ch_fls_png <- paste0(ch_dir_outdata, "figure03.png")
png(ch_fls_png, width = 20.5, height = 18, units = "cm", res = 500)
plot.new()

print(p_mk_comb, newpage = FALSE)

# add key caption
grid.text(bquote(bold("Kendall's " ~ tau)), x = 0.5, y = 1, 
          just = c("centre", "top"), gp = gpar(font = 2, cex = .85))

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

dev.off()


################################################################################
## statistics
################################################################################

## new product labels
labels <- list(expression("NDVI"["3g"]),    
               expression("NDVI"["Terra-C5"]), expression("NDVI"["Aqua-C5"]), 
               expression("NDVI"["Terra-C6"]), expression("NDVI"["Aqua-C6"]))


## statistics (p < 0.001)
fls_mk001 <- list.files(ch_dir_outdata, pattern = "0312_tau001.tif$", 
                       full.names = TRUE)[c(1, 2, 4, 3, 5)]
rst_mk001 <- lapply(fls_mk001, raster)

# trend differences
df_mk_stats <- foreach(i = rst_mk001, .combine = "rbind") %do% mkStats(i)
df_mk_stats <- cbind(product = unlist(labels), df_mk_stats)

# reformat
df_mk_stats[, 1] <- as.character(df_mk_stats[, 1])
rownames(df_mk_stats) <- df_mk_stats[, 1]
df_mk_stats <- df_mk_stats[, -c(1, 4, 6)]
names(df_mk_stats) <- c("Trend pixels", "Trends (%)", "Greening (%)", "Browning (%)")
rownames(df_mk_stats)[1] <- "GIMMS NDVI3g"

stargazer(df_mk_stats, summary = FALSE)

## statistics (p < 0.05)
fls_mk05 <- list.files(ch_dir_outdata, pattern = "0312_tau05.tif$", 
                       full.names = TRUE)[c(1, 2, 4, 3, 5)]
rst_mk05 <- lapply(fls_mk05, raster)

# trend differences
df_mk_stats <- foreach(i = rst_mk05, .combine = "rbind") %do% mkStats(i)
df_mk_stats <- cbind(product = products, df_mk_stats)

# reformat
df_mk_stats[, 1] <- as.character(df_mk_stats[, 1])
rownames(df_mk_stats) <- df_mk_stats[, 1]
df_mk_stats <- df_mk_stats[, -c(1, 4, 6)]
names(df_mk_stats) <- c("Trend pixels", "Trends (%)", "Greening (%)", "Browning (%)")
rownames(df_mk_stats)[1] <- "GIMMS NDVI3g"

stargazer(df_mk_stats, summary = FALSE)

# mean difference
## calculate ioa
dat_meandiff <- foreach(i = 2:length(rst_mk001), .combine = "rbind") %do% {
  
  foreach(j = 2:length(rst_mk001), .combine = "rbind") %do% {
    
    # if stacks are identical, ioa equals 1        
    if (i == j) {
      val_meandiff <- 0
    } else {
      
      rst1 <- rst_mk001[[i]]
      
      # resample modis
      if (i == 1 & j != 1) {
        rst2 <- resample(rst_mk001[[j]], rst1)
      } else {
        rst2 <- rst_mk001[[j]]
        
        if (i != 1 & j == 1) 
          rst1 <- resample(rst1, rst2)
      }
      
      # extract values
      val1 <- rst1[]
      val2 <- rst2[]
      
      # calculate mean difference
      val_meandiff <- Orcs::meanDifference(val1, val2)
    }
    
    data.frame(ref1 = products[i], ref2 = products[j], 
               md = val_meandiff)  
  }
}
save(dat_meandiff, file = "data/meandiff.RData")
load("data/meandiff.RData")

products <- products[2:5]
mat_meandiff <- matrix(ncol = length(products), nrow = length(products))
rst_meandiff <- raster(ncols = length(products), nrows = length(products), 
                       xmn = 0, xmx = 4, ymn = 0, ymx = 4)
for (i in 1:length(products)) {
  sub <- subset(dat_meandiff, ref1 == rev(products)[i])
  mat_meandiff[i, ] <- sub$md
}

rst_meandiff <- raster(mat_meandiff, xmn = 0, xmx = 4, ymn = 0, ymx = 4)
rst_meandiff[seq(4, ncell(rst_meandiff)-1, 3)] <- NA

txt <- c(expression("NDVI"['Terra-C5']), expression("NDVI"['Aqua-C5']), 
         expression("NDVI"['Terra-C6']), expression("NDVI"['Aqua-C6']))
cols <- colorRampPalette(brewer.pal(11, "RdBu"))
p_meandiff <- 
  spplot(rst_meandiff, col.regions = cols(100), at = seq(-.05, .05, .01), 
         scales = list(draw = TRUE, at = seq(.5, 4.5, 1), labels = txt, 
                       cex = .8, x = list(rot = 45)), 
         colorkey = list(space = "top", labels = list(cex = .8), width = .7)) + 
  latticeExtra::layer(sp.polygons(rasterToPolygons(rst_meandiff)))

png(paste0(ch_dir_outdata, "figure02.png"), width = 10.5, height = 12, 
    units = "cm", res = 500)
# main figure
grid.newpage()
vp0 <- viewport(x = 0, y = 0, width = 1, height = .9, 
                just = c("left", "bottom"), name = "vp_figure")
pushViewport(vp0)
print(p_meandiff, newpage = FALSE)

# key caption
downViewport(trellis.vpname("figure"))
grid.text(expression(bold("MD"[tau])), x = 0.5, y = 1.3, just = c("centre", "bottom"), 
          gp = gpar(font = 2, cex = .9))
dev.off()

tiff(paste0(ch_dir_outdata, "figure02.tiff"), width = 10.5, height = 12, 
     units = "cm", res = 500, compression = "lzw")
# main figure
grid.newpage()
vp0 <- viewport(x = 0, y = 0, width = 1, height = .9, 
                just = c("left", "bottom"), name = "vp_figure")
pushViewport(vp0)
print(p_meandiff, newpage = FALSE)

# key caption
downViewport(trellis.vpname("figure"))
grid.text(expression(bold("MD"[tau])), x = 0.5, y = 1.3, just = c("centre", "bottom"), 
          gp = gpar(font = 2, cex = .9))
dev.off()

## deregister parallel backend
stopCluster(cl)
