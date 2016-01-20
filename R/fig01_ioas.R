### environmental stuff

## clear workspace
rm(list = ls(all = TRUE))

## packages
lib <- c("doParallel", "Rsenal", "rgdal", "stargazer", "RColorBrewer", "grid", 
         "lattice", "maptools", "gimms")
Orcs::loadPkgs(lib)

## functions
source("R/insideNP.R")
source("R/visDEM.R")

## parallelization
supcl <- makeCluster(3)
registerDoParallel(supcl)

## folders
ch_dir_data <- "/media/fdetsch/XChange/kilimanjaro/ndvi_comparison/data/"
ch_dir_extdata <- "/media/fdetsch/XChange/kilimanjaro/ndvi_comparison/data/rst/"
ch_dir_outdata <- "/media/fdetsch/XChange/kilimanjaro/ndvi_comparison/out/"


### data processing

## reference extent
fls_ref <- list.files(paste0(ch_dir_extdata, "MOD13Q1.005"), full.names = TRUE,
                      pattern = "^SCL_AGGMAX", recursive = TRUE)[1]
rst_ref <- raster(fls_ref)

## national park boundary
spy_np <- readOGR(dsn = paste0(ch_dir_data, "/shp"), 
                  layer = "fdetsch-kilimanjaro-new_np-1420532792846_epsg21037")
spy_np <- spy_np[1, ]

## start and end year
st_year <- "2003"
nd_year <- "2012"

## avl products and corresponding file patterns
products <- c("GIMMS3g", 
              "MOD13Q1.005", "MYD13Q1.005", 
              "MOD13Q1.006", "MYD13Q1.006", 
              "MOD13C2.005", "MYD13C2.005", 
              "MOD13C2.006", "MYD13C2.006")

pattern <- paste(c("^MVC_WHT", "^SCL_AGGMAX_WHT", "^SCL_AGGMAX_WHT", "^MVC", 
                   "^MVC", rep("^AGG_WHT", 4)), ".tif$", sep = ".*")

## import data
ls_rst_ndvi <- foreach(i = products, j = pattern, 
                       .packages = c("raster", "gimms")) %dopar% {
    
  # list avl files  
  fls_ndvi <- if (i == "GIMMS3g") {
    rearrangeFiles(dsn = paste0(ch_dir_extdata, i), pattern = j, 
                   pos = c(4, 6, 11) + 23, full.names = TRUE, 
                   recursive = TRUE)
  } else {
    list.files(paste0(ch_dir_extdata, i), 
               pattern = j, full.names = TRUE, recursive = TRUE)
  }
  
  # import temporal subset
  st <- grep(ifelse(i == "GIMMS3g", "03jan", st_year), fls_ndvi)[1]
  nd <- grep(ifelse(i == "GIMMS3g", "12dec", nd_year), fls_ndvi)
  nd <- nd[length(nd)]
  
  fls_ndvi <- fls_ndvi[st:nd]
  rst_ndvi <- stack(fls_ndvi)
  
  if (i %in% c("MOD13Q1.006", "MYD13Q1.006"))
    rst_ndvi <- crop(rst_ndvi, rst_ref)
  
  return(rst_ndvi)
}


### long-term means ------------------------------------------------------------

# ## calculate long-term means
# dir_ltm <- paste0(ch_dir_outdata, "ltm/")
# lst_ltm <- foreach(i = ls_rst_ndvi, j = products, .packages = "raster") %dopar% {
#   fls_ltm <- paste0(dir_ltm, "LTM_", j, ".tif")
#   
#   if (file.exists(fls_ltm))
#     raster(fls_ltm)
#   else
#     overlay(i, fun = function(...) mean(..., na.rm = TRUE), 
#             filename = fls_ltm, format = "GTiff", overwrite = TRUE)
# }

## create modis mask for non-vegetated surfaces (tile products, ndvi < 0.15; 
## see also wittich and hansing 1995)
# rst_ltm_modis <- stack(lst_ltm[2:5])
# rst_ltm_valid <- overlay(rst_ltm_modis, fun = function(w, x, y, z) {
#   w_num <- w[]; x_num <- x[]; y_num <- y[]; z_num <- z[]
#   
#   out <- rep(1, length(w_num))
#   out[w_num < .15 | x_num < .15 | y_num < .15 | z_num < .15] <- 0
#   return(out)
# })
# 
# modis_nonveg <- which(rst_ltm_valid[] == 0)
# saveRDS(modis_nonveg, file = "data/modis_nonvegetated.rds")
modis_nonveg <- readRDS("data/modis_nonvegetated.rds")

# create modis mask for non-vegetated surfaces (cmg products, ndvi < 0.15)
# rst_ltm_modis <- stack(lst_ltm[6:9])
# rst_ltm_valid <- overlay(rst_ltm_modis, fun = function(w, x, y, z) {
#   w_num <- w[]; x_num <- x[]; y_num <- y[]; z_num <- z[]
#   
#   out <- rep(1, length(w_num))
#   out[w_num < .15 | x_num < .15 | y_num < .15 | z_num < .15] <- 0
#   return(out)
# })
# 
# modis_nonveg_cmg <- which(rst_ltm_valid[] == 0)
# saveRDS(modis_nonveg_cmg, file = "data/modis_nonvegetated_cmg.rds")
modis_nonveg_cmg <- readRDS("data/modis_nonvegetated_cmg.rds")

## cells inside national park
# gimms_inside <- insideNP(ls_rst_ndvi[[1]][[1]], spy_np, id = TRUE)
# saveRDS(gimms_inside, file = "data/gimms_inside_np.rds")
# gimms_inside <- readRDS("data/gimms_inside_np.rds")

# modis_inside <- insideNP(ls_rst_ndvi[[2]][[1]], spy_np, id = TRUE)
# saveRDS(modis_inside, file = "data/modis_inside_np.rds")
# modis_inside <- readRDS("data/modis_inside_np.rds")

## cells above 4,000 m
ch_fls_dem <- paste0(ch_dir_extdata, "../dem/DEM_ARC1960_30m_Hemp.tif")
rst_dem <- raster(ch_fls_dem)
# rst_dem[rst_dem[] < 4000] <- NA
# spy_dem <- rasterToPolygons(rst_dem)
# spy_dem <- unionSpatialPolygons(spy_dem, IDs = rep(1, length(spy_dem)))
# 
# gimms_summit <- insideNP(ls_rst_ndvi[[1]][[1]], spy_dem, limit = .2, id = TRUE)
# saveRDS(gimms_summit, file = "data/gimms_near_summit.rds")
# gimms_summit <- readRDS("data/gimms_near_summit.rds")
# 
# modis_summit <- insideNP(ls_rst_ndvi[[2]][[1]], spy_dem, limit = .2, id = TRUE)
# saveRDS(modis_summit, file = "data/modis_near_summit.rds")
# modis_summit <- readRDS("data/modis_near_summit.rds")

# calculate ioa
# dat_ioa <- foreach(i = 1:5, .combine = "rbind") %do% {
#   
#   # status message
#   cat("Processing list entry no. ", i, "...\n", sep = "")
#   
#   foreach(j = 1:5, .combine = "rbind", 
#           .packages = c("raster", "Rsenal")) %dopar% {
#     
#     # if stacks are identical, ioa equals 1        
#     if (i == j) {
#       val_ioa <- 1
#     } else {
#       
#       rst1 <- ls_rst_ndvi[[i]]
#       
#       # # reject cells located above 4,000 m a.s.l.
#       # rst1[if (i == 1) gimms_summit else modis_summit] <- NA
#       
#       # discard non-vegetated pixels
#       if (i != 1) {
#         rst1[modis_nonveg] <- NA
#       } else {
#         rst1[modis_nonveg_cmg] <- NA
#       }
# 
#       # resample modis
#       if (i == 1 & j != 1) {
#         # rst2 <- resample(ls_rst_ndvi[[j]], rst1)
#         rst2 <- ls_rst_ndvi[[j+4]]
#         rst2[modis_nonveg_cmg] <- NA
#         
#       } else if (i > 1 & j > 1) {
#         rst2 <- ls_rst_ndvi[[j]]
#         rst2[modis_nonveg] <- NA
# 
#       } else if (i != 1 & j == 1) {
#         rst1 <- ls_rst_ndvi[[i+4]]
#         rst1[modis_nonveg_cmg] <- NA
#         rst2 <- ls_rst_ndvi[[j]]
#         rst2[modis_nonveg_cmg] <- NA
#       }
#       
#       # extract values
#       val1 <- rst1[]
#       val2 <- rst2[]
#       
#       # calculate ioa
#       val_ioa <- mean(sapply(1:nrow(val1), function(k) {
#         Rsenal:::ioaC(val1[k, ], val2[k, ])
#       }), na.rm = TRUE)
#     }
#             
#     data.frame(ref1 = products[i], ref2 = products[j], ioa = val_ioa)  
#   }
# }
# 
# save(dat_ioa, file = "data/ioa.RData")
load("data/ioa.RData")

## reformat table
out <- matrix(ncol = 5, nrow = 5)
for (i in 1:5) {
  sub <- subset(dat_ioa, ref1 == products[i])
  out[i, ] <- sub$ioa
}

out <- round(out, 3)
out <- data.frame(out)
products[1] <- "GIMMS NDVI3g"
rownames(out) <- products[1:5]
colnames(out) <- products[1:5]

# # csv output
# write.csv(out, "/media/fdetsch/XChange/kilimanjaro/ndvi_comparison/pub/data/ioas.csv", 
#           quote = FALSE)

# latex output
stargazer(out, summary = FALSE)

# ################################################################################
# ## figure (similar to fig. 2 on mean difference)
# ################################################################################
# 
# ## raster template
# products[1] <- "GIMMS3g"
# mat_ioa <- matrix(ncol = length(products), nrow = length(products))
# for (i in 1:length(products)) {
#   sub <- subset(dat_ioa, ref1 == rev(products)[i])
#   mat_ioa[i, ] <- sub$ioa
# }
# mat_ioa[mat_ioa == 1] <- NA
# rst_ioa <- raster(mat_ioa, xmn = 0, xmx = 5, ymn = 0, ymx = 5)
# 
# col_ioa <- colorRampPalette(brewer.pal(7, "YlOrRd"))
# txt_ioa <- c(expression(bold("NDVI"['3g'])), 
#              expression(bold("NDVI"['Terra-C5'])), expression(bold("NDVI"['Aqua-C5'])), 
#              expression(bold("NDVI"['Terra-C6'])), expression(bold("NDVI"['Aqua-C6'])))
# 
# ## create figure
# 
# # labels
# lbl <- seq(.86, .94, .02)
# for (i in seq(1, length(lbl)*2-2, 2))
#   lbl <- append(lbl, "", i)
# 
# p_ioa <- 
#   spplot(rst_ioa, col.regions = col_ioa, at = seq(.75, .87, .01), 
#          scales = list(draw = TRUE, at = seq(.5, 4.5, 1), labels = txt_ioa, 
#                        cex = .8, x = list(rot = 45)), 
#          colorkey = list(space = "top", labels = list(cex = .8, labels = lbl), 
#                          width = .7, at = seq(.75, .87, .01))) + 
#   latticeExtra::layer(sp.polygons(rasterToPolygons(rst_ioa)))
# 
# ## write to file
# png(paste0(ch_dir_outdata, "figure00.png"), width = 10.5, height = 12, 
#     units = "cm", res = 500)
# # main figure
# grid.newpage()
# vp0 <- viewport(x = 0, y = 0, width = 1, height = .9, 
#                 just = c("left", "bottom"), name = "vp_figure")
# pushViewport(vp0)
# print(p_ioa, newpage = FALSE)
# 
# # key caption
# downViewport(trellis.vpname("figure"))
# grid.text(expression(bold("IOA"["s"])), x = 0.5, y = 1.3, 
#           just = c("centre", "bottom"), gp = gpar(font = 2, cex = .9))
# dev.off()
# 
# tiff(paste0(ch_dir_outdata, "figure00.tiff"), width = 10.5, height = 12, 
#      units = "cm", res = 500, compression = "lzw")
# # main figure
# grid.newpage()
# vp0 <- viewport(x = 0, y = 0, width = 1, height = .9, 
#                 just = c("left", "bottom"), name = "vp_figure")
# pushViewport(vp0)
# print(p_ioa, newpage = FALSE)
# 
# # key caption
# downViewport(trellis.vpname("figure"))
# grid.text(expression(bold("IOA"["s"])), x = 0.5, y = 1.3, 
#           just = c("centre", "bottom"), gp = gpar(font = 2, cex = .9))
# dev.off()


################################################################################
### gimms grid
################################################################################

## reference extent
fls_ndvi <- paste0(ch_dir_extdata, "MOD13Q1.006/whittaker_dsn/DSN_SCL_MVC_200301.tif")

rst_ndvi <- raster(fls_ndvi[1])
rst_ndvi <- projectRaster(rst_ndvi, crs = "+init=epsg:4326")
rst_ndvi <- trim(rst_ndvi)

num_xmin <- xmin(rst_ndvi)
num_xmax <- xmax(rst_ndvi)
num_ymin <- ymin(rst_ndvi)
num_ymax <- ymax(rst_ndvi)

## digital elevation model (dem)
rst_dem <- aggregate(rst_dem, fact = 10)
rst_dem <- projectRaster(rst_dem, crs = "+init=epsg:4326")
p_dem <- visDEM(rst_dem, labcex = .6, cex = 1.6, col = "black")

## gimms grid
rst_gimms <- raster(paste0(ch_dir_outdata, "/GIMMS3g_mk_0312_tau.tif"))
rst_gimms <- projectRaster(rst_gimms, crs = "+init=epsg:4326")
spy_gimms <- rasterToPolygons(rst_gimms)

## study area
rst_kili <- kiliAerial(upperLeft = c(num_ymax, num_xmin), 
                       lowerRight = c(num_ymin, num_xmax),
                       minNumTiles = 12L, projection = "+init=epsg:4326", 
                       type = "google")

# # create figure (bing)
# p_bing <- spplot(rst_kili[[1]], col.regions = NA, colorkey = FALSE, 
#                  sp.layout = list(rgb2spLayout(rst_kili, quantiles = c(.005, .9775)), 
#                                   list("sp.text", loc = c(37.02, -2.86), 
#                                        txt = "a)", font = 2, cex = .6, 
#                                        adj = c(.1, 1), col = "black"), 
#                                   list("sp.text", loc = c(37.6, -3.4), 
#                                        txt = "\uA9 OpenStreetMap contributors", 
#                                        font = 2, cex = .4, col = "grey90")),
#                  maxpixels = ncell(rst_kili), xlim = c(num_xmin, num_xmax), 
#                  ylim = c(num_ymin, num_ymax), 
#                  scales = list(draw = TRUE, cex = .5, 
#                                y = list(at = seq(-2.9, -3.3, -.2))))

# create figure (google)
scale <- list("SpatialPolygonsRescale", layout.scale.bar(), scale = 0.08998623, 
              offset = c(37.05, -3.38), fill = c("transparent", "black"))
text1 = list("sp.text", c(37.05, -3.36), "0", cex = .5, font = 2)
text2 = list("sp.text", c(37.1575, -3.36), "10 km", cex = .5, font = 2)

arrow <- list("SpatialPolygonsRescale", layout.north.arrow(type = 1), 
             offset = c(37, -3.41), scale = .075)

p_bing <- spplot(rst_kili[[1]], col.regions = NA, colorkey = FALSE, 
                 sp.layout = list(rgb2spLayout(rst_kili, quantiles = c(0, .999)), 
                                  list("sp.text", loc = c(37.02, -2.86), 
                                       txt = "a)", font = 2, cex = .6, 
                                       adj = c(.1, 1), col = "black"), 
                                  scale, arrow, text1, text2), 
#                                   list("sp.text", loc = c(37.65, -3.425), 
#                                        txt = "\uA9 Google Maps", 
#                                        font = 2, cex = .4, col = "black")),
                 maxpixels = ncell(rst_kili), xlim = c(num_xmin, num_xmax), 
                 ylim = c(num_ymin, num_ymax), 
                 scales = list(draw = TRUE, cex = .5, 
                               y = list(at = seq(-2.9, -3.3, -.2))))

## insert values
## gimms stack
rst1 <- ls_rst_ndvi[[1]]
rst1[modis_nonveg_cmg] <- NA
val1 <- rst1[]

## loop over modis stacks
dat_ioa <- foreach(j = 6:length(ls_rst_ndvi), .combine = "cbind",
                   .packages = c("raster", "Rsenal")) %dopar% {
                     
                     # modis stack
                     rst2 <- ls_rst_ndvi[[j]]
                     rst2[modis_nonveg_cmg] <- NA
                     val2 <- rst2[]
                     
                     # calculate ioa
                     val_ioa <- sapply(1:nrow(val1), function(k) {
                       if (all(is.na(val1[k, ])) | all(is.na(val2[k, ])))
                         return(NA)
                       else
                         return(Rsenal:::ioaC(val1[k, ], val2[k, ]))
                     })
                     
                     matrix(val_ioa, ncol = 1)
                   }

# val_ioa_gimms <- apply(dat_ioa, 1, mean)
val_ioa_gimms <- dat_ioa[, 2] # aqua-modis, collection5 -> best agreement
rst_ioa_gimms <- rst1[[1]]
rst_ioa_gimms[] <- val_ioa_gimms
rst_ioa_gimms <- projectRaster(rst_ioa_gimms, crs = "+init=epsg:4326", 
                               method = "ngb")
rst_ioa_gimms <- trim(rst_ioa_gimms)

spg_ioa_gimms <- as(rst_ioa_gimms, "SpatialGridDataFrame")
spg_ioa_gimms@data <- data.frame("ID" = formatC(1:63, width = 2, flag = "0"))
p_lbl <- spplot(spg_ioa_gimms, "ID", colorkey = FALSE, col.regions = "transparent",
                panel = function(...) {
                  panel.gridplot(..., border = "transparent")
                  panel.text(..., labels = spg_ioa_gimms@data$ID, 
                             col = "white", font = 2, adj = c(-.7, 2.4), 
                             cex = .5)
                })

p_bing <- p_bing + 
  latticeExtra::as.layer(p_lbl)

## create figure
p_ioa_gimms <- spplot(rst_ioa_gimms, col.regions = envinmrPalette(500), 
                      xlim = c(num_xmin, num_xmax), 
                      ylim = c(num_ymin, num_ymax), 
                      at = seq(.425, .975, .005), 
                      colorkey = list(space = "top", labels = list(cex = .8), width = .7),
                      scales = list(draw = TRUE), alpha.regions = 1) + 
  latticeExtra::layer(sp.text("b)", loc = c(37.02, -2.86), font = 2, cex = .6, 
                              adj = c(.1, 1), col = "black"))

## combine figures, add gimms grid and contour lines
p_ioa_comb <- latticeCombineGrid(list(p_bing, p_ioa_gimms), 
                                 layout = c(2, 1))
p_ioa_comb <- p_ioa_comb + 
  latticeExtra::layer(sp.polygons(rasterToPolygons(rst_ioa_gimms), 
                                  lwd = 1, lty = 3, col = "white")) + 
  latticeExtra::as.layer(p_dem)

install.packages("inst/extdata/ggplot2_1.0.1.tar.gz", repos = NULL)
source("R/fig01_homogeneities.R")
g_cellts <- ggplot2::ggplotGrob(p_cellts)
for (i in seq(3, 27, 4))
  g_cellts$heights[[i]] = unit(.1,"in")

## write to file
png(paste0(ch_dir_outdata, "figure01.png"), width = 20, height = 20, 
    units = "cm", res = 500)
plot.new()
vp0 <- viewport(x = 0, y = .55, width = .9, height = .5, 
                just = c("left", "bottom"), name = "vp_fig")
pushViewport(vp0)
print(p_ioa_comb, newpage = FALSE)

upViewport()
vp1 <- viewport(x = .865, y = .5, width = .1, height = .5, 
                just = c("left", "bottom"), name = "vp_key")
pushViewport(vp1)
draw.colorkey(key = list(col = envinmrPalette(500), 
                         width = .6, height = .425,
                         at = seq(.35, .95, .005), 
                         space = "right"), draw = TRUE)

grid.text(expression("IOA"["s"]), x = 1.05, y = .5, rot = -90, 
          gp = gpar(font = 2, cex = .85))

upViewport(n = 0)
vp2 <- viewport(x = .01, y = 0, width = .99, height = .6, 
                just = c("left", "bottom"), name = "vp_fig2")
pushViewport(vp2)
# print(p_cellts, newpage = FALSE)
grid.draw(g_cellts)

grid.text("c)", x = .09, y = .985, gp = gpar(font = 2, cex = .6))
dev.off()

tiff(paste0(ch_dir_outdata, "figure01.tiff"), width = 20, height = 20, 
     units = "cm", res = 500, compression = "lzw")
plot.new()
vp0 <- viewport(x = 0, y = .55, width = .9, height = .5, 
                just = c("left", "bottom"), name = "vp_fig")
pushViewport(vp0)
print(p_ioa_comb, newpage = FALSE)

upViewport()
vp1 <- viewport(x = .865, y = .5, width = .1, height = .5, 
                just = c("left", "bottom"), name = "vp_key")
pushViewport(vp1)
draw.colorkey(key = list(col = envinmrPalette(500), 
                         width = .6, height = .425,
                         at = seq(.35, .95, .005), 
                         space = "right"), draw = TRUE)

grid.text(expression("IOA"["s"]), x = 1.05, y = .5, rot = -90, 
          gp = gpar(font = 2, cex = .85))

upViewport(n = 0)
vp2 <- viewport(x = .01, y = 0, width = .99, height = .6, 
                just = c("left", "bottom"), name = "vp_fig2")
pushViewport(vp2)
# print(p_cellts, newpage = FALSE)
grid.draw(g_cellts)

grid.text("c)", x = .09, y = .985, gp = gpar(font = 2, cex = .6))
dev.off()

## close parallel backend
stopCluster(supcl)

detach("package:ggplot2", unload=TRUE)
install.packages("inst/extdata/ggplot2_2.0.0.tar.gz", repos = NULL)
