### environmental stuff

## clear workspace
rm(list = ls(all = TRUE))

## packages
lib <- c("grid", "Rsenal", "foreach", "latticeExtra", "ggplot2", 
         "Orcs")
Orcs::loadPkgs(lib)

## folders
ch_dir_extdata <- "/media/fdetsch/XChange/kilimanjaro/ndvi_comparison/data/rst/"
ch_dir_outdata <- "/media/fdetsch/XChange/kilimanjaro/ndvi_comparison/out/"

### data processing

## gimms grid
rst_gimms <- raster(paste0(ch_dir_outdata, "/GIMMS3g_mk_0312_tau.tif"))

## mann-kendall trend tests (2003-2012; p < 0.05)
st_year <- "2003"
nd_year <- "2012"

products <- list("GIMMS3g", 
                 "MOD13Q1.005", "MYD13Q1.005", 
                 "MOD13Q1.006", "MYD13Q1.006")

## create and visualize mann-kendall trend layers
lst_p_mk <- lapply(c(.05, .001), function(p_value) {
  # status message
  cat("Processing p =", p_value, "\n")
  
  foreach(i = products, .packages = "gimms",
          .export = ls(envir = environment())) %dopar% {
            
            fls_ndvi <- list.files(paste0(ch_dir_extdata, i), recursive = TRUE,
                                   pattern = "^DSN_.*.tif$", full.names = TRUE)
            
            st <- grep(st_year, fls_ndvi)[1]
            nd <- grep(nd_year, fls_ndvi)[length(grep(nd_year, fls_ndvi))]
            
            fls_ndvi <- fls_ndvi[st:nd]
            rst_ndvi <- stack(fls_ndvi)
            
            if (i != "GIMMS3g")
              rst_ndvi <- resample(rst_ndvi, rst_gimms)
            
            fls_out <- paste0(ch_dir_outdata, "tmp/", i, "_", 
                              formatC(substr(p_value, 3, nchar(p_value)), width = 2, flag = 0)
                              , ".tif")
            
            gimms::significantTau(rst_ndvi, p = p_value, filename = fls_out, 
                                  overwrite = TRUE, format = "GTiff")
            
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
ch_fls_png <- paste0(ch_dir_outdata, "figure02.png")
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
ch_fls_tif <- paste0(ch_dir_outdata, "figure02.tiff")
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
ch_fls_png <- paste0(ch_dir_outdata, "figure04.png")
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
ch_fls_tif <- paste0(ch_dir_outdata, "figure04.tiff")
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

