## non-vegetated modis pixels
modis_nonveg_cmg <- readRDS("data/modis_nonvegetated_cmg.rds")

products <- list("MOD13C2.005", "MYD13C2.005", 
                 "MOD13C2.006", "MYD13C2.006")

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
  
  labels <- list(expression(bold("a) NDVI"["Terra-C5"])), 
                 expression(bold("b) NDVI"["Aqua-C5"])), 
                 expression(bold("c) NDVI"["Terra-C6"])), 
                 expression(bold("d) NDVI"["Aqua-C6"])))

  foreach(i = products, txt = labels, 
          .export = ls(envir = environment())) %do% {
            
            # list avl files  
            fls_ndvi <- list.files(paste0(ch_dir_extdata, i), full.names = TRUE, 
                                   pattern = "^DSN_.*.tif$", recursive = TRUE)

            # import temporal subset
            st <- grep(st_year, fls_ndvi)[1]
            nd <- grep(nd_year, fls_ndvi)
            nd <- nd[length(nd)]
            
            fls_ndvi <- fls_ndvi[st:nd]
            rst_ndvi <- stack(fls_ndvi)
            
            rst_ndvi[modis_nonveg_cmg] <- NA
            
            p <- visMannKendall(rst = rst_ndvi, 
                                xlab = "", ylab = "",
                                p_value = p_value, crs = "+init=epsg:4326",
                                filename = paste0(ch_dir_outdata, i, "_mk_0312.tif"), 
                                at = seq(-.55, .55, .01), 
                                format = "GTiff", 
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
## combination final figure, p < 0.001
################################################################################
p_mk_comb <- latticeCombineGrid(lst_p_mk[[2]], 
                                layout = c(2, 2))
