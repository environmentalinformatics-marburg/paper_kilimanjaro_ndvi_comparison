### environmental stuff

## clear workspace
rm(list = ls(all = TRUE))

## packages
lib <- c("doParallel", "RColorBrewer", "raster", "grid", "latticeExtra", "Orcs")
Orcs::loadPkgs(lib)

## parallelization
supcl <- makeCluster(2)
registerDoParallel(supcl)

## folders
ch_dir_outdata <- "/media/fdetsch/XChange/kilimanjaro/ndvi_comparison/out/"

## compute mean difference for p < 0.05 and p < 0.001 separately
invisible(
  foreach(p_value = c(.05, .001), .packages = lib) %dopar% {
    
    # products
    products <- c("GIMMS3g", 
                  "MOD13Q1.005", "MYD13Q1.005", 
                  "MOD13Q1.006", "MYD13Q1.006")
    
    # import mann-kendall layers (p < 0.05)
    pattern <- paste0("0312_tau", substr(p_value, 3, nchar(p_value)), ".tif$")
    
    fls_mk <- list.files(ch_dir_outdata, pattern = pattern, 
                         full.names = TRUE)
    
    fls_mk <- fls_mk[sapply(products, function(i) grep(i, fls_mk))]
    rst_mk <- lapply(fls_mk, raster)
    
    # loop over layers  
    dat_md <- foreach(i = 2:length(rst_mk), 
                      .combine = "rbind") %do% {
      
      foreach(j = 2:length(rst_mk), 
              .combine = "rbind") %do% {
        
        # if stacks are identical, mean difference equals 0
        if (i == j) {
          val_md <- 0
        } else {
          
          rst1 <- rst_mk[[i]]
          
          # resample modis
          if (i == 1 & j != 1) {
            rst2 <- rst1
            rst2[] <- NA
            
            spy1 <- rasterToPolygons(rst1)
            rst2[which(!is.na(rst1[]))] <- extract(rst_mk[[j]], spy1, 
                                                   fun = mean, na.rm = TRUE)
          } else {
            rst2 <- rst_mk[[j]]
            
            if (i != 1 & j == 1) {
              rst1 <- rst2
              rst1[] <- NA
              
              spy2 <- rasterToPolygons(rst2)
              rst1[which(!is.na(rst2[]))] <- extract(rst_mk[[i]], spy2, 
                                                     fun = mean, na.rm = TRUE)
            }
          }
          
          # extract values
          val1 <- rst1[]
          val2 <- rst2[]
          
          # calculate mean difference
          val_md <- Orcs::meanDifference(val1, val2)
        }
        
        data.frame(ref1 = products[i], ref2 = products[j], 
                   md = val_md)  
      }
    }
    
    # write results to file
    file_out <- paste0("data/md_tau", substr(p_value, 3, nchar(p_value)), ".RData")
    save(dat_md, file = file_out)

    # remove gimms from 'products'
    products <- products[-1]
    
    # insert values into raster template
    mat_md <- matrix(ncol = length(products), nrow = length(products))
    for (i in 1:length(products)) {
      sub <- subset(dat_md, ref1 == rev(products)[i])
      mat_md[i, ] <- sub$md
    }
    mat_md[mat_md == 0] <- NA
    
    rst_md <- raster(mat_md, xmn = 0, xmx = length(products), 
                     ymn = 0, ymx = length(products))
    
    ## labels
    #     if (p_value == 0.05) {
    #       lbl <- c(expression(bold("NDVI"["3g"])), 
    #                expression(bold("NDVI"["Terra-C5"])), 
    #                expression(bold("NDVI"["Aqua-C5"])), 
    #                expression(bold("NDVI"["Terra-C6"])), 
    #                expression(bold("NDVI"["Aqua-C6"])))
    #     } else {
    lbl <- c(expression(bold("NDVI"["Terra-C5"])), 
               expression(bold("NDVI"["Aqua-C5"])), 
               expression(bold("NDVI"["Terra-C6"])), 
               expression(bold("NDVI"["Aqua-C6"])))
    # }
    
    # colors
    cols <- colorRampPalette(brewer.pal(11, "RdBu"))
    
    # create figure
    p_md <- spplot(rst_md, col.regions = cols(100), at = seq(-.05, .05, .01), 
                   scales = list(draw = TRUE, at = seq(.5, xmax(rst_md)-.5, 1), 
                                 cex = .8, x = list(rot = 45), labels = lbl), 
                   colorkey = list(space = "top", width = .7,
                                   labels = list(cex = .8, 
                                                 at = seq(-.09, .09, .03)))) + 
      latticeExtra::layer(sp.polygons(rasterToPolygons(rst_md)), 
                          data = list(rst_md = rst_md))
    
    # in-text .png version
    file_out <- paste0(ch_dir_outdata, "figure0", ifelse(p_value == .05, "3", "5"))
    png(paste0(file_out, ".png"), width = 10.5, height = 12, 
        units = "cm", res = 500)
    # main figure
    grid.newpage()
    vp0 <- viewport(x = 0, y = 0, width = 1, height = .9, 
                    just = c("left", "bottom"), name = "vp_figure")
    pushViewport(vp0)
    print(p_md, newpage = FALSE)
    
    # key caption
    downViewport(trellis.vpname("figure"))
    grid.text(expression(bold("MD"[tau])), x = 0.5, y = 1.3, 
              just = c("centre", "bottom"), gp = gpar(font = 2, cex = .9))
    dev.off()
    
    # stand-alone .tiff version
    tiff(paste0(file_out, ".tiff"), width = 10.5, height = 12, units = "cm", 
         res = 500, compression = "lzw")
    # main figure
    grid.newpage()
    vp0 <- viewport(x = 0, y = 0, width = 1, height = .9, 
                    just = c("left", "bottom"), name = "vp_figure")
    pushViewport(vp0)
    print(p_md, newpage = FALSE)
    
    # key caption
    downViewport(trellis.vpname("figure"))
    grid.text(expression(bold("MD"[tau])), x = 0.5, y = 1.3, 
              just = c("centre", "bottom"), gp = gpar(font = 2, cex = .9))
    dev.off()
    
  }
)

## deregister parallel backend
stopImplicitCluster()
