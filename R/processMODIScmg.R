### environmental stuff --------------------------------------------------------

## clear workspace
rm(list = ls(all = TRUE))

## set working directory
library(Orcs)
setwdOS(path_lin = "/media/fdetsch/XChange/", path_win = "D:/", 
        path_ext = "kilimanjaro/ndvi_comparison")

## load packages
lib <- c("raster", "rgdal", "MODIS", "doParallel", "Kendall", "RColorBrewer", 
         "reshape2", "ggplot2", "zoo", "GSODTools", "remote", "Rsenal", "rgeos")
loadPkgs(lib)

## parallelize
supcl <- makeCluster(3)
registerDoParallel(supcl)

## functions
dir_repo <- "/media/permanent/repositories/paper_kilimanjaro_ndvi_comparison/"
source(paste0(dir_repo, "R/qcMCD13.R"))


### data import ----------------------------------------------------------------

# ## re-organize files
# orgStruc(move = TRUE)

## 'MODIS' global settings
MODISoptions(localArcPath = paste0(getwd(), "/data/MODIS_ARC/"), 
             outDirPath = paste0(getwd(), "/data/MODIS_ARC/PROCESSED/"), 
             MODISserverOrder = c("LPDAAC", "LAADS"), quiet = TRUE, 
             outProj = "+init=epsg:21037")

## download and extract monthly ndvi and associated quality flag
for (collection in c("005", "006")) {
  for (product in c("MOD13C2", "MYD13C2")) {
    runGdal(product = product, job = paste(product, collection, sep = "."), 
            collection = collection, SDSstring = "1000000000001")
  }
}

## reference extents

# kili
rst_kili <- kiliAerial(minNumTiles = 20L)


### area contributions of modis cmg to gimms -----------------------------------

## gimms reference grid
rst_gimms <- raster("data/rst/GIMMS3g/trm/ndvi/TRM_PRJ_CRP_geo00apr15a.n14.VI3g_ndvi.tif")
spy_gimms <- rasterToPolygons(rst_gimms)

## modis cmg reference grid
rst_cmg <- raster("data/rst/MOD13C2.006/crp/CRP_MOD13C2.A2000032.CMG_0.05_Deg_Monthly_NDVI.tif")
spy_qc <- rasterToPolygons(rst_cmg)

## calculate proportionate area contributions of modis cmg cell per gimms cell
lst_areas <- lapply(1:length(spy_gimms), function(i) {
  cutset <- gIntersects(spy_qc, spy_gimms[i, ], byid = TRUE)
  cutset <- as(cutset, "logical")
  spy_qc_sub <- spy_qc[cutset, ]
  
  spy_qc_areas <- sapply(1:length(spy_qc_sub), function(j) {
    gArea(intersect(spy_qc_sub[j, ], spy_gimms[i, ])) / 
      gArea(spy_gimms[j, ])
  })
  
  data.frame(cell = which(cutset), area = spy_qc_areas)
})


### process data ---------------------------------------------------------------

## NDVI data
lst_cmg_agg <- lapply(c("MOD13C2", "MYD13C2"), function(product) {
  
  ## quality check
  dir_base <- paste0("data/rst/", product, ".006/")
  #   rst_qc <- qcMCD13(product, ref_ext = rst_kili, doy = FALSE, type = "cmg",
  #                     inpath = paste0(options()$MODIS_outDirPath, product, ".005"), 
  #                     apply_crop = TRUE, apply_qc = TRUE, 
  #                     apply_tso = FALSE, apply_adj = FALSE,
  #                     dsn = paste0("data/rst/", product, ".005/"), 
  #                     cores = 2L)
  #   
  #   fls_qc <- list.files(paste0("data/rst/", product, ".006/qc"), 
  #                        pattern = "^QA_.*.tif$", full.names = TRUE)
  #   rst_qc <- stack(fls_qc)
  #   
  #   ## scale factor
  #   dir_scl <- paste0(dir_base, "scl/")
  #   if (!dir.exists(dir_scl)) dir.create(dir_scl)
  #   
  #   lst_scl <- foreach(i = unstack(rst_qc), .packages = "raster") %dopar%
  #     calc(i, fun = function(x) {x * 0.0001}, 
  #          filename = paste0(dir_scl, "SCL_", names(i), ".tif"), 
  #          format = "GTiff", overwrite = TRUE)
  #   
  #   rst_scl <- stack(lst_scl)
  #   
  #   ## whittaker smoother
  #   fls_raw <- list.files(paste0("data/MODIS_ARC/PROCESSED/", product, ".005"),
  #                         pattern = paste(product, "NDVI.tif$", sep = ".*"), 
  #                         full.names = TRUE, recursive = TRUE)
  #   
  #   dir_wht <- paste0(dir_base, "whittaker")
  #   if (!dir.exists(dir_wht)) dir.create(dir_wht)
  #   
  #   lst_wht <- whittaker.raster(vi = rst_scl, removeOutlier = TRUE, 
  #                               threshold = .2, lambda = 6000, nIter = 3, 
  #                               timeInfo = orgTime(fls_raw, pillow = 0), 
  #                               outDirPath = dir_wht, 
  #                               overwrite = TRUE, format = "raster")
  #   
  #   rst_wht <- stack(lst_wht)
  #   rst_wht <- writeRaster(rst_wht, filename = paste0(dir_wht, "/WHT"), 
  #                          bylayer = TRUE, suffix = names(rst_scl), 
  #                          format = "GTiff", overwrite = TRUE)
  #   
  #   ## spatial aggregation
  dir_agg <- paste0(dir_base, "agg/")
  #   if (!dir.exists(dir_agg)) dir.create(dir_agg)
  #   
  #   lst_agg <- foreach(i = unstack(rst_wht), .packages = "raster") %dopar% {
  #     
  #     num_val <- sapply(lst_areas, function(j) {
  #       sum(j$area * i[j$cell], na.rm = TRUE)
  #     })
  # 
  #     rst_out <- setValues(rst_gimms, num_val)
  #     writeRaster(rst_out, filename = paste0(dir_agg, "AGG_", names(i), ".tif"), 
  #                 format = "GTiff", overwrite = TRUE)
  #   }
  #   
  #   rst_agg <- stack(lst_agg)
  
  ## deseasoning
  fls_agg <- list.files(dir_agg, pattern = "^AGG.*.tif$", full.names = TRUE)
  
  # start and end month
  st <- grep("2003", fls_agg)[1]
  nd <- grep("2013", fls_agg); nd <- nd[length(nd)]
  
  # deseason
  rst_agg <- raster::stack(fls_agg[st:nd])
  rst_dsn <- remote::deseason(rst_agg, use.cpp = TRUE)
  
  # store
  dir_dsn <- paste0(dir_base, "whittaker_dsn/")
  if (!dir.exists(dir_dsn)) dir.create(dir_dsn)
  
  rst_dsn <- writeRaster(rst_dsn, filename = paste0(dir_dsn, "DSN"), 
                         bylayer = TRUE, suffix = names(rst_agg), 
                         format = "GTiff", overwrite = TRUE)

  return(rst_dsn)
})  

## deregister parallel backend
stopCluster(cl)

## re-install most recent 'MODIS' version
detach("package:MODIS")
install.packages(paste0(dir_repo, "inst/extdata/MODIS_0.10-33.tar.gz"), 
                 repos = NULL)
