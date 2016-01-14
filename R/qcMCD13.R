qcMCD13 <- function(product, ref_ext, type = c("tile", "cmg"), doy = TRUE, 
                    inpath = options()$MODIS_outDirPath, dsn = getwd(), 
                    apply_crop = TRUE, apply_qc = TRUE, 
                    apply_tso = TRUE, apply_adj = TRUE,
                    cores = 1L) {
  
  ## packages
  lib <- c("doParallel", "raster", "rgdal", "GSODTools")
  jnk <- sapply(lib, function(x) library(x, character.only = TRUE))

  ## parallelize  
  if (cores > 1L) {
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
  }
  
  type <- type[1]
  
  ## modis file patterns
  pttrn <- if (doy) {
    paste(product, c("NDVI.tif$", "pixel_reliability.tif$", 
                     "composite_day_of_the_year.tif$"), sep = ".*")
  } else {
    paste(product, c("NDVI.tif$", "pixel_reliability.tif$"), sep = ".*")
  }
  
  ## crop
  if (!dir.exists(dsn))
    dir.create(dsn)
    
  if (!missing(ref_ext) & apply_crop) {
    cat("Initializing 'crop' ...\n")
    
    suppressWarnings(
      ndvi.rst <- foreach(i = pttrn, .packages = lib, 
                          .export = ls(envir = globalenv())) %dopar% {                                      
                            
                            # available files
                            fls <- list.files(inpath, 
                                              pattern = i, full.names = TRUE) 
                            
                            # create folder
                            dsn_crp <- paste0(dsn, "/crp")
                            if (!dir.exists(dsn_crp))
                              dir.create(dsn_crp)
                            
                            # stack, crop and store
                            rst <- raster::stack(fls)
                            raster::crop(rst, ref_ext, snap = "out", 
                                         format = "GTiff", overwrite = TRUE,
                                         filename = paste0(dsn_crp, "/CRP"), 
                                         bylayer = TRUE, suffix = names(rst))
                          }
    )
    
  ## if no reference extent is supplied, proceed without cropping  
  } else if (missing(ref_ext) | !apply_crop) {
    suppressWarnings(
      ndvi.rst <- foreach(i = pttrn, .packages = lib, 
                          .export = ls(envir = globalenv())) %dopar% {
                            fls <- list.files(paste0(dsn, "/crp"), 
                                              pattern = paste0("^CRP.*", i),
                                              full.names = TRUE)
                            raster::stack(fls)
                          }
    )
  }

  ## quality assurance
  
  # overlay ndvi and qa layers
  if (apply_qc) {
    cat("Initializing 'overlay' ...\n")
    
    dsn_qa <- paste0(dsn, "/qa")
    if (!dir.exists(dsn_qa))
      dir.create(dsn_qa)
    
    suppressWarnings(
      lst_ndvi_qa <- foreach(i = unstack(ndvi.rst[[1]]), j = ndvi.rst[[2]], 
                             .packages = lib, 
                             .export = ls(envir = globalenv())) %dopar% {
                               raster::overlay(i, j, fun = function(x, y) {
                                 x[!y[] %in% c(0:2)] <- NA
                                 return(x)
                               }, filename = paste0(dsn_qa, "/QA_", names(i)), 
                               format = "GTiff", overwrite = TRUE)
                             }
    )
    
    ndvi.rst.qa <- raster::stack(lst_ndvi_qa)
    
  } else {
    fls_qa <- list.files(paste0(dsn, "/qa"), pattern = "^QA.*.tif$", 
                         full.names = TRUE)
    ndvi.rst.qa <- raster::stack(fls_qa)
  }  
  
  if (type == "cmg")
    return(ndvi.rst.qa)
  
  ## additional outlier check
  
  # calc, tsOutliers
  if (apply_tso) {
    cat("Initializing 'tsOutliers' ...\n")
    
    ndvi_mat_qa <- raster::as.matrix(ndvi.rst.qa)
    ndvi_rst_sd <- ndvi.rst.qa
    
    ndvi_lst_sd <-  
      foreach(i = 1:nrow(ndvi_mat_qa), .packages = lib, 
              .export = ls(envir = globalenv())) %dopar% {
                val <- ndvi_mat_qa[i, ]
                id <- GSODTools::tsOutliers(val, lower_quantile = .4, 
                                            upper_quantile = .9, index = TRUE)
                val[id] <- NA
                return(matrix(val, ncol = length(val), byrow = TRUE))
              }
    
    ndvi_mat_sd <- do.call("rbind", ndvi_lst_sd)
    
    ndvi_rst_sd <- raster::setValues(ndvi_rst_sd, ndvi_mat_sd)
    
    dsn_sd <- paste0(dsn, "/sd")
    if (!dir.exists(dsn_sd))
      dir.create(dsn_sd)
    
    ndvi_rst_sd <- writeRaster(ndvi_rst_sd, format = "GTiff",
                               filename = paste0(dsn_sd, "/SD"),  
                               bylayer = TRUE, suffix = names(ndvi.rst.qa), 
                               overwrite = TRUE)
    
  } else {
    fls_sd <- list.files(paste0(dsn, "/sd"), pattern = "^SD.*.tif$", 
                         full.names = TRUE)
    ndvi_rst_sd <- stack(fls_sd)
  }
  
  ## reject neighboring pixels
  
  # adjacent
  if (apply_adj) {
    cat("Initializing 'adjacent' ...\n")
    
    suppressWarnings(
      lst_ndvi_fc <- foreach(i = unstack(ndvi_rst_sd), .packages = lib, 
                             .export = ls(envir = globalenv())) %dopar% {
                               cells <- which(is.na(i[]))
                               id <- raster::adjacent(i, cells = cells, 
                                                      directions = 8, pairs = FALSE)
                               i[id] <- NA
                               return(i)
                             }
    )
    
    rst_ndvi_fc <- stack(lst_ndvi_fc)
    
    # store
    dsn_fc <- paste0(dsn, "/adj")
    if (!dir.exists(dsn_fc))
      dir.create(dsn_fc)
    
    rst_ndvi_fc <- writeRaster(rst_ndvi_fc, format = "GTiff",
                               filename = paste0(dsn_fc, "/ADJ"),  
                               bylayer = TRUE, suffix = names(ndvi_rst_sd), 
                               overwrite = TRUE)

  } else {
    fls_fc <- list.files(paste0(dsn, "/adj"), pattern = "^ADJ.*.tif$", 
                         full.names = TRUE)
    rst_ndvi_fc <- stack(fls_fc)
  }  
  
  ## deregister parallel backend
  if (cores > 1L)
    parallel::stopCluster(cl)
  
  ## return quality controlled images
  return(rst_ndvi_fc)
}
