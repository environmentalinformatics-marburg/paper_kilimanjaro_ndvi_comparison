aggregateGimms <- function(files,
                           files_out, 
                           start = 5,
                           stop = 11, 
                           fun = max,
                           nodes = 1L,
                           format = "%Y%j",
                           ...) {

  ## packages
  library(raster)

  ## parallelization
  library(doParallel)
  cl_agg <- cluster(nodes)
  registerDoParallel(cl_agg)

  ## monthly indexes
  dates <- substr(basename(files), start, stop)
  indices <- as.numeric(as.factor(dates))

  ## monthly aggregation
  rst <- stack(files)
  
  rst_agg <- if (missing(files_out)) {
    foreach(i = unique(indices), .packages = c("raster", "rgdal")) %dopar% {
      rst_tmp <- rst[[which(indices == i)]]
      overlay(rst_tmp, fun = fun)
    }
    
  } else {
    foreach(i = unique(indices), j = 1:length(unique(indices)), 
            .packages = c("raster", "rgdal")) %dopar% {
      rst_tmp <- rst[[which(indices == i)]]
      overlay(rst_tmp, fun = fun, filename = file_out[j], ...)
    }
  }
  
  #   rst <- stack(files)
  #   rst_agg <- stackApply(rst, indices = indices, fun = fun, ...)
  
  return(rst_agg)
}