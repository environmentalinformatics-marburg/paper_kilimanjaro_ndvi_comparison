aggregateGimms <- function(files_in,
                           files_out, 
                           start = 5,
                           stop = 11, 
                           fun = max,
                           nodes = 1L,
                           format = "%Y%j",
                           ...) {
  
  ## packages
  library(raster)
  
  ## import data
  rst <- raster::stack(files_in)
  
  ## monthly indexes
  dates <- substr(basename(files_in), start, stop)
  indices <- as.numeric(as.factor(dates))
  
  ## parallelization
  if (nodes > 1) {
    
    library(doParallel)
    cl_agg <- makeCluster(nodes)
    registerDoParallel(cl_agg)

    # monthly aggregation
    rst_agg <- if (missing(files_out)) {
      foreach(i = unique(indices), .packages = c("raster", "rgdal")) %dopar% {
        rst_tmp <- rst[[which(indices == i)]]
        overlay(rst_tmp, fun = fun)
      }
      
    } else {
      foreach(i = unique(indices), j = 1:length(unique(indices)), 
              .packages = c("raster", "rgdal")) %dopar% {
                rst_tmp <- rst[[which(indices == i)]]
                overlay(rst_tmp, fun = fun, filename = files_out[j], ...)
              }
    }
    
    # deregister parallel backend
    stopImplicitCluster()
    
    ## single-core processing  
  } else {
    
    # monthly aggregation
    rst_agg <- stackApply(rst, indices = indices, fun = fun, ...)
    
  }
  
  ## return mvc data
  return(rst_agg)
}