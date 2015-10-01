aggregateGimms <- function(files,
                           start = 5,
                           stop = 11, 
                           fun = max,
                           format = "%Y%j",
                           ...) {

  ## packages
  library(raster)
  library(zoo)
  
  ## monthly indexes
  dates <- substr(basename(files), start, stop)
  dates <- as.Date(dates, format = format)
  dates <- as.yearmon(dates)
  indices <- as.numeric(as.factor(dates))
  
  ## monthly aggregation
  rst <- stack(files)
  rst_agg <- stackApply(rst, indices = indices, fun = fun, ...)
  
  return(rst_agg)
}