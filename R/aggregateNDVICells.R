aggregateNDVICells <- function(rst, 
                               rst_doy, 
                               dates, 
                               cores = 1,
                               save_output = FALSE,
                               ...) {
  
  ### Environmental settings
  lib <- c("raster", "zoo", "foreach", "dplyr")
  jnk <- sapply(lib, function(x) library(x, character.only = TRUE))


  ### Data processing
  months <- unique(as.yearmon(as.character(dates), format = "%Y%j"))
  months_len <- length(months)
  
  years <- as.numeric(substr(dates, 1, 4))
  
  mat <- as.matrix(rst)
  mat_doy <- as.matrix(rst_doy)
  
  ls_agg <- foreach(h = 1:nrow(mat)) %do% { 
    
    # status message
    if ((h %% 1000) == 0)
      cat("Currently processing cell", h, "...\n")
      
    val <- mat[h, ]
    doy <- mat_doy[h, ]
    
    df <- data.frame(dates, years, doy, val)
    
    if (any(is.na(df[, 3])))
      df[is.na(df[, 3]), 3] <- as.numeric(substr(df[is.na(df[, 3]), 1], 5, 7))
    
    df[, 1] <- sapply(1:nrow(df), function(i) {
      old_year <- df[i, 2]
      if (i == 1) {
        new_year <- old_year
      } else if ((df[i, 3] < df[i-1, 3]) & (df[i, 2] == df[i-1, 2])) {
        new_year <- old_year+1
      } else {
        new_year <- old_year
      }
      
      new_doy <- formatC(df[i, 3], width = 3, flag = "0")
      new_doy <- paste0(new_year, new_doy)
      return(new_doy)
    })
    
    df %>%
      dplyr::mutate(Month = as.yearmon(dates, format = "%Y%j")) %>%
      dplyr::group_by(Month) %>%
      dplyr::summarise(max(val, na.rm = TRUE)) %>%
      data.frame() -> val_agg
    
    names(val_agg) <- c("Month", as.character(h))
    
    if (!identical(months, val_agg[, 1]))
      val_agg <- merge(months, val_agg, by = 1, all.x = TRUE)
    
    return(val_agg[, 2])
  }
  
  ls_agg <- do.call("rbind", ls_agg)
  
  rst_template <- rst[[1]]
  lst_agg <- foreach(h = 1:ncol(ls_agg)) %do% {
    raster::setValues(rst_template, as.numeric(ls_agg[, h]))
  }
  
  rst_agg <- raster::stack(lst_agg)
  
  if (save_output)
    rst_agg <- writeRaster(rst_agg, ...)
  
  return(rst_agg)
}
