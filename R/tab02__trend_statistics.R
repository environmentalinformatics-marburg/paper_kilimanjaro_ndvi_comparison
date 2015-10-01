rst_terra_005 <- raster(paste0(ch_dir_outdata, "MOD13Q1.005_mk05_0312.tif"))
rst_aqua_005 <- raster(paste0(ch_dir_outdata, "MYD13Q1.005_mk05_0312.tif"))

rst_terra_006 <- raster(paste0(ch_dir_outdata, "MOD13Q1.006_mk05_0312.tif"))
rst_aqua_006 <- raster(paste0(ch_dir_outdata, "MYD13Q1.006_mk05_0312.tif"))

## trends
foreach(i = list(rst_terra_005, rst_aqua_005, rst_terra_006, rst_aqua_006), 
        j = list("MOD13Q1.005", "MYD13Q1.005", "MOD13Q1.006", "MYD13Q1.006"), 
        .combine = "rbind") %do% {
  
  # trends
  num_tau <- i[]
  num_trd <- sum(!is.na(num_tau) / length(num_tau))
  num_trd_rel <- round(num_trd / length(num_tau), 3)
  
  # greening and browning
  num_tau <- num_tau[!is.na(num_tau)]
  num_grn <- round(sum(num_tau > 0) / length(num_tau), 3)
  num_brn <- round(sum(num_tau < 0) / length(num_tau), 3)
  
  # return values
  data.frame(product = j, trends_abs = num_trd, trends_rel = num_trd_rel, 
             greening = num_grn, browning = num_brn)
}