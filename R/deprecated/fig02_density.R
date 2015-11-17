source("R/visDensity.R")

## folders
ch_dir_outdata <- "/media/fdetsch/XChange/kilimanjaro/ndvi_comparison/out/"

foreach(p_value = c(.05, .001), file_out = c("/figure02", "/figure04")) %do% {
  
  # create figure
  p <- visDensity(p = p_value, dsn = ch_dir_outdata)

#   # write to file
#   png(paste0(ch_dir_outdata, file_out, ".png"), width = 9, height = 6, 
#       units = "cm", res = 500)
#   print(p)
#   dev.off()
}
