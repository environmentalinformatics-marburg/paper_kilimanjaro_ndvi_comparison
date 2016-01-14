# ## Environmental stuff
# 
# ## clear workspace
# rm(list = ls(all = TRUE))
# 
# ## Packages
# lib <- c("raster", "doParallel", "reshape2", "plyr", "dplyr", 
#          "Rsenal", "scales", "RColorBrewer", "latticeExtra")
# Orcs::loadPkgs(lib)
# 
# # Parallelization
# supcl <- makeCluster(3)
# registerDoParallel(supcl)
# 
# ## folders
# ch_dir_extdata <- "/media/fdetsch/XChange/kilimanjaro/ndvi_comparison/data/rst/"
# ch_dir_outdata <- "/media/fdetsch/XChange/kilimanjaro/ndvi_comparison/out/"
# 
# 
# ### data processing
# 
# ## reference extent
# fls_ref <- list.files(paste0(ch_dir_extdata, "MOD13Q1.005"), full.names = TRUE,
#                       pattern = "^SCL_AGGMAX", recursive = TRUE)[1]
# rst_ref <- raster(fls_ref)
# 
# ## start and end year
# st_year <- "2003"
# nd_year <- "2012"
# 
# ## avl products and corresponding file patterns
# products <- c("GIMMS3g", 
#               "MOD13Q1.005", "MYD13Q1.005", 
#               "MOD13Q1.006", "MYD13Q1.006")
# 
# pattern <- paste(c("^MVC_WHT", "^SCL_AGGMAX_WHT", "^SCL_AGGMAX_WHT", "^MVC", "^MVC"), 
#                  ".tif$", sep = ".*")
# 
# ## import files
# lst_rst_ndvi <- foreach(i = products, j = pattern) %dopar% {
#   
#   # list avl files              
#   fls_ndvi <- if (i == "GIMMS3g") {
#     rearrangeFiles(dsn = paste0(ch_dir_extdata, i), pattern = j, 
#                    pos = c(4, 6, 11) + 23, full.names = TRUE, 
#                    recursive = TRUE)
#   } else {
#     list.files(paste0(ch_dir_extdata, i), 
#                pattern = j, full.names = TRUE, recursive = TRUE)
#   }
#   
#   # import temporal subset
#   st <- grep(ifelse(i == "GIMMS3g", "03jan", st_year), fls_ndvi)[1]
#   nd <- grep(ifelse(i == "GIMMS3g", "12dec", nd_year), fls_ndvi)
#   nd <- nd[length(nd)]
#   
#   fls_ndvi <- fls_ndvi[st:nd]
#   rst_ndvi <- raster::stack(fls_ndvi)
#   
#   if (i %in% c("MOD13Q1.006", "MYD13Q1.006"))
#     rst_ndvi <- raster::crop(rst_ndvi, rst_ref)
#   
#   return(rst_ndvi)
# }
# 
# ## create gimms grid template
# rst_gimms <- lst_rst_ndvi[[1]]
# spy_gimms <- rasterToPolygons(rst_gimms[[1]])
# 
# ## extract timestamps
# dates <- names(rst_gimms)
# dates <- paste0(substr(dates, 27, 31), "01")
# 
# systime_locale <- Sys.getlocale(category = "LC_TIME")
# if (Sys.info()[["sysname"]] == "Windows") {
#   invisible(Sys.setlocale(category = "LC_TIME", locale = "C"))
# } else {
#   invisible(Sys.setlocale(category = "LC_TIME", locale = "en_US.UTF-8"))
# }
# 
# dates <- strptime(dates, format = "%y%b%d")
# dates <- strftime(dates, format = "%Y%m")
# 
# # revoke locale time adjustment
# Sys.setlocale(category = "LC_TIME", locale = systime_locale)
# 
# ## reformat gimms data
# mat_gimms <- getValues(lst_rst_ndvi[[1]])
# dat_gimms <- data.frame("cell" = 1:nrow(mat_gimms), mat_gimms)
# names(dat_gimms)[2:ncol(dat_gimms)] <- dates
# dat_gimms_mlt <- reshape2::melt(dat_gimms, id.vars = 1, 
#                                 variable.name = "month", value.name = "gimms")
# 
# ### MODIS
# 
# ## modis terra, collection 5
# rst_modis <- lst_rst_ndvi[[3]]
# mat_modis <- as.matrix(rst_modis)
# lst_cells <- cellFromPolygon(rst_modis[[3]], spy_gimms)
# 
# 
# ## calculate statistics per gimms cell
# lst_stats <- foreach(i = lst_cells) %do% {
#   
#   # extract values
#   mat <- mat_modis[i, ]
# 
#   # calculate means
#   num_mean <- colMeans(mat, na.rm = TRUE)
#   
#   # calculate quantiles
#   mat_quan <- apply(mat, 2, FUN = function(x) {
#     quantile(x, probs = c(.1, .5, .9), na.rm = TRUE)
#   })
# 
#   # return data
#   data.frame("mean" = num_mean, 
#              "quan10" = mat_quan[1, ], 
#              "quan50" = mat_quan[2, ], 
#              "quan90" = mat_quan[3, ])
# }
# 
# ## reformat data to match up with the gimms dataset
# lst_mlt <- foreach(i = 1:4, j = c("mean", "quan10", "quan50", "quan90")) %do% {
#   mat <- t(sapply(lst_stats, "[[", i))
#   dat <- data.frame("cell" = 1:63, mat)
#   names(dat)[2:ncol(dat)] <- dates
#   dat_mlt <- reshape2::melt(dat, id.vars = 1, variable.name = "month", 
#                             value.name = j)
# }
# 
# ## merge all datasets, including gimms
# dat_mlt <- Reduce(function(...) merge(..., by = 1:2, sort = FALSE), 
#                   append(list(dat_gimms_mlt), lst_mlt))
# 
# dat_mlt$cell <- factor(dat_mlt$cell)
# dat_mlt$date <- as.Date(paste0(dat_mlt$month, "01"), format = "%Y%m%d")
# 
# saveRDS(dat_mlt, file = "data/cellstats.rds")

## reimport data
dat_mlt <- readRDS("data/cellstats.rds")

## create figure 
library(ggplot2)
library(scales)

# p_cellts <- ggplot(aes(x = date), data = dat_mlt) + 
#   geom_ribbon(aes(ymin = quan10, ymax = quan90), fill = "grey85") + 
#   geom_line(aes(y = mean), color = "grey50", size = .6) + 
#   geom_line(aes(y = gimms), size = .5, color = "black", linetype = "solid") + 
#   # geom_text(aes(label = paste("IOA:", ioa)), data = df_ioa_mean, fontface = "bold", 
#   #           x = Inf, y = -Inf, hjust = 1.2, vjust = -.4, size = 2.5) +
#   facet_wrap(~ cell, ncol = 9) + 
#   scale_x_date(labels = date_format("%Y"), breaks = date_breaks("4 years"), 
#                limits = as.Date(c("2003-01-01", "2012-12-01"))) + 
#   scale_y_continuous(breaks = seq(0, .75, .25), 
#                      labels = c("0.0", "", "0.5", "")) + 
#   theme_bw() + 
#   labs(x = "Time (months)", y = "NDVI") +
#   theme(panel.grid = element_blank(), 
#         axis.title.x = element_text(size = 9, vjust = -.5), 
#         axis.title.y = element_text(size = 9, vjust = 1), 
#         axis.text = element_text(size = 7), 
#         strip.text = element_text(size = 5))

dat_mlt <- dat_mlt[, c(1, 2, 8, 3, 4, 5, 7)]
dat_mlt2 <- reshape2::melt(dat_mlt, id.vars = c(1:3, 6:7))

dat_mlt2[dat_mlt2$cell %in% c(23, 32, 33), c("quan10", "quan90", "value")] <- NA

cols <- brewer.pal(4, "PuOr")[c(1, 4)]
names(cols) <- c("mean", "gimms")
cols_ribbon <- brewer.pal(9, "PuOr")[3]
lwds <- c("mean" = .6, "gimms" = .5)
ltys <- c("mean" = "31", "gimms" = "solid")

p_cellts <- ggplot(aes(x = date, y = value, color = variable), 
                   data = subset(dat_mlt2, variable %in% c("mean", "gimms"))) + 
  geom_ribbon(aes(ymin = quan10, ymax = quan90), 
              fill = cols_ribbon, colour = "transparent", alpha = .75) + 
  geom_line(data = subset(dat_mlt2, variable %in% c("mean", "gimms"))) + 
  # geom_text(aes(label = paste("IOA:", ioa)), data = df_ioa_mean, fontface = "bold", 
  #           x = Inf, y = -Inf, hjust = 1.2, vjust = -.4, size = 2.5) +
  facet_wrap(~ cell, ncol = 9) + 
  scale_x_date(labels = date_format("%Y"), breaks = date_breaks("4 years"), 
               limits = as.Date(c("2003-01-01", "2012-12-01"))) + 
  scale_y_continuous(breaks = seq(0, .75, .25), 
                     labels = c("0.0", "0.25", "0.50", "0.75")) + 
  scale_colour_manual("", values = cols, 
                      labels = c(expression("NDVI"["3g"]), 
                                 expression("NDVI"["Aqua-C5"]))) + 
#   scale_linetype_manual("", values = ltys, 
#                         labels = c(expression("NDVI"["3g"]), 
#                                    expression("NDVI"["Aqua-C5"]))) +
#  scale_size_manual("", values = lwds) + 
  theme_bw() + 
  labs(x = "Time (months)", y = "NDVI") +
  theme(panel.grid = element_blank(), 
        axis.title.x = element_text(size = 9, vjust = -.5), 
        axis.title.y = element_text(size = 9, vjust = 1), 
        axis.text = element_text(size = 7), 
        strip.text = element_text(size = 5), 
        legend.position = c(.5, 1.06), legend.direction = "horizontal", 
        legend.key.size = unit(1, "cm"), 
        legend.key = element_rect(colour = "transparent", fill = NA), 
        legend.text = element_text(size = 7), 
        legend.background = element_rect(fill = NA))
  
# ## deregister parallel backend
# stopCluster(supcl)
