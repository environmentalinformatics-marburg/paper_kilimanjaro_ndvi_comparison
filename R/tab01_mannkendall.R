### environmental stuff

## clear workspace
rm(list = ls(all = TRUE))

## packages
lib <- c("raster", "stargazer")
Orcs::loadPkgs(lib)

## functions
source("R/mkStats.R")

## folders
ch_dir_outdata <- "/media/fdetsch/XChange/kilimanjaro/ndvi_comparison/out/"


################################################################################
## p < 0.05
################################################################################

## new product labels
labels <- c("NDVI\\textsubscript{3g}",    
            "NDVI\\textsubscript{Terra-C5}", "NDVI\\textsubscript{Aqua-C5}", 
            "NDVI\\textsubscript{Terra-C6}", "NDVI\\textsubscript{Aqua-C6}")

## statistics (p < 0.05)
fls_mk05 <- list.files(ch_dir_outdata, pattern = "0312_tau05.tif$", 
                        full.names = TRUE)[c(1, 2, 4, 3, 5)]
rst_mk05 <- lapply(fls_mk05, raster)

# trend differences
df_mk_stats <- foreach(i = rst_mk05, .combine = "rbind") %do% mkStats(i)
df_mk_stats <- cbind(product = labels, df_mk_stats)

# reformat
df_mk_stats[, 1] <- as.character(df_mk_stats[, 1])
rownames(df_mk_stats) <- df_mk_stats[, 1]
df_mk_stats <- df_mk_stats[, -c(1, 4, 6)]
names(df_mk_stats) <- c("Trend pixels", "Trends (%)", "Greening (%)", "Browning (%)")

stargazer(df_mk_stats, summary = FALSE)


################################################################################
## p < 0.001
################################################################################

## statistics (p < 0.001)
fls_mk001 <- list.files(ch_dir_outdata, pattern = "0312_tau001.tif$", 
                       full.names = TRUE)[c(1, 2, 4, 3, 5)]
rst_mk001 <- lapply(fls_mk001, raster)

# trend differences
df_mk_stats <- foreach(i = rst_mk001, .combine = "rbind") %do% mkStats(i)
df_mk_stats <- cbind(product = labels, df_mk_stats)

# reformat
df_mk_stats[, 1] <- as.character(df_mk_stats[, 1])
rownames(df_mk_stats) <- df_mk_stats[, 1]
df_mk_stats <- df_mk_stats[, -c(1, 4, 6)]
names(df_mk_stats) <- c("Trend pixels", "Trends (%)", "Greening (%)", "Browning (%)")

stargazer(df_mk_stats, summary = FALSE)
