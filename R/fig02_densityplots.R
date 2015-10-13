## folders
ch_dir_extdata <- "/media/fdetsch/XChange/kilimanjaro/ndvi_comparison/data/rst/"
ch_dir_outdata <- "/media/fdetsch/XChange/kilimanjaro/ndvi_comparison/out/"

products <- c("GIMMS", 
              "MOD13Q1.005", "MYD13Q1.005", 
              "MOD13Q1.006", "MYD13Q1.006")

## trend values
fls_mk05 <- list.files(ch_dir_outdata, pattern = "mk05_0312.tif$", 
                       full.names = TRUE)[c(1, 2, 4, 3, 5)]
rst_mk05 <- lapply(fls_mk05, raster)
dat_mk05 <- foreach(i = 1:length(rst_mk05), .combine = "rbind") %do% {
  val_mk05 <- getValues(rst_mk05[[i]])
  val_mk05 <- as.numeric(na.omit(val_mk05))
  data.frame(product = products[i], value = val_mk05)
}

## add 'collection' column
dat_mk05$collection <- NA
dat_mk05$collection[grep("005", dat_mk05$product)] <- "005"
dat_mk05$collection[grep("006", dat_mk05$product)] <- "006"

## reformat 'product' column
dat_mk05$product <- as.character(dat_mk05$product)
dat_mk05$product <- sapply(strsplit(dat_mk05$product, "\\."), "[[", 1)

## subset with gimms data only
dat_mk05_gimms <- subset(dat_mk05, product == "GIMMS")
dat_mk05_gimms <- rbind(cbind(dat_mk05_gimms[, 1:2], collection = "005"), 
                        cbind(dat_mk05_gimms[, 1:2], collection = "006"))

## density plot

# collection 005
cols <- brewer.pal(9, "RdBu")[c(3, 7, 9)]
cols <- envinmrPalette(5)[c(3, 2, 5)]
ltys <- c("solid", "solid", "longdash")
names(ltys) <- names(cols) <- unique(dat_mk05$product)

p_dens <- ggplot(aes(x = value, group = product, colour = product, linetype = product), 
                 data = subset(dat_mk05, product != "GIMMS")) + 
  geom_hline(yintercept = 0, colour = "grey70", size = .2) + 
  geom_line(stat = "density", size = .8) + 
  facet_wrap(~ collection) + 
  geom_line(stat = "density", data = dat_mk05_gimms, size = .8) + 
  scale_linetype_manual("Product", values = ltys) + 
  scale_colour_manual("Product", values = cols) +
  labs(x = expression("Kendall's " * tau), y = "Density") + 
  xlim(-.9, .9) + 
  theme_bw() + 
  theme(text = element_text(size = 9), 
        panel.grid = element_blank(),
        legend.key.height = unit(.3, "cm"),
        legend.key.size = unit(1, "cm"), 
        legend.key = element_rect(colour = "transparent"),
        legend.text = element_text(size = 7),
        legend.position = c(.5, .75), legend.justification = c("center", "center"), 
        legend.background = element_rect(fill = "white", colour = "grey70"),
        plot.margin = unit(rep(0, 4), units = "mm"), 
        panel.border = element_rect(colour = "black"), 
        axis.title.y = element_text(vjust = 1), 
        axis.title.x = element_text(vjust = 0))

# in-text png version
png(paste0(ch_dir_outdata, "figure02.png"), width = 14, height = 7, units = "cm", 
    res = 500)
print(p_dens)
dev.off()

# standalone tiff version
tiff(paste0(ch_dir_outdata, "figure02.tiff"), width = 14, height = 7, 
     units = "cm", res = 500, compression = "lzw")
print(p_dens)
dev.off()