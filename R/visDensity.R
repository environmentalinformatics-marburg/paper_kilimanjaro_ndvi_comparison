visDensity <- function(dsn = getwd(), p = 0.001, combined = TRUE) {
  
  ## packages
  lib <- c("Rsenal", "rgdal", "foreach", "ggplot2")
  jnk <- sapply(lib, function(x) library(x, character.only = TRUE))
  
  ## products
  products <- c("MOD13Q1.005", "MYD13Q1.005",
                "MOD13Q1.006", "MYD13Q1.006")
  
  ## trend values depending on 'p'
  fls_mk <- list.files(dsn, 
                       pattern = paste0("Q1.*mk_0312_tau", substr(p, 3, nchar(p))),
                       full.names = TRUE)[c(1, 3, 2, 4)]
  rst_mk <- lapply(fls_mk, raster)
  dat_mk <- foreach(i = 1:length(rst_mk), .combine = "rbind") %do% {
    val_mk <- getValues(rst_mk[[i]])
    val_mk <- as.numeric(na.omit(val_mk))
    data.frame(product = products[i], value = val_mk)
  }
  
  ## add 'collection' information
  dat_mk$collection <- NA
  dat_mk$collection[grep("005", dat_mk$product)] <- "C5"
  dat_mk$collection[grep("006", dat_mk$product)] <- "C6"
  
  ## reformat 'product' column
  dat_mk$product <- as.character(dat_mk$product)
  dat_mk$product <- sapply(strsplit(dat_mk$product, "\\."), "[[", 1)
  
  # ## subset with gimms data only
  # dat_mk_gimms <- subset(dat_mk, product == "GIMMS")
  # dat_mk_gimms <- rbind(cbind(dat_mk_gimms[, 1:2], collection = "005"),
  #                         cbind(dat_mk_gimms[, 1:2], collection = "006"))
  
  ## create figure
  
  # colors
  cols <- c("grey50", "black")
  
  # linetypes
  ltys <- c("solid", "41")
  names(ltys) <- names(cols) <- unique(dat_mk$product)
  
  # labels
  lbls <- c("MOD13Q1" = expression("NDVI"["Terra"]), 
            "MYD13Q1" = expression("NDVI"["Aqua"]))
  
  # ggplot
  if (combined) {
    obj <- list(p_value = p)
    p_dens <- ggplot(aes(x = value, group = product, colour = product, linetype = product),
                     data = subset(dat_mk, product != "GIMMS")) +
      geom_hline(yintercept = 0, colour = "grey70", size = .2) +
      geom_line(stat = "density", size = 1) +
      facet_wrap(~ collection) +
      # geom_line(stat = "density", data = dat_mk_gimms, size = .8) +
      scale_linetype_manual("Product", values = ltys, labels = lbls) +
      scale_colour_manual("Product", values = cols, labels = lbls) +
      labs(x = bquote("Kendall's " * tau ~ (p < .(obj$p_value))), y = "Density") +
      xlim(-.7, .7) +
      theme_bw() +
      theme(text = element_text(size = 9),
            panel.grid = element_blank(),
            legend.key.height = unit(.3, "cm"),
            legend.key.size = unit(1.2, "cm"),
            legend.key = element_rect(colour = "transparent"),
            legend.text = element_text(size = 8),
            legend.position = c(.5, .75), legend.justification = c("center", "center"),
            legend.background = element_rect(fill = "white", colour = "grey70"),
            plot.margin = unit(rep(0, 4), units = "mm"),
            panel.border = element_rect(colour = "black"),
            axis.title.y = element_text(vjust = 1),
            axis.title.x = element_text(vjust = 0), 
            strip.text = element_text(size = 9))
    
    return(p_dens)

  } else {
    
    lst_p_dens <- lapply(c("C5", "C6"), function(i) {
      ggplot(aes(x = value, group = product, colour = product, 
                 linetype = product),
             data = subset(dat_mk, product != "GIMMS" & collection == i)) +
      geom_hline(yintercept = 0, colour = "grey70", size = .2) +
      geom_line(stat = "density", size = 1) +
      scale_linetype_manual("", values = ltys, labels = lbls) +
      scale_colour_manual("", values = cols, labels = lbls) +
      labs(x = NULL, y = NULL) +
      xlim(-.6, .6) + ylim(0, ifelse(p == .05, 16, 8)) + 
      theme_bw() +
      theme(text = element_text(size = 9),
            panel.grid = element_blank(),
            legend.key.height = unit(.3, "cm"),
            legend.key.size = unit(1.2, "cm"),
            legend.key = element_rect(colour = "transparent"),
            legend.text = element_text(size = 8),
            legend.position = c(ifelse(i == "C5", .5, 5), .875), 
            legend.justification = c("center", "center"),
            legend.background = element_rect(fill = NA, colour = "transparent"),
            plot.margin = unit(rep(1, 4), units = "mm"),
            panel.border = element_rect(colour = "black", size = .5),
            axis.ticks = element_line(size = .5),
            axis.title.y = element_text(vjust = 1),
            axis.text.y = element_blank(),
            axis.title.x = element_text(vjust = 0), 
            strip.text = element_text(size = 9))
    })
    
    return(lst_p_dens)
    
  }    
}
