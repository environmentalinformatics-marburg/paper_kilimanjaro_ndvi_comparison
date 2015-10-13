## trend values
val_mk <- lapply(rst_mk05, function(i) i[])

## density plot
cols <- c("black", "black", "grey50")
ltys <- c("longdash", "solid", "solid")
names(ltys) <- names(cols) <- c("GIMMS", "Terra", "Aqua")

p_dens <- ggplot() + 
  geom_line(aes(x = val_mk[[2]], y = ..count../sum(..count..), 
                linetype = "Terra", colour = "Terra"), 
            stat = "density", lwd = .8) + 
  geom_line(aes(x = val_mk[[4]], y = ..count../sum(..count..), 
                linetype = "Aqua", colour = "Aqua"), 
            stat = "density", lwd = .8) + 
  geom_line(aes(x = val_mk[[1]], y = ..count../sum(..count..), 
                linetype = "GIMMS", colour = "GIMMS"), 
            stat = "density", lwd = .8) + 
  geom_text(aes(x = -.675, y = .016), label = "d)", fontface = "bold", size = 4) + 
  scale_linetype_manual("", breaks = c("Terra", "Aqua", "GIMMS"), values = ltys) + 
  scale_colour_manual("", breaks = c("Terra", "Aqua", "GIMMS"), values = cols) +
  labs(x = expression("Kendall's " * tau), y = "Density") + 
  theme_bw() + 
  theme(text = element_text(size = 9.6), 
        panel.grid = element_blank(),
        legend.key.height = unit(.5, "cm"),
        legend.key.size = unit(1, "cm"), 
        legend.key = element_rect(colour = "transparent"),
        legend.text = element_text(size = 8),
        legend.position = c(.75, .7), legend.justification = c("center", "center"), 
        legend.background = element_rect(fill = "transparent"),
        plot.margin = unit(rep(0, 4), units = "mm"), 
        panel.border = element_rect(colour = "black"), 
        axis.title.y = element_text(vjust = 1), 
        axis.title.x = element_text(vjust = -.2))
