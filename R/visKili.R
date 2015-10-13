visKili <- function(lwd = 2, col = "black", fill = "white", cex = 2) {
  
  library(rworldmap)
  library(Rsenal)
  
  data("countriesCoarse")
  
  ## africa
  spy_africa <- subset(countriesCoarse, REGION == "Africa")
  
  ## tanzania
  spy_tanzania <- subset(countriesCoarse, POSTAL == "TZ")
  
  ## kilimanjaro
  rst_kili <- kiliAerial(projection = proj4string(spy_africa), rasterize = TRUE)
  
  spt_kili <- data.frame(x = (xmin(rst_kili) + xmax(rst_kili)) / 2, 
                         y = (ymin(rst_kili) + ymax(rst_kili)) / 2)
  coordinates(spt_kili) <- ~ x + y
  
  ## visualization
  ylim <- c(ymin(spy_africa) - 2, ymax(spy_africa) + 2)
  
  p <- spplot(spy_africa, col = "transparent", "ADMIN", colorkey = FALSE, col.regions = "grey85", ylim = ylim,
              par.settings = list(panel.background = list(col = "white"))) 
  
  p +
    latticeExtra::layer(panel.refline(h = 0, col = "grey60", lty = 2, lwd = .5)) + 
    latticeExtra::layer(sp.polygons(spy_tanzania, fill = "grey50", col = "transparent"), 
                        data = list(spy_tanzania = spy_tanzania)) + 
    latticeExtra::layer(sp.points(spt_kili, col = "black", pch = 24, cex = .5, fill = "white"), 
                        data = list(spt_kili = spt_kili))
}