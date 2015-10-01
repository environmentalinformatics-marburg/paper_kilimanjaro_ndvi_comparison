visKili <- function(lwd = 2, col = "red", cex = 2) {
  
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
  
  spplot(spy_africa, "ADMIN", colorkey = FALSE, col.regions = "grey85", ylim = ylim,
         sp.layout = list(list("sp.lines", as(spy_africa, "SpatialLines"), col = "grey50"), 
                          list("sp.lines", as(spy_tanzania, "SpatialLines"), 
                               lwd = lwd), 
                          list("sp.points", spt_kili, col = col, pch = 20, cex = cex)), 
         par.settings = list(panel.background = list(col = "white")))
}