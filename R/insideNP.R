insideNP <- function(x, y, limit = .2, id = FALSE) {
  
  ## cells intersecting polygon
  cells_inside <- raster::extract(x, y, cellnumbers = TRUE)
  cells_inside <- cells_inside[[1]][, 1] 
  
  ## cells intersecting polygon boundary
  cells_spl <- raster::extract(x, as(y, "SpatialLines"), cellnumbers = TRUE)
  cells_spl <- cells_spl[[1]][, 1]
  
  ## subset boundary cells
  spy_x <- raster::rasterToPolygons(x)
  spy_x <- spy_x[cells_spl, ]

  ## loop over single boundary cells
  id_outside <- sapply(1:length(spy_x), function(i) {
    
    # subset with current polygon
    spy_x_tmp <- spy_x[i, ]
    
    # intersect subset with national park area
    suppressWarnings(  
      spy_intersect <- rgeos::gIntersection(spy_x_tmp, y)
    )
    
    # if objects do not intersect...
    if (is.null(spy_intersect)) {
      return(TRUE)
      
    # ...else compute percent coverage in current polygon  
    } else {
      num_intersect <- rgeos::gArea(spy_intersect) / rgeos::gArea(spy_x_tmp)
      
      if (num_intersect <= limit) 
        return(TRUE) 
      else 
        return(FALSE)
    }
  })
  
  ## remove cells with small portion of national park coverage from mask
  cells_outside <- cells_spl[id_outside]
  if (length(which(cells_inside %in% cells_outside)) > 0)
    cells_inside <- cells_inside[!cells_inside %in% cells_outside]

  ## reject cells inside the national park
  x[cells_inside] <- NA
  
  if (id) 
    return(cells_inside) 
  else 
    return(x)
}