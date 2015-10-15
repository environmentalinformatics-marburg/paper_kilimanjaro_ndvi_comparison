### environmental stuff

## packages
lib <- c("doParallel", "Rsenal", "stargazer")
jnk <- sapply(lib, function(x) library(x, character.only = TRUE))

## parallelization
cl <- makeCluster(3)
registerDoParallel(cl)

## folders
ch_dir_extdata <- "/media/fdetsch/XChange/kilimanjaro/ndvi_comparison/data/rst/"
ch_dir_outdata <- "/media/fdetsch/XChange/kilimanjaro/ndvi_comparison/out/"


### data processing

## reference extent
fls_ref <- list.files(paste0(ch_dir_extdata, "MOD13Q1.005"), full.names = TRUE,
                      pattern = "^SCL_AGGMAX", recursive = TRUE)[1]
rst_ref <- raster(fls_ref)

## index of association (ioa; 2003-2012)
st_year <- "2003"
nd_year <- "2012"

## avl products and corresponding file patterns
products <- c("GIMMS3g", 
              "MOD13Q1.005", "MYD13Q1.005", 
              "MOD13Q1.006", "MYD13Q1.006")

pattern <- paste(c("^MVC_WHT", "^SCL_AGGMAX_WHT", "^SCL_AGGMAX_WHT", "^MVC", "^MVC"), 
                 ".tif$", sep = ".*")

## import data
ls_rst_ndvi <- foreach(i = products, j = pattern, .packages = "raster") %dopar% {
    
  # list avl files              
  fls_ndvi <- list.files(paste0(ch_dir_extdata, i), 
                         pattern = j, full.names = TRUE, recursive = TRUE)
  
  # create temporal subset
  st <- grep(st_year, fls_ndvi)[1]
  nd <- grep(nd_year, fls_ndvi)[length(grep(nd_year, fls_ndvi))]
  
  # import data
  fls_ndvi <- fls_ndvi[st:nd]
  rst_ndvi <- stack(fls_ndvi)
  
  if (i == "GIMMS3g") {
    spy <- rasterToPolygons(rst_ndvi[[1]])
    rst_ndvi <- crop(rst_ndvi, spy)
  }
  
  if (i %in% c("MOD13Q1.006", "MYD13Q1.006"))
    rst_ndvi <- crop(rst_ndvi, rst_ref)
  
  return(rst_ndvi)
}

## calculate ioa
dat_ioa <- foreach(i = 1:length(ls_rst_ndvi), .combine = "rbind") %do% {
  
  # status message
  cat("Processing list entry no. ", i, "...\n", sep = "")
  
  foreach(j = 1:length(ls_rst_ndvi), .combine = "rbind", 
          .packages = c("raster", "Rsenal")) %dopar% {
    
    # if stacks are identical, ioa equals 1        
    if (i == j) {
      val_ioa <- 1
    } else {
      
      rst1 <- ls_rst_ndvi[[i]]
      
      # resample modis
      if (i == 1 & j != 1) {
        rst2 <- resample(ls_rst_ndvi[[j]], rst1)
      } else {
        rst2 <- ls_rst_ndvi[[j]]
        
        if (i != 1 & j == 1) 
          rst1 <- resample(rst1, rst2)
      }
      
      # extract values
      val1 <- rst1[]
      val2 <- rst2[]
      
      # calculate ioa
      val_ioa <- mean(sapply(1:nrow(val1), function(k) {
        Rsenal:::ioaC(val1[k, ], val2[k, ])
      }))
    }
            
    data.frame(ref1 = products[i], ref2 = products[j], ioa = val_ioa)  
  }
}

save(dat_ioa, file = "data/ioa.RData")

## reformat table
out <- matrix(ncol = 5, nrow = length(products))
for (i in 1:length(products)) {
  sub <- subset(dat_ioa, ref1 == products[i])
  out[i, ] <- sub$ioa
}

out <- round(out, 3)
out <- data.frame(out)
products[1] <- "GIMMS NDVI3g"
rownames(out) <- products
colnames(out) <- products

# latex output
stargazer(out, summary = FALSE)

## close parallel backend
stopCluster(cl)