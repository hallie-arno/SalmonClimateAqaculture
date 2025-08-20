library(tidyverse)
library(raster)
library(data.table)
library(sdmpredictors)

lays <- c(
  "MS_sst05_5m",
  "MS_sst06_5m",
  "MS_sst07_5m",
  
  #"MS_biogeo13_sst_mean_5m",
  #"MS_biogeo14_sst_min_5m",
  #"MS_biogeo15_sst_max_5m",
  #"MS_biogeo16_sst_range_5m",
  #"MS_biogeo17_sst_variance_5m",
  
  
  "MS_sss05_5m",
  "MS_sss06_5m",
  'MS_sss07_5m',
  
  "BO_ph"#,
  #"BO2_templtmax_ss",
  #"BO2_templtmin_ss"
  
)


options(sdmpredictors_datadir= paste0(getwd(), "/Metadata"))
load_layers(lays, rasterstack = FALSE)


## List of lat longs to get data for
samplePoints <- fread("Metadata/CIGENE_220K_Metadata_All_Aug_15_2023.txt") %>%
  filter(Region == "NorthAm") %>% 
  dplyr::select(FID, Lat, Lon) %>%
  distinct() %>%
  mutate(Lat = as.numeric(Lat)) %>%
  mutate(Lon = as.numeric(Lon))

## Example of downloading a raster -- these can be bulk downloaded in list using c()
#options(sdmpredictors_datadir= getwd())
#load_layers("MS_sst05_5m")

## Standard euclidean distance function
euclidDist <- function(x1, x2, y1, y2) {
  csq <- (x1 - x2)^2 + (y1-y2)^2
  dist <- sqrt(csq)
  return(dist)
}

## Function to find closest point in layer from a given river/ sampled point
findClosest <- function(samplePtx, samplePty, x, y, val) {
  df <- data.frame(cbind(x, y, val))
  
  smallDf <- df %>%
    # Filter coarsely to speed up the process
    # Could improve to filter a circle radius but in this case it doesn't matter
    filter(y < samplePty + 1 & y > samplePty - 1) %>%
    filter(x < samplePtx + 1 & x > samplePtx - 1) %>%
    rowwise() %>%
    # Get euclidean distance from sample pt to each pt in raster
    mutate(dist = euclidDist(samplePtx, x, samplePty, y)) %>%
    # Sort lowest first
    arrange(dist)
  # Extract data from first (shortest) distance
  
  varToExtract <- smallDf[1,3]
  return(varToExtract)
  
}


getLayerData <- function(layerName) {
  # Read from hard drive method -- can also download directly from sdmpredictors
  r <- raster(paste0("Metadata/",layerName, "_lonlat.tif"))
  # Crop to Northwest Atlantic to save computational power
  e <- as(extent(-80,-50,  42,  60), 'SpatialPolygons')
  crs(e) <- "+proj=longlat +datum=WGS84 +no_defs"
  rc <- crop(r, e)
  
  # Convert cropped raster to dataframe
  m <- as.data.frame(rc,xy=TRUE); coordinates(m)=~x+y
  df <- na.omit(as.data.frame(m))
  samplePoints <- samplePoints %>%
    rowwise() %>%
    mutate({{layerName}} := findClosest(samplePtx = Lon, samplePty = Lat, x = df$x, y = df$y, val = df[,3])) %>%
    unnest()
  
  
  
  # Apply findClostest function from above
  # samplePoints[[layerName]] <- mapply(FUN = findClosest, samplePtx = samplePoints$Lon, samplePty = samplePoints$Lat, x = df[,1], y = df[,2], val = df[,3])
  #print(samplePoints2)
  return(samplePoints)
  
  
}

envDataframe <- getLayerData("MS_sst05_5m")


for (lay in lays) {
  samplePoints <- getLayerData(lay)
  print(paste(lay, "done"))
}

colnames(samplePoints) <- c("FID", "Lat", "Lon", lays)
samplePoints

write.csv(samplePoints, "Metadata/estuary_env.csv")
