###Function for connecting to Pleiades site and downloading .csv files
pleiadesFunction = function(file){
  con = gzcon(url(paste("http://atlantides.org/downloads/pleiades/dumps/", file, sep="")))
  txt = readLines(con)
  dat = read.csv(textConnection(txt))
}

###Function for processing .hdf tiles and merging
library(raster)
library(rgdal)
library(gdalUtils)

setwd("") # set directory

projection = "+init=EPSG:4326"   # set projection (WGS84)
crs = crs(projection) # create crs objection

dir = dir("folder containing tiles", ".hdf") # filter through directory for .hdf files

for(i in 1:length(dir)){  # for each file in directory, access desired raster, and write to new directory (example is for NPP)
  sds = get_subdatasets(dir[i])
  gdal_translate(sds[2], dst_dataset = "NPP.tif", overwrite = TRUE)
  r = raster("NPP.tif")
  r = writeRaster(r, paste0("tilesProj/NPP", i, ".tif"), overwrite = TRUE)
}

dir2 = list.files("tilesProj", ".tif")   # filter through directory for .tif files
rasterList = lapply(paste0("tilesProj/",dir2), raster)  #create list of rasters for each time 
rasterMerge =  do.call("merge",rasterList)   # merge all rasters to single raster
rasterMerge = projectRaster(rasterMerge, crs = crs)
rasterMerge = writeRaster(rasterMerge, "NPPmerge.tif", overwrite = TRUE) # export raster

###reclassify raster
rasterMerge = reclassify(rasterMerge, c(6, Inf, NA))  # remove values greater than 6



