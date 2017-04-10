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

dir = dir("tiles_raw/", ".hdf") # filter through directory for .hdf files

for(i in 1:length(dir)){  # for each file in directory, access NPP raster, and write to new directory
  sds = get_subdatasets(paste0(getwd(),"/tiles_raw/",dir[i]))
  gdal_translate(sds[2], dst_dataset = "NPP.tif", overwrite = TRUE)
  r = raster("NPP.tif")
  r = writeRaster(r, paste0("tiles_R/NPP", i, ".tif"), overwrite = TRUE)
}

dir2 = list.files("tiles_R", ".tif")   # filter through directory for .tif files
rasterList = lapply(paste0("tiles_R/",dir2), raster)  #create list of rasters for each time 
rasterMerge =  do.call("merge",rasterList)   # merge all rasters to single raster
rasterMerge = projectRaster(rasterMerge, crs = crs)
rasterMerge = writeRaster(rasterMerge, "NPPmerge.tif", overwrite = TRUE) # export raster

###reclassify raster
rasterMerge = reclassify(rasterMerge, c(6, Inf, NA))  # reclassifies values greater than 6 as NA

###export .pdf of iterative outputs
pdf("/Users/danplekhov/Desktop/maps.pdf")
for (i in 1:5){
plot(AFG,xlim=xy$x,ylim=xy$y, main = period)
plot(sites[sites[[period]] ==1,], pch =1, cex=0.75, add=TRUE)
} 
dev.off()

###over function
library(sp)

for (i in 10:32){
period = names(sites[,i])
table = as.data.frame(table(over(sites[sites[[period]] ==1,], watersheds)[,1]))
analysis[[period]][analysis$watershed %in% table$Var1] = table$Freq 
}  




