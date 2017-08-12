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

###graphs
#stacked bar
par(xpd=T, mar=par()$mar+c(0,0,0,6))
x = barplot(prop,
            col=cols,
            xaxt="n",
            width=1,
            main = "Percentages of Dining/Food Ceramic Vessels",
            bty = "L")
labs = paste(names(sara))
text(cex=1, x=x+0.25, y=-0.02, labels = labs, xpd=TRUE, srt=45, pos=2)
legend("right", inset = c(-0.15, 0), fill=cols, legend=rownames(data), cex = 0.5, xpd = TRUE, bty = "n",
       text.width = 0.25,
       x.intersp = 0.25)

dev.off()

#grouped bar
par(xpd=TRUE, mar=par()$mar+c(0,0,0,6))
x = barplot(prop,
            col=cols,
            xaxt="n",
            width=1,
            beside = T,
            main = "Percentages of Dining/Food Ceramic Vessels")
labs = paste(names(sara))
text(cex=0.75, x= colMeans(x)+0.25, y=-0.03, labels = labs, xpd=TRUE, srt=45, pos=2)
legend("right", inset = c(-0.1, 0), fill=cols, legend=rownames(data), cex = 0.5, xpd = TRUE, bty = "n",
       text.width = 0.25,
       x.intersp = 0.25)

#pie
x = pie$V2
labels = pie$V1
pie(x, labels, radius = 2.5, main = "Faunal Finds", col = cols, clockwise = T)

#exif data
library(exifr)

dir = dir("/Users/danplekhov/Desktop/photos/photos", full.names=T)
exifinfo = exifr(dir)[,81:83]

write.csv(exifinfo, "/Users/danplekhov/Desktop/photoData.csv")

#gpx data
library(rgdal)

list = ogrListLayers("/Users/danplekhov/Desktop/Waypoints_01-JUN-17.gpx")
gpx_waypoints = readOGR("/Users/danplekhov/Desktop/Waypoints_01-JUN-17.gpx", layer = list[1])
writeOGR(gpx_waypoints, dsn = "/Users/danplekhov/Desktop/", layer = "gpx_waypoints", driver = "ESRI Shapefile")

#download files
dest = "/Users/danplekhov/Desktop/Working Files/stonewalls/BU/2014 USDA NAIP Digital True Color Orthophotography"
temp = tempfile()

url = "http://www.rigis.org/data/img/2014usda/tif/BU"
page = read_html(url)
links = html_attr(html_nodes(page, "a"), "href")
zipfiles = grep(".zip", links) 
links = links[zipfiles]

for(i in 1:length(links)){
  download.file(paste0("http://www.rigis.org/", links[i]), destfile = paste0(temp,basename(links[i])), mode = "wb")
  unzip(paste0(temp,basename(links[i])), files = paste0(file_path_sans_ext(basename(links[i])), ".tif"), exdir = dest)
}




