##packages used
##link here -- http://geog.uoregon.edu/bartlein/courses/geog490/week07-RMaps.html
library(sp)
library(maptools)
library(rgdal)
library(rgeos)
library(raster)
library(rasterVis)
library(classInt)
library(RColorBrewer)

##create folder on desktop with source files from Natural Earth Data
setwd("~/Desktop/RMaps")
### Read the shape files, including global coastlines, land, large lakes, and a graticule
##Set file names:
shape_path <- "/Users/nishattasnim/Desktop/RMaps/"
coast_shapefile <- paste(shape_path, "ne_50m_coastline/ne_50m_coastline.shp", sep="")
ocean_shapefile <- paste(shape_path, "ne_50m_ocean/ne_50m_ocean.shp", sep="")
admin0_shapefile <- paste(shape_path, "ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp", sep="")
admin1_shapefile <- paste(shape_path, "ne_50m_admin_1_states_provinces_lakes/ne_50m_admin_1_states_provinces_lakes.shp", sep="")
lakes_shapefile <- paste(shape_path, "ne_50m_lakes/ne_50m_lakes.shp", sep="")
bb_shapefile <- paste(shape_path, "ne_50m_graticules_all/ne_50m_wgs84_bounding_box.shp", sep="")
grat30_shapefile <- paste(shape_path, "ne_50m_graticules_all/ne_50m_graticules_30.shp", sep="")

##Read and plot coastlines, ocean polygons, land polygons, large lakes 
layer <- ogrListLayers(coast_shapefile)
layer <- ogrListLayers(ocean_shapefile)
layer <- ogrListLayers(admin0_shapefile)
layer <- ogrListLayers(admin1_shapefile)
layer <- ogrListLayers(lakes_shapefile)
layer <- ogrListLayers(grat30_shapefile)
layer <- ogrListLayers(bb_shapefile)

coast_lines <- readOGR(coast_shapefile, layer=layer)
ocean_poly <- readOGR(ocean_shapefile, layer=layer)
admin0_poly <- readOGR(admin0_shapefile, layer=layer)
admin1_poly <- readOGR(admin1_shapefile, layer=layer)
lakes_poly <- readOGR(lakes_shapefile, layer=layer)
grat30_lines <- readOGR(grat30_shapefile, layer=layer)
bb_poly <- readOGR(bb_shapefile, layer=layer)

##convert polygon to lines
bb_lines <- as(bb_poly, "SpatialLines")

plot(coast_lines, col="black")
plot(admin0_poly, bor="gray", add=TRUE)
plot(admin1_poly, bor="pink", add=TRUE)
plot(lakes_poly, bor="lightblue", add=TRUE)
plot(grat30_lines, col="lightblue", add=TRUE)
plot(bb_lines, col="black", add=TRUE)
plot(coast_lines, col="purple", add=TRUE)

plot(ocean_poly, col="lightblue")

unproj_proj4string <- proj4string(coast_lines)
unproj_proj4string

#set to Robinson projection system
unproj_crs <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
unproj_crs
robin_crs <- CRS("+proj=robin +lon_0=0w")
robin_crs
## transform everything except ocean 
bb_poly_proj <- spTransform(bb_poly, robin_crs)
coast_lines_proj <- spTransform(coast_lines, robin_crs)
admin0_poly_proj <- spTransform(admin0_poly, robin_crs)
admin1_poly_proj <- spTransform(admin1_poly, robin_crs)
lakes_poly_proj <- spTransform(lakes_poly, robin_crs)
grat30_lines_proj <- spTransform(grat30_lines, robin_crs)

# convert polygons to spatial lines
admin0_lines_proj <- as(admin0_poly_proj, "SpatialLines")
admin1_lines_proj <- as(admin1_poly_proj, "SpatialLines")
lakes_lines_proj <- as(lakes_poly_proj, "SpatialLines")
bb_lines_proj <- as(bb_poly_proj, "SpatialLines")

##plot projected files. Note: ORDER MATTERS 
plot(bb_poly_proj, col="gray95")
plot(coast_lines_proj, col="green", add=TRUE)
plot(admin0_lines_proj, col="lightblue", add=TRUE)
plot(admin1_lines_proj, col="lightblue", add=TRUE)
plot(lakes_lines_proj, col="blue", add=TRUE)
plot(grat30_lines_proj, col="gray", add=TRUE)
plot(coast_lines_proj, col="black", add=TRUE)
plot(bb_lines_proj, col="black", add=TRUE)

##write out projected shapefiles. First make the glrob_50m folder
outpath <- "/Users/nishattasnim/Desktop/RMaps/derived/glrob_50m"
outshape <- coast_lines_proj
outfile <- "glRob_50m_coast_lines"
outshapefile <- paste(outpath,outfile,sep="")
spdf <- data.frame(as.numeric(row.names(outshape)))
row.names(spdf) <- row.names(outshape)
outshape <- SpatialLinesDataFrame(outshape, spdf)
writeLinesShape(outshape, outshapefile, factor2char=TRUE) ##ignore warning 

outshape <- bb_poly_proj
outfile <- "glRob_50m_bb_poly"
outshapefile <- paste(outpath,outfile,sep="")
spdf <- data.frame(as.numeric(row.names(outshape)))
row.names(spdf) <- row.names(outshape)
outshape <- SpatialPolygonsDataFrame(outshape, spdf)
writePolyShape(outshape, outshapefile, factor2char=TRUE)

outshape <- bb_lines_proj
outfile <- "glRob_50m_bb_lines"
outshapefile <- paste(outpath,outfile,sep="")
spdf <- data.frame(as.numeric(row.names(outshape)))
row.names(spdf) <- row.names(outshape)
outshape <- SpatialLinesDataFrame(outshape, spdf)
writeLinesShape(outshape, outshapefile, factor2char=TRUE)

outshape <- admin0_poly_proj
outfile <- "glRob_50m_admin0_poly"
outshapefile <- paste(outpath,outfile,sep="")
spdf <- data.frame(as.numeric(row.names(outshape)))
row.names(spdf) <- row.names(outshape)
outshape <- SpatialPolygonsDataFrame(outshape, spdf)
writePolyShape(outshape, outshapefile, factor2char=TRUE)

outshape <- admin0_lines_proj
outfile <- "glRob_50m_admin0_lines"
outshapefile <- paste(outpath,outfile,sep="")
spdf <- data.frame(as.numeric(row.names(outshape)))
row.names(spdf) <- row.names(outshape)
outshape <- SpatialLinesDataFrame(outshape, spdf)
writeLinesShape(outshape, outshapefile, factor2char=TRUE)

outshape <- admin1_poly_proj
outfile <- "glRob_50m_admin1_poly"
outshapefile <- paste(outpath,outfile,sep="")
spdf <- data.frame(as.numeric(row.names(outshape)))
row.names(spdf) <- row.names(outshape)
outshape <- SpatialPolygonsDataFrame(outshape, spdf)
writePolyShape(outshape, outshapefile, factor2char=TRUE)

outshape <- admin1_lines_proj
outfile <- "glRob_50m_admin1_lines"
outshapefile <- paste(outpath,outfile,sep="")
spdf <- data.frame(as.numeric(row.names(outshape)))
row.names(spdf) <- row.names(outshape)
outshape <- SpatialLinesDataFrame(outshape, spdf)
writeLinesShape(outshape, outshapefile, factor2char=TRUE)

outshape <- lakes_poly_proj
outfile <- "glRob_50m_lakes_poly"
outshapefile <- paste(outpath,outfile,sep="")
spdf <- data.frame(as.numeric(row.names(outshape)))
row.names(spdf) <- row.names(outshape)
outshape <- SpatialPolygonsDataFrame(outshape, spdf)
writePolyShape(outshape, outshapefile, factor2char=TRUE)

outshape <- lakes_lines_proj
outfile <- "glRob_50m_lakes_lines"
outshapefile <- paste(outpath,outfile,sep="")
spdf <- data.frame(as.numeric(row.names(outshape)))
row.names(spdf) <- row.names(outshape)
outshape <- SpatialLinesDataFrame(outshape, spdf)
writeLinesShape(outshape, outshapefile, factor2char=TRUE)

outshape <- grat30_lines_proj
outfile <- "glRob_50m_grat30_lines"
outshapefile <- paste(outpath,outfile,sep="")
spdf <- data.frame(as.numeric(row.names(outshape)))
row.names(spdf) <- row.names(outshape)
outshape <- SpatialLinesDataFrame(outshape, spdf)
writeLinesShape(outshape, outshapefile, factor2char=TRUE)

# read QIITA coordinates
csvpath <- "/Users/nishattasnim/Desktop/RMaps/"
csvname <- "QIITACoordinates.csv"
qmap <- read.csv(paste(csvpath, csvname, sep=""))
qmap_pts <- data.frame(cbind(qmap$Lon,qmap$Lat,qmap$n))
names(qmap_pts) <- c("lon","lat","samples")
coordinates(qmap_pts) <- ~lon+lat
proj4string(qmap_pts) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
##set size and colors for QIITA coordinates
nsamp_cutpts <- c(5,50,500)
nsamp_colors <- c("green3","deepskyblue2","purple")
nsamp_cex <- c(0.5,0.5,0.5)
nsamp_num <- findInterval(qmap_pts$samples, nsamp_cutpts)+1
##project QIITA data
robin.crs <- CRS("+proj=robin +lon_0=0w")
qmap_pts.proj <- spTransform(qmap_pts, robin.crs)

#Read Soil Data
csvpath <- "/Users/nishattasnim/Desktop/RMaps/"
csvname <- "SoilDataQIITA.csv"
smap <- read.csv(paste(csvpath, csvname, sep=""))
deduped.data <- unique(smap)
smap_pts <- data.frame(cbind(deduped.data$Lon,deduped.data$Lat)) ##thought I had some duplicates
names(smap_pts) <- c("lon","lat")
coordinates(smap_pts) <- ~lon+lat
proj4string(smap_pts) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
smap_pts.proj <- spTransform(smap_pts, robin.crs)
samp_colors <- c("red")
samp_cex <- c(0.5)

###PLOT MAP
pdffile <- "qmap_nsamp.pdf"
pdf(paste(pdffile,sep=""), paper="letter", width=8, height=8)

plot(bb_poly_proj, col="gray90", bor="black", lwd=0.1)
plot(grat30_lines_proj, col="black", lwd=0.3, add=TRUE)
plot(admin0_poly_proj, col="white", bor="gray50", lwd=0.4, add=TRUE)
plot(coast_lines_proj, col="black", lwd=0.5, add=TRUE)
plot(bb_lines_proj, col="black", lwd=1.0, add=TRUE)
plot(qmap_pts.proj, pch=2, col="white", cex=nsamp_cex[nsamp_num], lwd=1.5, add=TRUE)
plot(qmap_pts.proj, pch=2, col=nsamp_colors[nsamp_num], cex=nsamp_cex[nsamp_num], lwd=0.6, add=TRUE)
plot(smap_pts.proj, pch=2, col="white", cex=nsamp_cex[nsamp_num], lwd=1.5, add=TRUE)
plot(smap_pts.proj, pch=2, col=samp_colors, cex = samp_cex,lwd=0.6, add=TRUE)

text(-17000000, 9100000, pos=4, cex=0.8, "Global Soil and Gut Studies in QIITA ")
legend(-17000000, -5000000, legend=c("< 5","5 - 50","> 500"), bg="white",
       title="Number of Gut Samples", pch=2, pt.lwd=0.6, col=nsamp_colors, cex=0.5)
dev.off()
