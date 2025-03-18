#master map for all coiba stuff
#library(rethinking)
library(mapview)
library(rgdal)
library(sf)
library(raster)
##### load camera traps
cam_ids <- read.csv(file="coiba_camtrap_ids_gps.csv") #generated with this script file https://github.com/bjbarrett/camtrap_coiba/blob/main/coiba_camtrap_gps_visualization.R
all_cams <- st_as_sf(cam_ids , coords = c("longitude", "latitude"), crs = 4326) #do it again if rading csv
all_cams$toolusers <- ifelse(all_cams$tool_site == 1, "Tool-using", "Non-tool-using")
all_cams_map <- mapview(all_cams , zcol = "toolusers"  )
all_cams_map
#filelist <- list.files(path="~/Dropbox/camtrap_coiba/map/gpx/bjb_raw_garmin_folders" , full.names=TRUE)
filelist <- list.files(path="map/gpx/bjb_raw_garmin_folders" , full.names=TRUE)

bjb_cm_se_sw_thru_20170324  <- readOGR(dsn =filelist[7] , layer="tracks")
mapview(bjb_cm_se_sw_thru_20170324)

for(i in 1:length(filelist)){
  x  <- readOGR(dsn =filelist[i] , layer="tracks")
  print(mapview(x))
}
filelist

##### load stream gps points
#####here is how we turn streams to lines and trim

#quebrada catarata

s1 <- readOGR(dsn ="map/gpx/bjb_raw_garmin_folders/Track_CAARATA QUEB2353.gpx" , layer="track_points")
mapview(s1)
s2 <- s1[3:6,]
coord<-as.data.frame(coordinates(s2))
queb_cat <- Line(as.data.frame(coordinates(s2)))
queb_cat <-Lines(list(queb_cat),ID="queb_catarata")

###discovery/watersource
s1 <- readOGR(dsn ="map/gpx/bjb_raw_garmin_folders/Track_DISCO QUEB 121828.gpx" , layer="track_points")
mapview(s1)
coord<-as.data.frame(coordinates(s1))
queb_watersource <- Line(as.data.frame(coordinates(s1)))
queb_watersource <-Lines(list(queb_watersource),ID="queb_watersource")

###bbc/alicia
s1 <- readOGR(dsn ="map/gpx/bjb_raw_garmin_folders/Track_QUEB BBC28 145539.gpx", layer="track_points")
mapview(s1)
s2 <- s1[4:108,]
coord<-as.data.frame(coordinates(s2))
queb_alicia <- Line(as.data.frame(coordinates(s2)))
queb_alicia <-Lines(list(queb_alicia),ID="queb_alicia")

###posa
s1 <- readOGR(dsn ="map/gpx/bjb_raw_garmin_folders/Track_QUEB POSA8 135357.gpx" , layer="track_points")
mapview(s1[4:30,])
s2 <- s1[4:30,]
coord<-as.data.frame(coordinates(s2))
queb_posa <- Line(as.data.frame(coordinates(s2)))
queb_posa <-Lines(list(queb_posa),ID="queb_posa")

streams <- SpatialLines(list(queb_cat,queb_watersource,queb_alicia,queb_posa) , proj4string = s2@proj4string )
mapview(streams , color="blue" , lw=3)

all_streams_map <- mapview(streams, color="blue" , lw=3)

#####lets add in tool sites
tools_w201707201901 <- readOGR(dsn = "map/gpx/cleaned/Tool Use Sites July 2017 Jan 2018.GPX", layer="waypoints")
tools_w201803 <- readOGR(dsn = "map/gpx/cleaned/Tool Use Sites Mar 2018.GPX", layer="waypoints")
tools_w201807 <- readOGR(dsn = "map/gpx/cleaned/Tool Use Site July 2018.GPX", layer="waypoints")
tools_w201901 <- readOGR(dsn = "map/gpx/cleaned/Tool Use Sites Jan 2019.GPX", layer="waypoints")
tools_w201903 <- readOGR(dsn = "map/gpx/cleaned/Tool Sites Mar 2019.GPX", layer="waypoints")

all_tools_map <- mapview(tools_w201707201901 , col.region="red" , cex=3) + mapview(tools_w201803 , col.region="orange" , cex=3) + mapview(tools_w201807 , col.region="yellow" , cex=3) + mapview(tools_w201901 , col.region="violet" , cex=3) + mapview(tools_w201903 , col.region="blue" , cex=3)

all_tools_map + all_streams_map + all_cams_map


### add in almendras
#note, i did not get to finish punta ursula with camera trapsl
almendras <- readOGR(dsn = "map/gpx/cleaned/T_catappa_points.GPX", layer="waypoints")
almendras$name
#mccir <- readOGR(dsn ="~/Dropbox/Capuchin Monkeys, Coiba National Park_July2018/GPS Points/MCCinreachJan2020.GPX" , layer="waypoints")##inreaches need a space on xml code in first line removed to code
most_almendras_map <- mapview(almendras , col.regions="green" , alpha.regions = 0.2 , cex=1.2)

all_tools_map + all_streams_map + all_cams_map + most_almendras_map #map of all trees, tools, camerasm etc

all_tools_map + all_streams_map +  most_almendras_map #map of all trees, tools, camerasm etc


##grids? 

##########new way from Kate#######
#maybe use gridshapes

grid_shapes <- rbind(tools_w201707201901,tools_w201803,tools_w201807) # define grid area using a subset of tools w/ same number columns
str(grid_shapes)
e <- as(raster::extent(min(grid_shapes@coords[,1]), max(grid_shapes@coords[,1]), min(grid_shapes@coords[,2]), max(grid_shapes@coords[,2])), "SpatialPolygons")
proj4string(e) <- crs(grid_shapes)
e3 <- spTransform(e, CRS("+init=EPSG:32616"))
ebuf <- buffer(e3, width = 1500) #add 500 m buffer
e2 <- st_as_sf(e)
e2b <- st_as_sf(ebuf)
grid_1000m <- st_make_grid(e2b, square = T, cellsize = c(1000, 1000) ) %>% 
  st_sf() #1000m grid
grid_500m <- st_make_grid(e2b, square = T, cellsize = c(500, 500) ) %>% 
  st_sf() #500m grid
grid_250m <- st_make_grid(e2b, square = T, cellsize = c(250, 250)) %>% # the grid, covering bounding box
  st_sf() # not really required, but makes the grid nicer to work with later
mapview(grid_250m , col.regions="white")  + all_tools_map + all_streams_map + all_cams_map + most_almendras_map
mapview(grid_500m , col.regions="white") + all_tools_map + all_streams_map + all_cams_map + most_almendras_map
mapview(grid_1000m , col.regions="white") + all_tools_map + all_streams_map + all_cams_map + most_almendras_map

#####brendan needs to pick up here before publishing
pois <- readOGR(dsn = "map/gpx/cleaned/Coiba POI.GPX", layer="waypoints")
pois2 <- readOGR(dsn = "map/gpx/cleaned/Points of Interest.GPX", layer="waypoints")
mapview(pois)
mapview(pois2)
pois$name
pois2$name

pois_4_map_and_gps <- pois2[pois2$name!="SE Landing Spo",]
pois_4_map_and_gps$name

mono_encounters <-
  rbind(
  pois[str_detect(pois$name, "TROPA"), ] ,
  pois[str_detect(pois$name, "Cebus track s"), ] ,
  pois[str_detect(pois$name, "MONO"), ] , 
  pois[str_detect(pois$name, "MONKEY"), ] , 
  pois[str_detect(pois$name, "SOUTHERN GROUP"), ]
  )

mapview(mono_encounters)

####tracks we walked

#bind in raster package can deal with combining
#2018.03.28
trk_20180328<- readOGR(dsn =filelist[2] , layer="tracks")
trk_20180328$name
print(mapview(tracks_4_site))
#2017.03.24 bjb and claudio initial cross
trk_20170324  <- readOGR(dsn =filelist[7] , layer="tracks")
#2017.07.29 BJB PC cross
trk_20170729  <- readOGR(dsn =filelist[8] , layer="tracks")
#2017.12.20 cmm cross
trk_20171220   <- readOGR(dsn =filelist[9] , layer="tracks")

mapview(trk_20180328) + mapview(trk_20170324) + mapview(trk_20170729) + mapview(trk_20171220 )

mccir_tr <- readOGR(dsn ="map/gpx/gps_dumps/MCCinreachJan2020.GPX" , layer="tracks")##inreaches need a space on xml code in first line removed to code
mapview(mccir_tr , zcol="name") + mapview(trk_20180328) + mapview(trk_20170324) + mapview(trk_20170729) + mapview(trk_20171220 )
mccir_tr$name
all_tracks <- bind(mccir_tr ,trk_20180328 ,trk_20170324,trk_20170729 ,trk_20171220)
mapview(all_tracks , zcol="name") + mapview(grid_250m , col.regions="white") + all_tools_map + all_streams_map + mapview(mono_encounters)
####tracks
# mccir_wp <- readOGR(dsn ="map/gpx/gps_dumps/MCCinreachJan2020.GPX" , layer="waypoints")##inreaches need a space on xml code in first line removed to code
# mccir_tr <- readOGR(dsn ="map/gpx/gps_dumps/MCCinreachJan2020.GPX" , layer="tracks")##inreaches need a space on xml code in first line removed to code
# mapview(mccir_wp) + mapview(mccir_tr)
# 
# mccir_wp$name
# almendras$name
# mccir_wp@data$name
# mccir_wp@coords[mccir_wp@data$name=="TC 107",]
#useful resource https://cmerow.github.io/RDataScience/04_Spatial.html


#gpx to write-- camera traps, all tracks, almendras, pois, tools, 

all_cams
writeOGR(pois_4_map_and_gps, dsn = "map/gpx/coiba_pois_4_map_and_gps.GPX", dataset_options="GPX_USE_EXTENSIONS=yes",layer="waypoints",driver="GPX", overwrite_layer = T)
writeOGR(all_cams, dsn = "map/gpx/all_cameratraps_2020.GPX", dataset_options="GPX_USE_EXTENSIONS=yes",layer="waypoints",driver="GPX", overwrite_layer = T)
#fix above to sp

