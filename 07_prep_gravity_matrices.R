
## clear environment
rm(list=ls())

## Set repo location and indicator group
user               <- Sys.info()['user']
core_repo          <- sprintf('FILEPATH/',user)
vax_repo           <- sprintf('FILEPATH/',user)
measles_repo       <- sprintf('FILEPATH/',user)
remote             <- 'origin'
branch             <- 'develop'
pullgit            <- FALSE

## sort some directory stuff
commondir      <- sprintf('FILEPATH')
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0('FILEPATH/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

# Custom load indicator-specific functions
source(paste0('FILEPATH/misc_vaccine_functions.R'))
source(paste0('FILEPATH/gravity_functions.R'))


## load iso3 and shapefile version
#load_from_parallelize()
shp_version <- '2020_05_21'
iso3 <- 'ETH'
processing_date <- '2023_01_03'

if(iso3 == 'ETH') ctry_name <- 'Ethiopia'
######

ad0_shp <- sf::st_read(paste0("FILEPATH",shp_version,"/lbd_standard_admin_0.shp"))
ad0_shp <- subset(ad0_shp, ADM0_NAME == ctry_name)  

# Load an admin 2 shape and grab centroids
ad2_shp <- sf::st_read(paste0("FILEPATH",shp_version,"/lbd_standard_admin_2.shp"))
ad2_shp <- subset(ad2_shp, ADM0_NAME == ctry_name)



##### test pop weighted centroids
library(data.table)
library(spatialEco, lib.loc = 'FILEPATH')

library(raster)
library(rgdal)
dat <- raster('FILEPATH/worldpop_total_1y_2019_00_00.tif')
data_points <- rasterToPoints(dat, spatial=TRUE)

# intersect with polygons
point.in.poly <- function(x, y, sp = TRUE, duplicate = TRUE, ...) {
  if(!any(class(x)[1] == c("SpatialPoints", "SpatialPointsDataFrame", "sf"))) {
    stop("x is not a suitable point feature object class") }
  if(!any(class(y)[1] == c("SpatialPolygons", "SpatialPolygonsDataFrame", "sf"))) {
    stop("y is not a suitable polygon feature object class") }
  if(any(class(x) == "sfc")) { x <- sf::st_sf(x) }
  if(duplicate == FALSE) {
    if(!any(class(x)[1] == c("SpatialPoints", "SpatialPointsDataFrame"))) {
      x <- methods::as(x, "Spatial")
      if(dim(x@data)[2] == 0) stop("There are no attributes associated with points")
    }
    if(!any(class(y)[1] == c("SpatialPolygons", "SpatialPolygonsDataFrame"))) {
      y <- methods::as(y, "Spatial")
      if(dim(y@data)[2] == 0) stop("There are no attributes associated with polygons")
    }
    o <- sp::over(x, y, returnList = TRUE)
    m <- max(unlist(lapply(o, nrow)))
    ids <- row.names(y)
    xy <- data.frame(t(sapply(1:length(o),
                              function(i) c(ids[i], c(o[[i]][,1], rep(NA, m))[1:m])
    )))
    colnames(xy) <- c("p",paste0("pid", 1:m))
    x@data <- data.frame(x@data, xy)
    if( sp == FALSE ) sf::st_as_sf(x)
    return( x )
  } else {
    if(any( class(x) == c("SpatialPoints", "SpatialPointsDataFrame"))) {
      x <- sf::st_as_sf(x)
    }
    if(any(class(y) == "sfc")) { x <- sf::st_sf(y) }
    if(any( class(y) == c("SpatialPolygons", "SpatialPolygonsDataFrame"))) {
      y <- sf::st_as_sf(y)
    }
    if(dim(x)[2] == 1) x$pt.ids <- 1:nrow(x)
    if(dim(y)[2] == 1) y$poly.ids <- 1:nrow(y)
    o <- sf::st_join(x, y, ...)
    if( sp ) o <- methods::as(o, "Spatial")
    # if( sp ) o <- sf::as_Spatial(o)
    return( o )
  }
}
sf_use_s2(FALSE)

grid_centroids <- point.in.poly(data_points, ad2_shp)

# calculate weighted centroids
grid_centroids <- as.data.frame(grid_centroids)
w.centroids <- setDT(grid_centroids)[, lapply(.SD, weighted.mean, w=worldpop_total_1y_2019_00_00), by=ADM2_CODE, .SDcols=c('coords.x1','coords.x2')]
w.centroids
w.centroids.sf <- as(w.centroids, 'Spatial')
# Generate centroids
ad2_shp3 <- subset(ad2_shp, ADM2_CODE %in% c(setdiff(ad2_shp$ADM2_CODE,unique(w.centroids$ADM2_CODE))))

cents_special <- sf::st_centroid(ad2_shp)
cents_special$X
centroids = cents_special

point_locations <- data.table(sf::st_coordinates(centroids))
names(point_locations) <- c("X_COORD", "Y_COORD")

shape = ad0_shp

# Adapted from https://medium.com/@abertozz/mapping-travel-times-with-malariaatlas-and-friction-surfaces-f4960f584f08

# Pull friction surface
message("Retrieving friction surface")

friction <- raster('FILEPATH/2020_motorized_friction_surface2.geotiff')

masked <- mask(friction, shape)
cropped <- crop(masked, extent(shape))

friction <- copy(cropped)

message("Generating transition matrix")
T <- gdistance::transition(friction, function(x) 1/mean(x), 8) 
T.GC <- gdistance::geoCorrection(T)

# Convert to SpatialPointsDataFrame for this code to run
point_locations$poly_id <- 1:nrow(point_locations)
coordinates(point_locations) <- ~ X_COORD + Y_COORD

message("Looping over all points to generate travel time matrix")

# Find travel times between all sets of points

matrix_rows <- pbapply::pblapply(1:nrow(point_locations), function(i) {
  
  # Use point (i) as the target point
  target_point <- as.matrix(point_locations[i, ]@coords)
  
  # Generate an access raster using accumulated cost algorithm & friction surface
  access_raster <- gdistance::accCost(T.GC, target_point)
  
  # Extract the values for all other coordinates
  travel_times <- raster::extract(access_raster, point_locations)
  
  return(travel_times)
})

tt_matrix <- do.call(rbind, matrix_rows) %>% as.matrix

# Save file
saveRDS(tt_matrix, paste0('FILEPATH/',iso3,'_big6.RDS'))

# ##### now let's try the gravity matrix part......
# library(viridis)
# 
# ## gravity matrix formula:
# 
# ## M_ij = k ((P_i * P_j) / D_ij)
# ## P_i = population of i
# ## P_j = population of j
# ## D_ij = "distance" or travel time from i to j
# ## k = some constant ... 

load(paste0('FILEPATH/', iso3,'/',processing_date,'/', iso3,'_processed_population_by_epiweek_measles_model_age_bins_by_adm2.RData'))
pop_array <- pop_array_red

P <- pop_array[2094,,] ## taking the "2019 population"
P_district <- colSums(P) ## summing over all ages

G <- matrix(NA, length(P_district), length(P_district))
for(i in 1:length(P_district)){
  for(j in 1:length(P_district)){
    G[i,j] = (P_district[i] * P_district[j]) / tt_matrix[i,j]
  }
}
G <- ifelse(G == "Inf", 0, G)

save(G, file=paste0('FILEPATH/', iso3,'/',processing_date,'/', iso3,'_processed_gravity_matrix.RData'))

