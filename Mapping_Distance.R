##Data import and visualization
library(sp)
library(shapefiles)
library(terra)
library(sf)
library(raster)
library(plyr)
library(RColorBrewer)
library(geodata)
library(geosphere) 

# Get map data for countries of interest (and DRC at the province level), then assign centroids to
# center of each location polygon and calculate a matrix of distances.

# country_codes("Africa")
# drc <- gadm(country = "Democratic Republic of the Congo", level = 1, resolution = 2,
#            path = "../data/maps/")
# 
# ang <- gadm(country = "Angola", level = 0, resolution = 2,
#             path = "../data/maps/")
# 
# bur <- gadm(country = "Burundi", level = 0, resolution = 2,
#             path = "../data/maps/")
# 
# ken <- gadm(country = "KEN", level = 0, resolution = 2,
#             path = "../data/maps/")
# 
# rwa <- gadm(country = "RWA", level = 0, resolution = 2,
#             path = "../data/maps/")
# 
# zam <- gadm(country = "ZMB", level = 0, resolution = 2,
#             path = "../data/maps/")
# 
# zim <- gadm(country = "ZWE", level = 0, resolution = 2,
#             path = "../data/maps/")
# 
# car <- gadm(country = "Central African Republic", level = 0, resolution = 2,
#             path = "../data/maps/")
# 
# sud <- gadm(country = "Sudan", level = 0, resolution = 2,
#             path = "../data/maps/")
# 
# ssud <- gadm(country = "South Sudan", level = 0, resolution = 2,
#             path = "../data/maps/")
# 
# plot(drc, lwd = 2)
# plot(bur, add=TRUE, border='green')
# plot(ang, add=TRUE, border='green')
# plot(ken, add=TRUE, border='green')
# plot(rwa, add=TRUE, border='green')
# plot(zam, add=TRUE, border='green')
# plot(zim, add=TRUE, border='green')
# plot(car, add=TRUE, border='green')
# plot(sud, add=TRUE, border='green')
# plot(ssud, add=TRUE, border='green')
# 
# all <- rbind(drc,ang,bur,ken,rwa,zam,zim,car,sud,ssud)
# 
# plot(all, lwd=2)
# plot(centroids, add=TRUE)
# 
# 
# centroids <- centroids(all)
# centroids.df <- as.data.frame(centroids)
# 
# distance <- distance(centroids, unit="km")


# Create map objects for each country of interest. DRC is at the province level and
# other countries are all at the national level.

country_codes("Africa")
drc <- gadm(country = "Democratic Republic of the Congo", level = 1, resolution = 2, path = "../data/maps/")
ang <- gadm(country = "Angola", level = 0, resolution = 2, path = "../data/maps/")
bur <- gadm(country = "Burundi", level = 0, resolution = 2, path = "../data/maps/")
ken <- gadm(country = "KEN", level = 0, resolution = 2, path = "../data/maps/")
rwa <- gadm(country = "RWA", level = 0, resolution = 2, path = "../data/maps/")
zam <- gadm(country = "ZMB", level = 0, resolution = 2, path = "../data/maps/")
zim <- gadm(country = "ZWE", level = 0, resolution = 2, path = "../data/maps/")
car <- gadm(country = "Central African Republic", level = 0, resolution = 2, path = "../data/maps/")
sud <- gadm(country = "Sudan", level = 0, resolution = 2, path = "../data/maps/")
ssud <- gadm(country = "South Sudan", level = 0, resolution = 2, path = "../data/maps/")
uga <- gadm(country = "Uganda", level = 0, resolution = 2, path = "../data/maps/")

# Convert to sf objects
drc_sf <- st_as_sf(drc)
ang_sf <- st_as_sf(ang)
bur_sf <- st_as_sf(bur)
ken_sf <- st_as_sf(ken)
rwa_sf <- st_as_sf(rwa)
zam_sf <- st_as_sf(zam)
zim_sf <- st_as_sf(zim)
car_sf <- st_as_sf(car)
sud_sf <- st_as_sf(sud)
ssud_sf <- st_as_sf(ssud)
uga_sf <- st_as_sf(uga)

# Combine all country sf objects into one
all_sf <- rbind(ang_sf, bur_sf, ken_sf, rwa_sf, zam_sf, zim_sf, car_sf, sud_sf, ssud_sf, uga_sf)


# To merge DRC data to the rest of the country data, some changes need to be made.
# Add province data element to all country data, and assign NA.
all_sf$Province <- NA

# For the DRC data, create Province data element with Province names (from Name_1)
drc_sf$Province <- drc_sf$NAME_1

# Limit DRC data to only important variables (also present in all country data)
drc_sf2 <- drc_sf[,c("GID_0","COUNTRY","Province","geometry")]

# Bind together all country data and DRC data
all_sf2 <- rbind(all_sf, drc_sf2)

# Compute centroids for each location, so we can then calculate distances
centroids <- st_centroid(all_sf2)

# Extract centroid coordinates
centroid_coords <- st_coordinates(centroids)

# Compute pairwise distances in meters
distance_matrix <- distm(centroid_coords, fun = distHaversine)


# To organize matrix- assign names from DRC provinces and country names
drc_province_names <- unique(drc_sf$NAME_1)
country_names <- unique(all_sf$COUNTRY)

location_names <- c(country_names,drc_province_names)

# Assign names to matrix
rownames(distance_matrix) <- location_names
colnames(distance_matrix) <- location_names

# Print distance matrix
print(distance_matrix)

write.csv(distance_matrix, file="Data/distance_matrix.csv")

# Using Kramer et al: crossing core borders: 
# 1/ 1+ e^(B0 + B1*(dij/(pi*pj)^B2)
# If crossing a core border: B3 * (1/ 1+ e^(B0 + B1*(dij/(pi*pj)^B2))
# Where: B0= 5.792, B1 = 105.7, B2=0.186, B3= 0.15
# Can be found in Kramer et al Supplement
# 

# Using Kramer et al: crossing country borders: 
# 1/ 1+ e^(B0 + B1*(dij/(pi*pj)^B2)
# If crossing a core border: B3 * (1/ 1+ e^(B0 + B1*(dij/(pi*pj)^B2))
# Where: B0= 5.166, B1 = 157.1, B2=0.189, B3= 0.507
# Can be found in Kramer et al Supplement
# 
