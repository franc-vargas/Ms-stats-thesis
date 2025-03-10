library(sp)
library(sf)
library(stars)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(leaflet)
library(igraph)
library(gdistance)

# Transition ------------------------------------------------------------
library(raster)
library(sp)

columnas = 6000
filas = 0.61*columnas
r_aysen2 <- raster::raster(ncol = columnas, nrow = filas)
r_aysen2@extent <- raster::extent(c(-75.5, -72, -46.5, -43.5))
r_aysen2 <- fasterize::fasterize(chilemapas::generar_regiones() |> dplyr::filter(codigo_region == 11), r_aysen2)


r_aysen2[is.na(r_aysen2)] <- -999 # this turns all landmass to missing
r_aysen2[r_aysen2>-999] <- NA # assign unit cost to all grid cells in water
r_aysen2[r_aysen2==-999] <- 1 # make landmass impassable
r_aysen2
rm(rs2)
raster::plot(r_aysen2)
raster::plot(rs2)

load(file = "shape_regions.RData")

readr::write_rds(rs.transition2, "1_data/aysen_transition.matrix.11.rds")

## CONTINUACIÓN DEL PROCESO

rs.transition2 <- readr::read_rds("1_data/aysen_transition.matrix.11.rds")
rs.transition2 <- gdistance::transition(r_aysen2, mean, directions = 16)

rs.transition2@crs <- 
  sp::CRS("+ellps=WGS84 +proj=longlat +datum=WGS84 +no_defs +units=m")

rs.transition2 <- gdistance::geoCorrection(rs.transition2, "c")
readr::write_rds(rs.transition2, "1_data/aysen_transition.matrix.11.rds")

sites_coords <- shape_aysen |> dplyr::select(N_CODIGOCE, coordinates)
sites_coords$lon <- sites_coords$coordinates[,1]
sites_coords$lat <- sites_coords$coordinates[,2]
sites_coords <- sites_coords |> dplyr::select(lon, lat, N_CODIGOCE)
sites_coords <- sites_coords|> as.matrix()

spatpoints <- sp::SpatialPoints(
  sites_coords[,1:2], 
  proj4string=CRS("+ellps=WGS84 +proj=longlat +datum=WGS84 +no_defs +units=m")
  )

## REVISIÓN DE PUNTOS QUE PODRÍAN CAUSAR PROBLEMAS POR ESTAR MUY CERCANOS A TIERRA
plot(r_aysen2)
points(spatpoints[spatpoints@coords[,1] <= -74.3 &spatpoints@coords[,1] > -74.5
                  & spatpoints@coords[,2] <= -45.5], pch = 19, cex = 0.1)
rownames(sites_coords) <- sites_coords[,3]
spatpoints2 <- spatpoints[
  spatpoints@coords[,1] <= -74.3 &spatpoints@coords[,1] > -74.5
                          & spatpoints@coords[,2] <= -45.5]

dist.m.aysen <- gdistance::shortestPath(rs.transition2, 
                                        sites_coords[,1:2], 
                                        sites_coords[,1:2])


plot(r_aysen2)
points(spatpoints2[rownames(spatpoints2@coords) == 110953,], pch = 19, cex = 0.1)
moved_site <- spatpoints2[rownames(spatpoints2@coords) == 110953,]
moved_site2 <- spatpoints2[rownames(spatpoints2@coords) == 110953,]
moved_site2@coords[1] <- moved_site2@coords[1]+ 0.11

raster::pointDistance(moved_site, moved_site2, lonlat = T)


pointDistance(c(0, 0), a, lonlat=TRUE)
a <- cbind(c(1,5,55,31),c(3,7,20,22))

sf::st_distance(moved_site, moved_site2, lonlat = T)
moved_site[2,] <-spatpoints2[rownames(spatpoints2@coords) == 110953,] +0.11


dist.m1 <- readr::read_rds("1_data/aysen_distance.matrix_raw.rds") |> as.data.frame()

dim(dist.m.aysen)

dist.aysen.df <- dist.m.aysen |> as.data.frame()

inf_cols <- dist.aysen.df |> 
  dplyr::select(where(~sum(is.infinite(.x)) == nrow(dist.aysen.df)-1))
dim(inf_cols)
inf_cols.sites <- colnames(inf_cols)
inf_cols.sites <- gsub("X", "", inf_cols.sites) |> as.numeric()
max(dist.m.aysen)
is.matrix(inf_cols.sites)
shape_aysen$lon <- shape_aysen$coordinates[,1]
shape_aysen$lat <- shape_aysen$coordinates[,2]
shape_aysen |> glimpse()

## SE REMUEVE UN CODIGO DUPLICADO
shape_aysen2 <- shape_aysen |> filter(N_CODIGOCE == 110005)
shape_aysen <- shape_aysen |> filter(N_CODIGOCE != 110005)
shape_aysen2 <- shape_aysen2[1,]
shape_aysen <- rbind(shape_aysen, shape_aysen2)

sites_coords <- shape_aysen |> dplyr::select(N_CODIGOCE, coordinates)


sites_coords$lon <- sites_coords$coordinates[,1]
sites_coords$lat <- sites_coords$coordinates[,2]
sites_coords <- sites_coords |> dplyr::select(lon, lat, N_CODIGOCE)

sites_coord_inf <- sites_coords |> dplyr::filter(N_CODIGOCE %in% inf_cols.sites)
sites_coords |> filter(N_CODIGOCE == 110005)
sites_coord_inf <- sites_coord_inf |> select(N_CODIGOCE, lon, lat)
## MODIFICACIÓN DE SITIOS PROBLEMA QUE ESTAN MUY CERCANOS A TIERRA
{
sites_coord_inf2 <- sites_coord_inf
sites_coord_inf2[1,1:2] <- c(sites_coord_inf2[1,1]-0.02, sites_coord_inf2[1,2]+ 0.04)
sites_coord_inf2[2,1:2] <- c(sites_coord_inf2[2,1]-0.008, sites_coord_inf2[2,2]+ 0.01)
sites_coord_inf2[3,1:2] <- c(sites_coord_inf2[3,1]+0.11, sites_coord_inf2[3,2])
sites_coord_inf2[4,1:2] <- c(sites_coord_inf2[4,1], sites_coord_inf2[4,2]+0.029)
sites_coord_inf2[5,1:2] <- c(sites_coord_inf2[5,1]+0.007, sites_coord_inf2[5,2])
sites_coord_inf2[6,1:2] <- c(sites_coord_inf2[6,1], sites_coord_inf2[6,2]-0.007)
sites_coord_inf2[7,1:2] <- c(sites_coord_inf2[7,1], sites_coord_inf2[7,2]-0.007)
sites_coord_inf2[8,1:2] <- c(sites_coord_inf2[8,1]-0.0052, sites_coord_inf2[8,2])
sites_coord_inf2[9,1:2] <- c(sites_coord_inf2[9,1], sites_coord_inf2[9,2]+0.006)
sites_coord_inf2[10,1:2] <- c(sites_coord_inf2[10,1], sites_coord_inf2[10,2]+0.006)

sites_coord_inf3 <- sites_coord_inf2[c(1:9, 11),]
spatpoints <- sp::SpatialPoints(
  sites_coord_inf3[,1:2], 
  proj4string=CRS("+ellps=WGS84 +proj=longlat +datum=WGS84 +no_defs +units=m")
  )

}

## REVISIÓN Y VALIDACIÓN DE MODIFICACIONES
# Punto 1
plot(rs2, xlim = c(-73, -72 ), ylim = c(-45.7, -45.1))
points(spatpoints[1,], pch = 19, cex = 0.1)

# Punto 2
plot(r_aysen2, xlim = c(-75, -74 ), ylim = c(-46, -45.5))
points(spatpoints[2,], pch = 19, cex = 0.1)

# Punto 3
plot(r_aysen2, xlim = c(-75, -72 ), ylim = c(-46, -43.5))
points(spatpoints[3,], pch = 19, cex = 0.1)

# Punto 4
plot(r_aysen2, xlim = c(-74.5, -73.5 ), ylim = c(-45.5, -45.3))
points(spatpoints[4,], pch = 19, cex = 0.1)

# Punto 5
plot(r_aysen2, xlim = c(-74, -73.5 ), ylim = c(-44.8, -44.5))
points(spatpoints[5,], pch = 19, cex = 0.1)

# Punto 6
plot(r_aysen2, 
     xlim = c(-73.5, -73 ), ylim = c(-44, -43.5))
points(spatpoints[6,], pch = 19, cex = 0.1)

# Punto 7
plot(r_aysen2,
     xlim = c(-73.5, -72.5 ), ylim = c(-45.5, -45.1))
points(spatpoints[7,], pch = 19, cex = 0.1)

# Punto 8
plot(r_aysen2, 
     xlim = c(spatpoints@coords[8,1]-0.1, spatpoints@coords[8,1]+0.1 ), 
     ylim = c(spatpoints@coords[8,2]-0.1, spatpoints@coords[8,2]+0.1))
points(spatpoints[8,], pch = 19, cex = 0.1)

# Punto 9
index <- 9
plot(r_aysen2, 
     xlim = c(spatpoints@coords[index,1]-0.1, spatpoints@coords[index,1]+0.1 ), 
     ylim = c(spatpoints@coords[index,2]-0.1, spatpoints@coords[index,2]+0.1))
points(spatpoints[index,], pch = 19, cex = 0.1)

# Punto 10
index <- 10
plot(r_aysen2, 
     xlim = c(spatpoints@coords[index,1]-0.1, spatpoints@coords[index,1]+0.1 ), 
     ylim = c(spatpoints@coords[index,2]-0.1, spatpoints@coords[index,2]+0.1))
points(spatpoints[index,], pch = 19, cex = 0.1)


sites_coord_inf2.df <- sites_coord_inf2 |> as.data.frame()
sites_coords
sites_coords |> pull(N_CODIGOCE) |> unique() |> length()
sites_coords2 <- sites_coords
sites_coords2 <- sites_coords2 |> 
  filter(!(N_CODIGOCE %in% sites_coord_inf2.df$N_CODIGOCE))
333+153
sites_coords2 <- rbind(sites_coords2, sites_coord_inf3)
sites_coords2 |> nrow()
sites_coord_inf3 |> nrow()


# Aysen Distance Calculations -----------------------------------------------
sites_coords.matrix <- as.matrix(sites_coords2)

spatpoints[indice,1]
sites_coords.matrix[467,] 

plot(rs2, xlim = c(-75, -72 ), ylim = c(-46, -43.5))
points(sites_coords.matrix[467,1] ,
       sites_coords.matrix[467,2],
       pch = 19, cex = 0.1)
rownames(sites_coords.matrix) <- sites_coords.matrix[,3]
dist.m.aysen.corrected

readr::write_rds(sites_coords.matrix, "1_data/sites_coords.matrix.Fixed.rds")

dist.m.aysen.corrected <- 
  gdistance::shortestPath(rs.transition2, 
                          sites_coords.matrix[,1:2],
                          sites_coords.matrix[,1:2])
readr::write_rds(dist.m.aysen.corrected, 
                 "1_data/Aysen_DistanceMatrix.rds")

dist.aysen.df.2 <- dist.m.aysen.corrected |> as.data.frame()

## REVISIÓN / VALIDACIÓN SI ALGUNA COLUMNA QUEDA COMO INFINITO
inf_cols2 <- dist.aysen.df.2 |> 
  dplyr::select(where(~sum(is.infinite(.x)) == nrow(dist.m.aysen.corrected)-1))
dim(inf_cols2)
inf_cols2.sites <- colnames(inf_cols2)
inf_cols2.sites <- gsub("X", "", inf_cols2.sites) |> as.numeric()


