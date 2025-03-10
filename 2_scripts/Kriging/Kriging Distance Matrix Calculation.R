# Libraries
library(sp)
library(sf)
library(stars)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(leaflet)
library(igraph)
library(gdistance)


# Object - Loading --------------------------------------------------------



##### MOVIMIENTO LEVE DE PUNTOS ######
# Para sacarlos fuera del polígono
aysen_polygon <- sf::as_Spatial(chilemapas::generar_regiones() |> 
                                  dplyr::filter(codigo_region == 11))

# Aysen Calculations
library(raster)
rdr <- raster::raster(system.file("rs_aysen_v2.tif", package="raster"))
load(file = "shape_regions.RData")
rm(shape_lagos)

shape_aysen2 <- shape_aysen |> filter(N_CODIGOCE == 110005)
shape_aysen <- shape_aysen |> filter(N_CODIGOCE != 110005)
shape_aysen2 <- shape_aysen2[1,]
shape_aysen <- rbind(shape_aysen, shape_aysen2)

sites_coords <- shape_aysen |> dplyr::select(N_CODIGOCE, coordinates)


sites_coords$lon <- sites_coords$coordinates[,1]
sites_coords$lat <- sites_coords$coordinates[,2]
sites_coords <- sites_coords |> dplyr::select(lon, lat, N_CODIGOCE)
pts <- sites_coords
coordinates(pts) <- ~lon+lat
pts_modified <- SpatialPointsDataFrame(pts, data = data.frame(id = pts@data$N_CODIGOCE))
crs(pts_modified) <- crs(aysen_polygon)
over_points_sites <- extract(aysen_polygon, pts_modified)
pts_modified@data$in_polygon <- over_points_sites$codigo_region
pts_modified@data$in_polygon <- ifelse(!is.na(pts_modified@data$in_polygon),1, 0)
pts_modified@data$colour <- ifelse(pts_modified@data$in_polygon == 1,"darkred", "darkgreen")
plot(pts_modified, col = pts_modified@data$colour)
{
sites_clear <- subset(pts_modified, in_polygon == 0)
sites_tochange <- subset(pts_modified, in_polygon == 1)
sites_coords_modify <- sites_coords |> filter(N_CODIGOCE %in% sites_tochange@data$id)
sites_coords_clear <- sites_coords |> filter(!N_CODIGOCE %in% sites_tochange@data$id)
plot(aysen_polygon, xlim = c(-74.45376, -73), ylim = c(-46, -45))
points(sites_tochange[1,], pch = 2, col = "darkred")
sites_coords_modify[
  sites_coords_modify$N_CODIGOCE == 110804,
  c(1:2)] <- c(-74.45376 - 0.04,-45.73357)
points(sites_coords_modify$lon[1],sites_coords_modify$lat[1])

plot(aysen_polygon, xlim = c(-74.45376, -73), ylim = c(-46, -45))
points(sites_tochange[2,], pch = 2, col = "darkred")
sites_coords_modify[
  sites_coords_modify$N_CODIGOCE == 110953,
  c(1:2)] <- c(-74.33687,-45.51854 + 0.08)
points(sites_coords_modify$lon[2],sites_coords_modify$lat[2])

plot(aysen_polygon, xlim = c(-72.45376, -73.5), ylim = c(-44, -43))
points(sites_tochange[3,], pch = 2, col = "darkred")
sites_coords_modify[
  sites_coords_modify$N_CODIGOCE == 110402,
  c(1:2)] <- c(-73.00848 - 0.01,-43.74188 - 0.01)
points(sites_coords_modify$lon[3],sites_coords_modify$lat[3])

plot(aysen_polygon, xlim = c(-73, -72), ylim = c(-46, -45))
points(sites_tochange[4,], pch = 2, col = "darkred")
sites_coords_modify[
  sites_coords_modify$N_CODIGOCE == 110125,
  c(1:2)] <- c(-72.85284-0.05,-45.42931 - 0.01)
points(sites_coords_modify$lon[4],sites_coords_modify$lat[4])

plot(aysen_polygon, xlim = c(-74, -73), ylim = c(-45, -44))
points(sites_tochange[5,], pch = 2, col = "darkred")
sites_coords_modify[
  sites_coords_modify$N_CODIGOCE == 110594,
  c(1:2)] <- c(-73.66930 + 0.05,-44.15644)
points(sites_coords_modify$lon[5],sites_coords_modify$lat[5])

plot(aysen_polygon, xlim = c(-74, -73), ylim = c(-45, -44))
points(sites_tochange[6,], pch = 2, col = "darkred")
sites_coords_modify[
  sites_coords_modify$N_CODIGOCE == 110705,
  c(1:2)] <- c(-73.93083,-44.20206 + 0.01)
points(sites_coords_modify$lon[6],sites_coords_modify$lat[6])

plot(aysen_polygon, xlim = c(-73, -72), ylim = c(-46, -45))
points(sites_tochange[7,], pch = 2, col = "darkred")
sites_coords_modify[
  sites_coords_modify$N_CODIGOCE == 110453,
  c(1:2)] <- c(-73.21510,-45.29663 + 0.01)
points(sites_coords_modify$lon[7],sites_coords_modify$lat[7])
}
validation <- extract(aysen_polygon, SpatialPoints(coords = sites_coords_modify[,c(1:2)]))
# ALL not in polygons
sites_coords_clean <- rbind(sites_coords_clear, sites_coords_modify)

pts_new <- sites_coords_clean
coordinates(pts_new) <- ~lon+lat
pts_new <- SpatialPointsDataFrame(pts_new, data = data.frame(id = pts_new@data$N_CODIGOCE))
crs(pts_new) <- crs(aysen_polygon)
over_points_sites <- extract(aysen_polygon, pts_new)
pts_new@data$in_polygon <- over_points_sites$codigo_region
pts_new@data$in_polygon <- ifelse(!is.na(pts_new@data$in_polygon),1, 0)
pts_new@data$colour <- ifelse(pts_new@data$in_polygon == 1,"darkred", "darkgreen")
plot(pts_new, col = pts_new@data$colour)

#write_rds(sites_coords_clean, "1_data/jittered_aysen_sites.rds")
aysen_polygon <- sf::as_Spatial(chilemapas::generar_regiones() |> 
                                  dplyr::filter(codigo_region == 11),
                                cast = F)
class(mapa_base)
plot(aysen_polygon)
aysen_polygon@bbox <- matrix(c(-74.7, -74, -46.5, -43.7), nrow = 2,
                             byrow = T)
plot(aysen_polygon)
c(-75.5, -72, -46.5, -43.5)
s0=expand.grid(x=seq(-75 ,to = -72, by = 0.02), y=seq(-46.5,to = -43.5, by = 0.02))

dim(s0)
s0_pts <- SpatialPointsDataFrame(coords = s0, data = data.frame(id=1:nrow(s0)))
crs(s0_pts) <- crs(aysen_polygon)

r_aysen <- raster::raster(ncol = 4500, nrow = 4500)
r_aysen@extent <- raster::extent(c(-74.7, -72.4, -47, -43.7))
r_aysen <- fasterize::fasterize(
  chilemapas::generar_regiones() |> 
    dplyr::filter(codigo_region == 11), r_aysen)


over_points_s0 <- extract(r_aysen, s0_pts , method = "simple",
                          fun = mean, na.rm = F)
beepr::beep(0)
s0_pts@data$in_polygon <- over_points_s0
s0_pts@data$in_polygon <- ifelse(!is.na(s0_pts@data$in_polygon),1, 0)
s0_pts@data$colour <- ifelse(s0_pts@data$in_polygon == 1,"darkred", "darkgreen")

s0_sea <- subset(s0_pts, in_polygon == 0)

write_rds(s0_sea, file= "1_data/Kriging_points_Extra.rds")
s0

mapa_base <- chilemapas::generar_regiones() |> 
  dplyr::filter(codigo_region == 11)
ggplot(mapa_base)+
    geom_raster(data = s0_sea@data,
                aes(x = s0_sea@coords[,1] , y = s0_sea@coords[,2]),
                alpha = 0.65, interpolate = T, fill = "darkred",
                show.legend =T)+
    geom_sf()+
  geom_point(aes(x = s0_sea@coords[1,1],y = s0_sea@coords[1,2]), data = s0_matrix)
    coord_sf(xlim = c(-74.62, -72.64745), ylim = c(-46.25, -44))

##### Maximum Summaries #####

sites <-
  readr::read_rds("1_data/joined_data.rds") |>
  janitor::clean_names()
sites <-
  sites |>
  dplyr::mutate(
    adultos_totales = prom_ho + prom_am,
    prop.juv = prom_juv/(prom_juv + adultos_totales)) |>
  dplyr::select(codigo, especie, nombre_region,  fecha_primer_muestreo, fecha_declaracion,
                biomasa_centro, peso_promedio, numero_total_peces,
                prom_juv, prom_ho, prom_am, adultos_totales,
                prop.juv, temperatura_promedio,
                salinidad_promedio) |>
  dplyr::mutate(
    ddiff =  if_else(
      is.na(lag(fecha_declaracion)),
      interval(fecha_declaracion-days(7), ymd(fecha_declaracion)),
      interval(ymd(lag(fecha_declaracion)), ymd(fecha_declaracion))
    ),
    ddiff = ddiff/weeks(1),
    .by = c(codigo, especie)) |>
  dplyr::mutate(new_cycle = if_else(ddiff > 7 | fecha_declaracion == ymd(min(fecha_declaracion)), 1, 0), .by = c(codigo, especie)) |>
  dplyr::mutate(cycle = cumsum(new_cycle), .by = c(codigo, especie)) |>
  dplyr::mutate(
    month = month(fecha_declaracion),
    year = year(fecha_declaracion),
    carga_total = (adultos_totales + prom_juv)*40,
    estacion = case_when(
      month == 1|month ==2 |month == 12~ "summer",
      month == 3|month ==4 |month == 5~ "autumn",
      month == 6|month ==7 |month == 8~ "winter",
      month == 9|month ==10 |month == 11~ "spring",
      
    )
  )

meanMax.Cycles <-
  sites |>
  dplyr::filter(
    especie %in% c("SALMON DEL ATLANTICO", "TRUCHA ARCOIRIS")) |>
  dplyr::filter(cycle > 1 & cycle < max(cycle),
                .by = c(codigo, especie)) |>
  dplyr::mutate(adultos_totales = adultos_totales * 40) |> 
  summarize(
    max.adultos_totales = max(adultos_totales, na.rm = T),
    max.carga_total = max(carga_total, na.rm = T),
    .by = c(codigo, nombre_region,
            especie,cycle)) |>
  filter(max.carga_total >1) |>
  left_join(
    sites[,c(1:2,6:15,16:22)],
    by = c("codigo","especie",
           "cycle",
           "max.carga_total" = "carga_total"))
max.summary <-
  meanMax.Cycles |>
  summarize(
    max.adultos_totales = mean(max.adultos_totales, na.rm = T),
    max.carga_total = mean(max.carga_total, na.rm = T),
    log_max.adultos_totales = mean(log1p(max.adultos_totales), na.rm = T),
    log_max.carga_total = mean(log1p(max.carga_total), na.rm = T),
    peso_promedio = mean(peso_promedio, na.rm = T),
    numero_total_peces = mean(numero_total_peces, na.rm = T),
    prom_juv = mean(prom_juv, na.rm = T),
    prop.juv = mean(prop.juv, na.rm = T),
    temperatura_promedio = mean(temperatura_promedio, na.rm = T),
    salinidad_promedio = mean(salinidad_promedio, na.rm = T),
    .by=c(codigo, nombre_region,especie, cycle)) |>
  summarize(
    max.carga_total = mean(max.carga_total, na.rm = T),
    max.adultos_totales = mean(max.adultos_totales, na.rm = T),
    log_max.carga_total = mean(log_max.carga_total, na.rm = T),
    log_max.adultos_totales = mean(log_max.adultos_totales, na.rm = T),
    peso_promedio = mean(peso_promedio, na.rm = T),
    numero_total_peces = mean(numero_total_peces, na.rm = T),
    prom_juv = mean(prom_juv, na.rm = T),
    prop.juv = mean(prop.juv, na.rm = T),
    temperatura_promedio = mean(temperatura_promedio, na.rm = T),
    salinidad_promedio = mean(salinidad_promedio, na.rm = T),
    .by=c(codigo,nombre_region, especie)) |>
  left_join(sites_coords,
            by = c("codigo" = "N_CODIGOCE"))


##### CALCULO MATRIZ ADYACENCIA #####

sites_coords <- readr::read_rds("1_data/jittered_aysen_sites.rds")
s0_pts <- readr::read_rds("1_data/Kriging_points.rds")

###### Rasters ######
1433*12/60/24
columnas = 4500
filas = 4500
r_aysen <- raster::raster(ncol = columnas, nrow = filas)
r_aysen@extent <- raster::extent(c(-75.5, -72, -46.7, -43.5))
r_aysen <- fasterize::fasterize(
  chilemapas::generar_regiones() |> 
    dplyr::filter(codigo_region == 11), r_aysen)
plot(r_aysen)
points(s0_sea)
r_aysen[is.na(r_aysen)] <- -999 # this turns all landmass to missing
r_aysen[r_aysen>-999] <- NA # assign unit cost to all grid cells in water
r_aysen[r_aysen==-999] <- 1 # make landmass impassable
gc()
rs.transition <- gdistance::transition(r_aysen, mean, directions = 16)
gc()

rs.transition@crs <- sp::CRS("+ellps=WGS84 +proj=longlat +datum=WGS84 +no_defs +units=m")
#readr::write_rds(rs.transition, "1_data/Aysen_Transision_Kriging.rds")
rs.transition <- readr::read_rds("1_data/Aysen_Transision_Kriging.rds")
rs.transition <- gdistance::geoCorrection(rs.transition, "r")
gc()
readr::write_rds(rs.transition, "1_data/Aysen_Transision_Kriging_GeoCorrected.rds")


###### Distance Calculations ###### 

rs.transition <- readr::read_rds("1_data/Aysen_Transision_Kriging_GeoCorrected.rds")
sites_coords <- readr::read_rds("1_data/jittered_aysen_sites.rds")
s0_pts <- readr::read_rds("1_data/Kriging_points_Extra.rds")
max.summary <- readr::read_rds("1_data/max_summary_AysenLagos.rds")
aysen.summary <- max.summary |> filter(nombre_region == "XI REGION")
aysen_coords <- sites_coords |> filter(N_CODIGOCE %in% unique(aysen.summary$codigo))

matrix_coords <- as.matrix(aysen_coords[,c(1:2)])
rownames(matrix_coords) <- aysen_coords$N_CODIGOCE
s0_matrix <- coordinates(s0_pts)
s0_matrix |> head()
lista_distancias <- list()
# rs.transition@extent <- raster::extent(c(min(s0_matrix[,1]), max(s0_matrix[,1]),
#                                          min(s0_matrix[,2]), max(s0_matrix[,2])))

full_distances <- gdistance::costDistance(rs.transition,
                                matrix_coords,
                                s0_matrix
)
beepr::beep(8)
dim(full_distances)

write_rds(full_distances, "kriging_distances_extra_points.rds")

sum(apply(full_distances,2,function(x){
  sum(is.infinite(x))
})>0)


# Extra Points fixing -----------------------------------------------------
s0_distances <- readr::read_rds("kriging_distances_extra_points.rds")
coordinates(s0)
origin <- st_as_sf(aysen_coords[,c(1:2)], coords = c("lon", "lat"),
                   crs = st_crs("WGS84"))

colnames(s0_matrix) <- c("lon", "lat")
s0_matrix
dest_s0 <- st_as_sf(data.frame(s0_matrix), coords = c("lon", "lat"),
                    crs = st_crs("WGS84"))
distances_krig_euc <- st_distance(origin, dest_s0, by_element = F )
distances_krig_euc <-units::set_units(distances_krig_euc, "km")
distances_krig_euc <- matrix(distances_krig_euc, ncol = dim(distances_krig_euc)[2])
rownames(distances_krig_euc) <- aysen_coords$codigo
s0_distances_aux <- s0_distances


for(j in ncol(s0_distances_aux):1){
  aux_sum <- sum(is.infinite(s0_distances_aux[,j]))
  if(aux_sum > 0){
    dist_diff_euc <- distances_krig_euc[,j] - distances_krig_euc[,j+1]
    s0_distances_aux[,j] <- s0_distances_aux[,j+1] + dist_diff_euc
  }
  print(j)
}
dist_diff_euc <- distances_krig_euc[,10769] - distances_krig_euc[,10769-1]

s0_distances_aux[,10769] <- s0_distances_aux[,10769-1] + dist_diff_euc
inf_cols <- apply(s0_distances_aux, 2, function(x) sum(is.infinite(x)))
s0_matrix[10769,]
s0_matrix_aux <- s0_matrix
s0_matrix_aux[1:length(inf_cols),3] <- inf_cols

write_rds(s0_distances_aux, "1_data/kriging_distance_S0_Fixed.rds")


s0_matrix_aux <- s0_matrix_aux |> data.frame()

s0_matrix_aux$inf_cols <- inf_cols
mapa_base <- chilemapas::generar_regiones() |> dplyr::filter(codigo_region == 11)

# ggplot(mapa_base)+
#   geom_sf()+
#   geom_point(data = s0_matrix_aux, 
#              aes(x =  lon,
#                  y = lat,
#                  colour = inf_cols ),
#              show.legend =F , alpha = 0.4)+
#   coord_sf(xlim = c(-74.62, -72.64745), ylim = c(-46.25, -44))
  
