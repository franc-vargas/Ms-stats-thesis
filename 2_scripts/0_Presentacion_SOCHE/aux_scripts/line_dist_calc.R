library(gdistance)
library(raster)
rs.transition <- readr::read_rds("1_data/Aysen_Transision_Kriging_GeoCorrected.rds")
sites_coords <- readr::read_rds("1_data/jittered_aysen_sites.rds")
aysen_coords <- sites_coords |> filter(N_CODIGOCE %in% unique(aysen.summary$codigo))
matrix_coords <- as.matrix(aysen_coords[,c(1:2)])

AtoB <- list() # create a list

AtoB[[1]] <- gdistance::shortestPath(rs.transition, matrix_coords[16,],
                          matrix_coords[19,],
                          output = "SpatialLines")

mapa_base <- chilemapas::generar_regiones() |> dplyr::filter(codigo_region == 11)
plot(mapa_base,xlim = extent(rs.transition)[1:2], ylim = extent(rs.transition)[3:4])
lines(AtoB[[1]])
lineas <- AtoB[[1]]
crs(AtoB[[1]]) <- crs(mapa_base)
extent(AtoB[[1]] )
asd |> as.vector()
ggplot(mapa_base) + 
  geom_sf() +
  geom_point(
    fill = "black",
    data = aysen.summary[c(16,19),],
    alpha = 0.7,
    aes(x =  lon,
        y = lat)) + 
  geom_path(data =data.frame(lineas_coords),
            aes(x=x, y=y,
                colour = "Non-Euclidian"),
            alpha = 0.7) + 
  geom_segment(data = data.frame(matrix_coords), 
            aes(x = matrix_coords[16,1], 
                y = matrix_coords[16,2],
                xend = matrix_coords[19,1],
                yend = matrix_coords[19,2],
                colour = "Euclidian"),
            alpha = 0.7) + 
  theme_bw() + 
  labs(colour = "", x = "longitude", y = "latitude") +
  scale_colour_manual(values = c("Non-Euclidian" = "darkred", 
                                 "Euclidian" = "darkgreen")) + 
  coord_sf(xlim = c(-74.62, -72.64745), ylim = c(-45.73, -44.4))
  

lineas@lines[["coords"]]
readr::write_rds(lineas_coords, "2_scripts/0_Presentación/lineas_noneuc.rds")
lineas_coords <- (coordinates(lineas))[[1]][[1]]
lineas_coords
