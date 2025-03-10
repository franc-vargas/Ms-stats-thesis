##### Carga Datos #####
library(ggplot2)
library(sf)
library(tidyverse)
library(sp)
mapa_base <- chilemapas::generar_regiones() |> dplyr::filter(codigo_region == 11)

predicted_points <- readr::read_rds("1_data/Kriging_points.rds")
max.summary <- readr::read_rds("1_data/max_summary_AysenLagos.rds")
aysen.summary <- max.summary |> filter(nombre_region == "XI REGION")
aysen.summary$re_scaled_total_adults <- aysen.summary$max.adultos_totales / 40
aysen_coords <- 
  aysen.summary |>
  select(lon, lat, codigo)

vg_gstat <- gstat::variogram(max.adultos_totales ~ 1, aysen.sf)
plot(vg_gstat)

aysen.sf <- aysen.summary
coordinates(aysen.sf) <- ~lon+lat

##### Distance Matrices #####
## Non Euclidian
aysen.distance <- readr::read_rds("1_data/Aysen_DistanceMatrix.rds")
aysen.distance <- aysen.distance/1000
ids <- unique(aysen.summary$codigo)
aysen.dist.subset <- aysen.distance[rownames(aysen.distance)%in%ids, 
                                    colnames(aysen.distance) %in% ids]
aysen.dist.subset <- aysen.dist.subset[order(match(rownames(aysen.dist.subset),
                                                   ids)), 
                                       order(match(colnames(aysen.dist.subset), ids))]
S0 <- readr::read_rds("kriging_distances.rds")
S0_Non_Euc <- S0[order(match(rownames(S0),
                             ids)),]/1000

## Euclidian
origin <- st_as_sf(aysen_coords[,c(1:2)], coords = c("lon", "lat"),
                   crs = st_crs("WGS84"))
dest <- st_as_sf(aysen_coords[,c(1:2)], coords = c("lon", "lat"),
                 crs = st_crs("WGS84"))
distances <- st_distance(origin, dest, by_element = F )


distances <-units::set_units(distances, "km")
distances <- as.dist(distances) |> as.matrix()
rownames(distances) <- aysen_coords$codigo
colnames(distances) <- aysen_coords$codigo

predicted_points_euc <- coordinates(predicted_points)
colnames(predicted_points_euc) <- c("lon", "lat")
predicted_points_euc <- data.frame(predicted_points_euc)

dest_s0 <- st_as_sf(predicted_points_euc, coords = c("lon", "lat"),
                    crs = st_crs("WGS84"))
distances_krig_euc <- st_distance(origin, dest_s0, by_element = F )
distances_krig_euc <-units::set_units(distances_krig_euc, "km")
distances_krig_euc <- matrix(distances_krig_euc, ncol = dim(distances_krig_euc)[2])
rownames(distances_krig_euc) <- aysen_coords$codigo

##### Empirical Semivariogram #####

vg_custom <- cs.variog("max.adultos_totales", "codigo", aysen.distance, 
                      aysen.summary, breaks = 15, 
                      max.dist = round(max(aysen.distance)/4,0), plot = F)

vg_custom_euc <- cs.variog("max.adultos_totales", "codigo", distances, 
                           aysen.summary, breaks = 15, 
                           max.dist = round(max(distances)/4,0), plot = F)

vg_custom_scaled <- cs.variog("re_scaled_total_adults", "codigo", aysen.distance, 
                       aysen.summary, breaks = 15, 
                       max.dist = round(max(aysen.distance)/4,0), plot = F)

vg_custom_euc_scaled <- cs.variog("re_scaled_total_adults", "codigo", distances, 
                           aysen.summary, breaks = 15, 
                           max.dist = round(max(distances)/4,0), plot = F)


par(mfrow = c(1,2))

plot(vg_custom_euc$plot.data$lag, 
     vg_custom_euc$plot.data$semivariog 
     ,ylim = c(1e+05, 9e+05),
     main = "Euclidian Distance")

plot(vg_custom$plot.data$lag, 
     vg_custom$plot.data$semivariog 
     ,ylim = c(1e+05, 9e+05),
     main = "Least-Cost Distance")

plot(vg_custom_euc_scaled$plot.data$lag, 
     vg_custom_euc_scaled$plot.data$semivariog 
     ,
     main = "Euclidian Distance")

plot(vg_custom_scaled$plot.data$lag, 
     vg_custom_scaled$plot.data$semivariog 
     ,
     main = "Least-Cost Distance")

par(mfrow = c(1,1))


mu <- mean(aysen.summary$max.adultos_totales/ 40)
nugget <- vg_custom$plot.data$semivariog[1]
sigma2 <- vg_custom$plot.data$semivariog[9]
phi <- 21

h_new <-  seq(1,round(max(aysen.distance)/4,0), by=0.1)
plot(vg_custom)
lines(h_new, 
      variogram(h_new, nugget=nugget, sigma2=sigma2, phi=phi), type="l",
      col = "darkred")

plot(vg_custom_euc)
lines(h_new, 
      variogram(h_new, nugget=nugget, sigma2=sigma2, phi=20), type="l",
      col = "darkred")


##### Matrix Generation #####

# Based on Empirical Semivariogram
param.OK <- data.frame(mu = mu,sigma2 = sigma2, phi = phi, nugget = nugget)
# Distance Covariance
Sigma_Non_Euc <- geoR::cov.spatial(aysen.dist.subset, cov.model="matern", 
                           cov.pars=c(
                             sigma2, 
                             phi)) + nugget*diag(nrow(aysen.summary))

Sigma_Euc <- 
  geoR::cov.spatial(distances, cov.model="matern", 
                    cov.pars=c(
                      sigma2, 
                      phi)) + nugget*diag(nrow(aysen.summary))

c0 <- geoR::cov.spatial(S0_Non_Euc, 
                        cov.model="matern", 
                        cov.pars=c(sigma2, 
                                   nugget)
)
c0[h0==0] <- sigma2 + nugget

##### Kriging Ordinario #####

OK_Non_Euc <- krig.manual(data=aysen.summary$max.adultos_totales,
                          Sigma=Sigma_Non_Euc, 
                          param = param.OK,
                          cov.model="matern", 
                          h0=S0_Non_Euc,
                          type = "ord")
OK_Euclid <- krig.manual(data=aysen.summary$max.adultos_totales,
                         Sigma=Sigma_Euc, 
                         param = param.OK,
                         cov.model="matern", 
                         h0=distances_krig_euc,
                         type = "ord")

hist(OK_Non_Euc$pred)
dim(OK_Non_Euc$var.pred)

predicted_points@data$predictions_OK <- OK_Non_Euc$pred
predicted_points@data$var_pred_OK <- OK_Non_Euc$var.pred
predicted_points@data$predictions_OK_Euc <- OK_Euclid$pred
predicted_points@data$var_pred_OK_Euc <- OK_Euclid$var.pred

ggplot()+
  geom_histogram(aes(x = OK_Non_Euc$var.pred))
Sig_NE_Eigen <- eigen(Sigma_Non_Euc)
Sig_NE_Eigen$values

Kriging.Mean.Plot <- 
  ggplot(mapa_base)+
  geom_raster(data = predicted_points@data,
            aes(x = predicted_points@coords[,1] , y = predicted_points@coords[,2], 
                fill = predictions_OK/40),
            alpha = 0.65, interpolate = T,
            show.legend =T)+
  geom_sf()+
  geom_point(data = aysen.summary, 
             aes(x =  lon,
                 y = lat,
                 colour = max.adultos_totales/40),
             show.legend =F )+
  geom_point(shape = 1,
             colour = "black",
             data = aysen.summary,
             alpha = 0.7,
             aes(x =  lon,
                 y = lat))+
  coord_sf(xlim = c(-74.62, -72.64745), ylim = c(-46.25, -44))+
  scale_fill_gradientn(colors = rev(rainbow(5)),
                       limits = c(0, max(aysen.summary$max.adultos_totales)/40))+
  scale_colour_gradientn(colors = rev(rainbow(5)),
                         limits = c(0, max(aysen.summary$max.adultos_totales)/40))+
  theme_bw()+
  theme(legend.position = "bottom")+
  labs(x = "longitude",
       y = "latitude", 
       fill = "Carga Adultos Krig",
       colour = "Carga Adultos Centro",
       title = "Non-Euclidian Distance")

Kriging.Mean.Plot
Kriging_var_plot_noneuc <- 
  ggplot(mapa_base)+
  geom_raster(data = predicted_points@data,
              aes(x = predicted_points@coords[,1], y = predicted_points@coords[,2], 
                  fill = var_pred_OK),
              alpha = 0.65, interpolate = T,
              show.legend =T)+
  geom_sf()+
  geom_point(data = aysen.summary, 
             aes(x =  lon,
                 y = lat,
                 colour = max.adultos_totales),
             show.legend =F, alpha = 0.7 )+
  geom_point(shape = 1,
             colour = "black",
             data = aysen.summary,
             alpha = 0.7,
             aes(x =  lon,
                 y = lat))+
  coord_sf(xlim = c(-74.62, -72.64745), ylim = c(-46.25, -44))+
  scale_fill_gradientn(colors = rev(rainbow(5)),
                       )+
  scale_colour_gradientn(colors = rev(rainbow(5)),
                       limits = c(0, max(aysen.summary$max.adultos_totales)))+
  theme_bw()+
  theme(legend.position = "bottom")+
  labs(x = "longitude",
       y = "latitude", 
       fill = "Krig. Varianza",
       colour = "Carga Adultos Centro",
       title = "Non-Euclidian Distance")



# Ord. Kriging Euclidian -------------------------------------------------



hist(OK_Euclid$var.pred)
min(OK_Euclid$var.pred)
Kriging.Mean.Plot_euc <- 
  ggplot(mapa_base)+
  geom_raster(data = predicted_points@data,
              aes(x = predicted_points@coords[,1] , y = predicted_points@coords[,2], 
                  fill = predictions_OK_Euc/40),
              alpha = 0.65, interpolate = T,
              show.legend =T)+
  geom_sf()+
  geom_point(data = aysen.summary, 
             aes(x =  lon,
                 y = lat,
                 colour = max.adultos_totales/40),
             show.legend =F )+
  geom_point(shape = 1,
             colour = "black",
             data = aysen.summary,
             alpha = 0.7,
             aes(x =  lon,
                 y = lat))+
  coord_sf(xlim = c(-74.62, -72.64745), ylim = c(-46.25, -44))+
  scale_fill_gradientn(colors = rev(rainbow(5)),
                       limits = c(0, max(aysen.summary$max.adultos_totales)/40))+
  scale_colour_gradientn(colors = rev(rainbow(5)),
                         limits = c(0, max(aysen.summary$max.adultos_totales)/40))+
  theme_bw()+
  theme(legend.position = "bottom")+
  labs(x = "longitude",
       y = "latitude", 
       fill = "Carga Adultos Krig",
       colour = "Carga Adultos Centro",
       title = "Euclidian Distance")
library(patchwork)

Kriging_var_plot_euc <- 
  ggplot(mapa_base)+
  geom_raster(data = predicted_points@data,
              aes(x = predicted_points@coords[,1] , y = predicted_points@coords[,2], 
                  fill = var_pred_OK_Euc),
              alpha = 0.65, interpolate = T,
              show.legend =T)+
  geom_sf()+
  geom_point(data = aysen.summary, 
             aes(x =  lon,
                 y = lat,
                 colour = max.adultos_totales/40),
             show.legend =F, alpha = 0.7 )+
  geom_point(shape = 1,
             colour = "black",
             data = aysen.summary,
             alpha = 0.7,
             aes(x =  lon,
                 y = lat))+
  coord_sf(xlim = c(-74.62, -72.64745), ylim = c(-46.25, -44))+
  scale_fill_gradientn(colors = rev(rainbow(5)),
  )+
  scale_colour_gradientn(colors = rev(rainbow(5)),
                         limits = c(0, max(aysen.summary$max.adultos_totales)/40))+
  theme_bw()+
  theme(legend.position = "bottom")+
  labs(x = "longitude",
       y = "latitude", 
       fill = "Krig. Varianza",
       colour = "Carga Adultos Centro",
       title = "Euclidian Distance")


Kriging.Mean.Plot + Kriging.Mean.Plot_euc
Kriging_var_plot_noneuc + Kriging_var_plot_euc


ggplot(mapa_base)+
  geom_raster(data = predicted_points@data,
              aes(x = predicted_points@coords[,1] , y = predicted_points@coords[,2], 
                  fill = (predictions_OK - predictions_OK_Euc)),
              alpha = 0.65, interpolate = T,
              show.legend =T)+
  geom_sf()+
  
  geom_point(shape = 1,
             colour = "black",
             data = aysen.summary,
             alpha = 0.7,
             aes(x =  lon,
                 y = lat))+
  coord_sf(xlim = c(-74.62, -72.64745), ylim = c(-46.25, -44))+
  scale_fill_gradientn(colors = rev(rainbow(5)),
  )+
  scale_colour_gradientn(colors = rev(rainbow(5)))+
  theme_bw()+
  theme(legend.position = "bottom")+
  labs(x = "longitude",
       y = "latitude", 
       fill = "Diff. Krig. Prediction")

ggplot(mapa_base)+
  geom_raster(data = predicted_points@data,
              aes(x = predicted_points@coords[,1] , y = predicted_points@coords[,2], 
                  fill = (var_pred_OK/var_pred_OK_Euc)),
              alpha = 0.65, interpolate = T,
              show.legend =T)+
  geom_sf()+
  geom_point(shape = 1,
             colour = "black",
             data = aysen.summary,
             alpha = 0.7,
             aes(x =  lon,
                 y = lat))+
  coord_sf(xlim = c(-74.62, -72.64745), ylim = c(-46.25, -44))+
  scale_fill_gradientn(colors = rev(rainbow(5)),
                        limits = c(0,2.2))+
  theme_bw()+
  theme(legend.position = "bottom")+
  labs(x = "longitude",
       y = "latitude", 
       fill = "Diff. Krig. Varianza")

hist(OK_Non_Euc$var.pred)
hist(OK_Euclid$var.pred)
sum(ks1$var.pred < 0)
sum(OK_Euclid$var.pred < 0)
