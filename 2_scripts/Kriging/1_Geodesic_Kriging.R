library(sp)
library(sf)
library(stars)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(patchwork)

aysen.base <- chilemapas::generar_regiones() |> dplyr::filter(codigo_region == 11)

coordenadas <- max.summary.aysen@coords
info <- as.data.frame(max.summary.aysen@data)
plot.data <- cbind(info, coordenadas)


load(file = "shape_regions.RData")
shape_aysen <- shape_aysen |> janitor::clean_names() |> rename(codigo = n_codigoce)
shape_lagos <- shape_lagos |> janitor::clean_names() |> rename(codigo = n_codigoce)
full_shape <- rbind(shape_lagos, shape_aysen)

source("2_scripts/1_ImprovedSemivariogram.R")
source("2_scripts/DistObjects.R")

lagos.distance <- readr::read_rds("1_data/Lagos_DistanceMatrix.rds")
aysen.distance <- readr::read_rds("1_data/Aysen_DistanceMatrix.rds")

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

x_Max.Site <-
  sites |>
  dplyr::filter(
    especie %in% c("SALMON DEL ATLANTICO", "TRUCHA ARCOIRIS")) |>
  dplyr::filter(cycle > 1 & cycle < max(cycle),
                .by = c(codigo, especie)) |>
  filter(fecha_declaracion <= "2019-10-15" & fecha_declaracion >= "2018-10-15") |> 
  dplyr::mutate(adultos_totales = adultos_totales*40) |> 
  dplyr::summarize(
    adultos_totales = max(adultos_totales, na.rm = T),
    .by = c(codigo,nombre_region, especie)) |>
  dplyr::filter(adultos_totales >1) |>
  dplyr::left_join(
    sites[,c(1:2,5:14,16)],
    by =c("codigo","especie",  "adultos_totales")) |>
  left_join(full_shape[,c(3,22)], by = c("codigo" = "N_CODIGOCE")) |> 
  filter(nombre_region == "XI REGION")

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
  left_join(full_shape[,c(3,22)],
            by = c("codigo" = "N_CODIGOCE"))

# max.summary.lagos <-
#   max.summary |>
#   filter(codigo != 101558 & nombre_region == "X REGION") |>
#   arrange(codigo)
# max(max.summary.lagos$max.adultos_totales)
# coordinates(max.summary.lagos) <- ~coordinates
# proj4string(max.summary.lagos) <- "+proj=longlat +ellps=WGS84"
# ggplot(max.summary.lagos@data, aes(x = coordinates[,1],
#                                    y = coordinates[,2],
#                                    size = residuos.log))+
#   geom_point()
# exp(2)
# nrow(max.summary.lagos)
# nrow(max.summary.aysen)
# lagos.subset <-
#   lagos.distance[rownames(lagos.distance)%in% max.summary.lagos@data$codigo ,
#                  colnames(lagos.distance)%in% max.summary.lagos@data$codigo]
# 
# dist.list.Lagos <- DistObjects(max.summary.lagos)
# dist.list.Lagos$LeastCost <- lagos.subset/1000
# 
# dist.list.Lagos <- lapply(dist.list.Lagos, function(x){
#   rownames(x) <- max.summary.lagos@data$codigo
#   colnames(x) <- max.summary.lagos@data$codigo
#   x
# })
# 

# Aysen -------------------------------------------------------------------

max.summary.aysen <-
  max.summary |>
  filter(nombre_region == "XI REGION") |>
  arrange(codigo)

coordinates(max.summary.aysen) <- ~coordinates
proj4string(max.summary.aysen) <- "+proj=longlat +ellps=WGS84"

aysen.subset <-
  aysen.distance[rownames(aysen.distance)%in% max.summary.aysen@data$codigo ,
                 colnames(aysen.distance)%in% max.summary.aysen@data$codigo]

dist.list.Aysen <- DistObjects(max.summary.aysen)
dist.list.Aysen$LeastCost <- aysen.subset/1000

dist.list.Aysen <- lapply(dist.list.Aysen, function(x){
  rownames(x) <- max.summary.aysen@data$codigo
  colnames(x) <- max.summary.aysen@data$codigo
  x
})

m.aysen.raw <- lm(I((max.adultos_totales)) ~ temperatura_promedio+
                    I((salinidad_promedio)),max.summary.aysen@data)

m.aysen.log <- lm(I((log_max.adultos_totales)) ~ peso_promedio + temperatura_promedio+
                    I((salinidad_promedio)),max.summary.aysen@data)

max.summary.aysen@data$residuos <- resid(m.aysen.raw)
max.summary.aysen@data$residuos.log <- resid(m.aysen.log)

ggplot(max.summary.aysen@data, aes(x = coordinates[,1],
                                   y = coordinates[,2],
                                   size = residuos))+
  geom_point()


# Preparacion data --------------------------------------------------------
library(GeoModels)
coordenadas <- max.summary.aysen@coords
info <- as.data.frame(max.summary.aysen@data)
rownames(coordenadas) <- info$codigo 
mean <- mean(info$max.adultos_totales)
# mean1 <- median(info$max.adultos_totales)

info <- info |> dplyr::select(max.adultos_totales, temperatura_promedio, salinidad_promedio, codigo)
info2 <- info |> tidyr::pivot_longer(-codigo) |> tidyr::pivot_wider(names_from = codigo, values_from = value)
info.matrix <- as.matrix(info2)



# GeoModels Modeling ------------------------------------------------------
# info <- load("1_data/Region_Maximums.RData")

variograma <- 
  cs.variog("log_max.adultos_totales", id_col = "codigo",
          data = max.summary.aysen, 
          max.dist = 50,
          breaks = 25,
          dist_object = dist.list.Aysen$Geodesic,
          plot = F)
plot(variograma, robust = F)
plot(variograma, robust = T)

max(dist.list.Aysen$Geodesic)
max(dist.list.Aysen$Geodesic)
sum(as.dist(dist.list.Aysen$LeastCost) < 40)

dim(as.dist(dist.list.Aysen$LeastCost))
length(dist.list.Aysen$LeastCost)

# dist.list.Aysen$Geodesic |> view()
# dist.list.Aysen$LeastCost |> view()

min(variograma$plot.data$semivariog)/mean(variograma$plot.data$semivariog)
sill <- 40
nugget <- 0.4
scale <- 268/3
smooth=0.5
I=Inf
lower<-list(mean=0,
            scale=0
            ,sill=0
            )
upper<-list(mean=I,
            scale=I
            ,sill=I 
            )
 fixed<-list(nugget=nugget,smooth=smooth)
 start<-list(mean=mean(info$log_max.adultos_totales),
             scale=scale
             ,sill=sill
             )

fit <- 
  GeoModels::GeoFit(
  data=info$log_max.adultos_totales,
  coordx=coordenadas,
  corrmodel="Matern",
  optimizer="nlminb",
  upper=upper,
  lower=lower,
  likelihood="Marginal",
  type="Pairwise",
  model="Gaussian",
  # neighb = 3,
  maxdist = 268/3,
  start=start,
  fixed=fixed,
  distance = "Geod")

x <- seq(-74.7,to = -72.4, 0.01)
y <- seq(-46.341917,to = -43.7, 0.01)
d1 <- expand.grid(x = x, y = y)

Ord.kriging <- GeoModels::GeoKrig(fit, data=info$log_max.adultos_totales, 
                           coordx=coordenadas,corrmodel="Matern",
                           distance = "Geod",loc = d1)

Ord.prediction_data <- cbind(d1, exp(Ord.kriging$pred))
Ord.prediction_data <- as.data.frame(prediction_data)
colnames(Ord.prediction_data) <- c("x", "y", "Kriging")
Ord.kriging$param
# quilt.plot(coordenadas,info$max.adultos_totales,col=colour)
# simple kriging map prediction
# imagePlot(x, y, z=matrix(kriging$pred,nrow=length(x)))

Kriging.Mean.Plot <- 
  ggplot(aysen.base)+
  geom_tile(data = Ord.prediction_data,
            aes(x = x, y = y, fill = Kriging),
            alpha = 0.65,
            show.legend =F)+
  geom_sf()+
  geom_point(data = plot.data, 
             aes(x =  coordinates.coords.x1,
                 y = coordinates.coords.x2,
                 colour = max.adultos_totales),
             show.legend =F )+
  geom_point(shape = 1,colour = "black",
             data = plot.data, 
             aes(x =  coordinates.coords.x1,
                 y = coordinates.coords.x2))+
  scale_fill_gradientn(colors = rev(rainbow(5)))+
  coord_sf(xlim = c(-74.62, -72.64745), ylim = c(-46.25, -44))+
  labs(x = "longitud",
       y = "latitud")+
  scale_colour_gradientn(colors = rev(rainbow(5)))+
  theme_bw()+
  labs(subtitle = "Aysén - Media constante",
       colour = "Maximo de Adultos Totales")+
  theme(plot.subtitle = element_text(hjust = 0.5))
Kriging.Mean.Plot
# With Covariates ---------------------------------------------------------
mean <- mean(info$max.adultos_totales)
# mean1 <- median(info$max.adultos_totales)

sill <- 40
nugget <- 0.4
scale <- 268/3
smooth=0.5
I=Inf
mean <- mean(info$log_max.adultos_totales)
mean1 <- mean(info$temperatura_promedio)
mean2 <- mean(log(info$salinidad_promedio))
lower2<-list(mean=0,
             mean1 = 0,
             mean2 = 0,
            scale=0
            ,sill=0
)
upper2<-list(mean=I,
             mean1=I,
             mean2=I,
            scale=I
            ,sill=I 
)
fixed2<-list(nugget=nugget,
            smooth=smooth)
start2<-list(mean=mean,
             mean1=log(mean1),
             mean2=log(mean2),
            scale=scale
            ,sill=sill
)


# distancia <- geosphere::distm(coordenadas, fun = geosphere::distGeo)|> as.matrix()
# max(distancia)/3
covariables <- model.matrix(log_max.adultos_totales ~ 
                               I(log(temperatura_promedio))+
                              I(log(salinidad_promedio)),
                            info)

fit2 <- 
  GeoModels::GeoFit(
    data=info$log_max.adultos_totales,
    coordx=coordenadas,
    corrmodel="Matern",
    optimizer="nlminb",
    upper=upper2,
    lower=lower2,
    likelihood="Marginal",
    type="Pairwise",
    model="Gaussian",
    # neighb = 3,
    maxdist = 268/3,
    start=start2,
    fixed=fixed2,
    distance = "Geod",
    X = covariables)

Univ.pred.cov <- matrix(nrow = nrow(d1), ncol = 3, NA)
Univ.pred.cov[,1] <- rep(1, nrow(d1))
Univ.pred.cov[,2] <- rep(log(8), nrow(d1))
Univ.pred.cov[,3] <- rep(log(10), nrow(d1))

Univ.kriging <- GeoModels::GeoKrig(fit2, 
                              data=info$log_max.adultos_totales, 
                              coordx=coordenadas,
                              corrmodel="Matern",
                              distance = "Geod",
                              loc = d1,
                              X = covariables,
                              Xloc = Univ.pred.cov
                              )
Univ.prediction_data <- cbind(d1, Univ.kriging$pred)

Univ.prediction_data <- as.data.frame(Univ.prediction_data)
colnames(Univ.prediction_data) <- c("x", "y", "Kriging")
# quilt.plot(coordenadas,info$max.adultos_totales,col=colour)
# simple kriging map prediction
# imagePlot(x, y, z=matrix(kriging$pred,nrow=length(x)))
round(Univ.kriging$param,3)
Kriging.Cov.Plot <- 
  ggplot(aysen.base)+
  geom_tile(data = Univ.prediction_data,
            aes(x = x, y = y, fill = Kriging),
            alpha = 0.5)+
  geom_sf()+
  geom_point(data = plot.data, 
             aes(x =  coordinates.coords.x1,
                 y = coordinates.coords.x2,
                 colour = log_max.adultos_totales),
             fill = "black")+
  geom_point(shape = 1,colour = "black",
             data = plot.data, 
             aes(x =  coordinates.coords.x1,
                 y = coordinates.coords.x2))+
  scale_fill_gradientn(colors = rev(rainbow(5)))+
  coord_sf(xlim = c(-74.62, -72.64745), ylim = c(-46.25, -44))+
  labs(x = "longitud",
       y = "latitud")+
  scale_colour_gradientn(colors = rev(rainbow(5)))+
  theme_bw()+
  labs(subtitle = "Aysén - Kriging con Covariables",
       colour = "Maximo de log(Adultos Totales)")+
  theme(plot.subtitle = element_text(hjust = 0.5))

Kriging.Cov.Plot

# Solo Gráficos -----------------------------------------------------------
Kriging.Mean.Plot+Kriging.Cov.Plot
info |> select(max.adultos_totales, salinidad_promedio, temperatura_promedio) |> plot()
plot(log(info$temperatura_promedio), log(info$max.adultos_totales))
plot(1/(log(info$salinidad_promedio)), log(info$max.adultos_totales))
plot((log(info$salinidad_promedio)), log(info$temperatura_promedio))

fit2
