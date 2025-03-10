# Library & Data ----------------------------------------------------------
library(sp)
library(sf)
library(stars)
library(ggplot2)
library(dplyr)
library(tidyverse)

load(file = "shape_regions.RData")
full_shape <- rbind(shape_lagos, shape_aysen)
source("2_scripts/custom_semivariog.R")
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
  dplyr::summarize(
    adultos_totales = max(adultos_totales, na.rm = T),
    .by = c(codigo,nombre_region, especie)) |> 
  dplyr::filter(adultos_totales >1) |> 
  dplyr::left_join(
    sites[,c(1:2,5:14,16)],
    by =c("codigo","especie",  "adultos_totales")) |> 
  left_join(full_shape[,c(3,22)], by = c("codigo" = "N_CODIGOCE"))

meanMax.Cycles <- 
  sites |> 
  dplyr::filter(
      especie %in% c("SALMON DEL ATLANTICO", "TRUCHA ARCOIRIS")) |> 
  dplyr::filter(cycle > 1 & cycle < max(cycle), 
                .by = c(codigo, especie)) |> 
  summarize(
    carga_total = max(carga_total, na.rm = T),
    .by = c(codigo, nombre_region,
            especie, year,cycle)) |> 
  filter(carga_total >1) |> 
  left_join(
    sites[,c(1:2,6:15,16:22)],
    by = c("codigo","especie", "year",
            "cycle",
           "carga_total" = "carga_total"))

max.summary <- 
  meanMax.Cycles |> 
  summarize(
    carga_total = mean(carga_total, na.rm = T),
    peso_promedio = mean(peso_promedio),
    numero_total_peces = mean(numero_total_peces),
    prom_juv = mean(prom_juv),
    prop.juv = mean(prop.juv),
    temperatura_promedio = mean(temperatura_promedio), 
    salinidad_promedio = mean(salinidad_promedio),
    .by=c(codigo, nombre_region,especie, year, cycle)) |> 
  summarize(carga_total = mean(carga_total, na.rm = T),
            peso_promedio = mean(peso_promedio),
            numero_total_peces = mean(numero_total_peces),
            prom_juv = mean(prom_juv),
            prop.juv = mean(prop.juv),
            temperatura_promedio = mean(temperatura_promedio), 
            salinidad_promedio = mean(salinidad_promedio),
            .by=c(codigo,nombre_region, especie)) |> 
  left_join(full_shape[,c(3,22)], 
            by = c("codigo" = "N_CODIGOCE"))


# Custom variogram ----------------------------------------------------------
# Robert H. Shumway
# Kauffman Clustering
max.summary.lagos <-
  max.summary |> 
  filter(codigo != 101558 & nombre_region == "X REGION") |> 
  arrange(codigo) 

coordinates(max.summary.lagos) <- ~coordinates
proj4string(max.summary.lagos) <- "+proj=longlat +ellps=WGS84"

lagos.subset <- 
  lagos.distance[rownames(lagos.distance)%in% max.summary.lagos@data$codigo ,
                               colnames(lagos.distance)%in% max.summary.lagos@data$codigo]

dist.list.Lagos <- DistObjects(max.summary.lagos)
dist.list.Lagos$LeastCost <- lagos.subset/1000

dist.list.Lagos <- lapply(dist.list.Lagos, function(x){
  rownames(x) <- max.summary.lagos@data$codigo
  colnames(x) <- max.summary.lagos@data$codigo
  x
})

lay.matrix <- matrix(c(1,1,2,2,3,3,4,4, 5,5, 6, 6),3, 4, byrow = TRUE)
par(mai = c(0.38, 0.4 , 0.3 ,0.1), mgp = c(1.3, 0.4,0))
layout(lay.matrix,
       widths=c(1,1), heights=c(1,1))

semivariog.list.Lagos$Cosine$
semivariog.list.Lagos <- lapply(dist.list.Lagos, 
                          FUN = semivariog, 
                          response_variable = "carga_total", 
                          id_col = "codigo", 
                          data = max.summary.lagos,
                          max.dist = NULL,
                          plot = F,
                          breaks = 40 
                         )



lapply(names(semivariog.list.Lagos),function(x){plot(semivariog.list.Lagos[[x]], main = (x))})

test.matrix <- matrix(c(0,1,1,2,
                        1,0,2,1,
                        1,2,0,1,
                        2,1,1,0),4,4, byrow = T)
cov.test <- round(20*exp((-lagos.distance^2)/4),2)
eigen.dist <- eigen(cov.test, only.values = T)
max(eigen.dist$values)


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

lay.matrix <- matrix(c(1,1,2,2,3,3,4,4, 5,5, 6, 6),3, 4, byrow = TRUE)
par(mai = c(0.38, 0.4 , 0.3 ,0.1), mgp = c(1.3, 0.4,0))
layout(lay.matrix,
       widths=c(1,1), heights=c(1,1))


semivariog.list.Aysen <- lapply(dist.list.Aysen, 
                          FUN = semivariog, 
                          response_variable = "carga_total", 
                          id_col = "codigo", 
                          data = max.summary.aysen,
                          max.dist = NULL,
                          plot = F,
                          breaks = 60 
)

dim(max.summary.lagos)

lapply(names(semivariog.list.Aysen),function(x){plot(semivariog.list.Aysen[[x]], main = (x))})
dist.list.Aysen$LeastCost[1:nrow(dist.list.Aysen$LeastCost) <10 ,1:ncol(dist.list.Aysen$LeastCost) <10 ]

test <- dist.list.Aysen$LeastCost[lower.tri(dist.list.Aysen$LeastCost)]
diag(test) <- NA
sum(test < 5, na.rm = T)
length(test)
min(test, na.rm = T)
