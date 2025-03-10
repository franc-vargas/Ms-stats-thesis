library(sp)
library(sf)
library(stars)
library(ggplot2)
library(dplyr)
library(tidyverse)

load(file = "shape_regions.RData")
full_shape <- rbind(shape_lagos, shape_aysen)
source("1_ImprovedSemivariogram.R")
 source("DistObjects.R")

 lagos.distance <- readr::read_rds("Lagos_DistanceMatrix.rds")
 aysen.distance <- readr::read_rds("Aysen_DistanceMatrix.rds")

 sites <-
   readr::read_rds("joined_data.rds") |>
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

max.summary.lagos <-
   max.summary |>
   filter(codigo != 101558 & nombre_region == "X REGION") |>
   arrange(codigo)
max(max.summary.lagos$max.adultos_totales)
 coordinates(max.summary.lagos) <- ~coordinates
 proj4string(max.summary.lagos) <- "+proj=longlat +ellps=WGS84"
ggplot(max.summary.lagos@data, aes(x = coordinates[,1],
           y = coordinates[,2],
           size = residuos.log))+
  geom_point()
exp(2)
nrow(max.summary.lagos)
nrow(max.summary.aysen)
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

 m.lagos.raw <- lm(I((max.adultos_totales)) ~ temperatura_promedio+
            I((salinidad_promedio)),max.summary.lagos@data)
 m.aysen.raw <- lm(I((max.adultos_totales)) ~ temperatura_promedio+
                 I((salinidad_promedio)),max.summary.aysen@data)

 m.lagos.log <- lm(I((log_max.adultos_totales)) ~ temperatura_promedio+
                     I((salinidad_promedio)),max.summary.lagos@data)
 m.aysen.log <- lm(I((log_max.adultos_totales)) ~ peso_promedio + temperatura_promedio+
                     I((salinidad_promedio)),max.summary.aysen@data)

 max.summary.lagos@data$residuos <- resid(m.lagos.raw)
 max.summary.lagos@data$residuos.log <- resid(m.lagos.log)
 max.summary.aysen@data$residuos <- resid(m.aysen.raw)
 max.summary.aysen@data$residuos.log <- resid(m.aysen.log)

save(max.summary.lagos, max.summary.aysen,dist.list.Lagos, dist.list.Aysen, file = "WorkData.RData")

 # Define server logic required to draw a histogram
function(input, output, session) {
load("WorkData.RData")

semivariog.list.Lagos <- reactive({lapply(dist.list.Lagos, 
                                    FUN = cs.variog, 
                                    response_variable = input$lagos.Response, 
                                    id_col = "codigo", 
                                    data = max.summary.lagos,
                                    max.dist = input$lagos.distmax,
                                    plot = F,
                                    breaks = input$lagos.breaks 
)})

semivariog.list.Aysen <- reactive({lapply(dist.list.Aysen, 
                                          FUN = cs.variog, 
                                          response_variable = input$aysen.Response, 
                                          id_col = "codigo", 
                                          data = max.summary.aysen,
                                          max.dist = input$aysen.distmax,
                                          plot = F,
                                          breaks = input$aysen.breaks 
)})

output$LC_Semivariogram.Lagos <- renderPlot({
  plot(x = semivariog.list.Lagos()$LeastCost$plot.data$lag,
       y = semivariog.list.Lagos()$LeastCost$plot.data$semivariog,
       xlab = "Distance", ylab = "Semivariance",
       main = "Least-Cost Standard Semivariogram")
})

output$LC_Robust_Semivariogram.Lagos <- renderPlot({
  plot(x = semivariog.list.Lagos()$LeastCost$plot.data$lag,
       y = semivariog.list.Lagos()$LeastCost$plot.data$robust,
       xlab = "Distance", ylab = "Semivariance",
       main = "Least-Cost Robust Semivariogram")
})
output$Geodesic.Semivariogram.Lagos <- renderPlot({
  plot(x = semivariog.list.Lagos()$Geodesic$plot.data$lag,
       y = semivariog.list.Lagos()$Geodesic$plot.data$semivariog,
       xlab = "Distance", ylab = "Semivariance",
       main = "Geodesic Standard Semivariogram")
})
output$Geodesic.Robust.Lagos <- renderPlot({
  plot(x = semivariog.list.Lagos()$Geodesic$plot.data$lag,
       y = semivariog.list.Lagos()$Geodesic$plot.data$robust,
       xlab = "Distance", ylab = "Semivariance",
       main = "Geodesic Robust Semivariogram")
})

output$LC_Semivariogram.Aysen <- renderPlot({
  plot(x = semivariog.list.Aysen()$LeastCost$plot.data$lag,
       y = semivariog.list.Aysen()$LeastCost$plot.data$semivariog,
       xlab = "Distance", ylab = "Semivariance",
       main = "Least-Cost Standard Semivariogram")
})

output$LC_Robust_Semivariogram.Aysen <- renderPlot({
  plot(x = semivariog.list.Aysen()$LeastCost$plot.data$lag,
       y = semivariog.list.Aysen()$LeastCost$plot.data$robust,
       xlab = "Distance", ylab = "Semivariance",
       main = "Least-Cost Robust Semivariogram")
})
output$Geodesic.Semivariogram.Aysen <- renderPlot({
  plot(x = semivariog.list.Aysen()$Geodesic$plot.data$lag,
       y = semivariog.list.Aysen()$Geodesic$plot.data$semivariog,
       xlab = "Distance", ylab = "Semivariance",
       main = "Geodesic Standard Semivariogram")
})
output$Geodesic.Robust.Aysen <- renderPlot({
  plot(x = semivariog.list.Aysen()$Geodesic$plot.data$lag,
       y = semivariog.list.Aysen()$Geodesic$plot.data$robust,
       xlab = "Distance", ylab = "Semivariance",
       main = "Geodesic Robust Semivariogram")
})

map <- leaflet() %>% addTiles %>% setView(-74, -45.5, zoom = 8)

output$leaflet.map <- renderLeaflet(map)
    
}
