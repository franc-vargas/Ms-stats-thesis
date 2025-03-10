# Libraries & Data --------------------------------------------------------

library(sp)
library(sf)
library(stars)
library(tidyverse)
library(leaflet)
library(patchwork)
library(igraph)
library(gdistance)

source("2_scripts/1_ImprovedSemivariogram.R")


# Calculo de maximos por ciclo --------------------------------------------

sites <- 
  readr::read_rds("1_data/joined_data.rds") |> janitor::clean_names()
sites <- 
  sites |> 
  dplyr::mutate(
    adultos_totales = prom_ho + prom_am,
    prop.juv = prom_juv/(prom_juv + adultos_totales))

# Seleccionamos solo campos relevantes para nuestro resumen.
# Consideraremos solo salmo salar y o. mykiss para el analisis

x.region <- 
  sites |> 
  dplyr::filter(
    nombre_region == "X REGION" & 
      especie %in% c("SALMON DEL ATLANTICO", "TRUCHA ARCOIRIS")) |> 
  dplyr::select(codigo, especie, fecha_primer_muestreo, fecha_declaracion,
         biomasa_centro, peso_promedio, numero_total_peces,
         prom_juv, prom_ho, prom_am, adultos_totales,  
         prop.juv, temperatura_promedio, 
         salinidad_promedio) |> 
  mutate(
    ddiff =  if_else(
      is.na(lag(fecha_declaracion)), 
      interval(fecha_declaracion-days(7), ymd(fecha_declaracion)),
      interval(ymd(lag(fecha_declaracion)), ymd(fecha_declaracion))
    ),
    ddiff = ddiff/weeks(1),
    .by = c(codigo, especie)) |> 
  mutate(new_cycle = 
           if_else(
             ddiff > 7 | fecha_declaracion == ymd(min(fecha_declaracion)),
             1, 
             0),
         .by = c(codigo, especie)) |> 
  mutate(cycle = cumsum(new_cycle), .by = c(codigo, especie))

x.region.f <- 
  x.region |>
  dplyr::filter(cycle > 1 & cycle < max(cycle), 
         .by = c(codigo, especie))

x.region.f <- x.region.f |> mutate(
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
  x.region.f |> 
  summarize(
    adultos_totales = max(adultos_totales, na.rm = T),
    .by = c(codigo, especie)) |> 
  filter(adultos_totales >1)
# vector_0 <- x_Max.Site |> filter(adultos_totales == 0) |> pull(codigo)

x_Max.Site <- 
  x_Max.Site |> 
  left_join(
    x.region.f[,c(1:2,5:14,16)],
    by =c("codigo","especie",  "adultos_totales"))

load(file = "shape_regions.RData")
rm(shape_aysen)

x_Max.Site <- 
  x_Max.Site |> 
  left_join(shape_lagos[,c(3,22)], by = c("codigo" = "N_CODIGOCE"))

x_Max.Site$lon <- x_Max.Site$coordinates[,1]
x_Max.Site$lat <- x_Max.Site$coordinates[,2]

max.cycle.sites <- 
  x.region.f |> 
  summarize(
    carga_total = max(carga_total, na.rm = T),
    .by = c(codigo, especie, year, estacion,cycle)) |> 
  filter(carga_total >1)

max.cycle.sites <- 
  max.cycle.sites |> 
  left_join(
    x.region.f[,c(1:2,6:14,16:21)],
    by =c("codigo","especie", "year", "estacion", "cycle", "carga_total" = "carga_total"))

max.summary <- 
  max.cycle.sites |> 
  summarize(
    carga_total = mean(carga_total, na.rm = T),
    peso_promedio = mean(peso_promedio),
    numero_total_peces = mean(numero_total_peces),
    prom_juv = mean(prom_juv),
    prop.juv = mean(prop.juv),
    temperatura_promedio = mean(temperatura_promedio), 
    salinidad_promedio = mean(salinidad_promedio),
    .by=c(codigo, especie, year, cycle)) |> 
  summarize(carga_total = mean(carga_total, na.rm = T),
            peso_promedio = mean(peso_promedio),
            numero_total_peces = mean(numero_total_peces),
            prom_juv = mean(prom_juv),
            prop.juv = mean(prop.juv),
            temperatura_promedio = mean(temperatura_promedio), 
            salinidad_promedio = mean(salinidad_promedio),
            .by=c(codigo, especie))



# RECUERDA DESPUES BUSCAR MOMENTOS SIMILARES EN EL TIEMPO PARA LOS CENTROS

lagos.base <- chilemapas::generar_regiones() |> dplyr::filter(codigo_region == 10)

  ggplot(lagos.base)+
  geom_sf()+
  geom_point(data = x_Max.Site, aes(x = lon,
                                       y = lat,
                                       colour = prom_ho,
                                    shape = especie
  ),
  alpha = 0.8)+
    coord_sf(ylim = c(-44, -41.3), xlim = c(-74.8, -71.5))+
    labs(x = "longitud",
       y = "latitud")+
  scale_colour_gradientn(colours = c("darkgreen","yellow", "darkred"))+
  theme_bw()+
  labs(subtitle = "Aysén",
       colour = "Maximo de Adultos Totales")+
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.position = "right",
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 9))


vec_cols <-
  x_Max.Site |> select(biomasa_centro, prop.juv, salinidad_promedio, temperatura_promedio, numero_total_peces,
                                peso_promedio) |> colnames()
x_Max.Site |> 
  pivot_longer(cols = vec_cols, 
               names_to = "covariates", 
               values_to = "values") |> 
  ggplot(aes(x = values, y = adultos_totales, colour = especie))+
  geom_point()+
  facet_wrap(~covariates, scales = "free_x")


# Geostats ----------------------------------------------------------------
library(sp)
library(sf)
library(stars)

m1 <- 
  spmodel::splm(
  I(log(adultos_totales))~biomasa_centro+prop.juv+salinidad_promedio:especie+temperatura_promedio, 
                    x_Max.Site,
  xcoord = lon,
  ycoord = lat,
  spcov_type = "gaussian")
summary(m1)
plot(x_Max.Site[,2:9])

xmaxnew <- x_Max.Site |> select(-codigo, -coordinates)

full_formula <- paste(
  "adultos_totales", 
  paste(colnames(xmaxnew[,c(1,4:5,9, 11:10)]), 
        collapse = "+"), 
  sep = "~")

ffull <- as.formula(full_formula)
m2 <- 
  spmodel::splm(
    formula = ffull, 
    data = xmaxnew,
    xcoord = lon,
    ycoord = lat,
    spcov_type = "gaussian")

summary(m2) 
ggplot()+
  geom_density(aes(x = m2$fitted$response), fill = "darkred", alpha = 0.5)+
  geom_density(aes((x_Max.Site$adultos_totales)), fill = "darkblue", alpha = 0.5)

library(spmodel)

single_site <- x.region.f |> filter(codigo == 100619) |> arrange(fecha_declaracion)

# geoR Tests --------------------------------------------------------------
library(geoR)
library(gstat)

max2 <- x_Max.Site |> select(-coordinates)
geomax <- as.geodata(max2, coords.col = 13:14) 
unique(max2$codigo) |> length()
max2 |> count(codigo) |> arrange(desc(n))
variograma <- variog(geomax, bin.cloud =  T,
                     max.dist = 2,
                     estimator.type = "modulus",
                     tolerance = pi/2)
plot(variograma)


max.cycle.sites <- 
  x.region.f |> 
  summarize(
    carga_total = max(carga_total, na.rm = T),
    .by = c(codigo, especie, year, estacion,cycle)) |> 
  filter(carga_total >1)

max.cycle.sites <- 
  max.cycle.sites |> 
  left_join(
    x.region.f[,c(1:2,6:14,16:21)],
    by =c("codigo","especie", "year", "estacion", "cycle", "carga_total" = "carga_total"))

max.summary <- 
  max.cycle.sites |> 
    summarize(
      carga_total = mean(carga_total, na.rm = T),
      peso_promedio = mean(peso_promedio),
      numero_total_peces = mean(numero_total_peces),
      prom_juv = mean(prom_juv),
      prop.juv = mean(prop.juv),
      temperatura_promedio = mean(temperatura_promedio), 
      salinidad_promedio = mean(salinidad_promedio),
      .by=c(codigo, especie, year, cycle)) |> 
  summarize(carga_total = mean(carga_total, na.rm = T),
            peso_promedio = mean(peso_promedio),
            numero_total_peces = mean(numero_total_peces),
            prom_juv = mean(prom_juv),
            prop.juv = mean(prop.juv),
            temperatura_promedio = mean(temperatura_promedio), 
            salinidad_promedio = mean(salinidad_promedio),
            .by=c(codigo, especie))

load(file = "shape_regions.RData")
rm(shape_aysen)

max.summary <- 
  max.summary |> 
  left_join(shape_lagos[,c(3,22)], by = c("codigo" = "N_CODIGOCE"))

max.summary$lon <- max.summary$coordinates[,1]
max.summary$lat <- max.summary$coordinates[,2]

geomax <- 
  geoR::as.geodata(
  max.summary, coords.col = 11:12) 
variograma <- geoR::variog(geomax, bin.cloud =  T)
plot(variograma)
sp::coordinates(max.summary) <- ~lon+lat
GeoModels::GeoVariogram()
g <- gstat::gstat(id = "ln.carga_total", 
                  formula = 
                    (carga_total)^(1/3)~
                    especie+prop.juv+
                    numero_total_peces+
                    peso_promedio+
                    temperatura_promedio+
                    salinidad_promedio, data = max.summary,
                  vdist = F)

# lagos.distance <- readr::read_rds("1_data/Lagos_DistanceMatrix.rds")
# ids <- max.summary[["codigo"]]
# lg.dist.subset <- lagos.distance[rownames(lagos.distance)%in%ids, colnames(lagos.distance) %in% ids]
g$data$ln.carga_total$data@proj4string <- sp::CRS("+ellps=WGS84 +proj=longlat +datum=WGS84 +no_defs +units=m")
vgram <- gstat::variogram(g,
                          cressie = F,
                          width = 2)

test <- gstat::variogramLine(vgm(0.4, "Sph", 50, 0.4), maxdist = 70)
par(mfrow = c(1,2))

plot(x = vgram$dist, y =vgram$gamma)
plot(test, type = "l", add = T, ylim = c(0, 1.2))

ggplot(max.summary@data)+geom_density(aes(x = carga_total))
ggplot(max.summary@data)+geom_density(aes(x = (carga_total)^(1/3)))


# Prueba en la misma fecha ------------------------------------------------

x.region.f |> glimpse()
same.date <- x.region.f2 |> filter(fecha_declaracion == date("2019-07-15"))
same.date <- same.date |> left_join(shape_lagos[,c(3,22)], by = c("codigo" = "N_CODIGOCE")) |> 
  filter(carga_total > 1)
same.date$lon <- same.date$coordinates[,1]
same.date$lat <- same.date$coordinates[,2]
sp::coordinates(same.date) <- ~lon+lat

g <- gstat::gstat(id = "ln.carga_total", 
                  formula = 
                    log(carga_total)~1,
                    # especie+
                    # numero_total_peces+
                    # peso_promedio+
                    # temperatura_promedio+
                    # salinidad_promedio, 
                    data = same.date,
                  vdist = F)
# lagos.distance <- readr::read_rds("1_data/Lagos_DistanceMatrix.rds")
# ids <- max.summary[["codigo"]]
# lg.dist.subset <- lagos.distance[rownames(lagos.distance)%in%ids, colnames(lagos.distance) %in% ids]
g$data$ln.carga_total$data@proj4string <- sp::CRS("+ellps=WGS84 +proj=longlat +datum=WGS84 +no_defs +units=m")
vgram <- gstat::variogram(g, width = 3)

# test <- gstat::variogramLine(vgm(0.4, "Sph", 50, 0.4), maxdist = 70)
# par(mfrow = c(1,1))
plot(vgram)
plot(x = vgram$dist, y =vgram$gamma)
plot(test, type = "l", ylim = c(0, 1.2))

lagos.base <- chilemapas::generar_regiones() |> dplyr::filter(codigo_region == 10)

ggplot(lagos.base)+
  geom_sf()+
  geom_point(data = as.data.frame(g$data$ln.carga_total$data), aes(x = lon,
                                    y = lat,
                                    colour = adultos_totales,
                                    shape = especie
  ),
  alpha = 0.8)+
  coord_sf(ylim = c(-44, -41.3), xlim = c(-74.8, -71.5))+
  labs(x = "longitud",
       y = "latitud")+
  scale_colour_gradientn(colours = c("darkgreen","yellow", "darkred"))+
  theme_bw()+
  labs(subtitle = "Aysén",
       colour = "Maximo de Adultos Totales")+
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.position = "right",
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 9))




# Custom variogram --------------------------------------------------------

max.summary <- 
  max.cycle.sites |> 
  summarize(across(where(is.numeric), ~mean(.x, na.rm = T)), 
            .by=c(codigo, especie, year, cycle)) |> 
  summarize(across(where(is.numeric), ~mean(.x, na.rm = T)), 
            .by=c(codigo, especie)) 

lagos.distance <- readr::read_rds("1_data/Lagos_DistanceMatrix.rds")


max.summary <- 
  max.summary |> 
  left_join(shape_lagos[,c(3,22)], by = c("codigo" = "N_CODIGOCE"))

max.summary$lon <- max.summary$coordinates[,1]
max.summary$lat <- max.summary$coordinates[,2]


ids <- max.summary[["codigo"]]
lg.dist.subset <- lagos.distance[rownames(lagos.distance)%in%ids, colnames(lagos.distance) %in% ids]
names(lagos.distance[, colnames(lagos.distance) == 100124])
lg.dist.subset <- lg.dist.subset[sort(rownames(lg.dist.subset)), sort(colnames(lg.dist.subset))]
distance <- as.dist(lg.dist.subset)/1000

max.summary.unique <- max.summary[c(1:37, c(39:nrow(max.summary))),]
euc.coords <- as.matrix(max.summary.unique[,c(18,19)])
row.names(euc.coords) <- max.summary.unique$codigo
new <- sf::as_Spatial()
sf::st_transform(euc.coords,"+ellps=WGS84 +proj=longlat +datum=WGS84 +no_defs +units=m")

euc.distance<- dist(euc.coords, method = "euclidean")
manh.distance<- dist(euc.coords, method = "manhattan")
canb.distance<- dist(euc.coords, method = "canberra")
mink.distance<- dist(euc.coords, method = "minkowski")
max.distance<- dist(euc.coords, method = "maximum")

names(distance)
names(euc.distance) <- max.summary[["codigo"]]
names(euc.distance)

ggplot()+geom_point(aes(x = euc.distance, y = distance),alpha = .3)
ggplot()+geom_point(aes(x = manh.distance, y = distance),alpha = .3)
ggplot()+geom_point(aes(x = canb.distance, y = distance),alpha = .3)
ggplot()+geom_point(aes(x = mink.distance, y = distance),alpha = .3)
ggplot()+geom_point(aes(x = max.distance, y = distance),alpha = .3)

(names(distance)) == names(euc.distance)

sort(names(distance))


plot(x = euc.distance, y = distance)

# Revisar funcion variogramLine


# nlme Tests --------------------------------------------------------------

library(nlme)
lagos.distance <- readr::read_rds("1_data/Lagos_DistanceMatrix.rds")
lg.dist.subset <- lg.dist.subset[sort(rownames(lg.dist.subset)), sort(colnames(lg.dist.subset))]

distance <- as.dist(lg.dist.subset)/1000


nlme::Variogram()

S3method(Variogram, lme)
utils::getS3method("Variogram", "lme")

# GeoModels Modiffication -------------------------------------------------
library(GeoModels)

lagos.distance <- readr::read_rds("1_data/Lagos_DistanceMatrix.rds")
max.cycle.sites$year   <- factor(max.cycle.sites$year)
max.cycle.sites$month   <- factor(max.cycle.sites$month)
max.cycle.sites$cycle   <- factor(max.cycle.sites$cycle)

max.summary <- 
  max.cycle.sites |> 
  summarize(across(where(is.numeric), ~mean(.x, na.rm = T)), 
            .by=c(codigo, especie, year, cycle)) |> 
  summarize(across(where(is.numeric), ~mean(.x, na.rm = T)), 
            .by=c(codigo, especie)) 



max.summary <- 
  max.summary |> 
  left_join(shape_lagos[,c(3,22)], by = c("codigo" = "N_CODIGOCE")) |> 
  filter(codigo != 101558)

max.summary$lon <- max.summary$coordinates[,1]
max.summary$lat <- max.summary$coordinates[,2]


ids <- max.summary[["codigo"]]
lg.dist.subset <- lagos.distance[rownames(lagos.distance)%in%ids, colnames(lagos.distance) %in% ids]
dist.subset <- as.dist(lg.dist.subset)/1000

lagos.distance[rownames(lagos.distance) == 101964,]


lagos.distance |> view()
rownames(lagos.distance)
lagos.distance |> dim()
site_101558r <- lagos.distance[rownames(lagos.distance) == 101558, ]
lagos.distance <- rbind(lagos.distance, site_101558r)

site_101558c <- lagos.distance[,colnames(lagos.distance) == 101558]

lagos.distance <- cbind(lagos.distance, site_101558c)
max.summary |> count(codigo) |> arrange(desc(n))
lg.dist.subset <- lg.dist.subset[sort(rownames(lg.dist.subset)), sort(colnames(lg.dist.subset))]

test <- cs.variog("adultos_totales","codigo", lagos.distance/1000, max.summary, breaks = 10,
                   max.dist = 250, plot = F)

plot(max.summary[,-c(1,3,4,9:10,15:18)])
glimpse(max.summary)

plot((adultos_totales) ~ I((prop.juv)^2), max.summary)
plot((adultos_totales) ~ I((1/(salinidad_promedio))), max.summary)
plot((numero_total_peces) ~ (peso_promedio), max.summary)


m1 <- lm(I((adultos_totales)) ~ peso_promedio + I((prom_juv))+ 
           I((salinidad_promedio)),max.summary)

summary(m1)

max.summary$residuos <- ls.diag(m1)$std.res
car::vif(m1, type="predictor")

variograma.residuos <- 
  cs.variog("residuos.abs","codigo", lagos.distance/1000, max.summary, breaks = 20, plot =F,
            max.dist = 50)
plot(variograma.residuos, robust = T)
plot(variograma.residuos, robust = F)

plot(variograma.residuos$data$distance,variograma.residuos$data$semivariog)
plot(variograma.residuos$data$distance,variograma.residuos$data$robust^4)

variograma.residuos$data[variograma.residuos$data$semivariog == max(variograma.residuos$data$semivariog),]
max.summary |> filter(codigo %in% c(100413, 100680)) |> glimpse()
max.summary$residuos.sq <- max.summary$residuos^2
max.summary$residuos.abs <- abs(max.summary$residuos)

max(max.summary$residuos^2)
