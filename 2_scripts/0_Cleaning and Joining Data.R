library(sp)
library(sf)
library(stars)
library(tidyverse)
library(httr)
library(xml2)
library(XML)
library(rvest)

lista <- list()

for(i in seq_along(test$description)){
x <- read_html(test$description[i])
x <- x %>% html_element("table")
x <- x %>% html_table()
x <- x[3:nrow(x),1:2] %>% 
  rename(category = X1, value = X2) %>% 
  pivot_wider(names_from = category, values_from = value)
lista[[i]] <- x
}
new_list <- do.call("rbind", lista)
new_list2 <- new_list %>% filter(N_CODIGOCE != "0" & (T_GRUPOESP == "PECES" | T_GRUPOESP == "SALMONES"))
new_list$T_GRUPOESP %>% unique()
getwd()
setwd("/home/franc/Magister/thesis_fvega/1_data")

write_csv2(new_list, "concesiones_completas.csv")
write_csv2(new_list2, "concesiones_peces.csv")

cg <- readxl::read_excel("caligus_fom.xlsx")
new_list2 %>% colnames()
new_list3 <- new_list2 %>% 
  filter(TIPOPORC == "AGUA Y FONDO") %>%  
  select(N_CODIGOCE, COORDENADA) %>% 
  mutate(N_CODIGOCE = as.double(N_CODIGOCE))

cg %>% left_join(new_list3, by = c("Codigo" = "N_CODIGOCE"), relationship = "many-to-many")

repeated <- new_list3 %>% count(N_CODIGOCE) %>% arrange(desc(n)) %>% filter(n>1) %>% pull(N_CODIGOCE)

REPEATED <- new_list2 %>% filter(!N_CODIGOCE %in% repeated) %>% view()
REPEATED %>% mutate(COORDENADA = st_cast(COORDENADA))
shape <- read_sf(dsn = "./shapefiles/ConcesionesdeAcuiculturaporestadodetrámite_F.shp")

shapes_notduplicated <- shape %>% filter(!N_CODIGOCE %in% repeated) %>%
  filter(N_CODIGOCE != "0" & (T_GRUPOESP == "PECES" | T_GRUPOESP == "SALMONES")) %>% 
  filter(TIPOPORC == "AGUA Y FONDO") %>%  
  select(N_CODIGOCE) %>% mutate(N_CODIGOCE = as.double(N_CODIGOCE)) 

shapes_notduplicated %>% count(N_CODIGOCE) %>% arrange(desc(n))

shapes_duplicated <- shape %>% filter(N_CODIGOCE %in% repeated)
shapes_duplicated <- shapes_duplicated[1:2,] %>% select(N_CODIGOCE)
shapes_new <- rbind(shapes_notduplicated, shapes_duplicated)
joined_database <- cg %>% left_join(shapes_new, by = c("Codigo" = "N_CODIGOCE"), relationship = "many-to-one")
getwd()

write_rds(joined_database, "joined_data.rds")

# 2_ Calculating Centroids ------------------------------------------------
library(sp)
library(sf)
library(stars)
library(tidyverse)
library(leaflet)

sites <- read_rds("1_data/joined_data.rds") %>% pull(Codigo) %>% unique()

shape <- read_sf(dsn = "1_data/shapefiles/ConcesionesdeAcuiculturaporestadodetrámite_F.shp")
shape <- shape %>% filter(N_CODIGOCE %in% sites)

shape <- as_Spatial(shape)
shape <- spTransform(shape,  CRS("+ellps=WGS84 +proj=longlat +datum=WGS84 +no_defs"))

l <- leaflet() %>% 
  addProviderTiles("CartoDB.Positron") %>%
  setView(lat = '-42.659464', lng = '-73.042132', zoom = 7)


l %>% 
  addCircles(data = shape, color = "red", weight = 2,
             popup = shape@data$TITULAR) 
  
coordenadas <-coordinates(shape)

shape <- SpatialPointsDataFrame(coordenadas, data = shape@data)
write_rds(shape, "1_data/sites_coordinates.rds")
shape <- read_rds("1_data/sites_coordinates.rds")
sites <- read_rds("1_data/joined_data.rds")
sites <- 
  sites %>% mutate(
  year = year(fechaDeclaracion),
  month = month(fechaDeclaracion),
  total_adults = promHo+promAm)
sites_avg <-
  sites %>% 
  group_by(Codigo, nombreRegion, ACS, Especie, year, month) %>% 
    summarize(mean_cg = sum(total_adults))
sites_avg <- 
  sites_avg %>% mutate(
  period = paste0(year,"-", month, "-01"),
  period = as.Date(period)
)

shape@data %>% glimpse()
sites_avg %>% left_join(shape[c(1,3)], by = c("N_CODIGOCE" = "Codigo"))
shape2 <- shape
shape2@data <- shape2@data[c(1,3)]

colnames(shape2@data)[2] <- "Codigo"
shape2@data$coordenadas <- coordinates(shape2)

sites_avg2 <- sites_avg %>% left_join(shape2@data, by = "Codigo", relationship = "many-to-many")
sites_avg2 <- sites_avg2 %>% filter(!is.na(coordenadas[1]))

sum(duplicated(sites_avg2))
sites_avg3 <- SpatialPointsDataFrame(coords = sites_avg2$coordenadas, data = sites_avg2[c(1,2,3,6,7,8)], match.ID = F)

shape %>% glimpse()

l <- leaflet() %>% 
  addProviderTiles("CartoDB.Positron") %>%
  setView(lat = '-42.659464', lng = '-73.042132', zoom = 7)

sites_avg3@data %>% glimpse()
sites_avg3@proj4string <- CRS("+ellps=WGS84 +proj=longlat +datum=WGS84 +no_defs")
los_lagos <- subset(sites_avg3, nombreRegion == "X REGION")
aysen <- subset(sites_avg3, nombreRegion == "XI REGION")

Neighbour_18c <- subset(aysen, ACS == "ACS 18 C")
aysen@data %>% glimpse()
l %>% 
  addCircles(data = sites_avg3, color = "red", weight = 2,
             popup = sites_avg3@data$mean_cg) %>% 
  addTimeslider(data = sites_avg3@data$period)

library(plotly)
install.packages("plotly")
library(plotly)

install.packages("stpp")
los_lagos@data %>% 
plotly::plot_ly(lon = los_lagos@coords[,1], lat = los_lagos@coords[,2],
                 colours = ~I(colores),
                 frame = ~period,
                 mode = "scattergeo")

  los_lagos@data %>% glimpse()
  

los_lagos@data$colores <- cut(los_lagos@data$mean_cg,
                              breaks = c(0,20, 50, 100, Inf),
                              labels = c("red","orange", "yellow", "green"))
  quantile(los_lagos@data$mean_cg, 0.99999)
  
  
aysenPlot <- 
  plot_geo(color = I("red")) %>% 
add_markers(
  data = Neighbour_18c@data, 
  x = ~Neighbour_18c@coords[,1],
  y = ~Neighbour_18c@coords[,2],
  text = ~paste(Codigo, ": ",mean_cg),
  size = ~I(mean_cg), hoverinfo = "text", alpha = 0.5,
  frame = ~Neighbour_18c@data$period) %>% 
  plotly::animation_slider()


aysenPlot
sites_avg3@data %>% glimpse()

p <- 
  sites_avg3@data %>% ungroup() %>% summarize(
  sum_load = sum(mean_cg), .by = c(nombreRegion, ACS, period)) %>% 
  ggplot(aes(y = ACS, x = sum_load, fill = nombreRegion))+
  geom_col()+
  facet_wrap(~nombreRegion, scales = "free_y")
p2 <- ggplotly(p)

library("gganimate")
gganimate::animate()