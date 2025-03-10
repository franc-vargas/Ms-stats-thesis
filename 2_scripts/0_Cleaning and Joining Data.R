library(sp)
library(sf)
library(stars)
library(tidyverse)



concesiones_peces <- read.csv2("1_data/concesiones_peces.csv")

cg <- readxl::read_excel("caligus_fom.xlsx")
concesiones_peces <- concesiones_peces %>% 
  filter(TIPOPORC == "AGUA Y FONDO") %>%  
  select(N_CODIGOCE, COORDENADA) %>% 
  mutate(N_CODIGOCE = as.double(N_CODIGOCE))

cg %>% left_join(concesiones_peces, by = c("Codigo" = "N_CODIGOCE"),
                 relationship = "many-to-many")

repeated <- 
  concesiones_peces %>% 
  count(N_CODIGOCE) %>% 
  arrange(desc(n)) %>% 
  filter(n>1) %>% 
  pull(N_CODIGOCE)

REPEATED <- new_list2 %>% filter(!N_CODIGOCE %in% repeated) %>% view()
REPEATED %>% mutate(COORDENADA = st_cast(COORDENADA))
shape <- read_sf(dsn = "./shapefiles/ConcesionesdeAcuiculturaporestadodetrámite_F.shp")

shapes_notduplicated <- shape %>% filter(!N_CODIGOCE %in% repeated) %>%
  filter(N_CODIGOCE != "0",
         (T_GRUPOESP == "PECES" | T_GRUPOESP == "SALMONES")) %>% 
  filter(TIPOPORC == "AGUA Y FONDO") %>%  
  select(N_CODIGOCE) %>% mutate(N_CODIGOCE = as.double(N_CODIGOCE)) 

shapes_notduplicated %>% count(N_CODIGOCE) %>% arrange(desc(n))

shapes_duplicated <- shape %>% filter(N_CODIGOCE %in% repeated)
shapes_duplicated <- shapes_duplicated[1:2,] %>% select(N_CODIGOCE)
shapes_new <- rbind(shapes_notduplicated, shapes_duplicated)
joined_database <- cg %>% left_join(shapes_new, by = c("Codigo" = "N_CODIGOCE"), 
                                    relationship = "many-to-one")

write_rds(joined_database, "joined_data.rds")

# 2_ Calculating Centroids ------------------------------------------------
library(sp)
library(sf)
library(stars)
library(tidyverse)
library(leaflet)
library(plotly)

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

