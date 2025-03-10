DistObjects <- function(x, kilometers = T){
  if(!"SpatialPointsDataFrame" %in% class(x)){
    stop("Needs SpatialPointsDataFrame")
  }
  lista.obj <- list()
  lista.obj$Cosine <- geosphere::distm(x, fun = geosphere::distCosine) |> as.matrix()
  lista.obj$Geodesic <- geosphere::distm(x, fun = geosphere::distGeo)|> as.matrix()
  lista.obj$Haversine <- geosphere::distm(x, fun = geosphere::distHaversine)|> as.matrix()
  lista.obj$VincentyEllipsoid <- geosphere::distm(x, fun = geosphere::distVincentyEllipsoid)|> as.matrix()
  lista.obj$Rhumb <- geosphere::distm(x, fun = geosphere::distRhumb)|> as.matrix()
  if(kilometers){
    lista.obj <- lapply(lista.obj, function(x) x/1000)
  }
  return(lista.obj)
}