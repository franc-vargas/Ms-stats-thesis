DistObjects <- function(x, kilometers = T){
  if(!"SpatialPointsDataFrame" %in% class(x)){
    stop("Needs SpatialPointsDataFrame")
  }
  lista.obj <- list()
  lista.obj$Geodesic <- geosphere::distm(x, fun = geosphere::distGeo)|> as.matrix()
  if(kilometers){
    lista.obj <- lapply(lista.obj, function(x) x/1000)
  }
  return(lista.obj)
}