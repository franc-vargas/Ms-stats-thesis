cs.variog <- function(response_variable, id_col, dist_object, 
                       data, breaks = 10, max.dist = NULL, plot = T) {
  
  if("SpatialPointsDataFrame" %in% class(data)){
    newdata <- data.frame(data@data)
  }else{
    newdata <- data
  }
  id_vector <- newdata[[id_col]]
  
  if(is.null(max.dist)){
    max.dist <- max(dist_object[rownames(dist_object)%in% id_vector,
                                colnames(dist_object)%in% id_vector])/3
  }
  
    
  value_vector <- as.vector(newdata[[response_variable]])
  names(value_vector) <- id_vector
  
  differences <- outer(value_vector, value_vector, "-")
  split.diff <- asplit(differences,1)
  split.diff <- sapply(split.diff, function(x) {
    idx <- zoo::index(x)[x == 0][1]+1
    x <- as.matrix(x[idx:length(x)])
    })[1:(length(split.diff)-1)]
  
  identifiers <- names(split.diff)
  
  split.diff <- 
    lapply(names(split.diff), function(x){
      y <- split.diff[[x]]
      colnames(y) <- x
      x <- y
      x
    })

  names(split.diff) <- identifiers
  split.diff <- sapply(split.diff, function(x){
    match.row <- match(rownames(x), rownames(dist_object))
    match.col <- match(colnames(x), colnames(dist_object))
    
    x <- cbind(x, as.matrix(dist_object[match.row, match.col]))
    colnames(x) <- c("value", "distance")
    x
  })

  split.diff <- as.list(split.diff)
  
  split.diff <- lapply(split.diff, function(x){
    x <- as.data.frame(x)
    x$id <- rownames(x)
    x
  })

  results <- do.call("rbind", split.diff)
  results$semivariog <- results$value^2
  results$robust <- sqrt(abs(results$value))

  if(is.null(max.dist)){
    rango <- diff(range(results$distance))/breaks
  }else{
    results <- results[results$distance <= max.dist,]
    rango <- diff(range(results$distance))/breaks
  }
  
  aux.labels <- c(1:breaks)  
  aux.labels <- vapply(aux.labels, function(x)x*rango, 1)
  aux.breaks<- c(0,aux.labels)
  results$lag <- cut(results$distance, 
                       aux.breaks, 
                       round(aux.labels,3))
  summarized <- 
    aggregate(semivariog ~ lag,
              data = results, 
              FUN = function(x) mean(x)/2)
  
  
  summarized.robust <- 
    aggregate(robust ~ lag,
              data = results, 
              FUN = function(x) {((sum(x)/length(x))^4)/(0.914 + (0.988)/length(x))}
              )
  summarized.length <- 
    aggregate(robust ~ lag,
              data = results, 
              FUN = function(x) {length(x)}
    ) |> rename(n = robust)
  
  summarized <- left_join(summarized, summarized.robust, by = c("lag"))
  summarized <- left_join(summarized, summarized.length, by = c("lag"))
  
  summarized$lag <- round(as.numeric(summarized$lag)*rango,0)
  final.list <- list(data = results,
                  plot.data = summarized, 
                  nbreaks = breaks)
  class(final.list) <- c("varioplot", class(final.list))
  if(plot){
    plot(final.list)
    invisible(final.list)
  }else{final.list}
  
  
}


plot.varioplot <- function(list, main = NULL, robust = F){
  
  if(robust != T){
    if(is.null(main)){main <- "Standard Semivariogram"}
    plot(semivariog ~ lag, data = list$plot.data, xlab = "Lag distance", ylab = "Semivariance",
         main = main)
  }else{
    if(is.null(main)){main <- "Cressie's Robust Semivariogram"}
    plot(robust ~ lag, data = list$plot.data, xlab = "Lag distance", ylab = "Semivariance",
         main = main)}
}
