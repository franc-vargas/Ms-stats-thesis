semivariog <- function(response_variable, id_col, dist_object, 
                       data, breaks = 10, max.dist = NULL, plot = T) {
  
  if("SpatialPointsDataFrame" %in% class(data)){
    newdata <- data.frame(data@data)
  }else{
    newdata <- data
  }
  
  id_vector <- newdata[[id_col]]
  x <- newdata[, c(id_col, response_variable)]
  prev_ids <- c()
  list_sums <- list()
  if(is.null(max.dist)){
    max.dist <- max(dist_object[rownames(dist_object)%in% id_vector,colnames(dist_object)%in% id_vector])/3
  }
  dist.frame <- dist_object |> as.data.frame()
  dist.frame$codigo <- rownames(dist.frame)
  dist.frame$codigo <- as.numeric(gsub("X", "",dist.frame$codigo))
  
  for (i in id_vector) {
    prev_ids <- append(prev_ids, i)
    if(length(id_vector) == length(prev_ids)) break
    zi <- x[x[[id_col]] == i, ][[response_variable]]
    
    new_vec <- x[!x[[id_col]] %in% prev_ids, ][[id_col]]
    new_vec <- sort(new_vec)
    ydist <- dist.frame[dist.frame[[id_col]] %in% new_vec, 
                        colnames(dist.frame) == i| colnames(dist.frame) == id_col]
    colnames(ydist)[1] <- "distance"
    
    x <- x[order(x[[id_col]]), c(id_col, response_variable)]
    x <- x[x[[id_col]]%in% new_vec, c(id_col,response_variable)]
    difference <- x[[response_variable]]
    
    x$diff <- difference
    x$value <- (x$diff-zi)
    x$semivariog <- (x$value)^2
    x$robust <- sqrt(abs(x$value))
    
    y <- left_join(x, ydist, by = c("codigo"))
    y$site <- i
    
    site_code <- paste(i)
    list_sums[[site_code]] <-  y
    
  }
  list_sums <- do.call("rbind", list_sums)
  if(is.null(max.dist)){
    rango <- diff(range(list_sums$distance))/breaks
  }else{
    list_sums <- list_sums[list_sums$distance <= max.dist,]
    rango <- diff(range(list_sums$distance))/breaks
  }
  
  aux.labels <- c(1:breaks)  
  aux.labels <- vapply(aux.labels, function(x)x*rango, 1)
  aux.breaks<- c(0,aux.labels)
  list_sums$lag <- cut(list_sums$distance, 
                       aux.breaks, 
                       round(aux.labels,3))

    summarized <- 
    aggregate(semivariog ~ lag,
              data = list_sums, 
              FUN = function(x) mean(x)/2)
  
  
  summarized.robust <- 
    aggregate(robust ~ lag,
              data = list_sums, 
              FUN = function(x) ((sum(x)/length(x))^4)/(0.914 + (0.988)/length(x)))
  
  summarized <- left_join(summarized, summarized.robust, by = c("lag"))
  summarized$lag <- round(as.numeric(summarized$lag)*rango,0)
  results <- list(data = list_sums,
                  plot.data = summarized, 
                  nbreaks = breaks)
  
  class(results) <- c("varioplot", class(results))
  if(plot){
    plot(results)
    invisible(results)
  }else{results}
  
  
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
