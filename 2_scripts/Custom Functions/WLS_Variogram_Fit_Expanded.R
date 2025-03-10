source("2_scripts/Custom Functions/Variogram_Objective_Functions.R")

fit_WLS_variogram <- function(custom_variogram, 
                              cov.model = "matern",
                              robust = F,
                              cloud = F,
                              fixed.kappa = NULL){
  if(!cloud){
    newdata <- custom_variogram$plot.data 
  }else{
    newdata <- custom_variogram$data[,-6]
    newdata$n <- 1
    robust <- F
    colnames(newdata)[colnames(newdata) == "distance"] <- "lag"
    
  }
  
    
    if(!robust){
    aux_nugget <- newdata$semivariog[1]
    aux_sigma <- mean(newdata$semivariog)
    aux_sigma <- density(newdata$semivariog)
    aux_sigma <- aux_sigma$x[aux_sigma$y == max(aux_sigma$y)]
    max_sigma <- max(newdata$semivariog)
    
  }else{
    aux_nugget <- newdata$robust[1]
    aux_sigma <- density(newdata$robust)
    aux_sigma <- aux_sigma$x[aux_sigma$y == max(aux_sigma$y)]
    max_sigma <- max(newdata$robust)
  }
  
  aux_kappa <- 0.5
  aux_phi <- 15
  
  if(is.null(fixed.kappa)){
    theta <- c(
      nugget = aux_nugget,
      sigma = aux_sigma, 
      kappa = aux_kappa,
      phi = aux_phi
               )
  }else{
    theta <- c(
      nugget = aux_nugget,
      sigma = aux_sigma, 
      phi = aux_phi
               )  
  }
  vec_nombres <- names(theta)
  
  if(cov.model == "matern"){
    cat("Using matern semivariogram")
    objective_function <- objective_matern
    arg_list <- list(fixed.kappa = fixed.kappa,
                     robust = robust)
    lower <- rep(0, length(theta))
    upper <- rep(Inf, length(theta))
    names(upper) <- vec_nombres
    upper[["kappa"]] <- 2
    theta[["nugget"]] <- 0
    upper[["nugget"]] <- aux_sigma
    upper[["sigma"]] <- max_sigma
    }else if(cov.model == "spherical"){
      cat("Using spherical semivariogram")
      objective_function <- objective_spherical
      arg_list <- list(robust = robust)
      theta <- theta[names(theta) != "kappa"]
      vec_nombres <- vec_nombres[vec_nombres != "kappa"]
      lower <- rep(0, length(theta))
      upper <- rep(Inf, length(theta))
        
      names(upper) <- vec_nombres
      theta[["nugget"]] <- 0
      upper[["nugget"]] <- aux_sigma
      upper[["sigma"]] <- max_sigma
      
    }else if(cov.model == "gaussian"){
      cat("Using gaussian semivariogram")
      objective_function <- objective_gaussian
      arg_list <- list(robust = robust)
      theta <- theta[names(theta) != "kappa"]
      vec_nombres <- vec_nombres[vec_nombres != "kappa"]
      
      lower <- rep(0, length(theta))
      upper <- rep(Inf, length(theta))
      
      names(upper) <- vec_nombres
      
      theta[["nugget"]] <- 0
      upper[["nugget"]] <- aux_sigma
      upper[["sigma"]] <- max_sigma
    }else if (cov.model == "cauchy"){
      cat("Using cauchy semivariogram")
      objective_function <- objective_cauchy
      arg_list <- list(fixed.kappa = fixed.kappa, 
                       robust = robust,
                       nugget = aux_nugget)
      
      lower <- rep(0, length(theta))
      upper <- rep(Inf, length(theta))
      names(upper) <- vec_nombres
      theta[["nugget"]] <- 0
      upper[["nugget"]] <- aux_sigma
      upper[["sigma"]] <- max_sigma
  }
  
  op <- optim(par = theta, 
              fn = objective_function, 
              data = newdata,
              args = arg_list,
              method = "L-BFGS-B",
              lower = lower,
              upper = upper,
              control = list(maxit = 500)
              )
  names(op$par) <- vec_nombres
  op$model <- cov.model
  op
}
  
