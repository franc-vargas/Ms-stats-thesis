fit_WLS_variogram <- function(custom_variogram, 
                              cov.model = "matern",
                              robust = F,
                              fixed.kappa = NULL){
newdata <- custom_variogram$plot.data 

if(!robust){
  aux_nugget <- newdata$semivariog[1]
  aux_sigma <- mean(newdata$semivariog)
  aux_sigma <- density(newdata$semivariog)
  aux_sigma <- aux_sigma$x[aux_sigma$y == max(aux_sigma$y)]
  
  
}else{
  aux_nugget <- newdata$robust[1]
  aux_sigma <- density(newdata$robust)
  aux_sigma <- aux_sigma$x[aux_sigma$y == max(aux_sigma$y)]
  }

aux_kappa <- 0.5
aux_phi <- 5

if(is.null(fixed.kappa)){
  theta <- c(nugget = aux_nugget,
             sigma = aux_sigma, 
             kappa = aux_kappa,
             phi = aux_phi)
  }else{
  theta <- c(nugget = aux_nugget, 
             sigma = aux_sigma, 
             phi = aux_phi)  
  }
vec_nombres <- names(theta)

if(cov.model == "matern"){
  cat("Using matern semivariogram")
}else{
    stop("No othe variogram function implemented")
}

op <- optim(par = theta, 
            fn = objective_function, 
            data = newdata,
            fixed.kappa = fixed.kappa,
            robust = robust,
            method = "L-BFGS-B",
            lower = rep(0, length(theta)))
names(op$par) <- vec_nombres
op
}

objective_function <- function(data,
                               theta,
                               fixed.kappa = NULL,
                               robust = F){
  h <- data[,"lag"]
  h <- h - min(h)
  
  n <- data[,"n"]
  if(!robust) {
    empirical <- data[,"semivariog"]
  }else{
    empirical <- data[,"robust"]
    }
  sigma <- theta[["sigma"]]
  if(is.null(fixed.kappa)){
    kappa <- theta[["kappa"]]
    }else{
    kappa <- fixed.kappa
  }
  phi <- theta[["phi"]]
  nugget <- theta[["nugget"]]
  
  theoretical <- ifelse(h > 0,
                        nugget + 
    (sigma *
    (1 -
       ((2^(1-kappa))/gamma(kappa)) * 
       ((h/phi)^kappa) * 
       (besselK(h/phi,kappa)) * 
       (h/phi))), nugget)
 aux <- ( (empirical/theoretical) - 1 )^2
 sum(n * aux)
}
