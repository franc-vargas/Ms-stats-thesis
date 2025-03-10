objective_matern <- function(data,
                               theta,
                             args = list(fixed.kappa = NULL,
                                         #nugget = nugget,
                               robust = F)
                             ){
  h <- data[,"lag"]

  robust <- args[["robust"]]
  fixed.kappa <- args[["fixed.kappa"]]
  n <- data[,"n"]
  if(!robust) {
    empirical <- data[,"semivariog"]
  }else{
    empirical <- data[,"robust"]
  }
  if(is.null(fixed.kappa)){
    kappa <- theta[["kappa"]]
  }else{
    kappa <- fixed.kappa
  }
  phi <- theta[["phi"]]
  sigma <- theta[["sigma"]]
  nugget <- theta[["nugget"]]
  
  
  theoretical <- 
           sigma+nugget - 
             geoR::cov.spatial(h, 
                               cov.model="matern", 
                               cov.pars=c(sigma, phi),
                               kappa)
  
  aux <-  (empirical - theoretical)^2
  sum((n/empirical^2) * aux)
}

objective_spherical <- function(data,
                             theta,
                             args = list(robust = F,
                                         #nugget = nugget
                                         )
                             ){
  h <- data[,"lag"]

  n <- data[,"n"]
  
  robust <- args[["robust"]]
  if(!robust) {
    empirical <- data[,"semivariog"]
  }else{
    empirical <- data[,"robust"]
  }
  sigma <- theta[["sigma"]]
  phi <- theta[["phi"]]
  nugget <- theta[["nugget"]]
  
  theoretical <- 
           sigma+nugget - geoR::cov.spatial(h, 
                                            cov.model="spherical", 
                                            cov.pars=c(sigma, phi))

  aux <-  (empirical - theoretical)^2
  sum((n/empirical^2) * aux)
}


objective_gaussian <- function(data,
                                theta,
                               args = list(robust = F,
                                           #nugget = nugget
                                           )
                               ){
  h <- data[,"lag"]

  n <- data[,"n"]
  
  robust <- args[["robust"]]
  if(!robust) {
    empirical <- data[,"semivariog"]
  }else{
    empirical <- data[,"robust"]
  }
  sigma <- theta[["sigma"]]
  phi <- theta[["phi"]]
  nugget <- theta[["nugget"]]
  
  theoretical <- 
           sigma+nugget - geoR::cov.spatial(h, 
                                            cov.model="gaussian", 
                                            cov.pars=c(sigma, phi))
    
  aux <-  (empirical - theoretical)^2
  sum((n/empirical^2) * aux)
}

objective_cauchy <- function(data,
                             theta,
                             args = list(fixed.kappa = NULL,
                                         #nugget = nugget,
                                         robust = F)
){
  h <- data[,"lag"]

  robust <- args[["robust"]]
  fixed.kappa <- args[["fixed.kappa"]]
#  nugget <- args[["nugget"]]
  n <- data[,"n"]
  if(!robust) {
    empirical <- data[,"semivariog"]
  }else{
    empirical <- data[,"robust"]
  }
  if(is.null(fixed.kappa)){
    kappa <- theta[["kappa"]]
  }else{
    kappa <- fixed.kappa
  }
  phi <- theta[["phi"]]
  sigma <- theta[["sigma"]]
  nugget <- theta[["nugget"]]
  
  theoretical <- 
    sigma+nugget - geoR::cov.spatial(h, 
                                     cov.model="cauchy", 
                                     cov.pars=c(sigma, phi),
                                     kappa)
  
  aux <-  (empirical - theoretical)^2
 sum((n/empirical^2) * aux)
}


