krig.manual <- function(media0, media, data, s0=NULL,
                      coords=NULL, InvSigma=NULL, 
                      Sigma=NULL, c0=NULL, 
                      H=NULL, h0=NULL,
                      covmodel="matern", param=NULL,
                      type = "ord",...,
                      kappa = NULL){
  # Se calcula Sigma en caso de omitirse
  if(is.null(Sigma)){
    if(is.null(H)){
      H=as.matrix(dist(coords))
    }
    n=nrow(H)
    Sigma <- geoR::cov.spatial(H, cov.model=covmodel, 
                         cov.pars=c(
                           param[['sigma2']], 
                           param[['phi']]),..., 
                         kappa = ifelse(is.null(kappa),
                                        0.5, 
                                        kappa) ) + param[['nugget']]*diag(n) 
  }
  
  # Se calcula la inversa de Sigma en caso de omitirse
  if(is.null(InvSigma)){
    InvSigma=solve(Sigma) 
  }
  
  # Se calcula la covarianza en caso de omitirse
  if(is.null(c0)){
    if(is.null(h0)){ 
      h0=dist0(coords, s0)  
    }
    c0 <- geoR::cov.spatial(h0, 
                            cov.model=covmodel, 
                      cov.pars=c(param[['sigma2']], 
                                 param[['phi']])
                      )
    c0[h0==0] <- param[['sigma2']] + param[['nugget']]
  }
  
  # determinación de los pesos optimos
  
  sigma2 <- param[["sigma2"]]
  nugget <- param[["nugget"]]
  # prediccion y varianza de predicción
  if(type == "simp"){
      lambda=InvSigma %*% c0
      pred = media0 + t(lambda) %*% (data - media)
    
    var.pred <- 
      diag(t(lambda) %*% Sigma %*% lambda) -
      2 * diag( t(lambda) %*% c0 ) +
      param[['sigma2']] + param[['nugget']]
  }
  
  if(type == "ord"){
    C <- solve(Sigma)
    
    pred <- apply(c0,2,function(x){
      
      x0 <- c(x)
      aux_c0 <- sigma2 + nugget
      aux_sigma0 <- sigma2 + nugget
      aux_vec <- rep(1, nrow(C))
      lambda_OK <- drop(C %*% (x + aux_vec %*% 
                                 ( solve(t(aux_vec) %*% C %*% aux_vec) )%*%
                                 (1 - t(aux_vec) %*% C %*% x)))

      ok_pred <- drop(lambda_OK  %*% data)
      
      ok_pred.var <- 
        aux_sigma0 - (t(x0) %*% C %*% x0) + 
        t(1 - t(aux_vec) %*% C %*% x0) %*% 
        solve(t(aux_vec) %*% C %*% aux_vec) %*%
        (1 - t(aux_vec) %*% C %*% x0)
      
      x <- matrix(c(ok_pred, ok_pred.var), ncol = 2)
      x
    }, simplify = T)
   
    pred <- t(pred)
    colnames(pred) <- c("pred", "var.pred")
    # var.pred <- apply(c0,2,function(x){
    #   
    #   x0 <- c(x)
    #   aux_sigma0 <- sigma2 + nugget
    #   aux_vec <- rep(1, nrow(C))
    #   
    #   lambda_OK <- drop(C %*% (x + aux_vec %*% 
    #                             ( solve(t(aux_vec) %*% C %*% aux_vec) )%*%
    #                         (1 - t(aux_vec) %*% C %*% x)))
    #   
    #  ok_pred.var <- 
    #     aux_sigma0 - (t(x0) %*% C %*% x0) + 
    #     t(1 - t(aux_vec) %*% C %*% x0) %*% 
    #     solve(t(aux_vec) %*% C %*% aux_vec) %*%
    #     (1 - t(aux_vec) %*% C %*% x0)
    # 
    #   ok_pred.var <- drop(ok_pred.var)
    #   ok_pred.var 
    # })
    
  }
  
  if(type == "lognormal"){
    C <- solve(Sigma)
    
    optim_normal <- 
      function(x, 
               theta){
        mu <- theta[["mu"]]
        sigma <- theta[["sigma"]]
        n  <- length(data)
        logl<- -.5*n*log(2*pi) -.5*n*log(sigma) -
          (1/(2*sigma))*sum((x-mu)**2)
        
        return(-logl)
      }
    theta <- c(mu = mean(data), sigma = var(data))
    
    parametros <- 
      optim(theta,
          fn = optim_normal, 
          x = data,
          lower = c(-Inf, 0), 
          upper = c(Inf, Inf), method = "L-BFGS-B")
    print(parametros$par)
    pred <- apply(c0,2,function(x){
      
      x0 <- c(x)
      aux_c0 <- sigma2 + nugget
      aux_sigma0 <- sigma2 + nugget
      aux_vec <- rep(1, nrow(C))
      lambda_OK <- drop(C %*% (x + aux_vec %*% 
                                 ( solve(t(aux_vec) %*% C %*% aux_vec) )%*%
                                 (1 - t(aux_vec) %*% C %*% x)))
      
      ok_pred <- drop(lambda_OK  %*% data)
      
     
      
      ok_pred.var <- 
        aux_sigma0 - (t(x0) %*% C %*% x0) + 
        t(1 - t(aux_vec) %*% C %*% x0) %*% 
        solve(t(aux_vec) %*% C %*% aux_vec) %*%
        (1 - t(aux_vec) %*% C %*% x0)
    
      m <- (1 - t(aux_vec) %*% C %*% x0 ) / (t(aux_vec %*% C %*% aux_vec))
      mu_y <- solve(t(aux_vec %*% C %*% aux_vec)) %*% aux_vec %*% C %*% data
      lognormal_pred <- exp(ok_pred - ok_pred.var/2 - m)

      lognormal_var <- 
         exp(2 * parametros$par[1] + parametros$par[2]) * 
        exp(parametros$par[2]) * 
        (1 + (exp(-ok_pred.var + m)) * (exp(m) - 2) )
      
      x <- matrix(c(lognormal_pred, lognormal_var), ncol = 2)
      x
    }, 
    simplify = T)
    
    
    pred <- t(pred)
   
      colnames(pred) <- c("pred", "var.pred")
      
   
    
    
  }
  
  
  return(list(pred=pred[,1], var.pred=pred[,2]))
}


variogram <- function(h, nugget=0, sigma2, phi, covmodel="matern", ...){
  Ch <- geoR::cov.spatial(h, cov.model=covmodel, cov.pars=c(sigma2, phi),...)
  C0 = sigma2+nugget
  vh <- C0-Ch
  return(vh)
}
