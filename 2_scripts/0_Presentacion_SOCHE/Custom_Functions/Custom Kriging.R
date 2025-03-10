krig.manual <- function(media0, media, data, s0=NULL,
                      coords=NULL, InvSigma=NULL, 
                      Sigma=NULL, c0=NULL, 
                      H=NULL, h0=NULL,
                      covmodel="matern", param=NULL,
                      type = "ord",...){
  # Se calcula Sigma en caso de omitirse
  if(is.null(Sigma)){
    if(is.null(H)){
      H=as.matrix(dist(coords))
    }
    n=nrow(H)
    Sigma <- geoR::cov.spatial(H, cov.model=covmodel, 
                         cov.pars=c(
                           param[['sigma2']], 
                           param[['phi']]),... ) + param[['nugget']]*diag(n) 
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
    data_aux <- data
    names(data_aux) <- rownames(C)
    loo <- lapply(colnames(C), function(x){
      
      x0 <- drop(C[rownames(C) != x, colnames(C) == x])
      C_aux <- C[rownames(C) != x, colnames(C) != x]
      aux_c0 <- sigma2 + nugget
      aux_sigma0 <- sigma2 + nugget
      aux_vec <- rep(1, nrow(C_aux))
      loo_data <- data_aux[names(data_aux) != x]
      lambda_OK <- drop(C_aux %*% (x0 + aux_vec %*% 
                                 ( solve(t(aux_vec) %*% C_aux %*% aux_vec) )%*%
                                 (1 - t(aux_vec) %*% C_aux %*% x0)))
      
      x <- drop(lambda_OK  %*% loo_data)
      x
      
    })
    
    loo <- unlist(loo)
  }
  
  return(list(pred=pred[,1], var.pred=pred[,2], fitted = loo))
}


variogram <- function(h, nugget=0, sigma2, phi, covmodel="matern", ...){
  Ch <- geoR::cov.spatial(h, cov.model=covmodel, cov.pars=c(sigma2, phi),...)
  C0 = sigma2+nugget
  vh <- C0-Ch
  return(vh)
}
