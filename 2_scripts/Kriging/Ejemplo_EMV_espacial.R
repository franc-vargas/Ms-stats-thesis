## librerias
 library(MASS) # Contiene la normal Multivariada
 library(geoR) # Contiene los modelos de covarianza

 
## Gráfica de los Modelos de variograma 
 
 variogram <- function(h, nugget=0, sigma2, phi, covmodel="matern", ...){
   Ch <- geoR::cov.spatial(h, cov.model=covmodel, cov.pars=c(sigma2, phi),...)
   C0 = sigma2+nugget
   vh <- C0-Ch
   return(vh)
 }
 
 h=seq(0,1, by=0.01)
 
 plot(h, variogram(h, nugget=0, sigma2=1, phi=0.1), type="l", ylim=c(0,1.5))
 lines(h, variogram(h, nugget=0.1, sigma2=1, phi=0.1), col=2)
 lines(h, variogram(h, nugget=0.15, sigma2=1, phi=0.1, kappa=1.5), col=3)
## Simulación de coordenadas
 
 n.sim=200 #cantidad de coordenadas 
 coords.sim=data.frame(x=round(runif(n.sim),4), y=round(runif(n.sim),4))
 plot(coords.sim, pch=20)
 hdist=dist(coords.sim) # distancias
 H=as.matrix(hdist)
  
## Parámetros para la simulación 
 
 mu.sim=10
 sigma2.sim=2
 phi.sim=0.05
 tau2.sim=0.5
 
 param.sim=data.frame(mu=mu.sim, sigma2=sigma2.sim, phi=phi.sim, nugget=tau2.sim)

## Vector de medias y matriz de covarianzas para la simulación  
 media.sim <- rep(mu.sim, n.sim)
 Sigma.sim <- cov.spatial(H, cov.model= "matern", kappa=1.5, 
                          cov.pars=c(sigma2.sim, phi.sim) ) + tau2.sim*diag(n.sim)
 
 media.sim[1:5]
 Sigma.sim[1:5,1:5]

## Simulación de campos aleatorios Gaussianos
 
 m.sim <- 10 # cantidad de replicas
 data.sim <- mvrnorm(n=m.sim, mu=media.sim, Sigma=Sigma.sim) # cada fila es una realización
 
 hist(data.sim[1,]) # un histograma
 hist(data.sim[2,]) # un histograma
 hist(data.sim[3,]) # un histograma
 
## Funcion de log verosimilitud
 
 loglik <- function(param, data, 
                    coords=NULL, H=NULL, 
                    covmodel= "matern", media="cte", ...){
   n=length(data)
   if(media=="cte"){
     media=matrix(rep(param[['mu']],n),nc=1) 
   }
   
   if(is.null(coords) & is.null(H)){
     stop("debe ingregar las coordenadas o la matriz de distancias")
   }
   
   if(is.null(H)){ 
     H=as.matrix(dist(coords)) 
   }

   Sigma <- geoR::cov.spatial(H, cov.model=covmodel, 
                            cov.pars=c(
                              param[['sigma2']], 
                              param[['phi']]),... ) + param[['nugget']]*diag(n)
   logdet = determinant(Sigma)$modulus
   
   fcuad = t(data-media) %*% solve(Sigma) %*% (data-media)
   
   lik= -(n*log(2*pi)/n+logdet/2+fcuad/2)
   return(lik)
 }
 
 # Parametros iniciales desde los datos
 param0=data.frame(mu=mean(apply(data.sim,1,mean)), 
                   sigma2=mean(apply(data.sim,1,var))*0.9, 
                   phi=seq(0.01, 0.15, by=0.001),
                   nugget=mean(apply(data.sim,1,var))*0.1)
 
 loglik(param=param.sim, data=data.sim[1,], H=H )
 loglik(param=param0[1,], data=data.sim[1,], H=H )
 
 
 par(mfrow=c(2,2))
 
 for(j in 1:4){
 plot(seq(0.01, 0.15, by=0.001), type="l",
      sapply(1:length(seq(0.01, 0.15, by=0.001)), 
             function(i){ loglik(param0[i,],
                                 data=data.sim[j,], H=H )}))
 points(param.sim$phi, loglik(param=param.sim, 
                              data=data.sim[j,], 
                              H=H ), 
        col=2, pch=20, lwd=2)
 }

## Estimación de parámetros por Máxima Verosimilitud 
 res.sim=list()
 for(j in 1:m.sim){
   res.sim[[j]] <- optim(par=param.sim, fn=loglik, 
                         data=data.sim[j,], H=H, 
                         control=list(fnscale=-1),
                         method = "L-BFGS-B",
                         lower=c(-Inf,0,0,0), upper=Inf)
   }
 # ejemplo estimadores
 res.sim[[m.sim]]$par
 
 param.fit <- c()
 for(j in 1:m.sim){
   param.fit=rbind(param.fit, res.sim[[j]]$par) }
 
par(mfrow = c(1,1))
 boxplot(param.fit/rep(unlist(param.sim), each=m.sim))
abline(h=1, lty=2)

##----------------------------------------------------------------------------##
## Ejemplo kriging

ll=101
s0=expand.grid(x=seq(0,1, l=ll), y=seq(0,1, l=ll))
plot(s0, pch=20)

# calculo de distancia de s0 a s
dist0 <- function(coords, s0){
  n=nrow(coords)
  m=nrow(s0)
  D=matrix(NA, nr=n, nc=m)
  for(i in 1:m){
    D[,i] <- dist(rbind(s0[i,], coords))[1:n] 
    } 
  return(D)
}

h0=dist0(coords.sim, s0)
dim(h0)
head(h0[1:5,1:5])
c0=cov.spatial(h0, cov.model="matern", 
              cov.pars=c(param.sim[['sigma2']], param.sim[['phi']]))
sum(h0==0)
if(h0==0){ c0[h0==0] <- param.sim[['sigma2']] + param.sim[['nugget']] }
dim(c0)
head(c0[1:5,1:5])

## Funcion para kriging simple
kr.simple <- function(media0, media, data, s0=NULL,
                      coords=NULL, InvSigma=NULL, 
                      Sigma=NULL, c0=NULL, 
                      H=NULL, h0=NULL,
                      covmodel="matern", param=NULL,...){
  # Se calcula Sigma en caso de omitirse
  if(is.null(Sigma)){
    if(is.null(H)){
      H=as.matrix(dist(coords))
      }
    n=nrow(H)
    Sigma <- cov.spatial(H, cov.model=covmodel, 
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
    c0 <- cov.spatial(h0, cov.model=covmodel, 
                         cov.pars=c(param[['sigma2']], 
                                    param[['phi']]),...)
    c0[h0==0] <- param[['sigma2']] + param[['nugget']]
  }
  
  # determinación de los pesos optimos
  lambda=InvSigma %*% c0
  
  # prediccion y varianza de predicción
  pred = media0 + t(lambda) %*% (data - media)
  var.pred <- 
    diag(t(lambda) %*% Sigma %*% lambda)-
    2 * diag( t(lambda) %*% c0 ) +
    param[['sigma2']] + param[['nugget']]
  
  return(list(pred=pred, var.pred=var.pred))
}

# se evalua en los verdaderos parámetros con la realización j

j=4
ks1 <- kr.simple(media0=rep(param.sim$mu, nrow(s0)), 
                 media=rep(param.sim$mu, length(data.sim[j,])), 
                data=data.sim[j,], Sigma=Sigma.sim, 
                h0=h0, param=param.sim)
par(mfrow=c(1,2))
image(unique(s0$x), unique(s0$y), matrix(ks1$pred, nr=ll, nc=ll))
points(coords.sim, pch=20, col="gray90")
image(unique(s0$x), unique(s0$y), matrix(ks1$var.pred, nr=ll, nc=ll))
points(coords.sim, pch=20, col="gray90")

filled.contour(matrix(ks1$pred, nr=ll, nc=ll), plot.axes = {}, main="Kriging")
filled.contour(matrix(ks1$var.pred, nr=ll, nc=ll), plot.axes = {}, main="Var Kriging")

hist(data.sim[j,])
hist(ks1$pred)


save.image(file = "current_session.RData")
load("current_session.RData")

