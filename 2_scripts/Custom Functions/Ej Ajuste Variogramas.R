test <- fit_WLS_variogram(vg_custom_scaled,
                          robust =T,
                          cov.model = "matern")

test$par
ggplot(vg_custom_scaled$plot.data) +
  geom_point(aes(x = lag, y = robust)) +
  # geom_label(aes(x = lag, y = semivariog, label = n))+
  labs(x = "Distancia (km)", y = "Semivarianza",
       caption = expression(paste("Modelo Matérn Least-Cost",": ",sigma^2 ,"=",394, ", ", Phi,"=",20, ", ",
                                  nugget, "=", 106 , ", ",kappa==0.5) )) + 
  geom_function(fun = variogram, 
                args = list(nugget=test$par[1], sigma2=test$par[2], phi=test$par[4],
                            covmodel = "matern",kappa = test$par[3]),
                colour = "darkred",
                xlim = c(0,90)) +
  ylim(c(0, max(vg_custom_scaled$plot.data$semivariog)))

vg_custom_scaled$data |> glimpse()
colnames(vg_custom_scaled$data)[colnames(vg_custom_scaled$data) == "distance"] <- "lag"
test2 <- fit_WLS_variogram(vg_custom_scaled,
                           robust = T,
                           cov.model = "matern")
test2$par
ggplot(vg_custom_scaled$plot.data) +
  geom_point(aes(x = lag, y = robust)) +
  geom_label(aes(x = lag, y = robust, label = n))+
  labs(x = "Distancia (km)", y = "Semivarianza",
       caption = expression(paste("Modelo Matérn Least-Cost",": ",sigma^2 ,"=",394, ", ", Phi,"=",20, ", ",
                                  nugget, "=", 106 , ", ",kappa==0.5) )) + 
  geom_function(fun = variogram, 
                args = list(nugget=test2$par[["nugget"]], sigma2=test2$par[["sigma"]], phi=test2$par[["phi"]],
                            covmodel = "matern",kappa = test$par[["kappa"]]),
                colour = "darkred") 
