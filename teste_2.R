rkTest <- function(data = "SpatialPointsDataFrame", newdata = c("SpatialPointsDataFrame", "SpatialPixelsDataFrame"), formula = "formula"){
  fitted_data <- fit.gstatModel(data, formula, newdata, fit.family = gaussian())    #Must use covariates in SpatialPixelsDataFrame type
  interpolated_data <- predict.gstatModel(fitted_data, newdata)
  return(interpolated_data)
}

#Limpa a variavel de resultados RMSE
rmse_result <- NULL
quant_outliers <- NULL
assimetria_sem_outlier <- NULL
anisotropia_sem_outlier <- NULL
dados <- NULL
dados.grid <- NULL
teste.grid <- NULL
produtividade.grid <- NULL
produtividade <- NULL

#dados
dados <- nutrientes.dados
teste2 <- removeNA(removeOutlier(dados, variavel), variavel)
produtividade <- produtividade_trigo_2009utm.rg

#grid
dados.grid <- nutrientes.grid
produtividade.grid <- makeGrid(contornoutm.rg, cellsize = 10)
teste.grid <- produtividade.grid


variavel <- "HAL"
covariavel <- "DIST_NCIA_"
covariaveis_names <- names(produtividade)
covariaveis_names <- c("DIST_NCIA_", "SENTIDO_DE", "DURA__O_HR", "ALTITUDE_M",
                      "LARG__CARR", "FLUXO_DA_C", "UMIDADE___", "PRODU__O_D", "VOL__DE_PR", "PROD_DE_MA", "VOL__DE_P2", "N__DA_VOLT", "VELOCIDADE", "PROD__HA_H")
#covariaveis_names <- c('PH','HAL', 'AL', 'CA', 'MG', 'K', 'SB', 'CTC', 'P', 'C', 'MO', 'V', 'ALSATURACA')
#covariaveis <- covariaveis_names[which(covariaveis_names != variavel)]


#tentar interpolar grande numero de pontos para menor numero de pontos, ou coletar os pontos exatos, ou aumentar o numero de pontos do talhao
#filtrar
produtividade_limpo <- removeOutlier(produtividade, covariaveis_names)
for(i in covariaveis_names){
  produtividade_limpo <- removeNA(produtividade_limpo, i)
}
length(produtividade_limpo)
for(i in covariaveis_names){
  formula <- makeFormula(i)
  produtividade.grid@data[i] <- inverseDistanceWeighted(produtividade_limpo, produtividade.grid, formula)
  spplot(produtividade.grid, zcol = i)
}
spplot(produtividade.grid, zcol = "SENTIDO_DE")


#verifica covariaveis para cada atributo
for(i in c("PH", "HAL", "AL", "CA", "MG", "K", "SB", "CTC", "P", "C", "MO", "V", "ALSATURACA")){
  print(i)
  print(checkCovariates(teste2, produtividade.grid, main_attribute_column = i, column_names = covariaveis_names))
}
nutrientes_col <- c("HAL", "AL", "CA", "MG", "K", "SB", "CTC", "P", "C", "MO", "V", "ALSATURACA")
for(i in c("PH", "HAL", "AL", "CA", "MG", "K", "SB", "CTC", "P", "C", "MO", "V", "ALSATURACA")){
  print(i)
  print(checkCovariates(teste2, main_attribute_column = i, column_names = nutrientes_col[which(nutrientes_col != i)]))
}

#Faz formula
#covariates <- checkCovariatesCKO(data = teste2, covariate_data = produtividade, main_attribute_column = variavel, column_names = covariaveis_names)
covariates <- checkCovariatesCKO(teste2, main_attribute_column = variavel, column_names = nutrientes_col[which(nutrientes_col != variavel)])
print(covariates)
formula <- makeFormula(variavel, covariates)

#COKRIGING
teste.grid$cokriging <- coKriging(teste2, produtividade.grid, formula, teste2)
teste.grid$cokriging <- coKriging(teste2, produtividade.grid, HAL~ALTITUDE_M, produtividade.grid)
teste.grid$cokriging <- coKriging(teste2, produtividade.grid, HAL~V, teste2)
spplot(teste.grid, zcol = "cokriging")
print(rmse(teste2, HAL~ALTITUDE_M, coKriging, produtividade.grid))
print(rmse(teste2, HAL~V, coKriging, teste2))

#OK
teste.grid$ok <- ordinaryKriging(teste2, produtividade.grid, HAL~1)
spplot(teste.grid, zcol = "ok")
print(rmse(teste2, HAL~1, ordinaryKriging))


#gwr
#teste.grid$gwr <- gwrKriging(teste2, produtividade.grid, CA~PH+C)
#spplot(teste.grid, zcol = "gwr")
# rk_test <- rkTest(teste2, produtividade.grid, CA~SENTIDO_DE+ALTITUDE_M)
# plot(rk_test)
# teste.grid$rk <- regressionKriging(teste2, produtividade.grid, CA~SENTIDO_DE+ALTITUDE_M)
# spplot(teste.grid, zcol = "rk")
#teste.grid$gwrk <- 


