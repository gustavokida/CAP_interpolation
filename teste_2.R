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






#Testes para o artigo

#elemento P
#sem tratamento
cap_P <- cap(formula = P~1, data = dados, newdata = highres.grid)
cap_P$handle_anisotropy <- FALSE
cap_P$handle_assimetry <- FALSE

#com tratamento de outlier
cap_P_no_outlier <- removeOutlier(object = cap_P)

#com tratamento de outlier e anisotropia
cap_P_no_outlier_aniso <- cap_P_no_outlier
cap_P_no_outlier_aniso$handle_anisotropy <- TRUE

#com tratamento de anisotropia e assimetria
cap_P_no_aniso_assimetry <- cap_P
cap_P_no_aniso_assimetry$handle_anisotropy <- TRUE
cap_P_no_aniso_assimetry$handle_assimetry <- TRUE



#IDW
cap_P_IDW <- inverseDistanceWeighted(cap_P)
cap_P_IDW <- rmse(cap_P_IDW)

cap_P_no_outlier_IDW <- inverseDistanceWeighted(cap_P_no_outlier)
cap_P_no_outlier_IDW <- rmse(cap_P_no_outlier_IDW)

cap_P_no_outlier_aniso_IDW <- inverseDistanceWeighted(cap_P_no_outlier_aniso)
cap_P_no_outlier_aniso_IDW <- rmse(cap_P_no_outlier_aniso_IDW)

cap_P_no_aniso_assimetry_IDW <- inverseDistanceWeighted(cap_P_no_aniso_assimetry)
cap_P_no_aniso_assimetry_IDW <- rmse(cap_P_no_aniso_assimetry_IDW)

#KO
cap_P_KO <- ordinaryKriging(cap_P)
cap_P_KO <- rmse(cap_P_KO)

cap_P_no_outlier_KO <- ordinaryKriging(cap_P_no_outlier)
cap_P_no_outlier_KO <- rmse(cap_P_no_outlier_KO)

cap_P_no_outlier_aniso_KO <- ordinaryKriging(cap_P_no_outlier_aniso)
cap_P_no_outlier_aniso_KO <- rmse(cap_P_no_outlier_aniso_KO)

cap_P_no_aniso_assimetry_KO <- ordinaryKriging(cap_P_no_aniso_assimetry)
cap_P_no_aniso_assimetry_KO <- rmse(cap_P_no_aniso_assimetry_KO)

cap_P_no_aniso_assimetry_UK <- universalKriging(cap_P_no_aniso_assimetry)
cap_P_no_aniso_assimetry_UK <- rmse(cap_P_no_aniso_assimetry_UK)

#plot
plot(cap_P_IDW)
plot(cap_P_no_outlier_IDW)
plot(cap_P_no_outlier_aniso_IDW)
plot(cap_P_no_aniso_assimetry_IDW)

plot(cap_P_KO)
plot(cap_P_no_outlier_KO)
plot(cap_P_no_outlier_aniso_KO)
plot(cap_P_no_aniso_assimetry_KO)


spplot(cap_P_IDW$newdata, zcol = "IDW_P", main = list(font = 1, cex = 1, label = "IDQ"), sub = list(font = 1, cex = 1, label = "(b) Sem tratamentos"))
spplot(cap_P_no_outlier_IDW$newdata, zcol = "IDW_P", main = list(font = 1, cex = 1, label = "IDQ"), sub = list(font = 1, cex = 1, label = "(d) Com remoção de outliers"))
spplot(cap_P_no_aniso_assimetry_IDW$newdata, zcol = "IDW_P", main = list(font = 1, cex = 1, label = "IDQ"), sub = list(font = 1, cex = 1, label = "(e) Com remoção de anisotropia"))
spplot(cap_P_no_outlier_aniso_IDW$newdata, zcol = "IDW_P", main = list(font = 1, cex = 1, label = "IDQ"), sub = list(font = 1, cex = 1, label = "(g) Com rem. de outliers/anisotropia"))

spplot(cap_P_no_outlier_aniso_IDW$newdata, zcol = "IDW_P", main = list(font = 1, cex = 1, label = "IDQ"), sub = list(font = 1, cex = 1, label = "(c) Seguindo o modelo"))

spplot(cap_P_KO$newdata, zcol = "OK_P", main = list(font = 1, cex = 1, label = "KO"), sub = list(font = 1, cex = 1, label = "(a) Sem tratamentos"))
spplot(cap_P_no_outlier_KO$newdata, zcol = "OK_P", main = list(font = 1, cex = 1, label = "KO"), sub = list(font = 1, cex = 1, label = "(c) Com remoção de outliers"))
spplot(cap_P_no_outlier_aniso_KO$newdata, zcol = "OK_P", main = list(font = 1, cex = 1, label = "KO"), sub = list(font = 1, cex = 1, label = "(f) Com rem. de outliers/anisotropia"))




#elemento HAL
cap_HAL <- cap(formula = HAL~1, data = dados, newdata = highres.grid)
cap_HAL$handle_anisotropy <- FALSE
cap_HAL$handle_assimetry <- FALSE

#com tratamento de outlier
cap_HAL_no_outlier <- removeOutlier(object = cap_HAL)

#com tratamento de outlier e anisotropia
cap_HAL_no_outlier_aniso <- cap_HAL_no_outlier
cap_HAL_no_outlier_aniso$handle_anisotropy <- TRUE

#com tratamento de anisotropia
cap_HAL_no_aniso <- cap_HAL
cap_HAL_no_aniso$handle_anisotropy <- TRUE

#com tratamento de anisotropia e assimetria
cap_HAL_no_aniso_assimetry <- cap_HAL
cap_HAL_no_aniso_assimetry$handle_anisotropy <- TRUE
cap_HAL_no_aniso_assimetry$handle_assimetry <- TRUE

#IDW
cap_HAL_IDW <- inverseDistanceWeighted(cap_HAL)
cap_HAL_IDW <- rmse(cap_HAL_IDW)

cap_HAL_no_outlier_IDW <- inverseDistanceWeighted(cap_HAL_no_outlier)
cap_HAL_no_outlier_IDW <- rmse(cap_HAL_no_outlier_IDW)

cap_HAL_no_outlier_aniso_IDW <- inverseDistanceWeighted(cap_HAL_no_outlier_aniso)
cap_HAL_no_outlier_aniso_IDW <- rmse(cap_HAL_no_outlier_aniso_IDW)

cap_HAL_no_aniso_IDW <- inverseDistanceWeighted(cap_HAL_no_aniso)
cap_HAL_no_aniso_IDW <- rmse(cap_HAL_no_aniso_IDW)

cap_HAL_no_aniso_assimetry_IDW <- inverseDistanceWeighted(cap_HAL_no_aniso_assimetry)
cap_HAL_no_aniso_assimetry_IDW <- rmse(cap_HAL_no_aniso_assimetry_IDW)


#KO
cap_HAL_KO <- ordinaryKriging(cap_HAL)
cap_HAL_KO <- rmse(cap_HAL_KO)

cap_HAL_no_outlier_KO <- ordinaryKriging(cap_HAL_no_outlier)
cap_HAL_no_outlier_KO <- rmse(cap_HAL_no_outlier_KO)

cap_HAL_no_outlier_aniso_KO <- ordinaryKriging(cap_HAL_no_outlier_aniso)
cap_HAL_no_outlier_aniso_KO <- rmse(cap_HAL_no_outlier_aniso_KO)

cap_HAL_no_aniso_KO <- ordinaryKriging(cap_HAL_no_aniso)
cap_HAL_no_aniso_KO <- rmse(cap_HAL_no_aniso_KO)

cap_HAL_no_aniso_assimetry_KO <- ordinaryKriging(cap_HAL_no_aniso_assimetry)
cap_HAL_no_aniso_assimetry_KO <- rmse(cap_HAL_no_aniso_assimetry_KO)

#plot
plot(cap_HAL_IDW)
plot(cap_HAL_no_outlier_IDW)
plot(cap_HAL_no_outlier_aniso_IDW)
plot(cap_HAL_no_aniso_IDW)
plot(cap_HAL_no_aniso_assimetry_IDW)

plot(cap_HAL_KO)
plot(cap_HAL_no_outlier_KO)
plot(cap_HAL_no_outlier_aniso_KO)
plot(cap_HAL_no_aniso_KO)
plot(cap_HAL_no_aniso_assimetry_KO)

spplot(cap_HAL_IDW$newdata, zcol = "IDW_HAL", main = list(font = 1, cex = 1, label = "IDQ"), sub = list(font = 1, cex = 1, label = "(b) Sem tratamentos"))
spplot(cap_HAL_no_outlier_IDW$newdata, zcol = "IDW_HAL", main = list(font = 1, cex = 1, label = "IDQ"), sub = list(font = 1, cex = 1, label = "(d) Com remoção de outliers"))
spplot(cap_HAL_no_aniso_IDW$newdata, zcol = "IDW_HAL", main = list(font = 1, cex = 1, label = "IDQ"), sub = list(font = 1, cex = 1, label = "(f) Com remoção de anisotropia"))
spplot(cap_HAL_no_aniso_assimetry_IDW$newdata, zcol = "IDW_HAL", main = list(font = 1, cex = 1, label = "IDQ"), sub = list(font = 1, cex = 1, label = "(h) Com remoção de assimetria"))

spplot(cap_HAL_KO$newdata, zcol = "OK_HAL", main = list(font = 1, cex = 1, label = "KO"), sub = list(font = 1, cex = 1, label = "(a) Sem tratamentos"))
#spplot(cap_HAL_no_outlier_KO$newdata, zcol = "OK_HAL", main = list(font = 1, cex = 1, label = "KO"), sub = list(font = 1, cex = 1, label = "(c) Com remoção de outliers"))

spplot(cap_HAL_no_outlier_KO$newdata, zcol = "OK_HAL", main = list(font = 1, cex = 1, label = "KO"), sub = list(font = 1, cex = 1, label = "(c) Seguindo o modelo"))
plot(cap_HAL_no_outlier_KO$variogram$exp_var, cap_HAL_no_outlier_KO$variogram$var_model, sub = list(font = 1, cex = 1, label = "(b) Seguindo o modelo"))

spplot(cap_HAL_no_aniso_KO$newdata, zcol = "OK_HAL", main = list(font = 1, cex = 1, label = "KO"), sub = list(font = 1, cex = 1, label = "(e) Com remoção de anisotropia"))
spplot(cap_HAL_no_aniso_assimetry_KO$newdata, zcol = "OK_HAL", main = list(font = 1, cex = 1, label = "KO"), sub = list(font = 1, cex = 1, label = "(g) Com remoção de assimetria"))

plot(cap_HAL_no_aniso_assimetry_KO$variogram$exp_var, cap_HAL_no_aniso_assimetry_KO$variogram$var_model, sub = list(font = 1, cex = 1, label = "(d) Com remoção de assimetria"))






#teste artigo_2
cap_P <- cap(formula = P~1, data = dados, newdata = highres.grid)
cap_HAL <- cap(formula = HAL~1, data = dados, newdata = highres.grid)

cap_P_auto <- autoInterpolation(cap_P)
cap_P_auto <- rmse(cap_P_auto)

cap_P_IDW <- inverseDistanceWeighted(cap_P)
cap_P_IDW <- rmse(cap_P_IDW)

cap_P_OK <- ordinaryKriging(cap_P)
cap_P_OK <- rmse(cap_P_OK)



cap_HAL_auto <- autoInterpolation(cap_HAL)
cap_HAL_auto <- rmse(cap_HAL_auto)

cap_HAL_IDW <- inverseDistanceWeighted(cap_HAL)
cap_HAL_IDW <- rmse(cap_HAL_IDW)

cap_HAL_OK <- ordinaryKriging(cap_HAL)
cap_HAL_OK <- rmse(cap_HAL_OK)


plot(cap_P_auto)
plot(cap_P_IDW)
plot(cap_P_OK)

plot(cap_HAL_auto)
plot(cap_HAL_IDW)
plot(cap_HAL_OK)





cap_P_EIDW <- eidw(cap_P_no_outlier_aniso)
cap_P_EIDW <- rmse(cap_P_EIDW)
plot(cap_P_EIDW)

