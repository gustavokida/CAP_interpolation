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

cap_P_IDW <- inverseDistanceWeighted(cap_P_IDW)
cap_P_IDW <- removeOutlier(cap_P)
cap_P_IDW$handle_anisotropy <- TRUE
cap_P_IDW <- eidw(cap_P_IDW)
cap_P_IDW <- rmse(cap_P_IDW)

spplot(cap_P_IDW$newdata, zcol = "eidw_P", main = list(font = 1, cex = 1, label = "IDQ"), sub = list(font = 1, cex = 1, label = "(c) Seguindo o modelo"))

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







data_aniso_test <- removeOutlier(data=dados, column_names="P")
new_data_aniso_test <- highres.grid
p_aniso_test <- checkAnisotropy(data_aniso_test, P~1)
p_aniso_test <- p_aniso_test$direction
#data_aniso_test@coords <- rotateAnisotropicData(data_aniso_test, p_aniso_test)@coords
new_data_aniso_test$result <- eidw(data=data_aniso_test, newdata=new_data_aniso_test, formula= P~1, anisotropy=p_aniso_test)

data_aniso_test2 <-removeOutlier(data=dados, column_names="P")
new_data_aniso_test2 <- highres.grid
p_aniso_test2 <- c(p_aniso_test$direction, p_aniso_test$ratio)
# data_aniso_test2@coords <- coords.aniso(data_aniso_test2@coords, p_aniso_test2)
# colnames(data_aniso_test2@coords) <- c("x", "y")
# new_data_aniso_test2@coords <- coords.aniso(new_data_aniso_test2@coords, p_aniso_test2)
# colnames(new_data_aniso_test2@coords) <- c("x", "y")
new_data_aniso_test2$result <- idw_test(data = data_aniso_test2, newdata=new_data_aniso_test2, formula = P~1, handle_anisotropy = FALSE, handle_assimetry = FALSE, idp = 2, anisotropy=p_aniso_test2)
# new_data_aniso_test2 <- handleAnisotropy(data=new_data_aniso_test2, formula=P~1, anisotropy=p_aniso_test2, reverse=TRUE)
# data_aniso_teste2 <- handleAnisotropy(data=data_aniso_test2, formula=P~1, anisotropy=p_aniso_test2, reverse=TRUE)


spplot(new_data_aniso_test, zcol = "result", main = list(font = 1, cex = 1, label = "IDQ"), sub = list(font = 1, cex = 1, label = "teste1"))
spplot(new_data_aniso_test2, zcol = "result", main = list(font = 1, cex = 1, label = "IDQ"), sub = list(font = 1, cex = 1, label = "teste2"))

teste1_rmse <- rmse(data=data_aniso_test, formula=P~1, funcInterpolation=eidw)
teste2_rmse <- rmse(data=data_aniso_test2, formula=P~1, funcInterpolation=idw_test)

print(teste1_rmse)
print(teste2_rmse)

idw_test <- function(data = "SpatialPointsDataFrame", newdata = c("SpatialPointsDataFrame", "SpatialPixelsDataFrame"), formula = "formula", handle_anisotropy = FALSE,
                         handle_assimetry = FALSE, rmse_data = NULL, idp = 2, smooth = 0, n = 10, anisotropy = NULL){
  #anisotropy <- NULL
  if(is.null(anisotropy)){
    if(is.null(rmse_data)){
      anisotropy <- checkAnisotropy(data, formula)
      anisotropy <- c(anisotropy$direction, anisotropy$ratio)
    }
    else{
      anisotropy <- checkAnisotropy(rmse_data, formula)
      anisotropy <- c(anisotropy$direction, anisotropy$ratio)
    }
  }
  
  data@coords <- coords.aniso(data@coords, anisotropy)
  colnames(data@coords) <- c("x", "y")
  newdata@coords <- coords.aniso(newdata@coords, anisotropy)
  colnames(newdata@coords) <- c("x", "y")
  
  attribute <- formulaToVector(formula=formula, side="left")
  newdata_len <- length(newdata)
  data_len <- length(data)
  result <- vector(length=newdata_len)
  for (j in 1:newdata_len){
    weight <- 0
    sum_of_weights <- 0
    weighted_values_sum <- 0
    weight_df <- data.frame(weight=numeric(data_len), value=numeric(data_len))
    for(i in 1:data_len){
      delta_x = newdata@coords[j,1] - data@coords[i,1]
      delta_y = newdata@coords[j,2] - data@coords[i,2]
      weight_df$weight[i] <- sqrt(delta_x^2 + delta_y^2)^idp
      weight_df$value[i] <- data@data[i, attribute]
      #sum_of_weights <- sum_of_weights + weight
      #weighted_values_sum <- weighted_values_sum + weight*data@data[i, attribute]
    }
    sorted_df <- weight_df[order(weight_df$weight),]
    selected_n <- sorted_df[1:n,]
    sum_of_weights <- sum(1/selected_n$weight)
    weighted_values_sum <- sum(selected_n$value / selected_n$weight)
    point <- unname(weighted_values_sum / sum_of_weights)
    result[j] <- point
  }
  newdata$result <- result
  newdata <- handleAnisotropy(data=newdata, formula=result~1, anisotropy=anisotropy, reverse=TRUE)
  #data <- handleAnisotropy(data=data, formula=P~1, anisotropy=anisotropy, reverse=TRUE)
  
  return(result)
}





cap_P_no_outlier_aniso_KO_test <- ordinaryKriging(cap_P_no_outlier_aniso)
cap_P_no_outlier_aniso_KO_test <- rmse(cap_P_no_outlier_aniso_KO_test)

plot(cap_P_no_outlier_aniso_KO_test)

cap_P_no_outlier_aniso_KO_test2 <- ordinaryKriging_aniso_test(cap_P_no_outlier_aniso)
cap_P_no_outlier_aniso_KO_test2 <- rmse(cap_P_no_outlier_aniso_KO_test2)

plot(cap_P_no_outlier_aniso_KO)



plot(ordinaryKriging(cap_P_no_outlier_aniso$data, cap_P_no_outlier_aniso$newdata, cap_P_no_outlier_aniso$formula))
plot(ordinaryKriging_aniso_test.default(cap_P_no_outlier_aniso$data, cap_P_no_outlier_aniso$newdata, cap_P_no_outlier_aniso$formula))
autokrigetest <- autoKrige(formula = cap_P_no_outlier_aniso$formula, input_data = cap_P_no_outlier_aniso$data, new_data = cap_P_no_outlier_aniso$newdata)


#ordinary kriging
ordinaryKriging_aniso_test.default <- function(data = "SpatialPointsDataFrame", newdata = c("SpatialPointsDataFrame", "SpatialPixelsDataFrame"), formula = "formula",
                                    data_variogram = NULL, handle_anisotropy = TRUE, handle_assimetry = TRUE, rmse_data = NULL, cap_return = FALSE, idp = NULL){
  #check and remove assimetry
  if(isTRUE(handle_assimetry)){
    if(is.null(rmse_data)){
      normal_distribution <- normalDistribution(data = data, formula = formula)
      if(isFALSE(normal_distribution)){
        lambda <- boxCoxLambda(formula = formula, data = data)
        main_var <- formulaToVector(formula = formula, side = "left")
        data@data[, main_var] <- boxCoxTransform(formula = formula, data = data, lambda = lambda, reverseBoxCox = FALSE)
      }
    }
    else{
      normal_distribution <- normalDistribution(data = rmse_data, formula = formula)
      if(isFALSE(normal_distribution)){
        lambda <- boxCoxLambda(formula = formula, data = rmse_data)
        main_var <- formulaToVector(formula = formula, side = "left")
        data@data[, main_var] <- boxCoxTransform(formula = formula, data = data, lambda = lambda, reverseBoxCox = FALSE)
        rmse_data@data[, main_var] <- boxCoxTransform(formula = formula, data = rmse_data, lambda = lambda, reverseBoxCox = FALSE)
      }
    }
  }
  #check and remove anisotropy
  if(isTRUE(handle_anisotropy)){
    if(is.null(rmse_data)){
      anisotropy <- checkAnisotropy(data=data, formula=formula)
      data_variogram <- handleAnisotropy(data = data, formula = formula, anisotropy = anisotropy)
      data <- handleAnisotropy(data = data, formula = formula, anisotropy = anisotropy)
      newdata <- handleAnisotropy(data = newdata, formula = formula, anisotropy = anisotropy)
    }
    else{
      anisotropy <- checkAnisotropy(data=rmse_data, formula=formula)
      data_variogram <- handleAnisotropy(rmse_data, formula, anisotropy = anisotropy)
      data <- handleAnisotropy(data = data, formula = formula, anisotropy = anisotropy)
      newdata <- handleAnisotropy(data = newdata, formula = formula, anisotropy = anisotropy)
    }
  }
  else{
    if(is.null(rmse_data) && is.null(data_variogram)){
      data_variogram <- data
    }
    else if(!is.null(rmse_data) && is.null(data_variogram)){
      data_variogram <- rmse_data
    }
  }
  
  #do the kriging with the generated variogram
  interpolated_data <- autoKrige(formula, data, newdata, data_variogram = data_variogram)
  
  # #rotates coordinates back to anisotropic
  if(isTRUE(handle_anisotropy)){
    interpolated_data$krige_output <- handleAnisotropy(data = data$krige_output, formula = formula, anisotropy = anisotropy, reverse = TRUE)
    data <- handleAnisotropy(data = data, formula = formula, anisotropy = anisotropy, reverse = TRUE)
    newdata <- handleAnisotropy(data = newdata, formula = formula, anisotropy = anisotropy, reverse = TRUE)
  }
  
  #back-transforms data to assimetry
  if(isTRUE(handle_assimetry)){
    if(isFALSE(normal_distribution)){
      if(is.null(rmse_data)){
        data@data[, main_var] <- boxCoxTransform(formula = formula, data = data, lambda = lambda, reverseBoxCox = TRUE)
      }
      else{
        data@data[, main_var] <- boxCoxTransform(formula = formula, data = data, lambda = lambda, reverseBoxCox = TRUE)
        rmse_data@data[, main_var] <- boxCoxTransform(formula = formula, data = rmse_data, lambda = lambda, reverseBoxCox = TRUE)
      }
      interpolated_data$krige_output$var1.pred <- boxCoxTransform(formula = formula, data = interpolated_data, lambda = lambda, reverseBoxCox = TRUE)
    }
  }
  
  
  if(isTRUE(cap_return) && is.null(rmse_data)){
    return(interpolated_data)
  }
  else{
    return(interpolated_data$krige_output$var1.pred)
  }
}
