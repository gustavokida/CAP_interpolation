library(sp)
library(openxlsx)

ordinaryKrigingModified <- function(data = "SpatialPointsDataFrame", newdata = c("SpatialPointsDataFrame", "SpatialPixelsDataFrame"), formula = "formula"){
  interpolated_data <- autoKrige(formula, data, newdata)
  return(interpolated_data$krige_output$var1.pred)
}



logNormalOrdinaryKriging <- function(data = "SpatialPointsDataFrame", newdata = c("SpatialPointsDataFrame", "SpatialPixelsDataFrame"), formula = "formula"){
  left_side <- formulaToVector(formula, 'left')
  #data_normal <- logTransformation(data, formula))
  data_normal <- data
  data_normal@data[, left_side] <- boxCoxTransformation(formula, data)
  data_variogram <- handleAnisotropy(data_normal, formula)
  interpolated_data <- autoKrige(formula, data_normal, newdata, data_variogram = data_variogram)
  #backtransformed_data <- backTransformation(interpolated_data = interpolated_data)
  backtransformed_data <- interpolated_data
  backtransformed_data$krige_output$var1.pred <- boxCoxTransformation(formula, interpolated_data, BoxCox.lambda(data@data[, left_side], method = 'loglik'))
  #return(backtransformed_data)
  return(backtransformed_data$krige_output$var1.pred)
}



ordinaryKriging <- function(data = "SpatialPointsDataFrame", newdata = c("SpatialPointsDataFrame", "SpatialPixelsDataFrame"), formula = "formula"){
  data_variogram <- handleAnisotropy(data, formula)
  interpolated_data <- autoKrige(formula, data, newdata, data_variogram = data_variogram)
  return(interpolated_data$krige_output$var1.pred)
}



#try glm
ordinaryKrigingDetrended <- function(data = "SpatialPointsDataFrame", newdata = c("SpatialPointsDataFrame", "SpatialPixelsDataFrame"), formula = "formula"){
  modelo_lm <- NULL
  if(checkTrend(data, formula)){
    modelo_lm <- detrend(data, formula)
    name_residual <- paste(formulaToVector(formula, "left"), "residual", sep="_")
    data@data[name_residual] <- modelo_lm$residuals
    predicted_lm <- predict(modelo_lm, newdata)
    new_formula = as.formula(expr(!!as.symbol(name_residual)~1))
    interpolated_data <- autoKrige(new_formula, data, newdata)
    interpolated_data$krige_output$var1.pred <- interpolated_data$krige_output$var1.pred + predicted_lm
  }
  else{
    interpolated_data <- autoKrige(formula, data, newdata)
  }
  return(interpolated_data$krige_output$var1.pred)
}



plotOrdinaryKrigingModified <- function(data = "SpatialPointsDataFrame", newdata = c("SpatialPointsDataFrame", "SpatialPixelsDataFrame"), formula = "formula"){
  interpolated_data <- autoKrige(formula, data, newdata)
  return(interpolated_data)#plot(interpolated_data$exp_var, interpolated_data$var_model))
}



plotLogNormalOrdinaryKriging <- function(data = "SpatialPointsDataFrame", newdata = c("SpatialPointsDataFrame", "SpatialPixelsDataFrame"), formula = "formula"){
  left_side <- formulaToVector(formula, 'left')
  data_normal <- data
  #data_normal <- logTransformation(data, formula)
  data_normal@data[, left_side] <- boxCoxTransformation(formula, data)
  data_variogram <- handleAnisotropy(data_normal, formula)
  interpolated_data <- autoKrige(formula, data_normal, newdata, data_variogram = data_variogram)
  #backtransformed_data <- backTransformation(interpolated_data = interpolated_data)
  backtransformed_data <- interpolated_data
  backtransformed_data$krige_output$var1.pred <- boxCoxTransformation(formula, interpolated_data, BoxCox.lambda(data@data[, left_side], method = 'loglik'))
  return(backtransformed_data)#plot(interpolated_data$exp_var, interpolated_data$var_model))
}



plotOrdinaryKriging <- function(data = "SpatialPointsDataFrame", newdata = c("SpatialPointsDataFrame", "SpatialPixelsDataFrame"), formula = "formula"){
  data_variogram <- handleAnisotropy(data, formula)
  interpolated_data <- autoKrige(formula, data, newdata, data_variogram = data_variogram)
  return(interpolated_data)#plot(interpolated_data$exp_var, interpolated_data$var_model))
}



plotOrdinaryKrigingDetrended <- function(data = "SpatialPointsDataFrame", newdata = c("SpatialPointsDataFrame", "SpatialPixelsDataFrame"), formula = "formula"){
  modelo_lm <- NULL
  if(checkTrend(data, formula)){
    modelo_lm <- detrend(data, formula)
    name_residual <- paste(formulaToVector(formula, "left"), "residual", sep="_")
    data@data[name_residual] <- modelo_lm$residuals
    predicted_lm <- predict(modelo_lm, newdata)
    new_formula <- as.formula(expr(!!as.symbol(name_residual)~1))
    interpolated_data <- autoKrige(new_formula, data, newdata)
    interpolated_data$krige_output$var1.pred <- interpolated_data$krige_output$var1.pred + predicted_lm
  }
  else{
    interpolated_data <- autoKrige(formula, data, newdata)
  }
  return(interpolated_data)
}



plotUniversalKriging <- function(data = "SpatialPointsDataFrame", newdata = c("SpatialPointsDataFrame", "SpatialPixelsDataFrame"), formula = "formula"){
  new_formula <- detrendFormula(data, formula)
  interpolated_data <- autoKrige(new_formula, data, newdata)
  return(interpolated_data)
}



setwd("C:/Users/Gustavo/OneDrive/UFMT/CAP/2020-2021/Projeto/Algoritmos_R")
resistencia.rg <- read.xlsx("Dados de resis_limpo.xlsx")

#Nutrientes.rg arrumado, transformado na vari?vel dados
e1 <- data.frame(resistencia.rg$Norte)
e2 <- data.frame(resistencia.rg$Sul)
e3 <- cbind(e1,e2)
names(e3)[1] <- "x"
names(e3)[2] <- "y"
nutrientescoords <- e3
nutrientescrs <- nutrientesutm.rg@proj4string
nutrientesdata <- resistencia.rg[4:44]
resistencia_dados <- SpatialPointsDataFrame(coords = nutrientescoords, data = nutrientesdata)#, proj4string = nutrientescrs)
resistencia_dados_sem_outlier <- removeNA(removeOutlier(resistencia_dados, "V1"), "V1")

#spsample_test <- spsample(resistencia_dados, type = "regular")
#gridded(spsample_test) <- TRUE

#x_ori <- round(coordinates(resistencia_dados)[1, 1]/100) * 100
#y_ori <- round(coordinates(resistencia_dados)[1, 2]/100) * 100
x_cell <- 50
y_cell <- 75
cell_size <- 5
#resistencia_extent <- extent(x_ori, x_ori + (x_cell * cell_size), y_ori, y_ori + (y_cell * cell_size))
resistencia_extent <- extent(resistencia_dados@bbox)
#resistencia_extent <- extent(resistencia_dados)
ras <- raster(resistencia_extent)
res(ras) <- c(cell_size, cell_size)
ras[] <- 0
#projection(ras) <- nutrientescrs
resistencia.grid <- rasterToPoints(ras, spatial = TRUE)
gridded(resistencia.grid) <- TRUE
plot(resistencia.grid)
#resistencia.grid <- as(ras, "SpatialPixels")




#Variavel a utilizar  PH, HAL, AL, CA, MG, K, SB, CTC, P, C, MO, V, ALSATURACA
#variaveis com tendencia: CA, SB
nutrientes_attributes <- names(nutrientes.dados)

#variavel do dataset resistencia. atributos bons: R5.3
#mais promissores:R3, V1.1, V1.2, Ds.1, V2.1, Ds.2, V2.2, V4.2, R2.2, R3.2, R5.3.1, R5.3.2, R7.2.1, R7.2.2, Média.1, Média.2, Ds.3, V1.3, V4.3, R1.3, R2.3, R5.3.3, Média.3
resistencia_attributes <- names(resistencia_dados)

resistencia_attributes_promissores <- resistencia_dados@data[c("R3", "Ds.1", "V1.1", "V1.2", "V2.1", "Ds.2", "V2.2", "V4.2", "R2.2", "R3.2", "R5.3.1",
                                                               "R5.3.2", "R7.2.1", "R7.2.2", "Média.1", "Média.2", "Ds.3", "V1.3", "R1.3", "R2.3", "R5.3.3", "Média.3")]

resistencia_cor <- cor(resistencia_attributes_promissores, method = "pearson")

#SKEW: HAL, AL, ALSATURACA


#Limpa a variavel de resultados RMSE
rmse_result <- NULL
quant_outliers <- NULL
assimetria_sem_outlier <- NULL
anisotropia_sem_outlier <- NULL
dados <- NULL
dados.grid <- NULL
teste.grid <- NULL


variavel <- "HAL"
#nutrientes.dados, resistencia_dados
dados <- nutrientes.dados
#nutrientes.grid, resistencia.grid
dados.grid <- nutrientes.grid
teste.grid <- nutrientes.grid


#Faz formula
formula = makeFormula(variavel)
formula_xy = makeFormula(variavel, c("x", "y"))
formula_xy_2 = CA~poly(x,y,degree=2)
teste2 <- removeNA(removeOutlier(dados, variavel), variavel)


#quantidade outliers
quant_outliers <- length(dados) - length(teste2)


#Valores de assimetria
assimetria <- skew(dados@data[, variavel])
if(length(dados) != length(teste2)){
  assimetria_sem_outlier <- skew(teste2@data[, variavel])
}

#Valores de anisotropia
anisotropia <- checkAnisotropy(dados, formula)

if(length(dados) != length(teste2)){
  anisotropia_sem_outlier <- checkAnisotropy(teste2, formula)
}


#Interpolação OK
plot_ok <- plotOrdinaryKrigingModified(dados, dados.grid, formula)
plot(plot_ok$exp_var, plot_ok$var_model, sub = list(font = 1, cex = 1, label = "(a) Sem tratamentos"))
teste.grid$ok <- ordinaryKrigingModified(dados, dados.grid, formula)
spplot(teste.grid, zcol = "ok", main = list(font = 1, cex = 1, label = "KO"),  sub = list(font = 1, cex = 1, label = "(a) Sem tratamentos"))


#interpolação ok sem outliers
if(length(dados) != length(teste2)){
  plot_ok_no_outlier <- plotOrdinaryKrigingModified(teste2, dados.grid, formula)
  plot1 <- plot(plot_ok_no_outlier$exp_var, plot_ok_no_outlier$var_model, sub = list(font = 1, cex = 1, label = "(b) Com remoção de outliers"))
  teste.grid$ok_no_outlier <- ordinaryKrigingModified(teste2, dados.grid, formula)
  plot2 <- spplot(teste.grid, zcol = "ok_no_outlier", main = list(font = 1, cex = 1, label = "KO"),  sub = list(font = 1, cex = 1, label = "(c) Com remoção de outliers"))
  print(plot1)
  print(plot2)
}

#UNIVERSAL KRIGING ARRUMAR O DETREND
if(length(dados) != length(teste2)){
  plot_uk_no_outlier <- plotUniversalKriging(teste2, dados.grid, formula = formula)
  plot1 <- plot(plot_uk_no_outlier$exp_var, plot_uk_no_outlier$var_model, sub = list(font = 1, cex = 1, label = "(b) Com remoção de outliers"))
  teste.grid$uk_no_outlier <- universalKriging(teste2, dados.grid, formula)
  plot2 <- spplot(teste.grid, zcol = "uk_no_outlier", main = list(font = 1, cex = 1, label = "KO"),  sub = list(font = 1, cex = 1, label = "(c) Com remoção de outliers"))
  print(plot1)
  print(plot2)
}


#OK DETRENDED
if(length(dados) != length(teste2)){
  #Universal Kriging    ARRUMAR O RMSE PARA UK
  #plot(autoKrige(CA~poly(x, y, degree = 2), teste2, dados.grid))
  plot_ok_detrended_no_outlier <- plotOrdinaryKrigingDetrended(teste2, dados.grid, formula = formula)
  plot1 <- plot(plot_ok_detrended_no_outlier$exp_var, plot_ok_detrended_no_outlier$var_model, sub = list(font = 1, cex = 1, label = "(b) Com remoção de outliers"))
  teste.grid$ok_detrended_no_outlier <- ordinaryKrigingDetrended(teste2, dados.grid, formula)
  plot2 <- spplot(teste.grid, zcol = "ok_detrended_no_outlier", main = list(font = 1, cex = 1, label = "KO"),  sub = list(font = 1, cex = 1, label = "(c) Com remoção de outliers"))
  print(plot1)
  print(plot2)
}

#Para atributo com anisotropia e com outliers
if(!is.null(anisotropia)){
  plot_ok_isotropic <- plotOrdinaryKriging(dados, dados.grid, formula)
  plot1 <- plot(plot_ok_isotropic$exp_var, plot_ok_isotropic$var_model, sub = list(font = 1, cex = 1, label = "(c) Com remoção de anisotropia"))
  teste.grid$ok_isotropic <- ordinaryKriging(dados, dados.grid, formula)
  plot2 <- spplot(teste.grid, zcol = "ok_isotropic", main = list(font = 1, cex = 1, label = "KO"),  sub = list(font = 1, cex = 1, label = "(e) Com remoção de anisotropia"))
  print(plot1)
  print(plot2)
}

#Universal Kriging
#plot_ok_isotropic <- universalKriging(dados, dados.grid, formula_xy)
#plot1 <- plot(plot_ok_isotropic$exp_var, plot_ok_isotropic$var_model, sub = list(font = 1, cex = 1, label = "(c) Com remo??o de anisotropia"))


#Para atributo com anisotropia sem outliers
if(!is.null(anisotropia_sem_outlier)){
  plot_ok_isotropic_no_outlier <- plotOrdinaryKriging(teste2, dados.grid, formula)
  plot1 <- plot(plot_ok_isotropic_no_outlier$exp_var, plot_ok_isotropic_no_outlier$var_model, sub = list(font = 1, cex = 1, label = "(d) Com rem. de outliers/anisotropia"))
  teste.grid$ok_isotropic_no_outlier <- ordinaryKriging(teste2, dados.grid, formula)
  plot2 <- spplot(teste.grid, zcol = "ok_isotropic_no_outlier", main = list(font = 1, cex = 1, label = "KO"),  sub = list(font = 1, cex = 1, label = "(e) Com rem. de outliers/anisotropia"))
  print(plot1)
  print(plot2)
}


#Para atributos assim?tricos com outliers
if(!normalDistribution(dados, formula)){
  skew_dados <- dados
  #skew_normal_dist <- logTransformation(skew_dados, formula)
  #skew_assimetria_depois_transformacao <- skew(skew_normal_dist@data[, variavel])
  #anisotropia_skew_normal_dist <- checkAnisotropy(skew_normal_dist, formula)
  #teste.grid$ok_isotropic_no_outlier_normal_dist <- ordinaryKriging(teste2_normal_dist, dados.grid, formula, backtransform = TRUE)
  #plotLogNormalOrdinaryKriging(teste2, dados.grid, formula)
  plot_ok_isotropic_skew_normal_dist <- plotLogNormalOrdinaryKriging(skew_dados, dados.grid, formula)
  plot1 <- plot(plot_ok_isotropic_skew_normal_dist$exp_var, plot_ok_isotropic_skew_normal_dist$var_model, sub = list(font = 1, cex = 1, label = "(d) Com rem. de anisotropia/assimetria"))
  teste.grid$ok_isotropic_skew_normal_dist <- logNormalOrdinaryKriging(skew_dados, dados.grid, formula)
  plot2 <- spplot(teste.grid, zcol = "ok_isotropic_skew_normal_dist", main = list(font = 1, cex = 1, label = "KO"),  sub = list(font = 1, cex = 1, label = "(f) Com rem. de anisotropia/assimetria"))
  print(plot1)
  print(plot2)
}


#Para atributos assimétricos sem outliers
if(!normalDistribution(teste2, formula) && (length(dados) != length(teste2))){
  #teste2_normal_dist <- logTransformation(teste2, formula)
  #assimetria_depois_transformacao <- skew(teste2_normal_dist@data[, variavel])
  #anisotropia_sem_outlier_normal_dist <- checkAnisotropy(teste2, formula)
  #teste.grid$ok_isotropic_no_outlier_normal_dist <- ordinaryKriging(teste2_normal_dist, dados.grid, formula, backtransform = TRUE)
  #plotLogNormalOrdinaryKriging(teste2, dados.grid, formula)
  plot_ok_isotropic_no_outlier_normal_dist <- plotLogNormalOrdinaryKriging(teste2, dados.grid, formula)
  plot1 <- plot(plot_ok_isotropic_no_outlier_normal_dist$exp_var, plot_ok_isotropic_no_outlier_normal_dist$var_model, sub = list(font = 1, cex = 1, label = "(d) Com rem. de outliers/assimetria"))
  teste.grid$ok_isotropic_no_outlier_normal_dist <- logNormalOrdinaryKriging(teste2, dados.grid, formula)
  plot2 <- spplot(teste.grid, zcol = "ok_isotropic_no_outlier_normal_dist", main = list(font = 1, cex = 1, label = "KO"),  sub = list(font = 1, cex = 1, label = "(d) Com rem. de outliers/assimetria"))
  print(plot1)
  print(plot2)
}


#TESTAR IDW ISOTROPICO

  #Interpolação IDW
teste.grid$idw <- inverseDistanceWeighted(dados, dados.grid, formula)
spplot(teste.grid, zcol = "idw", main = list(font = 1, cex = 1, label = "IDQ"),  sub = list(font = 1, cex = 1, label = "(b) Sem tratamentos"))
if(length(dados) != length(teste2)){
  teste.grid$idw_no_outlier <- inverseDistanceWeighted(teste2, dados.grid, formula)
  spplot(teste.grid, zcol = "idw_no_outlier", main = list(font = 1, cex = 1, label = "IDQ"),  sub = list(font = 1, cex = 1, label = "(d) Com remoção de outliers"))
}



#Teste RMSE
rmse_result$ok <- rmse(dados, formula, ordinaryKrigingModified)
if(length(dados) != length(teste2)){
  rmse_result$ok_no_outlier <- rmse(teste2, formula, ordinaryKrigingModified)
  rmse_result$ok_detrended_no_outlier <- rmse(teste2, formula, ordinaryKrigingDetrended)
  rmse_result$uk_no_outlier <- rmse(teste2, formula, universalKriging)
}
if(!is.null(anisotropia)){
  rmse_result$ok_isotropic <- rmse(dados, formula, ordinaryKriging)
}
if(!is.null(anisotropia_sem_outlier)){
  rmse_result$ok_isotropic_no_outlier <- rmse(teste2, formula, ordinaryKriging)
}
if(!normalDistribution(skew_dados, formula)){
  rmse_result$ok_isotropic_skew_normal_dist <- rmse(skew_dados, formula, logNormalOrdinaryKriging)
}
if(!normalDistribution(teste2, formula)){
  rmse_result$ok_isotropic_no_outlier_normal_dist <- rmse(teste2, formula, logNormalOrdinaryKriging)
  }
rmse_result$idw <- rmse(dados, formula, inverseDistanceWeighted)
if(length(dados) != length(teste2)){
  rmse_result$idw_no_outlier <- rmse(teste2, formula, inverseDistanceWeighted)
}


#Mostra os valores das an?lises
print(variavel)
quant_outliers
assimetria
if(length(dados) != length(teste2)){
  assimetria_sem_outlier
}
if(!normalDistribution(dados, formula)){
  skew_assimetria_depois_transformacao
}
if(!normalDistribution(teste2, formula)){
  assimetria_depois_transformacao
}
anisotropia
if(length(dados) != length(teste2)){
  anisotropia_sem_outlier
}
if(!normalDistribution(dados, formula)){
  anisotropia_skew_normal_dist
}
if(!normalDistribution(teste2, formula)){
  anisotropia_sem_outlier_normal_dist
}
rmse_result
