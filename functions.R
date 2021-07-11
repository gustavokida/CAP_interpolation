
variavel_1 <- 'C'
dados_teste <- subset(dados, !variavel_1 %in% boxplot.stats(dados$C)$out)
teste123 <- dados[,variavel_1]
atributo_test <- dados %>% extract2(variavel_1)
zscore2 <- dados$MG[abs((dados$MG-mean(dados$MG))/sd(dados$MG)) <= 2]

library(dplyr)    #na_if
library(magrittr)   #pipeline
library(spatialEco)   #removeNA
library(psych)    #normal distribution
library(lmtest)   #bptest
library(gvlma)    #skew, kurtosis, heteroscedasticity
library(Hmisc)    #rcorr for pearson correlation
library(gstat)    #variogram
library(automap)    #autofitVariogram
library(ape)    #Moran's I, checks for spatial autocorrelation
library(spdep)  #Moran's I
library(geoR)   #trend spatial(ainda verificando), #coords.aniso
library(rlang)    #formulaToVector
library(intamap)    #estimateAnisotropy
library(car)    #vif
library(MASS)   #stepAIC, boxcox
library(forecast)    #boxcox
library(raster)   #spPixelsToRaster
library(EnvStats)   #Coefficient of variation

library(spatial)
library(nlme)   #gls
library(rms)    #gls
library(lme4)   #lmer



makeGrid <- function(contour = "SpatialPolygons", cellsize = "numeric"){
  grid <- spsample(contour, n = 0, cellsize = cellsize, "regular")
  d1 <- data.frame(grid@coords[,1])
  d2 <- data.frame(grid@coords[,2])
  d3 <- cbind(d1,d2)
  names(d3)[1] <- "x"
  names(d3)[2] <- "y"
  grid_spdf <- SpatialPixelsDataFrame(points = d3[c("x", "y")], data = d3 , proj4string = grid@proj4string)
  return(grid_spdf)
}


rasterToSpPixels <- function(data = c("RasterLayer", "Raster")){
  return(as(data, "SpatialPixelsDataFrame"))
}



spPixelsToRaster <- function(data = c("SpatialPointsDataFrame", "SpatialPixelsDataFrame")){
  return(raster(data))
}



#remover toda a linha(row)
removeOutlier <- function(data = "SpatialPointsDataFrame", column_names = c("character", "vector")){
  for(column in column_names){
    attribute <- data %>% extract2(column)
    outliers <- boxplot.stats(attribute)$out
    for(outlier in outliers){
      attribute <- na_if(attribute, outlier)
    }
    data@data[, column] <- attribute
  }
  return(data)
}



removeNA <- function(data = "SpatialPointsDataFrame", column_name = c("character", "vector")){
  if(anyNA(data@data[, column_name])){
    result <- sp.na.omit(data, col.name = column_name)
  }
  else{
    result <- data
  }
  return(result)
}



checkQuantity <- function(data = "SpatialPointsDataFrame", column_name = c("character", "vector")){
  row_number <- length(na.omit(data@data[, column_name]))
  if (row_number > 50){
    result <- TRUE
  }
  else{
    result <- FALSE
  }
  return(result)
}



logTransformation <- function(data = "SpatialPointsDataFrame", formula = "formula"){
  column_name <- formulaToVector(formula, "left")
  data@data[, column_name] <- log(data@data[, column_name])
  return(data)
}



#talvez arrumar (procurar por unbiased backtransformation)
backTransformation <- function(interpolated_data = "SpatialPointsDataFrame", backtransform = FALSE){
  backtransformed_var <- NULL
  if(is_true(backtransform)){
    backtransformed_var <- exp(interpolated_data$krige_output$var1.pred + (interpolated_data$krige_output$var1.var/2))
  }
  else{
    backtransformed_var <- interpolated_data
  }
  return(backtransformed_var)
  
}



boxCoxTransformation <- function(formula = "formula", data = c("SpatialPointsDataFrame", "autoKrige"), lambda = NULL){
  # lm_model <- lm(formula, data)
  # bc <- boxcox(lm_model)
  # best_lam <- bc$x[which(bc$y==max(bc$y))]
  bc <- NULL
  left_var <- formulaToVector(formula, "left")
  if (is.null(lambda)){
    lambda <- BoxCox.lambda(data@data[, left_var], method = 'loglik')
    bc <- BoxCox(data@data[, left_var], lambda)
  }
  else{
    if (is(data, 'autoKrige')){
      bc <- InvBoxCox(as.numeric(unlist(data$krige_output$var1.pred)), lambda)
    }
    else{
      bc <- InvBoxCox(data@data[, left_var], lambda)
    }
  }
  output <- as.numeric(unlist(bc))
  return(output)
}


#talvez remover o (column in column_names)
normalDistribution <- function(data = "SpatialPointsDataFrame", formula = "formula"){
  column_name <- formulaToVector(formula, "left")
  attribute <- data %>% extract2(column_name)
  skewness <- skew(attribute)
  if(skewness > 1){
    result <- FALSE
  }
  else{
    result <- TRUE
  }
  return(result)
}



makeFormula <- function(main_attribute_column = c("character", "vector"), column_names = "1"){
  formula <- reformulate(termlabels = column_names, response = main_attribute_column)
  return(formula)
}



formulaToVector <- function(formula = "formula", side = c("character", "vector")){
  if(side == "left"){
    return(all.vars(f_lhs(formula)))
  }
  else if(side == "right"){
    right_side = f_rhs(formula)
    if(is.numeric(right_side)){
      return(as.character(right_side))
    }
    else{
      return(all.vars(right_side))
    }
  }
  else if(side == "all"){
    column_names <- append(formulaToVector(formula, side = "left"), formulaToVector(formula, side = "right"))
    return(column_names)
  }
}



modelAnalysis <- function(model ="lm"){
  analysis_result <- gvlma(model)
  return(analysis_result)
}



#checar moran index monte carlo
checkSpatialStructure <- function(data = "SpatialPointsDataFrame", column_name = c("character", "vector")){
  moranI_coords.dists <- as.matrix(dist(cbind(data@coords[,2], data@coords[,1])))
  moranI_coords.dists.inv <- 1/moranI_coords.dists
  diag(moranI_coords.dists.inv) <- 0
  moranI_result <- Moran.I(data@data[, column_name], moranI_coords.dists.inv)
  zscore <- (moranI_result$observed - moranI_result$expected)/moranI_result$sd
  if(moranI_result$p.value <= 0.05 && zscore > 1.96){
    if(moranI_result$observed > 0.3){
      result <- TRUE
    }
    else{
      result <- FALSE
    }
  }
  else{
    result <- FALSE
  }
  return(result)
}

#checkSpatialStructure <- function(data = "SpatialPointsDataFrame", formula = "formula"){
  #variogram <- autofitVariogram(formula = formula, input_data = data)
  #result <- variogram$var_model$psill[1] / (variogram$var_model$psill[1] + variogram$var_model$psill[2])
  #return(result)
#}

#verificar ordem para rotacionar de volta os dados
#ajeitar as fun?oes de anisotropia
handleAnisotropy <- function(data = "SpatialPointsDataFrame",formula = "formula", reverse_data = NULL, reverse = FALSE){
  if(reverse == FALSE){
    anisotropy <- checkAnisotropy(data, formula)
    if(!is.null(anisotropy)){
      rotated_coords <- rotateAnisotropy(data, anisotropy)
      result <- data
      result@coords <- rotated_coords
      return(result)
    }
    else{
      return(data)
    }
  }
  else if (reverse == TRUE){
    result <- rotateBackAnisotropy(data, anisotropy)
    return(result)
  }
}


checkAnisotropy <- function(data = "SpatialPointsDataFrame", formula = "formula"){
  column_name <- formulaToVector(formula, "left")
  anisotropy <- estimateAnisotropy(data, depVar = column_name)
  if(anisotropy$doRotation == TRUE){
    return(c(anisotropy$direction, anisotropy$ratio))
  }
  else{
    return(NULL)  
  }
}

rotateAnisotropy <- function(data = "SpatialPointsDataFrame", anis_var = c("list", "vector")){
  result <- coords.aniso(data@coords, anis_var)
  return(result)
}

rotateBackAnisotropy <- function(data = "SpatialPointsDataFrame", anis_var = c("list", "vector")){
  result <- coords.aniso(data@coords, anis_var, reverse = TRUE)
  return(result)
}




#georob predict calcula drift(trend)
#package EcoGenetics discontinued
eco.detrend(data, coordinates(data))
surf.ls()   #package spatial
trend.spatial("1st", dados)


#ratio_1 maior ou menor que ratio_2?  #CHECAR QUAL RESIDUO TEM A MENOR MEDIA
checkVariogram <- function(formula_1 = "formula", formula_2 = NULL, data = "SpatialPointsDataFrame"){
  best_formula <- NULL
  if(is.null(formula_2)){
    best_formula <- formula_1
  }
  else{
    variogram_1 <- autofitVariogram(formula_1, data)
    variogram_2 <- autofitVariogram(formula_2, data)
    ratio_1 = variogram_1$var_model$psill[1]/sum(variogram_1$var_model$psill)
    ratio_2 = variogram_2$var_model$psill[1]/sum(variogram_2$var_model$psill)
    if(ratio_1 > ratio_2){  #ratio_1 > ratio_2
      best_formula <- formula_1
    }
    else{
      best_formula <- formula_2
    }
  }
  return(best_formula)
}


#dist(data@coords)


checkTrend <- function(data = "SpatialPointsDataFrame", formula = "formula"){
  auto_variogram <- autofitVariogram(formula, data)
  if(max(auto_variogram$exp_var$dist) > max(auto_variogram$var_model$range)){    #TALVEZ IMPLEMENTAR UM autofitvariogram PARA VERIFICAR A RANGE
    result <- FALSE
  }
  else{
    result <- TRUE
  }
  return(result)
}


detrendFormula <- function(data = "SpatialPointsDataFrame", formula = "formula"){
  degree <- 0
  new_formula <- formula
  best_formula <- formula
  left_formula <- as.symbol(formulaToVector(formula, "left"))
  if(checkTrend(data, new_formula)){
    best_formula <- NULL
    while(degree < 10){
      degree = degree + 1
      formulas <- c(expr(poly(x, degree = !!degree)), 
                    expr(poly(y, degree = !!degree)),
                    expr(poly(x, y, degree = !!degree)))
      for(right_formula in formulas){
        new_formula <- as.formula(expr(!!left_formula ~ !!right_formula))
        if(!checkTrend(data,new_formula)){
          best_formula <- checkVariogram(new_formula, best_formula, data)
        }
      }
      if(!is.null(best_formula)){
        break
      }
    }
  }
  return(best_formula)
}


#Talvez colocar um function regression_type para escolher o tipo de regressão
checkTrendLm <- function(data = "SpatialPointsDataFrame", formula = "formula"){ #regression_type = function
  lm_model <- lm(formula, data)
  data@data["residuals"] <- lm_model$residuals
  auto_variogram <- autofitVariogram(residuals~1, data)
  if(max(auto_variogram$exp_var$dist) > max(auto_variogram$var_model$range)){    #TALVEZ IMPLEMENTAR UM autofitvariogram PARA VERIFICAR A RANGE
    result <- FALSE
  }
  else{
    result <- TRUE
  }
  return(result)
}

detrend <- function(data = "SpatialPointsDataFrame", formula = "formula"){
  degree <- 0
  lm_model <- NULL
  new_formula <- formula
  left_formula <- as.symbol(formulaToVector(formula, "left"))
  while(checkTrendLm(data,new_formula)){
    degree <- degree + 1
    right_formula <- expr(poly(x, y, degree = !!degree))
    new_formula <- as.formula(expr(!!left_formula ~ !!right_formula))
    lm_model <- lm(new_formula, data = data)
  }
  return(lm_model)
}



checkCovariatesCKO <- function(data = "SpatialPointsDataFrame", covariate_data = NULL, main_attribute_column = c("character", "vector"), column_names = c("character", "vector")){
  if(!is.null(covariate_data)){
    for(i in column_names){
      formula <- makeFormula(i)
      data@data[i] <- inverseDistanceWeighted(covariate_data, data, formula)
    }
  }
  pearson_result <- pearsonCorrelation(data, main_attribute_column, column_names)
  return(pearson_result)
}



checkCovariates <- function(data = "SpatialPointsDataFrame", covariate_data = NULL, main_attribute_column = c("character", "vector"), column_names = c("character", "vector")){
  if(!is.null(covariate_data)){
    for(i in column_names){
      formula <- makeFormula(i)
      data@data[i] <- inverseDistanceWeighted(covariate_data, data, formula)
    }
  }
  pearson_result <- pearsonCorrelation(data, main_attribute_column, column_names)
  if(is.character(pearson_result) && pearson_result != 1 && length(pearson_result) > 1){
    vif_result <- handleVIF(lm(makeFormula(main_attribute_column, pearson_result), data))
    if(vif_result != 1){
      AIC_result <- handleAIC(lm(makeFormula(main_attribute_column, vif_result), data))
    }
    else{
      AIC_result = 1
    }
  }
  else{
    AIC_result = 1
  }
  return(AIC_result)
}



#ADICIONAR CONJUNTO DAS COVARIAVEIS
#Selec?o de covari?veis p-value <= 0,01, cor > 0.4 and cor < 0.8
pearsonCorrelation <- function(data = "SpatialPointsDataFrame", main_attribute_column = c("character", "vector"), column_names = c("character", "vector")){
  accepted_covariates <- c()
  for(column in column_names){
    correlation <- cor.test(data@data[, main_attribute_column], data@data[, column], method = "pearson")
    if((correlation$p.value <= 0.01) && (correlation$estimate > 0.1) && (correlation$estimate < 0.99)){
      accepted_covariates <- append(accepted_covariates, column)
    }
  }
  if(length(accepted_covariates) == 0){
    accepted_covariates = 1
  }
  return(accepted_covariates)
}


#Retirar VIF > 10 (remo??o de colinearidade)
handleVIF <- function(model = "lm"){
  vif_model <- vif(model)
  accepted_covariates <- names(vif_model)[which(vif_model <= 10)]
  if(length(accepted_covariates) == 0){
    accepted_covariates = 1
  }
  return(accepted_covariates)
}



#stepwise multiple regression
handleAIC <- function(model = "model"){
  selected_aic <- stepAIC(model,direction="both")
  accepted_covariates <- names(selected_aic$coefficients)[which(names(selected_aic$coefficients) != "(Intercept)")]
  if(length(accepted_covariates) == 0){
    accepted_covariates = 1
  }
  return(accepted_covariates)
}


#stepAIC(lm(C~1,data=dados),direction="both", scope=~CTC+HAL+x)   #scope start from 0 covariates and then add covariates step by step in the model



#kroenker bp  test para nao-estacionaridade(valor p < 0.05)
bpTest <-function(model = "lm"){
  bp_test <- bptest(model, studentize = FALSE)
  if(bp_test < 0.05){
    result <- TRUE
  }
  else{
    result <- FALSE
  }
  return(result)
}



checkAbruptChange <- function{
  
}


#moran's I
checkCluster <- function{
  
}



checkNoise <- function{
  
}



checkArtefact <- function{
  
}












