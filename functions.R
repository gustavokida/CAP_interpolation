
variavel_1 <- 'C'
dados_teste <- subset(dados, !variavel_1 %in% boxplot.stats(dados$C)$out)
teste123 <- dados[,variavel_1]
atributo_test <- dados %>% extract2(variavel_1)
zscore2 <- dados$MG[abs((dados$MG-mean(dados$MG))/sd(dados$MG)) <= 2]

library(methods)
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




#Transform SpatialPolygonsDataFrame to SpatialPixelsDataFrame
makeGrid <- function(contour = "SpatialPolygons", cellsize = "numeric"){
  #makes a SpatialPoints from SpatialPolygonsDataFrame contour
  grid <- spsample(contour, n = 0, cellsize = cellsize, "regular")
  #get coordinates
  d1 <- data.frame(grid@coords[,1])
  d2 <- data.frame(grid@coords[,2])
  d3 <- cbind(d1,d2)
  #name the coordinates
  names(d3)[1] <- "x"
  names(d3)[2] <- "y"
  #transform the SpatialPoints to SpatialPixelsDataFrame
  grid_spdf <- SpatialPixelsDataFrame(points = d3[c("x", "y")], data = d3 , proj4string = grid@proj4string)
  return(grid_spdf)
}


#Transform  raster to SpatialPixelsDataframe
rasterToSpPixels <- function(data = c("RasterLayer", "Raster")){
  return(as(data, "SpatialPixelsDataFrame"))
}


#Transform SpatialPixelsDataFrame  to raster
spPixelsToRaster <- function(data = c("SpatialPointsDataFrame", "SpatialPixelsDataFrame")){
  return(raster(data))
}



#Removes all outliers, replacing outliers with NA values
removeOutlier.default <- function(data = "SpatialPointsDataFrame", column_names = c("character", "vector")){
  #Transforms all attributes in column_names
  for(column in column_names){
    #get column
    attribute <- data %>% extract2(column)
    #get outliers
    outliers <- boxplot.stats(attribute)$out
    #replace any value in outliers to NA
    for(outlier in outliers){
      attribute <- na_if(attribute, outlier)
    }
    data@data[, column] <- attribute
  }
  removeNA(data, column_names)
  return(data)
}


#Removes all rows containing NA values
removeNA <- function(data = "SpatialPointsDataFrame", column_names = c("character", "vector")){
  for(column in column_names){
    if(anyNA(data@data[, column])){
      #remove entire row, even if only one value is NA
      data <- sp.na.omit(data, col.name = column)
    }
  }
  return(data)
}


#check if there is enough points in dataframe to do kriging (more than 50)
checkQuantity.default <- function(data = "SpatialPointsDataFrame", column_name = c("character", "vector")){
  row_number <- length(na.omit(data@data[, column_name]))
  #more than 50 points kriging is acceptable, making it TRUE
  if (row_number > 50){
    result <- TRUE
  }
  #less or equal 50 points kriging is not acceptable, making it FALSE
  else{
    result <- FALSE
  }
  return(result)
}


#make a Log Transformation on the dependent variable in the formula
logTransformation <- function(data = "SpatialPointsDataFrame", formula = "formula"){
  column_name <- formulaToVector(formula, "left")
  data@data[, column_name] <- log(data@data[, column_name])
  return(data)
}



#talvez arrumar (procurar por unbiased backtransformation)
#backtransform the log transformed data
backTransformation <- function(interpolated_data = "SpatialPointsDataFrame", backtransform = FALSE){
  backtransformed_var <- NULL
  if(is_true(backtransform)){
    #Backtransform the transformed data using exponential and then sum it with the half of the variance
    backtransformed_var <- exp(interpolated_data$krige_output$var1.pred + (interpolated_data$krige_output$var1.var/2))
  }
  else{
    backtransformed_var <- interpolated_data
  }
  return(backtransformed_var)
  
}


#Apply BoxCox Transformation using method loglik
boxCoxTransform.default <- function(formula = "formula", data = c("SpatialPointsDataFrame", "autoKrige"), lambda = "numeric", reverseBoxCox = "logical"){
  # lm_model <- lm(formula, data)
  # bc <- boxcox(lm_model)
  # best_lam <- bc$x[which(bc$y==max(bc$y))]
  bc <- NULL
  left_var <- formulaToVector(formula, "left")
  if (isFALSE(reverseBoxCox)){
    #Calculate the lambda
    #lambda <- BoxCox.lambda(data@data[, left_var], method = 'loglik')
    #Apply BoxCox with the lambda value
    bc <- BoxCox(data@data[, left_var], lambda)
  }
  #if lambda is supplied in the function, it will backtransform the data
  else{
    #data with class autoKrige
    if (is(data, 'autoKrige')){
      bc <- InvBoxCox(as.numeric(unlist(data$krige_output$var1.pred)), lambda)
    }
    #data with class SpatialPoints
    else{
      bc <- InvBoxCox(data@data[, left_var], lambda)
    }
  }
  output <- as.numeric(unlist(bc))
  return(output)
}


#checks the normal distribution, if skewness > 1, the data is not normally distributed
normalDistribution.default <- function(data = "SpatialPointsDataFrame", formula = "formula"){
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


#Makes formula by inserting a main attribute
#if there is no independent variable, it automatically makes a formula like x~1
makeFormula <- function(main_attribute_column = c("character", "vector"), column_names = "1"){
  formula <- reformulate(termlabels = column_names, response = main_attribute_column)
  return(formula)
}


#transform a formula to a vector of characters
#choose the "left", "right" or "all" sides
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


#Currently not in use
modelAnalysis <- function(model ="lm"){
  analysis_result <- gvlma(model)
  return(analysis_result)
}



#checar moran index monte carlo
#check Spatial Structure or spatial autocorrelation using Moran's Index
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
#Treats anisotropy on the attribute automatically
#if reverse is true, backtransform the data that was treated with isotropic transformation
handleAnisotropy.default <- function(data = "SpatialPointsDataFrame",formula = "formula", reverse_data = NULL, reverse = FALSE){
  if(reverse == FALSE){
    anisotropy <- checkAnisotropy(data, formula)
    if(!is.null(anisotropy)){
      rotated_coords <- coords.aniso(data@coords, anisotropy)
      colnames(rotated_coords) <- c("x", "y")
      result <- data
      result@coords <- rotated_coords
      return(result)
    }
    else{
      return(data)
    }
  }
  else if (reverse == TRUE){
    rotated_coords <- coords.aniso(data@coords, anisotropy, reverse = TRUE)
    colnames(rotated_coords) <- c("x", "y")
    result <- data
    result@coords <- rotated_coords
    return(result)
  }
}



#checks anisotropy, if detected, returns anisotropy direction and ratio
checkAnisotropy.default <- function(data = "SpatialPointsDataFrame", formula = "formula"){
  column_name <- formulaToVector(formula, "left")
  anisotropy <- estimateAnisotropy(data, depVar = column_name)
  if(anisotropy$doRotation == TRUE){
    return(c(anisotropy$direction, anisotropy$ratio))
  }
  else{
    return(NULL)  
  }
}



#georob predict calcula drift(trend)
#package EcoGenetics discontinued
eco.detrend(data, coordinates(data))
surf.ls()   #package spatial
trend.spatial("1st", dados)




#ratio_1 maior ou menor que ratio_2?  #CHECAR QUAL RESIDUO TEM A MENOR MEDIA

#compares 2 formulas using spatial structure of the variogram using nugget to sill ratio
#best formula is determined by a larger ratio
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


#in test
#check trend in variogram by analising if the distance of the farthest lag is bigger than the range of the variogram
#if the range is bigger than the distance, it indicates a trend
checkTrend <- function(data = "SpatialPointsDataFrame", formula = "formula"){
  #generates an automatic variogram
  auto_variogram <- autofitVariogram(formula, data)
  if(max(auto_variogram$exp_var$dist) > max(auto_variogram$var_model$range)){    #TALVEZ IMPLEMENTAR UM autofitvariogram PARA VERIFICAR A RANGE
    result <- FALSE
  }
  else{
    result <- TRUE
  }
  return(result)
}



#uses checkTrend to find a trend, and then detrend  by searching the best formula using only the coordinates
#it will search the best polynomial formula
detrendFormula.default <- function(data = "SpatialPointsDataFrame", formula = "formula"){
  degree <- 0
  new_formula <- formula
  best_formula <- formula
  left_formula <- as.symbol(formulaToVector(formula, "left"))
  if(checkTrend(data, new_formula)){
    best_formula <- NULL
    #checks polynomials from 1 to 10, e.g (x+y), (x+y)^2, (x+y)^3, ...
    while(degree < 10){
      degree = degree + 1
      #test formulas with x, y and x+y
      formulas <- c(expr(poly(x, degree = !!degree)), 
                    expr(poly(y, degree = !!degree)),
                    expr(poly(x, y, degree = !!degree)))
      for(right_formula in formulas){
        new_formula <- as.formula(expr(!!left_formula ~ !!right_formula))
      #choose the detrended formula with the lowest polynomial
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
#check trend in residual variogram, same way as the function checkTrend
checkResidualTrend <- function(data = "SpatialPointsDataFrame", formula = "formula"){ #regression_type = function
  lm_model <- lm(formula, data)
  data@data["residuals"] <- lm_model$residuals
  #makes a variogram with the residuals from the linear regression
  auto_variogram <- autofitVariogram(residuals~1, data)
  if(max(auto_variogram$exp_var$dist) > max(auto_variogram$var_model$range)){    #TALVEZ IMPLEMENTAR UM autofitvariogram PARA VERIFICAR A RANGE
    result <- FALSE
  }
  else{
    result <- TRUE
  }
  return(result)
}

#detrend a formula with x+y, selecting the lowest polynomial without trend
detrendResidual <- function(data = "SpatialPointsDataFrame", formula = "formula"){
  degree <- 0
  lm_model <- NULL
  new_formula <- formula
  left_formula <- as.symbol(formulaToVector(formula, "left"))
  #check for trend until it finds a detrented formula and then return the model
  while(checkTrendLm(data,new_formula)){
    degree <- degree + 1
    right_formula <- expr(poly(x, y, degree = !!degree))
    new_formula <- as.formula(expr(!!left_formula ~ !!right_formula))
    lm_model <- lm(new_formula, data = data)
  }
  return(lm_model)
}



#choose the best combination of covariates for Ordinary Co-kriging
#it can search covariates in the same dataset
#if covariate_data is specified, it will search for covariates in the auxiliary datasets
checkCovariatesCKO.default <- function(data = "SpatialPointsDataFrame", covariate_data = NULL, main_attribute_column = c("character", "vector"), column_names = c("character", "vector")){
  if(!is.null(covariate_data)){
    #this is used to reduce the number of points from larger datasets to allow tests like pearsonCorrelation
    for(i in column_names){
      formula <- makeFormula(i)
      data@data[i] <- inverseDistanceWeighted(covariate_data, data, formula)
    }
  }
  pearson_result <- pearsonCorrelation(data, main_attribute_column, column_names)
  return(pearson_result)
}



#choose the best combination of covariates for Geographically Weighted Regression
#search covariates in the same dataset or in a auxiliary dataset
checkCovariates.default <- function(data = "SpatialPointsDataFrame", covariate_data = NULL, main_attribute_column = c("character", "vector"), column_names = c("character", "vector")){
  if(!is.null(covariate_data)){
    #use this to reduce the number of larger auxiliary dataset to allow tests like pearsonCorrelation
    for(i in column_names){
      formula <- makeFormula(i)
      data@data[i] <- inverseDistanceWeighted(covariate_data, data, formula)
    }
  }
  #check correlation
  pearson_result <- pearsonCorrelation(data, main_attribute_column, column_names)
  if(is.character(pearson_result) && pearson_result != 1 && length(pearson_result) > 1){
    #check for multicollinearity
    vif_result <- handleVIF(lm(makeFormula(main_attribute_column, pearson_result), data))
    if(vif_result != 1){
      #check for best model
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
#choose covariables using pearson correlation with p-value <= 0,01 and correlation between cor > 0.4 and cor < 0.8
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



#Remove multicollinearity by selecting the covariables with the Variance Inflation Factor(VIF) <= 10
handleVIF <- function(model = "lm"){
  vif_model <- vif(model)
  accepted_covariates <- names(vif_model)[which(vif_model <= 10)]
  if(length(accepted_covariates) == 0){
    accepted_covariates = 1
  }
  return(accepted_covariates)
}



#Select the best model using Stepwise Multiple Regression
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












