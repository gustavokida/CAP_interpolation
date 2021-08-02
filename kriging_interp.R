
library(GSIF)   #Regression Kriging
library(gstat)    #cokriging
library(automap)    #autoKrige
library(GWmodel)


#ordinary kriging
ordinaryKriging.default <- function(data = "SpatialPointsDataFrame", newdata = c("SpatialPointsDataFrame", "SpatialPixelsDataFrame"), formula = "formula", rmse_data = NULL){
  if(is.null(rmse_data)){
    #check and remove anisotropy
    data_variogram <- handleAnisotropy(data, formula)
    #do the kriging with the generated variogram
    interpolated_data <- autoKrige(formula, data, newdata, data_variogram = data_variogram)
    return(interpolated_data$krige_output$var1.pred)
  }
  else{
    #check and remove anisotropy
    data_variogram <- handleAnisotropy(rmse_data, formula)
    #do the kriging with the generated variogram
    interpolated_data <- autoKrige(formula, data, newdata, data_variogram = data_variogram)
    return(interpolated_data$krige_output$var1.pred)
  }
}



universalKriging.default <- function(data = "SpatialPointsDataFrame", newdata = c("SpatialPointsDataFrame", "SpatialPixelsDataFrame"), formula = "formula", rmse_data = NULL){
  #select the best detrended formula
  new_formula <- detrendFormula(data, formula)
  #performs Universal Kriging with the selected formula
  interpolated_data <- autoKrige(new_formula, data, newdata)
  return(interpolated_data$krige_output$var1.pred)
}



#automatically do Universal Kriging if a trend is detected, otherwise performs ordinary kriging
autoKriging <- function(data = "SpatialPointsDataFrame", newdata = c("SpatialPointsDataFrame", "SpatialPixelsDataFrame"), formula = "formula"){
  #check trend
  if(checkTrend(data, formula)){
    #finds the best formula
    formula <- detrendFormula(data, formula)
    #generate a variogram
    data_variogram <- autofitVariogram(formula, data)
  }
  else{
    #if there is no trend, check and remove anisotropy
    data_variogram <- handleAnisotropy(data, formula)
  }
  #perform kriging
  interpolated_data <- autoKrige(formula, data, newdata, data_variogram = data_variogram)
}



#arrumar
#gstat

#Cokriging function, it is necessary to inform the covariate data
#if the covariate data is in the primary dataset, put the same dataset in the covariate_data e.g. data = dataset, covariate_data = dataset
coKriging.default <- function(data = "SpatialPointsDataFrame", newdata = c("SpatialPointsDataFrame", "SpatialPixelsDataFrame"),
                              formula = "formula", covariate_data = "SpatialPointsDataFrame", rmse_data = NULL){
  #separate the primary variable from the auxiliary variables
  left_var <- formulaToVector(formula = formula, side = "left")
  right_cov <- formulaToVector(formula = formula, side = "right")
  gstat_object <- NULL
  #makes a formula with the primary variable
  this_formula <- makeFormula(left_var)
  #generate a gstat object with the primary variable
  gstat_object <- gstat(gstat_object, id=left_var, formula=this_formula, data=data, nmax = 10)
  #generate others gstat objects inside the same one above using the auxiliary variables
  for(variable in right_cov){
    #makes a formula for each variable that is in the right side of the formula
    this_formula <- makeFormula(variable)
    gstat_object <- gstat(gstat_object, id=variable, formula=this_formula, data=covariate_data, nmax = 10)
  }
  #makes a variogram of the gstat object
  variogram_model <- variogram(gstat_object)#, cutoff=1000)
  #makes a vgm model
  model <- vgmModel(data = data, formula = makeFormula(left_var))
  #insert the model in the gstat object
  gstat_object <- gstat(gstat_object, fill.all=T, model=model)
  #fit the variogram in the gstat object
  gstat_object.fit <- fit.lmc(v = variogram_model,g = gstat_object)
  plot(variogram_model, model = gstat_object.fit)
  gstat_object_cv <- gstat.cv(gstat_object.fit)
  #try to interpolate, if the cross-variogram is bad, the predict function will throw an error
  #if successful, return the interpolated data
  interpolated_data <- tryCatch(
                        {
                          interpolated_data <- predict(gstat_object.fit, newdata = newdata)
                        },
                      error=function(cond){
                        return(NULL)
                      }
                    )
  if(is.null(interpolated_data)){
    return(NULL)
  }
  #return(interpolated_data@data[,1])
  return(c(gstat_object_cv, interpolated_data))
}


#fit a variogram automatically
vgmModel <- function(data = "SpatialPointsDataFrame", formula = "formula"){
  return(autofitVariogram(formula = formula, input_data = data)$var_model)
}



#newdata em fit.gstatModel pode ser diferente do newdata em predict.gstatModel    #newdata em fit.gstatmodel recebe dataframe de covariaveis

#Regression Kriging
regressionKriging.default <- function(data = "SpatialPointsDataFrame", newdata = c("SpatialPointsDataFrame", "SpatialPixelsDataFrame"),
                                      formula = "formula", covariate_data = "SpatialPointsDataFrame", rmse_data = NULL){
  fitted_data <- fit.gstatModel(data, formula, covariate_data, fit.family = gaussian())    #Must use covariates in SpatialPixelsDataFrame type
  interpolated_data <- predict.gstatModel(fitted_data, newdata)
  return(interpolated_data@predicted$var1.pred)
}




krigingED <- function(){
  
}


#work in progress
gwr <- function(data = "SpatialPointsDataFrame", newdata = c("SpatialPointsDataFrame", "SpatialPixelsDataFrame"), formula = "formula"){
  bw <- bw.gwr(formula, data = data, kernel = "bisquare", adaptive = TRUE)
  colin <- gwr.collin.diagno(formula, data = data, bw = bw, kernel = "bisquare", adaptive = TRUE)
  statio <- gwr.montecarlo(formula, data = data, bw = bw, kernel = "bisquare", adaptive = TRUE)
  model <- gwr.basic(formula, data = data, newdata, kernel = "bisquare", bw = bw, adaptive = TRUE)
  predict <- gwr.predict(formula, data = data, predictdata=newdata, bw = bw, kernel = "bisquare", adaptive = TRUE)
  return(predict)
}


gwrKriging <- function(data = "SpatialPointsDataFrame", newdata = c("SpatialPointsDataFrame", "SpatialPixelsDataFrame"), formula = "formula"){
  predict <- gwr(data, newdata, formula)
  data@data$gwr_residuals <- gwr.model$lm$residuals
  interpolated_data <- ordinaryKriging(data, newdata, gwr_residuals~1)
  return(cgwr.predict)
}



##interpola??o GWR
cgwr.bw <- bw.gwr(C~CTC+HAL, data = dados, kernel = "bisquare", adaptive = TRUE)
cgwr.bw <- bw.gwr(C~CTC+HAL, data = dados, kernel = "bisquare", adaptive = FALSE)
cgwr.bw <- bw.gwr.lcr(C~CTC+HAL, data = dados, kernel = "bisquare", adaptive = TRUE)
cgwr.colin <- gwr.collin.diagno(C~CTC+HAL, data = dados, bw = cgwr.bw, kernel = "bisquare", adaptive = TRUE)
cgwr.statio <- gwr.montecarlo(C~CTC+HAL, data = dados, bw = cgwr.bw, kernel = "bisquare", adaptive = TRUE)
cgwr.model <- gwr.basic(C~CTC+HAL, data = dados, dados.grid, kernel = "bisquare", bw = cgwr.bw, adaptive = TRUE)

#variograma com res?duos
residuos <- cgwr.model$lm$residuals
cgwr.residuos <- data.frame(residuos)
cgwr.residuos$fittedvalue <- cgwr.model$lm$fitted.values
cgwr.variogram <- SpatialPointsDataFrame(coords = nutrientescoords, data = cgwr.residuos, proj4string = nutrientescrs)
cgwr.vario<- variogram(residuos~1, data=cgwr.variogram)
#"Exp","Sph","Gau","Mat"
cgwr.fitvario <- fit.variogram(cgwr.vario, vgm("Exp"))
plot(cgwr.vario, cgwr.fitvario, xlab = "Distancia", ylab = "Semivariancia")

#interpola??o
cgwr.predict <- gwr.predict(C~CTC+HAL, data = dados, predictdata=dados.grid, bw = cgwr.bw, kernel = "bisquare", adaptive = TRUE)
#plot do mapa interpolado
dados.grid$C <- cgwr.predict$SDF$prediction
predc <- data.frame(dados.grid$C)
cgwr.predict_plot <- SpatialPixelsDataFrame(points = d3[c("x", "y")], data = predc , proj4string = interpolation_grid@proj4string)
plot(cgwr.predict_plot)



mgwrKriging <- function(){
  
}





