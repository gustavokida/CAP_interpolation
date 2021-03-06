library(spm)    #idwcv
library(utils)  #captureOutput


#RMSE   remover 1 ponto do observed e tentar predizer nas coordenadas do proprio observed(62)
rmse_test <- function(data = "SpatialPointsDataFrame", formula = "formula", funcInterpolation = "function"){
  column_name <- formulaToVector(formula, "left")
  right_side <- formulaToVector(formula, "right")
  all_points_predicted <- NULL
  for(i in 1:length(data)){
    new_formula = makeFormula(column_name, right_side)
    point_to_predict <- data[i,]
    data_without_point <- data[-i,]
    all_points_predicted[i] <- funcInterpolation(data_without_point, point_to_predict, new_formula)
  }
  return(sqrt(mean((data@data[, column_name] - all_points_predicted)^2)))
}


#ARRUMAR PARA ACEITAR COVARIAVEIS
#verificar colunas OK e UK
rmse <- function(data = "SpatialPointsDataFrame", formula = "formula", funcInterpolation = "function"){
  column_name <- formulaToVector(formula, "left")
  right_side <- formulaToVector(formula, "right")
  all_points_predicted <- NULL
  for(i in 1:length(data)){
    new_formula = makeFormula(column_name, right_side)
    point_to_predict <- data
    data_without_point <- data[-i,]
    points_predicted <- funcInterpolation(data_without_point, point_to_predict, new_formula)
    all_points_predicted[i] <- points_predicted[i]
  }
  return(sqrt(mean((data@data[, column_name] - all_points_predicted)^2)))
}


#use same model in cross validation.    reference: A tutorial guide to geostatistics computing and modelling variograms and kriging
rmse.default <- function(data = "SpatialPointsDataFrame", formula = "formula", funcInterpolation = "function",
                         covariate_data = NULL, handle_anisotropy = TRUE, handle_assimetry = TRUE, idp = 2){
  column_name <- formulaToVector(formula, "left")
  right_side <- formulaToVector(formula, "right")
  all_points_predicted <- NULL
  
  if(is.null(covariate_data)){
    for(i in 1:length(data)){
      new_formula <- makeFormula(column_name, right_side)
      point_to_predict <- data
      data_without_point <- data[-i,]
      capture.output(points_predicted <- funcInterpolation(data = data_without_point, newdata = point_to_predict,
                                                           formula = new_formula, handle_anisotropy = handle_anisotropy,
                                                           handle_assimetry = handle_assimetry, rmse_data = data, idp = idp))
      all_points_predicted[i] <- points_predicted[i]
    }
  }
  else{
    for(i in 1:length(data)){
      new_formula <- makeFormula(column_name, right_side)
      point_to_predict <- data
      data_without_point <- data[-i,]
      capture.output(points_predicted <- funcInterpolation(data = data_without_point, newdata = point_to_predict,
                                                           formula = new_formula, covariate_data = covariate_data,
                                                           handle_anisotropy = handle_anisotropy,
                                                           handle_assimetry = handle_assimetry,
                                                           rmse_data = data, idp = idp))
      all_points_predicted[i] <- points_predicted[i]
    }
  }
  return(sqrt(mean((data@data[, column_name] - all_points_predicted)^2)))
}


#IMPLEMENTAR CALCULO DE VARIANCIA DO MAPA INTERPOLADO POR KRIGAGEM



#CROSS VALIDATION
rmse.cap <- function(object){
  stopifnot(!is.null(object$interpolation_function))

  object$rmse <- rmse(data = object$data, formula = object$formula, funcInterpolation = object$interpolation_function, covariate_data = object$covariate_data,
       handle_anisotropy = object$handle_anisotropy, handle_assimetry = object$handle_assimetry, idp = object$idp)
  
  return(object)
}


