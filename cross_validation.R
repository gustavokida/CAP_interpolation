library(spm)    #idwcv


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



rmse.default <- function(data = "SpatialPointsDataFrame", formula = "formula", funcInterpolation = "function", covariate_data = NULL){
  column_name <- formulaToVector(formula, "left")
  right_side <- formulaToVector(formula, "right")
  all_points_predicted <- NULL
  if(is.null(covariate_data)){
    for(i in 1:length(data)){
      new_formula <- makeFormula(column_name, right_side)
      point_to_predict <- data
      data_without_point <- data[-i,]
      points_predicted <- funcInterpolation(data_without_point, point_to_predict, new_formula)
      all_points_predicted[i] <- points_predicted[i]
    }
  }
  else{
    for(i in 1:length(data)){
      new_formula <- makeFormula(column_name, right_side)
      point_to_predict <- data
      data_without_point <- data[-i,]
      points_predicted <- funcInterpolation(data_without_point, point_to_predict, new_formula, covariate_data)
      all_points_predicted[i] <- points_predicted[i]
    }
  }
  return(sqrt(mean((data@data[, column_name] - all_points_predicted)^2)))
}


#IMPLEMENTAR CALCULO DE VARIANCIA DO MAPA INTERPOLADO POR KRIGAGEM


