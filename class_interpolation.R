#kriging
ordinaryKriging.cap <- function(object){
  interpolation_function <- ordinaryKriging
  predicted_data <- ordinaryKriging(data = object$data, newdata = object$newdata, formula = object$formula, cap_return = TRUE, handle_anisotropy = object$handle_anisotropy)
  object$variogram$exp_var <- predicted_data$exp_var
  object$variogram$var_model <- predicted_data$var_model
  output_name <- paste("OK", left_side, sep = "_")
  object$newdata@data[, output_name] <- predicted_data$krige_output$var1.pred
  return(object)
}



coKriging.cap <- function(object){
  interpolation_function <- coKriging
  
  output_name <- paste("CK", left_side, sep = "_")
  object$newdata@data[, output_name] <- predicted_data
  return(object)
}




universalKriging.cap <- function(object){
  interpolation_function <- universalKriging
  predicted_data <- universalKriging(data = object$data, newdata = object$newdata, formula = object$formula, cap_return = TRUE)
  output_name <- paste("UK", left_side, sep = "_")
  object$newdata@data[, output_name] <- predicted_data$krige_output$var1.pred
  return(object)
}




regressionKriging.cap <- function(object){
 
  
  output_name <- paste("RK", left_side, sep = "_")
  object$newdata@data[, output_name] <- predicted_data
  return(object)
}



plot.cap <- function(object){
  if(is.null(object$variogram)){
    stop("No variogram found.")
  }
  name <- formulaToVector(object$formula, "left")
  plot(object$variogram$exp_var, object$variogram$var_model, sub = list(font = 1, cex = 1, label = name))
}



#deterministic

inverseDistanceWeighted.cap <- function(object){
  predicted_data <- inverseDistanceWeighted(data = object$data, newdata = object$newdata, formula = object$formula, handle_anisotropy = object$handle_anisotropy)
  output_name <- paste("IDW", left_side, sep = "_")
  object$newdata@data[, output_name] <- predicted_data 
}



nearestNeighbor.cap <- function(object){
  
  
  output_name <- paste("NN", left_side, sep = "_")
  object$newdata@data[, output_name] <- predicted_data
}



spline.cap <- function(object){
  
  
  output_name <- paste("Spline", left_side, sep = "_")
  object$newdata@data[, output_name] <- predicted_data
}



triangulation.cap <- function(object){
  
  
  output_name <- paste("Triangulation", left_side, sep = "_")
  object$newdata@data[, output_name] <- predicted_data
}



naturalNeighbor.cap <- function(object){
  
  
  output_name <- paste("NaturalNeighbor", left_side, sep = "_")
  object$newdata@data[, output_name] <- predicted_data
}




#CROSS VALIDATION
rmse.cap <- function(object){
  
  
}

