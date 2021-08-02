#kriging
ordinaryKriging.cap <- function(object){
  
  
  output_name <- paste("OK", left_side, sep = "_")
  object$newdata@data[, output_name] <- predicted_data
}



coKriging.cap <- function(object){
  
  
  output_name <- paste("CK", left_side, sep = "_")
  object$newdata@data[, output_name] <- predicted_data
}



universalKriging.cap <- function(object){
  
  
  output_name <- paste("UK", left_side, sep = "_")
  object$newdata@data[, output_name] <- predicted_data
}



regressionKriging.cap <- function(object){
 
  
  output_name <- paste("RK", left_side, sep = "_")
  object$newdata@data[, output_name] <- predicted_data 
}



#deterministic

inverseDistanceWeighted.cap <- function(object){
 
  
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


