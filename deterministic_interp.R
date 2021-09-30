
library(gstat)    #idw, nn
#valida??o cruzada C - IDW
library(spm) #valida??o cruzada - IDW
library(fields)   #tps
library(car)    #predict for tps
library(interp)   #TIN
library(raster)


#colocar tratamento anisotropia
#procurar por fun?ao que busca automaticamente o melhor idp
#IDW
inverseDistanceWeighted.default <- function(data = "SpatialPointsDataFrame", newdata = c("SpatialPointsDataFrame", "SpatialPixelsDataFrame"),
                                            formula = "formula", handle_anisotropy = TRUE, handle_assimetry = FALSE, rmse_data = NULL, idp = 2){
  # #check and remove assimetry
  # if(isTRUE(handle_assimetry)){
  #   if(is.null(rmse_data)){
  #     normal_distribution <- normalDistribution(data = data, formula = formula)
  #   }
  #   else{
  #     normal_distribution <- normalDistribution(data = rmse_data, formula = formula)
  #   }
  #   if(isFALSE(normal_distribution)){
  #     if(is.null(rmse_data)){
  #       lambda <- boxCoxLambda(formula = formula, data = data)
  #     }
  #     else{
  #       lambda <- boxCoxLambda(formula = formula, data = rmse_data)
  #     }
  #     main_var <- formulaToVector(formula = formula, side = "left")
  #     data@data[, main_var] <- boxCoxTransform(formula = formula, data = data, lambda = lambda, reverseBoxCox = FALSE)
  #   }
  # }
  
  #check and remove anisotropy
  # if(isTRUE(handle_anisotropy)){
  #   if(is.null(rmse_data)){
  #     anisotropy <- checkAnisotropy(data, formula)
  #   }
  #   else{
  #     anisotropy <- checkAnisotropy(rmse_data, formula)
  #   }
  #   data <- handleAnisotropy(data = data, formula = formula, anisotropy = anisotropy)
  #   newdata <- handleAnisotropy(data = newdata, formula = formula, anisotropy = anisotropy)
  # }
  
  #perform idw interpolation
  newdata@data$result <- idw(formula = formula, locations = data, newdata = newdata, idp = idp)$var1.pred
  
  #rotates coordinates back to anisotropic
  # if(isTRUE(handle_anisotropy)){
  #   data <- handleAnisotropy(data = data, formula = formula, anisotropy = anisotropy, reverse = TRUE)
  #   newdata <- handleAnisotropy(data = newdata, formula = formula, anisotropy = anisotropy, reverse = TRUE)
  # }
  
  # #back-transforms data to assimetry
  # if(isTRUE(handle_assimetry)){
  #   if(isFALSE(normal_distribution)){
  #     data@data[, main_var] <- boxCoxTransform(formula = formula, data = data, lambda = lambda, reverseBoxCox = TRUE)
  #     newdata@data$result <- boxCoxTransform(formula = result~1, data = newdata, lambda = lambda, reverseBoxCox = TRUE)
  #   }
  # }
  
  result <- newdata@data$result
  
  return(result)
}

#Maciej Tomczak 1998
#Comparison of different spatial interpolation methods for historical hydrographic data of the lowermost Mississippi River 2019
eidw.default <- function(data = "SpatialPointsDataFrame", newdata = c("SpatialPointsDataFrame", "SpatialPixelsDataFrame"), formula = "formula", handle_anisotropy = TRUE,
                         handle_assimetry = FALSE, rmse_data = NULL, idp = 2, smooth = 0, n = 10){
  
  if(isTRUE(handle_anisotropy)){
    if(is.null(rmse_data)){
      anisotropy <- checkAnisotropy(data, formula)
    }
    else{
      anisotropy <- checkAnisotropy(rmse_data, formula)
    }
  }
  #anisotropy$direction is angle, anisotropy$ratio is ratio
  #values for EIDW
  txx <- cos(anisotropy$direction)/anisotropy$ratio
  txy <- sin(anisotropy$direction)/anisotropy$ratio
  tyx <- -sin(anisotropy$direction)
  tyy <- cos(anisotropy$direction)
  axx <- txx^2 + tyx^2
  axy <- 2*(txx*txy + tyx*tyy)
  ayy <- tyy^2 + txy^2
  
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
      weight_df$weight[i] <- sqrt(axx*delta_x^2 + axy*delta_x*delta_y + ayy*delta_y^2)^idp
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
  return(result)
}

#verificar o nmax
nearestNeighbor.default <- function(data = "SpatialPointsDataFrame", newdata = c("SpatialPointsDataFrame", "SpatialPixelsDataFrame"), formula = "formula"){
  result <- idw(formula, data, data@coords, newdata = newdata, idp = 0, nmax=1)$var1.pred
  return(result)
}



naturalNeighbor.default <- function(){
  
}


#arrumar
triangulation.default <- function(data = "SpatialPointsDataFrame",  newdata = c("SpatialPointsDataFrame", "SpatialPixelsDataFrame"), formula = "formula"){
  newdata <- spPixelsToRaster(newdata)
  column_name <- formulaToVector(formula, "left")
  fit <- interp( # using {interp}
    x = data@coords[,1],           # the function actually accepts coordinate vectors
    y = data@coords[,2],
    z = data@data[, column_name],
    xo = newdata@coords[, 1],     # here we already define the target grid
    yo = newdata@coords[, 2]
    #output = "points"
  ) #%>% bind_cols()
  #interpolated_data <- raster::rasterFromXYZ(fit_TIN, crs = crs_raster_format)
  return(fit)
}



spline.default <- function(data = "SpatialPointsDataFrame", newdata = c("SpatialPointsDataFrame", "SpatialPixelsDataFrame"), formula = "formula", handle_anisotropy = TRUE,
                           handle_assimetry = FALSE, rmse_data = NULL, idp = NULL){
  newdata <- spPixelsToRaster(newdata)
  column_name <- formulaToVector(formula, "left")
  fit <- Tps(data@coords, data@data[, column_name], miles = FALSE, lon.lat = TRUE)
  interpolated_data <- interpolate(newdata, fit)
  return(rasterToSpPixels(interpolated_data)$layer)
}







