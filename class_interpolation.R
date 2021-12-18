

autoInterpolation.cap <- function(object){
  object <- removeOutlier(object)
  object <- checkQuantity(object)
  object <- checkAnisotropy(object)
  if(isFALSE(object$point_quantity)){
    if(isFALSE(object$anisotropy)){
      object <- inverseDistanceWeighted(object)
    }
    else{
      object <- eidw(object)
    }
  }
  else{
    if(isFALSE(object$trend)){
      if(is.null(object$covariate_data)){
        if(!isFALSE(object$anisotropy)){
          object <- checkQuantity(object, aniso = TRUE)
        }
        if(isFALSE(object$point_quantity)){
          object <- eidw(object)
        }
        else{
          object <- ordinaryKriging(object)
        }
      }
      else{
        object <- coKriging(object)
      }
    }
    else{
      if(is.null(object$covariate_data)){
        object <- universalKriging(object)
      }
      else{
        object <- regressionKriging(object)
      }
    }
  }
  return(object)
}


autoInterpolation.cap <- function(object){
  object <- removeOutlier(object)
  object <- checkQuantity(object)
  object <- checkAnisotropy(object)
  if(isFALSE(object$point_quantity)){
    if(isFALSE(object$anisotropy)){
      object <- inverseDistanceWeighted(object)
    }
    else{
      object <- eidw(object)
    }
  }
  else{
    if(isFALSE(object$trend)){
      if(!isFALSE(object$anisotropy)){
        object <- checkQuantity(object, aniso = TRUE)
      }
      if(isFALSE(object$point_quantity)){
        object <- eidw(object)
      }
      else{
        object <- ordinaryKriging(object)
      }
    }
    else{
      object <- universalKriging(object)
    }
  }
  return(object)
}



#kriging
ordinaryKriging.cap <- function(object){
  object$interpolation_function <- ordinaryKriging
  
  #name of the interpolated output
  main_var <- formulaToVector(formula = object$formula, side = "left")
  output_name <- paste("OK", main_var, sep = "_")
  object$formula_output = expr(!!as.symbol(output_name)~1)
  
  #check and remove assimetry
  if(isTRUE(object$handle_assimetry)){
    object <- normalDistribution(object = object)
    if(isFALSE(object$normal_distribution)){
      object <- boxCoxLambda(object = object)
      object <- boxCoxTransform(object = object)
    }
  }
  
  #check and remove anisotropy
  if(isTRUE(object$handle_anisotropy)){
    object <- checkAnisotropy(object = object)
    data_variogram <- handleAnisotropy(data = object$data, formula = object$formula, anisotropy = object$anisotropy)
  }
  else{
    data_variogram <- object$data
  }

  
  predicted_data <- ordinaryKriging(data = object$data, newdata = object$newdata, formula = object$formula, data_variogram = data_variogram, 
                                    handle_anisotropy = FALSE, handle_assimetry = FALSE, cap_return = TRUE)
  
  object$variogram$exp_var <- predicted_data$exp_var
  object$variogram$var_model <- predicted_data$var_model
  object$newdata@data[, output_name] <- predicted_data$krige_output$var1.pred
  
  #back-transforms data to assimetry
  if(isTRUE(object$handle_assimetry)){
    if(isFALSE(object$normal_distribution)){
      object <- boxCoxTransform(object = object, formula_output = TRUE)
      object <- boxCoxTransform(object = object)
    }
  }
  
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
  main_var <- formulaToVector(formula = object$formula, side = "left")
  output_name <- paste("UK", main_var, sep = "_")
  object$formula_output = expr(!!as.symbol(output_name)~1)

  predicted_data <- universalKriging(data = object$data, newdata = object$newdata, formula = object$formula, cap_return = TRUE)
  
  object$variogram$exp_var <- predicted_data$exp_var
  object$variogram$var_model <- predicted_data$var_model
  object$newdata@data[, output_name] <- predicted_data$krige_output$var1.pred
  
  return(object)
}




regressionKriging.cap <- function(object){
  object$interpolation_function <- regressionKriging
  
  output_name <- paste("RK", left_side, sep = "_")
  object$newdata@data[, output_name] <- predicted_data
  return(object)
}




#deterministic

inverseDistanceWeighted.cap <- function(object){
  object$interpolation_function <- inverseDistanceWeighted
  
  #name of the interpolated output
  main_var <- formulaToVector(formula = object$formula, side = "left")
  output_name <- paste("IDW", main_var, sep = "_")
  object$formula_output = expr(!!as.symbol(output_name)~1)
  
  # #check and remove assimetry
  # if(isTRUE(object$handle_assimetry)){
  #   object <- normalDistribution(object = object)
  #   if(isFALSE(object$normal_distribution)){
  #     object <- boxCoxLambda(object = object)
  #     object <- boxCoxTransform(object = object)
  #   }
  # }
  
  #check and remove anisotropy
  if(isTRUE(object$handle_anisotropy)){
    object <- checkAnisotropy(object = object)
    object <- handleAnisotropy(object = object)
    #object$newdata <- handleAnisotropy(data = object$newdata, formula = object$formula, anisotropy = object$anisotropy)
  }
  
  #object <- estimateIdp(object)
  
  #perform idw interpolation
  object$newdata@data[, output_name] <- inverseDistanceWeighted(data = object$data, newdata = object$newdata, formula = object$formula,
                                                                handle_anisotropy = FALSE, handle_assimetry = FALSE, idp=object$idp)
  #rotates coordinates back to anisotropic
  if(isTRUE(object$handle_anisotropy)){
    object <- handleAnisotropy(object = object)
    #object$newdata <- handleAnisotropy(data = object$newdata, formula = object$formula, anisotropy = object$anisotropy, reverse = TRUE)
  }

  # #back-transforms data to assimetry
  # if(isTRUE(object$handle_assimetry)){
  #   if(isFALSE(object$normal_distribution)){
  #     object <- boxCoxTransform(object = object, formula_output = TRUE)
  #     object <- boxCoxTransform(object = object)
  #   }
  # }
  
  return(object)
}



eidw.cap <- function(object){
  object$interpolation_function <- eidw
  main_var <- formulaToVector(formula = object$formula, side = "left")
  output_name <- paste("EIDW", main_var, sep = "_")
  object$formula_output = expr(!!as.symbol(output_name)~1)
  object <- checkAnisotropy(object)
  object <- estimateIdp(object)
  
  #Elliptical Inverse distance weighted
  object$newdata@data[, output_name] <- eidw(data = object$data, newdata = object$newdata, formula = object$formula, idp = object$idp)
  
  return(object)
}


nearestNeighbor.cap <- function(object){
  object$interpolation_function <- nearestNeighbor
  
  output_name <- paste("NN", left_side, sep = "_")
  object$newdata@data[, output_name] <- predicted_data
}



spline.cap <- function(object){
  object$interpolation_function <- spline
  main_var <- formulaToVector(formula = object$formula, side = "left")
  output_name <- paste("Spline", main_var, sep = "_")
  object$formula_output = expr(!!as.symbol(output_name)~1)
  
  #perform spline interpolation
  object$newdata@data[, output_name] <- spline(data = object$data, newdata = object$newdata, formula = object$formula)
  
  return(object)
}



triangulation.cap <- function(object){
  
  
  output_name <- paste("Triangulation", left_side, sep = "_")
  object$newdata@data[, output_name] <- predicted_data
}



naturalNeighbor.cap <- function(object){
  
  
  output_name <- paste("NaturalNeighbor", left_side, sep = "_")
  object$newdata@data[, output_name] <- predicted_data
}



