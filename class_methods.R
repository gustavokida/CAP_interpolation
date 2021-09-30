library(rlist)

#CLASS METHODS


removeOutlier.cap <- function(object){
  column_names <- formulaToVector(object$formula, side = "all")
  object$data <- removeOutlier(data = object$data, column_names = column_names)
  return(object)
}



checkQuantity.cap <- function(object){
  column_name <- formulaToVector(object$formula, side = "left")
  object$point_quantity <- checkQuantity(data = object$data, column_name = column_name)
  return(object)
}



boxCoxLambda.cap <- function(object){
  left_var <- formulaToVector(object$formula, "left")
  object$lambda <- BoxCox.lambda(object$data@data[, left_var], method = 'loglik')
  return(object)
}


#if formula_output equals TRUE, it will transform the main variable from formula and the main variable from formula_output
boxCoxTransform.cap <- function(object, formula_output = FALSE){
  stopifnot(!is.null(object$lambda))
  
  #transforms data from formula_output, which is the predicted data on the newdata
  if(isTRUE(formula_output)){
    output_name <- formulaToVector(object$formula_output, "left")
    object$newdata@data[, output_name] <- boxCoxTransformation(object$formula_output, object$newdata, object$lambda, object$reverse_boxcox)
  }
  else{
    #transforms the data from main variable
    output_name <- formulaToVector(object$formula, "left")
    object$data@data[, output_name] <- boxCoxTransformation(object$formula, object$data, object$lambda, object$reverse_boxcox)
  
    #after done the rotation, it will change the value to signalize if the data is rotated(TRUE) or not(FALSE).
    if(isFALSE(object$reverse_boxcox)){
      object$reverse_boxcox <- TRUE
    }
    else{
      object$reverse_boxcox <- FALSE
    }
  }
  return(object)
}

normalDistribution.cap <- function(object){
  object$normal_distribution <- normalDistribution(data = object$data, formula = object$formula)
  return(object)
}



estimateIdp.cap <- function(object){
  idp_range <- seq(0.5, 4, 0.5)
  rmse_values <- rep(NA, length(idp_range))
  for(i in seq(along=idp_range)){
    object$idp <- idp_range[i]
    rmse_values[i] <- rmse(object=object)$rmse
  }
  best <- which(rmse_values == min(rmse_values))[1]
  object$idp <- idp_range[best]
  return(object)
}



checkAnisotropy.cap <- function(object = cap()){
  object$anisotropy <- checkAnisotropy(object$data, object$formula)
  return(object)
}


#object$reverse_anisotropy = FALSE - original coords
#object$reverse_anisotropy = TRUE - rotated coords
handleAnisotropy.cap <- function(object){
  if(is.empty(object$anisotropy)){
    object <- checkAnisotropy(object = object)
  }
  if(isFALSE(object$anisotropy)){
    return(object)
  }
  else{
    
    #rotated_coords <- coords.aniso(object$data@coords, object$anisotropy, reverse = object$reverse_anisotropy)
    rotated_coords <- rotateAnisotropicData(object$data, object$anisotropy)
    #colnames(rotated_coords) <- c("x", "y")
    #object$data@coords <- rotated_coords
    object$data <- rotated_coords
    #rotated_newdata_coords <- coords.aniso(object$newdata@coords, object$anisotropy, reverse = object$reverse_anisotropy)
    #colnames(rotated_newdata_coords) <- c("x", "y")
    #$newdata@coords <- rotated_newdata_coords
    
    #object$reverse_anisotropy <- !object$reverse_anisotropy
    
    return(object)
  }
}



detrendFormula.cap <- function(object){
  object$formula <- detrendFormula(data = object$data, formula = object$formula)
  return(object)
}
  
  

checkCovariatesCKO.cap <- function(object){
  left_side <- formulaToVector(object$formula, "left")
  right_side <- formulaToVector(object$formula, "right")
  covariates <- checkCovariatesCKO(data = object$data, covariate_data = object$covariate_data, main_attribute_column = left_side, column_names = right_side)
  new_formula <- makeFormula(main_attribute_column = left_side, column_names = covariates)
  object$formula <- new_formula
  return(object)
  }



checkCovariates.cap <- function(object){
  left_side <- formulaToVector(object$formula, "left")
  right_side <- formulaToVector(object$formula, "right")
  covariates <- checkCovariates(data = object$data, covariate_data = object$covariate_data, main_attribute_column = left_side, column_names = right_side)
  new_formula <- makeFormula(main_attribute_column = left_side, column_names = covariates)
  object$formula <- new_formula
  return(object)
}



plot.cap <- function(object){
  return_list <- list()
  if(!is.null(object$formula_output)){
    name <- formulaToVector(object$formula_output, "left")
    map_plot <- spplot(object$newdata, zcol = name, sub = list(font = 1, cex = 1, label = name))
    return_list <- list.append(return_list, map_plot)
    
    if(!is.null(object$variogram)){
      vario_plot <- plot(object$variogram$exp_var, object$variogram$var_model, sub = list(font = 1, cex = 1, label = name))
      return_list <- list.append(return_list, vario_plot)
    }
    
  }
  if(!is.empty(object$rmse)){
    return_list <- list.append(return_list, object$rmse)
  }
  return(return_list)
}

  
  
  
  
  
  