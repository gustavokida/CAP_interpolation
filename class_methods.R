#CLASS METHODS


removeOutlier.cap <- function(object){
  column_names <- formulaToVector(object$formula, side = "all")
  object$data <- removeOutlier(data = object$data, column_names = column_names)
  return(object)
}



checkQuantity.cap <- function(object){
  column_name <- formulaToVector(object$formula, side = "left")
  object$point_quantity <- checkQuantity(data = object$data, column_name = column_name)
}



boxCoxLambda.cap <- function(object){
  left_var <- formulaToVector(object$formula, "left")
  object$lambda <- BoxCox.lambda(object$data@data[, left_var], method = 'loglik')
  return(object)
}



boxCoxTransform.cap <- function(object){
  left_var <- formulaToVector(object$formula, "left")
  object$data@data[, left_var] <- boxCoxTransformation(object$formula, object$data, object$lambda, object$reverse_boxcox)
  if(isFALSE(object$reverse_boxcox)){
    object$reverse_boxcox <- TRUE
  }
  else{
    object$reverse_boxcox <- FALSE
  }
  return(object)
}

normalDistribution.cap <- function(object){
  object$normal_distribution <- normalDistribution(data = object$data, formula = object$formula)
  return(object)
}



checkAnisotropy.cap <- function(object = cap()){
  object$anisotropy <- checkAnisotropy(object$data, object$formula)
  return(object)
}



handleAnisotropy.cap <- function(object){
  if(is.empty(object$anisotropy)){
    object <- checkAnisotropy(data, formula)
  }
  if(is.null(object$anisotropy)){
    return(object)
  }
  else{
    rotated_coords <- coords.aniso(object$data@coords, object$anisotropy, reverse = object$reverse_anisotropy)
    colnames(rotated_coords) <- c("x", "y")
    object$data@coords <- rotated_coords
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
  }



checkCovariates.cap <- function(object){
  left_side <- formulaToVector(object$formula, "left")
  right_side <- formulaToVector(object$formula, "right")
  covariates <- checkCovariates(data = object$data, covariate_data = object$covariate_data, main_attribute_column = left_side, column_names = right_side)
  new_formula <- makeFormula(main_attribute_column = left_side, column_names = covariates)
}


  
  
  
  
  
  