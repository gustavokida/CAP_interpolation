
#OBJECT S4
cap <- setClass("cap", slots = list(data = "SpatialPointsDataFrame", newdata = "SpatialPixelsDataFrame", formula = "formula",
                                    reverseBoxCox = "logical", trend = "logical", reverseAnisotropy = "logical", pointQuantity = "logical", lambda = "numeric",
                                    anisotropy = "numeric", rmse = c(NULL, "numeric")))


#OBJECT S3
cap <- function(data = "SpatialPointsDataFrame", newdata = "SpatialPixelsDataFrame", formula = "formula"){
  object <- list(data = data, newdata = newdata, formula = formula,
                  reverseBoxCox = FALSE, trend = FALSE, reverseAnisotropy = FALSE, pointQuantity = FALSE, lambda = numeric(0),
                  anisotropy = numeric(0), rmse = numeric(0))
  class(object) <- "cap"
  return(object)
}
  


#METHODS
boxCoxLambda <- function(.object){
  left_var <- formulaToVector(.object$formula, "left")
  .object$lambda <- BoxCox.lambda(.object$data@data[, left_var], method = 'loglik')
  return(.object)
}



boxCoxTransform <- function(.object){
  left_var <- formulaToVector(.object$formula, "left")
  .object$data@data[, left_var] <- boxCoxTransformation(.object$formula, .object$data, .object$lambda, .object$reverseBoxCox)
  if(isFALSE(.object$reverseBoxCox)){
    .object$reverseBoxCox <- TRUE
  }
  else{
    .object$reverseBoxCox <- FALSE
  }
  return(.object)
}



anisotropy <- function(.object){
  .object$anisotropy <- checkAnisotropy(.object$data, .object$formula)
  return(.object)
}
