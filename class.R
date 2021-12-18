
#OBJECT S3
cap <- function(object = NULL, formula = NULL, data = NULL, newdata = NULL, covariate_data = NULL){
  stopifnot(("cap" %in% is(object) | is.null(object)))
  stopifnot(("formula" %in% is(formula) | is.null(formula)))
  stopifnot(("SpatialPointsDataFrame" %in% is(data) | is.null(data)))
  stopifnot(("SpatialPixelsDataFrame" %in% is(newdata) | is.null(newdata)))
  stopifnot(("SpatialPointsDataFrame" %in% is(covariate_data) | is.null(covariate_data)))
  
  if(is.null(object)){
    object <- list(formula = formula, data = data, newdata = newdata, covariate_data = covariate_data,
                   formula_output = NULL, trend = FALSE, point_quantity = FALSE, handle_assimetry = TRUE,
                   normal_distribution = FALSE, lambda = numeric(0), reverse_boxcox = FALSE,
                   handle_anisotropy = TRUE, anisotropy = numeric(0), reverse_anisotropy = FALSE,
                   rmse = numeric(0), interpolation_function = NULL, variogram = NULL, idp = 2)
    class(object) <- "cap"
  }
  else{
    if(!is.null(formula)){
      object$formula <- formula
    }
    if(!is.null(data)){
      object$data <- data
    }
    if(!is.null(newdata)){
      object$newdata <- newdata
    }
    if(!is.null(covariate_data)){
      object$covariate_data <- covariate_data
    }
  }
  return(object)
}


