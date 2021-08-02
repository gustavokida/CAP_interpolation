
#OBJECT S4
cap <- setClass("cap", slots = list(data = "SpatialPointsDataFrame", newdata = "SpatialPixelsDataFrame", formula = "formula",
                                    reverse_boxcox = "logical", trend = "logical", reverse_anisotropy = "logical", point_quantity = "logical",
                                    normal_distribution = "logical", lambda = "numeric", anisotropy = "numeric", rmse = c(NULL, "numeric")))


#OBJECT S3
cap <- function(object = NULL, formula = NULL, data = NULL, newdata = NULL, covariate_data = NULL){
  stopifnot(("cap" %in% is(object) | is.null(object)))
  stopifnot(("formula" %in% is(formula) | is.null(formula)))
  stopifnot(("SpatialPointsDataFrame" %in% is(data) | is.null(data)))
  stopifnot(("SpatialPixelsDataFrame" %in% is(newdata) | is.null(newdata)))
  stopifnot(("SpatialPointsDataFrame" %in% is(covariate_data) | is.null(covariate_data)))
  
  if(is.null(object)){
    object <- list(formula = formula, data = data, newdata = newdata, covariate_data = covariate_data,
                   reverse_boxcox = FALSE, trend = FALSE, reverse_anisotropy = FALSE, point_quantity = FALSE, normal_distribution = FALSE,
                    lambda = numeric(0), anisotropy = numeric(0), rmse = numeric(0))
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


