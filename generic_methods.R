
#functions
removeOutlier <- function(object, ...) UseMethod("removeOutlier")

removeNA <- function(object, ...) UseMethod("removeNA")

checkQuantity<- function(object, ...) UseMethod("checkQuantity")

boxCoxLambda <- function(object, ...) UseMethod("boxCoxLambda")

boxCoxTransform <- function(object, ...) UseMethod("boxCoxTransform")

normalDistribution <- function(object, ...) UseMethod("normalDistribution")

checkAnisotropy <- function(object, ...) UseMethod("checkAnisotropy")

handleAnisotropy <- function(object, ...) UseMethod("handleAnisotropy")

detrendFormula <- function(object, ...) UseMethod("detrendFormula")

checkCovariatesCKO <- function(object, ...) UseMethod("checkCovariatesCKO")

checkCovariates <- function(object, ...) UseMethod("checkCovariates")

estimateIdp <- function(object, ...) UseMethod("estimateIdp")


#validation
rmse <- function(object, ...) UseMethod("rmse")


#interpolation
autoInterpolation <- function(object, ...) UseMethod("autoInterpolation")

#kriging
ordinaryKriging <- function(object, ...) UseMethod("ordinaryKriging")

coKriging <- function(object, ...) UseMethod("coKriging")

universalKriging <- function(object, ...) UseMethod("universalKriging")

regressionKriging <- function(object, ...) UseMethod("regressionKriging")

#deterministic

inverseDistanceWeighted <- function(object, ...) UseMethod("inverseDistanceWeighted")

eidw <- function(object, ...) UseMethod("eidw")

nearestNeighbor <- function(object, ...) UseMethod("nearestNeighbor")

spline <- function(object, ...) UseMethod("spline")

triangulation <- function(object, ...) UseMethod("triangulation")

naturalNeighbor <- function(object, ...) UseMethod("naturalNeighbor")


