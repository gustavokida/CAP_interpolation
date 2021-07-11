
library(gstat)    #idw, nn
#validação cruzada C - IDW
library(spm) #validação cruzada - IDW
library(fields)   #tps
library(car)    #predict for tps
library(interp)   #TIN
library(raster)


#colocar tratamento anisotropia
#procurar por funçao que busca automaticamente o melhor idp
#IDW
inverseDistanceWeighted <- function(data = "SpatialPointsDataFrame", newdata = c("SpatialPointsDataFrame", "SpatialPixelsDataFrame"), formula = "formula"){
  invisible(capture.output(result <- idw(formula, data, data@coords, newdata = newdata, idp = 2.0)$var1.pred))
  return(result)
}



#verificar o nmax
nearestNeighbor <- function(data = "SpatialPointsDataFrame", newdata = c("SpatialPointsDataFrame", "SpatialPixelsDataFrame"), formula = "formula"){
  result <- idw(formula, data, data@coords, newdata = newdata, idp = 0, nmax=1)$var1.pred
  return(result)
}



naturalNeighbor <- function(){
  
}


#arrumar
triangulation <- function(data = "SpatialPointsDataFrame",  newdata = c("SpatialPointsDataFrame", "SpatialPixelsDataFrame"), formula = "formula"){
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



spline <- function(data = "SpatialPointsDataFrame", newdata = c("SpatialPointsDataFrame", "SpatialPixelsDataFrame"), formula = "formula"){
  newdata <- spPixelsToRaster(newdata)
  column_name <- formulaToVector(formula, "left")
  fit <- Tps(data@coords, data@data[, column_name], miles = FALSE, lon.lat = TRUE)
  interpolated_data <- interpolate(newdata, fit)
  return(rasterToSpPixels(interpolated_data)$layer)
}







