rfillspgaps <- function(rasterLayer,
                      maskPol=NULL,
                      nmax =25,
                      ...
                      ){

  if(class(rasterLayer)=="RasterLayer") { r=rasterLayer} else{r= raster(rasterLayer) }
  crs= r@crs
  zcol= r@data@names
  
  if(!is.null(maskPol)){
    rr = extract(  r , maskPol, cellnumbers=T, df=T)
    rr[is.na(rr[ ,3]),3] <- c(-9999)
    r[rr$cell] <-  rr[ ,3]
  }else{
    r[is.na(r)] <- c(-9999)
  }
rr=NULL
  
  r=rasterToPoints(r, fun=NULL, spatial=TRUE)
  names(r)[1]= zcol

data <- r[(r@data[,1]!=c(-9999)),]
data@proj4string <- CRS(as.character(NA))
newdata=r[r@data[,1]==c(-9999),]
newdata@proj4string <- CRS(as.character(NA))

 gaps = idw(as.formula(paste(zcol,"~1",sep="")), data,newdata ,nmax =nmax,...)
 gaps=gaps[,1]
 names(gaps)= zcol
 r[(r@data[,zcol]==c(-9999)),] <- gaps@data[,1]
 
 r= rasterFromXYZ(as.data.frame(r[,zcol]) , crs=crs)
 
 if(class(rasterLayer)!="RasterLayer"){
   r= as(r,class(rasterLayer))
 }
  return(r)
}
  
