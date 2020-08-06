### Harvest data.frame of climate space
require(raster)
require(ncdf4)
require(here)


dpath1 <- '/ESS_EarthObs/CLIMATE_DATA/CRU/cru_ts4.00/'
tpath <- 'dataProcessing/climData/'
dir.create(path = tpath, showWarnings = F, recursive = T)

# years of interest
yrS <- 2008; yrE <- 2012 #

# input data
tmp <- stack(paste0(dpath1,'cru_ts4.00.1901.2015.tmp.dat.nc'))
ppt <- stack(paste0(dpath1,'cru_ts4.00.1901.2015.pre.dat.nc'))

# select in time 
tmp <- tmp[[(12*(yrS - 1901)+1):(12*(yrE - 1901)+12)]]
#tmp <- aggregate(tmp,fact=2)
ppt <- ppt[[(12*(yrS - 1901)+1):(12*(yrE - 1901)+12)]]
# ppt <- aggregate(ppt,fact=2)

# prepare empty stacks
TMPmean.yr <- stack()
TMPrange.yr <- stack()
PPTsum.yr <- stack()
AIndex.yr <- stack()
SNWsum.yr <- stack()

for(iT in 1:5){
  z1 <- raster::calc(tmp[[(1:12)+12*(iT-1)]],fun=mean,na.rm=T)
  TMPmean.yr <- addLayer(TMPmean.yr,z1)
  z2 <- raster::calc(tmp[[(1:12)+12*(iT-1)]],fun=function(x){(max(x)-min(x))})
  TMPrange.yr <- addLayer(TMPrange.yr,z2)
  z3 <- raster::calc(ppt[[(1:12)+12*(iT-1)]],fun=sum,na.rm=T)
  PPTsum.yr <- addLayer(PPTsum.yr,z3)
  z4 <- z3/(z1+30)
  AIndex.yr <- addLayer(AIndex.yr,z4)
  
  subZtemp <- tmp[[(1:12)+12*(iT-1)]]<0
  snowppt <- ppt[[(1:12)+12*(iT-1)]] * subZtemp
  z5 <- raster::calc(snowppt,fun=sum,na.rm=T)
  SNWsum.yr <- addLayer(SNWsum.yr,z5)
  
}


TMPmean.5yr <- raster::calc(TMPmean.yr,fun=median,na.rm=T)
TMPrange.5yr <- raster::calc(TMPrange.yr,fun=median,na.rm=T)
PPTsum.5yr <- raster::calc(PPTsum.yr,fun=median,na.rm=T)
AIndex.5yr <- raster::calc(AIndex.yr,fun=median,na.rm=T)
SNWsum.5yr <- raster::calc(SNWsum.yr,fun=median,na.rm=T)


dfClim <- data.frame(TMPmean=getValues(TMPmean.5yr),
                     TMPrange=getValues(TMPrange.5yr),
                     PPTsum=getValues(PPTsum.5yr),
                     AIndex=getValues(AIndex.5yr),
                     SNWsum=getValues(SNWsum.5yr))

dfClim <- dfClim[complete.cases(dfClim),]

dumLatLon <- rasterToPoints(TMPmean.5yr)

dfClim$Lat <- dumLatLon[,'y']
dfClim$Lon <- dumLatLon[,'x']

dfClim$Col <- colFromX(TMPmean.5yr, dfClim$Lon)
dfClim$Row <- rowFromY(TMPmean.5yr, dfClim$Lat)


save('dfClim', file=paste0(tpath,'df4ClimateSpace_05dd.RData'))


