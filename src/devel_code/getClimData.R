# getClimData.R
#
# script to harvest the necessary climate data and store it into an appropriate
# data.frame for ingestion in the gap-filling function
# 
# G. Duveiller
################################################################################


#### Initialize #### ----

library(raster)
library(ncdf4)
library(R.utils)
library(here)


#### Get input data #### ----

# define years of interest
yrS <- 2008 
yrE <- 2012 

# get the input_data path
source('data/input_data/___loadDataPaths___.R')    # dpath_cru

# decompress file
tmp.ncfile <- gunzip(paste0(dpath_cru,'/cru_ts4.04.1901.2019.tmp.dat.nc.gz'), remove = F)
pre.ncfile <- gunzip(paste0(dpath_cru,'/cru_ts4.04.1901.2019.pre.dat.nc.gz'), remove = F)

# make a raster stack
tmp <- stack(tmp.ncfile)
pre <- stack(pre.ncfile)

# select in time 
tmp_sub <- tmp[[(12*(yrS - 1901)+1):(12*(yrE - 1901)+12)]]
pre_sub <- pre[[(12*(yrS - 1901)+1):(12*(yrE - 1901)+12)]]

# aggregate from 0.5dd to 1dd to match S4T data
tmp_sub <- aggregate(tmp_sub, fact = 2)
pre_sub <- aggregate(pre_sub, fact = 2)

#### Calculate climatic indices #### ----

# prepare empty stacks
TMP_mean.yr <- stack()
TMP_range.yr <- stack()
PPT_sum.yr <- stack()
AridIndex.yr <- stack()
SNW_sum.yr <- stack()

# add each year to the stack for each index
for(iT in 1:5){
  
  z1 <- raster::calc(tmp_sub[[(1:12)+12*(iT-1)]], 
                     fun = mean, na.rm = T)
  TMP_mean.yr <- addLayer(TMP_mean.yr, z1)
  
  z2 <- raster::calc(tmp_sub[[(1:12)+12*(iT-1)]], 
                     fun = function(x){(max(x)-min(x))})
  TMP_range.yr <- addLayer(TMP_range.yr, z2)
  
  z3 <- raster::calc(pre_sub[[(1:12)+12*(iT-1)]],
                     fun = sum, na.rm = T)
  PPT_sum.yr <- addLayer(PPT_sum.yr, z3)
  
  z4 <- z3/(z1 + 30)
  AridIndex.yr <- addLayer(AridIndex.yr, z4)
  
  subZtemp <- tmp_sub[[(1:12)+12*(iT-1)]] < 0
  snowppt <- pre_sub[[(1:12)+12*(iT-1)]] * subZtemp
  z5 <- raster::calc(snowppt,fun = sum, na.rm = T)
  SNW_sum.yr <- addLayer(SNW_sum.yr, z5)
  
}

# calculate raster medians and store all in a stack
clim.5yr <- stack(
  raster::calc(TMP_mean.yr,  fun = median, na.rm = T),
  raster::calc(TMP_range.yr, fun = median, na.rm = T),
  raster::calc(PPT_sum.yr,   fun = median, na.rm = T),
  raster::calc(AridIndex.yr, fun = median, na.rm = T),
  raster::calc(SNW_sum.yr,   fun = median, na.rm = T))


# Convert to points in a data.frame
dfClim <- as.data.frame(rasterToPoints(clim.5yr))
colnames(dfClim) <- c('Lon', 'Lat', 'TMP_mean', 'TMP_range', 'PPT_sum', 'AridIndex', 'SNW_sum')

#### Export and clean-up #### ----

out_dpath <- 'data/inter_data/climData'
dir.create(path = out_dpath, showWarnings = F, recursive = T)

save('dfClim', file = paste0(out_dpath,'/df4ClimateSpace_1dd.RData'))

# delete the uncompressed files
file.remove(tmp.ncfile)
file.remove(pre.ncfile)
