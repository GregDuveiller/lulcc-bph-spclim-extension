


# gap filling 
source('codeProcessing/gapFill-with-RF.r')
gapFill(snowTag = 'withSNOW', outGrdSizDeg = c(0.5,0.5), inputTrans = 'FOR %=>% CaG')





gapFill <- function(snowTag = 'withSNOW', outGrdSizDeg = c(0.5,0.5),  inputTrans = 'FOR %=>% CaG'){
  
  require(raster)
  require(ncdf4)
  require(randomForest)
  require(dplyr)
  require(ggplot2)
  require(tidyr)
  require(RColorBrewer)
  require(scales)
  require(sf)
  require(here)
  
  # load vectors... to be defined after... 
  
  vpath <- '/ESS_Datasets/USERS/Duveiller/AncillaryDatasets/WorldVector/'
  
  wmap_df_land <- st_read(paste0(vpath,'ne_110m_land.shp'),quiet=T)
  wmap_df_ocean <- st_read(paste0(vpath,'ne_110m_ocean.shp'),quiet=T)
  

  # load CLIM data
  load('data/inter_data/climData/df4ClimateSpace_05dd.RData') # dfClim
  
  source('data/input_data/___loadDataPaths___.R')  # dpath_s4t
  
  
  
  
  # load and prepare a mask for permanent snow/ice & sparse vegetation

  minTHR <- 0.05
  
  LCfile2 <- paste0('data/input_data',
                    '/ESACCI-LC-L4-LCCS-Map-300m-P1Y-aggregated-1.000000Deg-2010-v2.0.7.nc')
  
  msk_FOR <- disaggregate(raster(LCfile2,varname='FOR'), fact = 2) > minTHR
  msk_SHR <- disaggregate(raster(LCfile2,varname='SHR'), fact = 2) > minTHR
  msk_CaG <- disaggregate(raster(LCfile2,varname='CaG'), fact = 2) > minTHR
  msk_SAV <- disaggregate(raster(LCfile2,varname='SAV'), fact = 2) > minTHR
  
  
  df_msk <- bind_rows(
    as.data.frame(msk_FOR|msk_SHR,xy=T) %>% mutate(Transition='FOR %=>% SHR'),
    as.data.frame(msk_FOR|msk_CaG,xy=T) %>% mutate(Transition='FOR %=>% CaG'),
    as.data.frame(msk_FOR|msk_SAV,xy=T) %>% mutate(Transition='FOR %=>% SAV'),
    as.data.frame(msk_SHR|msk_CaG,xy=T) %>% mutate(Transition='SHR %=>% CaG'),
    as.data.frame(msk_SHR|msk_SAV,xy=T) %>% mutate(Transition='SHR %=>% SAV'),
    as.data.frame(msk_CaG|msk_SAV,xy=T) %>% mutate(Transition='CaG %=>% SAV')) %>%
    rename(Lon=x,Lat=y,Mask=layer)
  
  df_msk$Mask[is.na(df_msk$Mask)] <- F  # because some had NAs due to msk sum
  
  
  

  tpath <- 'dataProcessing/df/'
  iVar <- 'albedo'
  print(paste('Working on',iVar))
  month <- 1
  type <- 'IGBPgen'
  
  #### Read ncdf and turn in into usable R dataframe ####
  
  rs <- brick(x = paste0(dpath_s4t, '/', iVar,'_',type,'.nc'), 
              varname = paste0('Delta_',iVar),
              level = month)
  
  # set a threshold on quality based on standard deviation within gridcell
  sd_thr <- 0.05
  sd <- brick(x = paste0(dpath_s4t, '/', iVar,'_',type,'.nc'), 
              varname = paste0('SD_Delta_',iVar),
              level = month)
  
  rs <- mask(x = rs, mask = sd > sd_thr, maskvalue = T, updatevalue = NA)
  
  # set a threshold on quality based on number of 0.05dd samples within gridcell
  nu_thr <- 5
  nu <- brick(x = paste0(dpath_s4t, '/', iVar,'_',type,'.nc'), 
              varname = paste0('N_Delta_',iVar),
              level = month)  
  rs <- mask(x = rs, mask = nu < nu_thr, maskvalue = T, updatevalue = NA)
  
  
  # pass the masked raster into a long dataframe
  df <- as.data.frame(x = rs, xy = T, long = T)
  colnames(df) <- c('Lon', 'Lat', 'Time', 'Delta_T') 
  

  # tweak to ensure no small decimal point uncertainties
  df <- df %>% 
    mutate(Lon = round(Lon,digits = 5), 
           Lat = round(Lat, digits = 5),
           Time = as.numeric(substring(Time,2)))

  # function to gap-fill with RF for a given transition
  estRF <- function(iTime){
    
    # prepare df
    dfdat <- df %>% 
      filter(Time == iTime) %>% 
      right_join(filter(df_msk, Transition == inputTrans), by = c('Lat', 'Lon')) %>% 
      filter(Mask == T) %>% 
      dplyr::select(Lon, Lat, Delta_T) %>% 
      left_join(dfClim, by = c('Lat', 'Lon'))
    
    
    nbins <- 25
    
    bin.TMPmean <- seq(floor(min(dfdat[,'TMPmean'],na.rm=T)), ceiling(max(dfdat[,'TMPmean'],na.rm=T)), length=nbins)
    bin.TMPrange <- seq(floor(min(dfdat[,'TMPrange'],na.rm=T)), ceiling(max(dfdat[,'TMPrange'],na.rm=T)), length=nbins)
    bin.PPTsum <- seq(floor(min(dfdat[,'PPTsum'],na.rm=T)), ceiling(max(dfdat[,'PPTsum'],na.rm=T)), length=nbins)
    bin.AIndex <- seq(floor(min(dfdat[,'AIndex'],na.rm=T)), ceiling(max(dfdat[,'AIndex'],na.rm=T)), length=nbins)
    bin.SNWsum <- seq(floor(min(dfdat[,'SNWsum'],na.rm=T)), ceiling(max(dfdat[,'SNWsum'],na.rm=T)), length=nbins)
    
    df.obs <- dfdat %>% 
      filter(is.na(Delta_T) == F) %>% # select only observations
      dplyr::select(TMPmean,TMPrange,PPTsum,AIndex,SNWsum)
    
    # find combined set of bins in which the given transition is present
    freq <-  as.data.frame(table(findInterval(as.vector(unlist(df.obs[,'TMPmean'])), bin.TMPmean),
                                 findInterval(as.vector(unlist(df.obs[,'TMPrange'])), bin.TMPrange),
                                 findInterval(as.vector(unlist(df.obs[,'PPTsum'])), bin.PPTsum),
                                 findInterval(as.vector(unlist(df.obs[,'AIndex'])), bin.AIndex),
                                 findInterval(as.vector(unlist(df.obs[,'SNWsum'])), bin.SNWsum)))
    
    freq0 <- filter(freq, Freq > 0)
    colnames(freq0) <- c(paste0('cat_', c('TMPmean','TMPrange','PPTsum','AIndex','SNWsum')),'Count')
    
    dfdat <- dfdat %>% 
      mutate(cat_TMPmean = factor(cut(TMPmean, breaks = bin.TMPmean, labels = F)),
             cat_TMPrange = factor(cut(TMPrange, breaks = bin.TMPrange, labels = F)),
             cat_PPTsum = factor(cut(PPTsum, breaks = bin.PPTsum, labels = F, include.lowest = T)),
             cat_AIndex = factor(cut(AIndex, breaks = bin.AIndex, labels = F, include.lowest = T)),
             cat_SNWsum = factor(cut(SNWsum, breaks = bin.SNWsum, labels = F, include.lowest = T))) %>%
      right_join(by=c(paste0('cat_', c('TMPmean','TMPrange','PPTsum','AIndex','SNWsum'))),freq0)
    
    # set seed
    set.seed(1982)
    
    # Make estimation of Delta_AST based on climate space
    #  fitRF <- randomForest(Delta_T ~ TMPmean + TMPrange + PPTsum + AIndex, data=filter(dfdat,is.na(Delta_T)==F))
    fitRF <- randomForest(Delta_T ~ TMPmean + TMPrange + PPTsum + SNWsum + AIndex, data=filter(dfdat,is.na(Delta_T)==F))
    dfdat$Delta_T_RF <- predict(fitRF,dfdat)
    
    # correct the bias on the residues (RF appears to underestimate extremes)
    fitBC <- lm(Delta_T ~ Delta_T_RF, data=filter(dfdat,is.na(Delta_T)==F))
    dfdat$Delta_T_RF <- predict(fitBC,dfdat)
    
    # add time in it...
    dfdat$Time <- iTime
    
    dir.create(tfpath <- paste0(tpath,'figuresRFfit___',iVar,'/'),showWarnings = F)
    
    sqshlims <- c(-0.3,0.3)
    gmap <- ggplot(gather(dfdat,key=Source,value=DeltaT,c('Delta_T','Delta_T_RF')))+
      geom_sf(data=wmap_df_land,fill='Grey50',colour='Grey50',size=0)+
      geom_raster(aes(x=Lon,y=Lat,fill=DeltaT))+
      geom_sf(data=wmap_df_ocean,fill='Grey20',colour='Grey20',size=0)+
      facet_wrap(~Source,nc=1)+
      scale_fill_gradientn(paste('Change\nin',iVar), colors = rev(brewer.pal(8,'RdBu')),limits=sqshlims,oob=squish)+
      coord_sf(ylim=c(-58,84),expand=F)+
      theme(panel.background=element_rect(fill='grey60'),
            panel.grid = element_blank())+
      ggtitle(paste0('Transition: ', gsub('%=>%','to',inputTrans)))
    
    ggsave(filename=paste0('Transition_', gsub(' %=>% ','2',inputTrans), '_time', iTime , '_RFmapCompare.png'),plot=gmap,path=tfpath,width=8,height=8)
    
    fulllims <- quantile(dfdat$Delta_T, probs = c(0.01,0.99), na.rm = T)
    gsca <- ggplot(dfdat)+
      geom_point(aes(x = Delta_T, y = Delta_T_RF),alpha=0.1,colour='cornflowerblue')+
      geom_abline()+
      coord_equal(xlim = fulllims, ylim = fulllims)+
      ggtitle(paste0('Transition_', gsub('%=>%','to',inputTrans),' | RMSE = ',round(sqrt(mean((dfdat$Delta_T_RF-dfdat$Delta_T)^2,na.rm=T)),digits=4)))
    
    ggsave(filename = paste0('Transition_', gsub(' %=>% ','2', inputTrans), '_time', iTime ,'_RFscatterplotCompare.png'),
           plot = gsca, path = tfpath, width = 6, height = 6)
    
    
    return(dfdat)
  }
  
  
  
  # run for all months
  df.RF <- NULL
  
  for(iTime in 1:12){
    df.RF <- bind_rows(estRF(iTime),df.RF)
    print(paste0('>>> Done for time ',iTime,'...'))
  }
  
  
  save('df.RF',file = paste0(tpath,'/RFest_df___',iVar,'_',snowTag,'.RData'))
  
  # get back to NetCDF
  
  
  varName <- 'Albedo'; units <- "";      
  if(snowTag != "withSNOW") { vName <- "snow free surface albedo"} else {vName <- "surface albedo"}
  
  var1 <- ncvar_def(name = paste0('Delta_',varName), units = units, dim = list(dimX, dimY, dimT), missval = NA, 
                    longname = paste('Difference in',vName,'for a given PFT transition'))
  
  nc <- nc_create(filename = paste0(opath, varName, '_', snowTag, '_spExt.nc'), vars = list(var1))
  
  for(iT in 1:12){
    
    df.RF.sub <- df.RF %>% filter(Time == iT) %>%
      select('Lon', 'Lat', 'Delta_T', 'Delta_T_RF')
    
    r <- rasterize(x = select(df.RF.sub, c('Lon', 'Lat')), 
                   y = raster(xmn = -180, xmx = 180, ymn = -90, ymx = 90, resolution = outGrdSizDeg),
                   field = df.RF.sub$Delta_T_RF)
    
    ncvar_put(nc = nc, varid = var1, vals = getValues(r), 
              start = c(1,1,iT), count = c(-1,-1,1), verbose = FALSE) 
    
    
  }
  
  nc_close(nc)
  
  print(paste('Finished with',iVar))
  
}