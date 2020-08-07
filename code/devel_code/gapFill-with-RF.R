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


# # gap filling 
# source('codeProcessing/gapFill-with-RF.r')
# gapFill( outGrdSizDeg = c(0.5,0.5), inputTrans = 'FOR %=>% CaG')


#### Initialisation #### ---- 


fig.path = 'data/inter_data/gapfilling_checkplots'

#### load data #### ---- 
source('data/input_data/___loadDataPaths___.R')  # dpath_s4t

# load vectors... to be defined after... 
vpath <- '../../AncillaryDatasets/WorldVector/'
wmap_df_land <- st_read(paste0(vpath,'ne_50m_land.shp'), quiet = T)
wmap_df_ocean <- st_read(paste0(vpath,'ne_50m_ocean.shp'), quiet = )



# load CLIM data
load('data/inter_data/climData/df4ClimateSpace_1dd.RData') # dfClim





#### Prepare masking to define extent #### ---- 

# >>> NEED TO SWITCH BETWEEN DATASET (gen vs det)

minTHR <- 0.05

LCfile2 <- paste0('data/input_data',
                  '/ESACCI-LC-L4-LCCS-Map-300m-P1Y-aggregated-1.000000Deg-2010-v2.0.7.nc')

msk_FOR <- raster(LCfile2,varname='FOR') > minTHR
msk_SHR <- raster(LCfile2,varname='SHR') > minTHR
msk_CaG <- raster(LCfile2,varname='CaG') > minTHR
msk_SAV <- raster(LCfile2,varname='SAV') > minTHR

df_msk <- bind_rows(
  as.data.frame(msk_FOR | msk_SHR, xy = T) %>% 
    mutate(Transition = 'FOR %=>% SHR', iTr = 1),
  as.data.frame(msk_FOR | msk_CaG, xy = T) %>%
    mutate(Transition ='FOR %=>% CaG', iTr = 2),
  as.data.frame(msk_FOR|msk_SAV, xy=T) %>%
    mutate(Transition='FOR %=>% SAV', iTr = 3),
  as.data.frame(msk_SHR|msk_CaG, xy=T) %>%
    mutate(Transition='SHR %=>% CaG', iTr = 4),
  as.data.frame(msk_SHR|msk_SAV, xy=T) %>%
    mutate(Transition='SHR %=>% SAV', iTr = 5),
  as.data.frame(msk_CaG|msk_SAV, xy=T) %>% 
    mutate(Transition='CaG %=>% SAV', iTr = 6)
) %>%
  rename(Lon = x, Lat = y, Mask = layer)

df_msk$Mask[is.na(df_msk$Mask)] <- F  # because some had NAs due to msk sum


# 


#### Define function to gap-fill a single layer #### ---- 

# function to gap-fill with RF for a given transition at a given time
gapfill_layer <- function(df_dat, checkplot = T){
  
  inputTrans <- unique(df_dat$Transition)
  iMonth <- unique(month)
  nbins <- 25
  
  bin.TMP_mean  <- seq(floor(min(df_dat[,'TMP_mean'], na.rm=T)), 
                       ceiling(max(df_dat[,'TMP_mean'], na.rm=T)), 
                       length=nbins)
  bin.TMP_range <- seq(floor(min(df_dat[,'TMP_range'], na.rm=T)), 
                       ceiling(max(df_dat[,'TMP_range'], na.rm=T)),
                       length=nbins)
  bin.PPT_sum   <- seq(floor(min(df_dat[,'PPT_sum'], na.rm=T)), 
                       ceiling(max(df_dat[,'PPT_sum'], na.rm=T)), 
                       length=nbins)
  bin.AridIndex <- seq(floor(min(df_dat[,'AridIndex'], na.rm=T)), 
                       ceiling(max(df_dat[,'AridIndex'], na.rm=T)), 
                       length=nbins)
  bin.SNW_sum   <- seq(floor(min(df_dat[,'SNW_sum'], na.rm=T)), 
                       ceiling(max(df_dat[,'SNW_sum'], na.rm=T)), 
                       length=nbins)
  
  df_obs <- df_dat %>% 
    filter(is.na(delta) == F) %>% # select only observations
    dplyr::select(TMP_mean, TMP_range, PPT_sum, AridIndex, SNW_sum)
  
  # find combined set of bins in which the given transition is present
  freq <-  as.data.frame(
    table(findInterval(as.vector(unlist(df_obs[,'TMP_mean'])), bin.TMP_mean),
          findInterval(as.vector(unlist(df_obs[,'TMP_range'])), bin.TMP_range),
          findInterval(as.vector(unlist(df_obs[,'PPT_sum'])), bin.PPT_sum),
          findInterval(as.vector(unlist(df_obs[,'AridIndex'])), bin.AridIndex),
          findInterval(as.vector(unlist(df_obs[,'SNW_sum'])), bin.SNW_sum)))
  
  freq0 <- filter(freq, Freq > 0)
  colnames(freq0) <- c(paste0('cat_', c('TMP_mean','TMP_range','PPT_sum','AridIndex','SNW_sum')),'Count')
  
  df_dat <- df_dat %>% 
    mutate(cat_TMP_mean = factor(cut(TMP_mean, breaks = bin.TMP_mean, labels = F)),
           cat_TMP_range = factor(cut(TMP_range, breaks = bin.TMP_range, labels = F)),
           cat_PPT_sum = factor(cut(PPT_sum, breaks = bin.PPT_sum, labels = F, include.lowest = T)),
           cat_AridIndex = factor(cut(AridIndex, breaks = bin.AridIndex, labels = F, include.lowest = T)),
           cat_SNW_sum = factor(cut(SNW_sum, breaks = bin.SNW_sum, labels = F, include.lowest = T))) %>%
    right_join(by=c(paste0('cat_', c('TMP_mean','TMP_range','PPT_sum','AridIndex','SNW_sum'))),freq0)
  
  # set seed
  set.seed(1982)
  
  # Make estimation of Delta_AST based on climate space
  #  fitRF <- randomForest(delta ~ TMPmean + TMPrange + PPTsum + AIndex, data=filter(df_dat,is.na(delta)==F))
  fitRF <- randomForest(delta ~ TMP_mean + TMP_range + PPT_sum + SNW_sum + AridIndex, 
                        data = filter(df_dat, is.na(delta) == F))
  df_dat$delta_RF0 <- predict(fitRF,df_dat)
  
  # correct the bias on the residues (RF appears to underestimate extremes)
  fitBC <- lm(delta ~ delta_RF0, data = filter(df_dat, is.na(delta) == F))
  df_dat$delta_RF1 <- predict(fitBC, df_dat)
  
  
  # Add FIG if requested ... 
  if(checkplot == T){
    
    require(ggplot2)
    require(grid)
    require(scales)
    
    monthTag <- ifelse(iMonth < 10, 
                       paste0('0',as.character(iMonth)), 
                       as.character(iMonth))
    
    sqshlims <- max(abs(quantile(df_dat$delta, probs = c(0.05, 0.95), na.rm = T))) * c(-1,1)
    
    gmap <- ggplot(gather(df_dat, key = Source, value = DeltaT, 
                          c('delta','delta_RF0','delta_RF1')))+
      geom_sf(data = wmap_df_land, fill='Grey50',colour='Grey50',size=0)+
      geom_raster(aes(x=Lon, y=Lat,fill=DeltaT))+
      geom_sf(data = wmap_df_ocean, fill='Grey20',colour='Grey20',size=0)+
      facet_wrap(~Source, nc=1)+
      scale_fill_gradientn(paste('Change in',iVar), colors = rev(brewer.pal(8,'RdBu')),
                           limits = sqshlims, oob = squish) +
      coord_sf(ylim = c(-58,84), expand=F)+
      theme(panel.background=element_rect(fill='grey60'),
            panel.grid = element_blank(),
            legend.position = 'bottom',
            legend.key.width = unit(2.4, "cm")) +
      guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5)) + 
      ggtitle(paste('Transition', gsub('%=>%','to',inputTrans),
                    'for', month.name[iMonth] ))
    
    
    fulllims <- quantile(df_dat$delta, probs = c(0.01, 0.99), na.rm = T)
    
    gsc0 <- ggplot(df_dat)+
      geom_point(aes(x = delta, y = delta_RF0), alpha = 0.1,
                 colour='cornflowerblue')+
      geom_abline()+
      coord_equal(xlim = fulllims, ylim = fulllims)+
      ggtitle(paste('Before bias corr. | RMSE =', 
                    round(sqrt(mean((df_dat$delta_RF0 - df_dat$delta)^2,
                                    na.rm = T)), digits = 4)))
    
    gsc1 <- ggplot(df_dat)+
      geom_point(aes(x = delta, y = delta_RF1), alpha = 0.1,colour='navyblue')+
      geom_abline()+
      coord_equal(xlim = fulllims, ylim = fulllims)+
      ggtitle(paste('After bias corr. | RMSE =', 
                    round(sqrt(mean((df_dat$delta_RF1 - df_dat$delta)^2,
                                    na.rm = T)), digits = 4)))    
    
    ## Printing the entire figure ----
    fig.name <-  paste(type, iVar, gsub(' %=>% ','2', inputTrans), monthTag ,'checkplot', sep = '_')
    fig.width <- 10; fig.height <- 8; fig.fmt <- 'png'
    fig.fullfname <- paste0(fig.path, '/', fig.name, '.', fig.fmt)
    
    dir.create(fig.path, showWarnings = F, recursive = T)
    
    if(fig.fmt == 'png'){png(fig.fullfname, width = fig.width, height = fig.height, units = "in", res= 150)}
    if(fig.fmt == 'pdf'){pdf(fig.fullfname, width = fig.width, height = fig.height)}
    
    print(gmap, vp = viewport(width = 0.6, height = 1.0, x = 0.0, y = 0.0, just = c(0,0)))
    print(gsc0, vp = viewport(width = 0.4, height = 0.5, x = 0.6, y = 0.5, just = c(0,0)))
    print(gsc1, vp = viewport(width = 0.4, height = 0.5, x = 0.6, y = 0.0, just = c(0,0)))
    
    dev.off()
    
  }
  
  df_out <- df_dat %>% select(Lon, Lat, delta_RF1) %>%
    rename(delta = delta_RF1)
  return(df_out)
}



#### Wrap around a file #### ---- 

iVar <- 'albedo' 

type <- 'IGBPgen'

nc_filename_in <- paste0(dpath_s4t, '/', iVar,'_',type,'.nc')
file.copy(from = nc_filename_in, 
          to = paste0('data/final_data/', iVar, '_', type, '_gapfilled.nc'))
nc_out <- paste0('data/final_data/', iVar, '_', type, '_gapfilled.nc')
nc <- nc_open(nc_filename_in)



iMonth <- 12
iTr <- 3


# get the data 
nc_dat <- ncvar_get(nc, varid = paste0('Delta_',iVar), 
                    start = c(1,1, iMonth, iTr), count = c(-1, -1, 1, 1))


# set a threshold on quality based on standard deviation within gridcell
sd_thr <- 0.05
sd_dat <- ncvar_get(nc, varid = paste0('SD_Delta_',iVar), 
                    start = c(1,1, iMonth, iTr), count = c(-1, -1, 1, 1))

# set a threshold on quality based on number of 0.05dd samples within gridcell
nu_thr <- 5
nu_dat <- ncvar_get(nc, varid = paste0('N_Delta_',iVar), 
                    start = c(1,1, iMonth, iTr), count = c(-1, -1, 1, 1))

nc_dat[sd_dat > sd_thr] <- NaN
nc_dat[nu_dat < nu_thr] <- NaN


df_dat <- dfClim %>% 
  mutate(delta = as.vector(nc_dat),
         month = iMonth) %>%
  left_join(df_msk %>% filter(iTr == 1), by = c('Lon', 'Lat')) %>%
  filter(Mask == T)

df_out <- gapfill_layer(df_dat)








# get back to NetCDF


var1 <- ncvar_def(name = paste0('Delta_',varName), units = units, dim = list(dimX, dimY, dimT), missval = NA, 
                  longname = paste('Difference in', vName,' for a given PFT transition'))

nc <- nc_create(filename = paste0(opath, varName, '_', snowTag, '_spExt.nc'), vars = list(var1))

for(iT in 1:12){
  
  df.RF.sub <- df.RF %>% filter(Time == iT) %>%
    select('Lon', 'Lat', 'delta', 'delta_RF')
  
  r <- rasterize(x = select(df.RF.sub, c('Lon', 'Lat')), 
                 y = raster(xmn = -180, xmx = 180, ymn = -90, ymx = 90, resolution = outGrdSizDeg),
                 field = df.RF.sub$delta_RF)
  
  ncvar_put(nc = nc, varid = var1, vals = getValues(r), 
            start = c(1,1,iT), count = c(-1,-1,1), verbose = FALSE) 
  
  
}

nc_close(nc)

print(paste('Finished with',iVar))

