spclimext_withRF <- function(dfClim, iVar, type, nu_thr = 10, minTHR = 0.05, 
                             comprLevel = 7, mzVal = NA, do_checkplot = T){
  
  
  require(raster)
  require(ncdf4)
  require(randomForest)
  require(dplyr)
  require(tidyr)
  
  
  
  #### Prepare masking to define extent #### ---- 
  if(type == 'IGBPgen'){
    
    LCfile <- paste0('results/cleaned_input_data/land_cover', '/CCI2', type,
                     '/ESACCI-LC-L4-LCCS-Map-300m-P1Y-aggregated-1.000000Deg-2010-v2.0.7.nc')
    
    
    msk_FOR <- raster(LCfile, varname = 'FOR') > minTHR
    msk_SHR <- raster(LCfile, varname = 'SHR') > minTHR
    msk_CaG <- raster(LCfile, varname = 'CaG') > minTHR
    msk_SAV <- raster(LCfile, varname = 'SAV') > minTHR
    
    df_msk <- bind_rows(
      as.data.frame(msk_FOR | msk_SHR, xy = T) %>% 
        mutate(Transition = 'FOR %=>% SHR', iTr = 1),
      as.data.frame(msk_FOR | msk_CaG, xy = T) %>%
        mutate(Transition = 'FOR %=>% CaG', iTr = 2),
      as.data.frame(msk_FOR | msk_SAV, xy = T) %>%
        mutate(Transition = 'FOR %=>% SAV', iTr = 3),
      as.data.frame(msk_SHR | msk_CaG, xy = T) %>%
        mutate(Transition = 'SHR %=>% CaG', iTr = 4),
      as.data.frame(msk_SHR | msk_SAV, xy = T) %>%
        mutate(Transition = 'SHR %=>% SAV', iTr = 5),
      as.data.frame(msk_CaG | msk_SAV, xy = T) %>% 
        mutate(Transition = 'CaG %=>% SAV', iTr = 6)
    ) %>%
      rename(Lon = x, Lat = y, Mask = layer)
    
    df_msk$Mask[is.na(df_msk$Mask)] <- F  # because some had NAs due to msk sum
    nTrans <- length(unique(df_msk$Transition))
    
  }
  if(type == 'IGBPdet'){
    
    LCfile <- paste0('results/cleaned_input_data/land_cover', '/CCI2', type,
                     '/ESACCI-LC-L4-LCCS-Map-300m-P1Y-aggregated-1.000000Deg-2010-v2.0.7.nc')
    
    
    
    vgt.grps <- c('EBF','DBF','ENF','DNF','MF','SAV','SHR','GRA','CRO','WET')
    
    for(iVG in vgt.grps){
      msk <- raster(LCfile, varname = iVG) > minTHR
      eval(parse(text = paste0('msk_', iVG, ' <- msk')))
    }
    
    n <- length(vgt.grps); nTr <- n*(n-1)/2
    Tr.From <- vector(length = nTr); Tr.To <- Tr.From; i <- 1
    
    df_msk <- NULL
    iTrNum <- 1
    for(iPFTfrom in vgt.grps[1:(n-1)]){
      for(iPFTto in vgt.grps[(1 + which(vgt.grps == iPFTfrom)):n]){
        
        eval(parse(text = paste0('df <- as.data.frame(', 'msk_', iPFTfrom, 
                                 ' | msk_', iPFTto, ', xy = T)')))
        df <- df %>% 
          mutate(Transition = paste(iPFTfrom, iPFTto, sep = ' %=>% '), 
                 iTr = iTrNum) %>%
          rename(Lon = x, Lat = y, Mask = layer)
        df_msk <- df_msk %>% bind_rows(df)
        iTrNum <- iTrNum + 1
      }}
    
    df_msk$Mask[is.na(df_msk$Mask)] <- F  # because some had NAs due to msk sum
    nTrans <- length(unique(df_msk$Transition))
  }
  
  #### Define function to spatially extend a single layer #### ---- 
  source('code/spextend_layer.R')
  
  #### Prepare run for a given NCDF file #### ---- 
  
  nc_filename_in <- paste0('data/bph-lulcc___S4Tdata/', iVar,'_',type,'.nc')
  
  # make a local copy to ensure connection to server is not lost whilst processing
  nc_filename_tmp <- paste0('scratch/', basename(nc_filename_in))
  file.copy(from = nc_filename_in, to = nc_filename_tmp, overwrite = T)
  nc <- nc_open(nc_filename_tmp)
  
  # preprare output file
  nc_filename_out <- paste0('results/final_products/', iVar, '_', type, '_ext.nc')
  
  print(paste('Starting with', iVar, 'in', type, 'classification...'))
  
  # make new nc
  
  
  # Prepare the 4 dimensions
  dimX <- ncdim_def(name = 'lon',
                    units = ncatt_get(nc, varid = 'lon', attname = 'units')$value, 
                    vals = ncvar_get(nc, varid = 'lon'), unlim = F,
                    longname = ncatt_get(nc, varid = 'lon', attname = 'long_name')$value)
  dimY <- ncdim_def(name = 'lat',
                    units = ncatt_get(nc, varid = 'lat', attname = 'units')$value, 
                    vals = ncvar_get(nc, varid = 'lat'), unlim = F,
                    longname = ncatt_get(nc, varid = 'lat', attname = 'long_name')$value)
  dimM <- ncdim_def(name = 'mon',
                    units = ncatt_get(nc, varid = 'mon', attname = 'units')$value, 
                    vals = ncvar_get(nc, varid = 'mon'), unlim = F,
                    longname = ncatt_get(nc, varid = 'mon', attname = 'long_name')$value)
  dimT <- ncdim_def(name = 'iTr',
                    units = 'transition', 
                    vals = ncvar_get(nc, varid = 'iTr'), unlim = F,
                    longname = ncatt_get(nc, varid = 'iTr', attname = 'long_name')$value)
  
  
  # prepare new (gap-filled) var
  var_longname <- ncatt_get(nc, varid = paste0('Delta_', iVar), 
                            attname = 'long_name')$value
  
  new_var_GFD <- ncvar_def(name = paste0('Delta_', iVar, '_ext'), 
                           dim = list(dimX, dimY, dimM, dimT), 
                           units = '', 
                           missval = mzVal, compression = comprLevel,
                           longname = paste0('Spatially extended ', 
                                             tolower(substr(var_longname, 0, 1)),
                                             substring(var_longname, 2)))
  
  varID <- paste0('Delta_', iVar)
  new_var_delta <- ncvar_def(name = varID, 
                             dim = list(dimX, dimY, dimM, dimT), 
                             units = '', 
                             missval = mzVal, compression = comprLevel,
                             longname =  ncatt_get(nc, varid = varID, 
                                                   attname = 'long_name')$value)
  
  varID <- paste0('SD_Delta_', iVar)
  new_var_sd <- ncvar_def(name = varID, 
                          dim = list(dimX, dimY, dimM, dimT), 
                          units = '', 
                          missval = mzVal, compression = comprLevel,
                          longname =  ncatt_get(nc, varid = varID, 
                                                attname = 'long_name')$value)
  
  varID <- paste0('N_Delta_', iVar)
  new_var_n <- ncvar_def(name = varID, 
                         dim = list(dimX, dimY, dimM, dimT), 
                         units = '', 
                         missval = mzVal, compression = comprLevel,
                         longname =  ncatt_get(nc, varid = varID, 
                                               attname = 'long_name')$value)
  
  
  
  # ensure target file does not yet exist
  file.remove(nc_filename_out)

  # create new nc (and thus overwrite)
  nc_new <- nc_create(filename = nc_filename_out,
                      vars = list(new_var_GFD, new_var_delta, new_var_sd, new_var_n), 
                      force_v4 = T, verbose = F)
  nc_close(nc_new)
  
  
  
  
  #### Wrap around a file #### ---- 
  
  nc_new <- nc_open(nc_filename_out, write = TRUE)
  
  for(iTrans in 1:nTrans){
    for(iMonth in 1:12){
      
      # get the original delta data ...
      nc_dat <- ncvar_get(nc, varid = paste0('Delta_',iVar), 
                          start = c(1,1, iMonth, iTrans), count = c(-1, -1, 1, 1))
      # ... and pass it to the new NetCDF
      ncvar_put(nc = nc_new, varid = paste0('Delta_',iVar), vals = nc_dat,
                start = c(1,1, iMonth, iTrans), count = c(-1, -1, 1, 1))
      
      
      # get the sd delta data ...
      sd_dat <- ncvar_get(nc, varid = paste0('SD_Delta_',iVar),
                          start = c(1,1, iMonth, iTrans), count = c(-1, -1, 1, 1))
      # ... and pass it to the new NetCDF
      ncvar_put(nc = nc_new, varid = paste0('SD_Delta_',iVar), vals = sd_dat,
                start = c(1,1, iMonth, iTrans), count = c(-1, -1, 1, 1))
      
      
      # get the N delta data
      nu_dat <- ncvar_get(nc, varid = paste0('N_Delta_',iVar), 
                          start = c(1,1, iMonth, iTrans), count = c(-1, -1, 1, 1))
      # ... and pass it to the new NetCDF
      ncvar_put(nc = nc_new, varid = paste0('N_Delta_',iVar), vals = nu_dat,
                start = c(1,1, iMonth, iTrans), count = c(-1, -1, 1, 1))    
      
      
      
      # filter out points that have few underlying samples
      nc_dat[nu_dat < nu_thr] <- NaN
      
      # prepare the data into a proper dataframe
      df_dat <- dfClim %>% 
        mutate(delta = as.vector(nc_dat),
               month = iMonth) %>%
        left_join(df_msk %>% filter(iTr == iTrans), by = c('Lon', 'Lat')) %>%
        filter(Mask == T)
      
      # spatially extend/gap-fill where possible
      df_dum <- spextend_layer(df_dat, checkplot = do_checkplot)
      
      # fill extra points with no data
      df_out <- dfClim %>% select(Lon, Lat) %>% 
        left_join(df_dum, by = c("Lon", "Lat"))
      
      # put the new variable inside the new NetCDF
      ncvar_put(nc = nc_new, 
                varid = paste0('Delta_', iVar, '_ext'), 
                vals = df_out$delta,
                start = c(1,1, iMonth, iTrans), count = c(-1, -1, 1, 1))
      
      
    }
    print(paste(' >>>', iVar, 'transition', iTrans, '/', nTrans, 'is done'))
  }
  
  nc_close(nc_new)
  nc_close(nc)
  file.remove(nc_filename_tmp)
  
  print(paste('Finished with',iVar))
  
  
}










