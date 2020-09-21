spextend_layer <- function(df_dat, checkplot = T){
  # function to spatially extend with RF for a given transition at a given time
  
  require(raster)
  require(ncdf4)
  require(randomForest)
  
  
  
  inputTrans <-  unique(df_dat$Transition)
  iMonth <- unique(df_dat$month)
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
    right_join(by = c(paste0('cat_', c('TMP_mean','TMP_range','PPT_sum','AridIndex','SNW_sum'))),freq0)
  
  # set seed
  set.seed(1982)
  
  # prepare data
  fltr.data <- filter(df_dat, is.na(delta) == F)
  
  # Make estimation of Delta_AST based on climate space
  #  fitRF <- randomForest(delta ~ TMPmean + TMPrange + PPTsum + AIndex, data=filter(df_dat,is.na(delta)==F))
  fitRF <- randomForest(delta ~ TMP_mean + TMP_range + PPT_sum + SNW_sum + AridIndex, 
                        data = fltr.data, ntree = 500, mtry = 3)
  df_dat$delta_RF0 <- predict(fitRF,df_dat)
  
  # correct the bias on the residues (RF appears to underestimate extremes)
  fitBC <- lm(delta ~ delta_RF0, data = filter(df_dat, is.na(delta) == F))
  df_dat$delta_RF1 <- predict(fitBC, df_dat)
  
  # Add FIG if requested ... 
  source('src/devel_code/make_checkplot.R')
  if(checkplot == T){ 
    ctrl.fig.path <- paste0('scratch/gapfilling_checkplots/', iVar, '_', type)
    make_checkplot(df_dat, iMonth, ctrl.fig.path) }
 
  # clean out df before returning it
  df_out <- df_dat %>% select(Lon, Lat, delta_RF1) %>%
    rename(delta = delta_RF1)
  return(df_out)
}
