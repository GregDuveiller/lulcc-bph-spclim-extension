make_checkplot <- function(df_dat, iMonth, ctrl.fig.path){
  
  require(ggplot2)
  require(grid)
  require(scales)
  require(sf)
  
  
  # load vectors... (for coasts on figures)
  wmap_df_land <- st_read('data/WorldVector/ne_50m_land.shp', quiet = T)
  wmap_df_ocean <- st_read('data/WorldVector/ne_50m_ocean.shp', quiet = T)
  
  
  
  monthTag <- ifelse(iMonth < 10, 
                     paste0('0', as.character(iMonth)), 
                     as.character(iMonth))
  
  inputTrans <-  unique(df_dat$Transition)
  
  sqshlims <- max(abs(quantile(df_dat$delta, probs = c(0.05, 0.95), na.rm = T))) * c(-1,1)
  
  gmap <- ggplot(gather(df_dat, key = Source, value = DeltaT, 
                        c('delta','delta_RF0','delta_RF1')))+
    geom_sf(data = wmap_df_land, fill='Grey50',colour='Grey50',size=0)+
    geom_raster(aes(x=Lon, y=Lat,fill=DeltaT))+
    geom_sf(data = wmap_df_ocean, fill='Grey20',colour='Grey20',size=0)+
    facet_wrap(~Source, nc = 1)+
    scale_fill_gradientn(paste('Change in',iVar), colors = rev(brewer.pal(8,'RdBu')),
                         limits = sqshlims, oob = squish) +
    coord_sf(ylim = c(-58,84), expand = F)+
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
               colour = 'cornflowerblue')+
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
  fig.fullfname <- paste0(ctrl.fig.path, '/', fig.name, '.', fig.fmt)
  
  dir.create(ctrl.fig.path, showWarnings = F, recursive = T)
  
  if(fig.fmt == 'png'){png(fig.fullfname, width = fig.width, height = fig.height, units = "in", res= 150)}
  if(fig.fmt == 'pdf'){pdf(fig.fullfname, width = fig.width, height = fig.height)}
  
  print(gmap, vp = viewport(width = 0.6, height = 1.0, x = 0.0, y = 0.0, just = c(0,0)))
  print(gsc0, vp = viewport(width = 0.4, height = 0.5, x = 0.6, y = 0.5, just = c(0,0)))
  print(gsc1, vp = viewport(width = 0.4, height = 0.5, x = 0.6, y = 0.0, just = c(0,0)))
  
  dev.off()
  
}
