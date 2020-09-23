vizchech_for_output <- function(iVar, type, iMon, iTrn, path){
  
  require(ncdf4)
  require(RColorBrewer)
  require(ggplot2)
  require(scales)
  require(patchwork)
  
  # # input params
  # iVar <- 'LSTnight' # 'HG'  # 'albedo'    # 'LSTnight'  # 'LSTday'  
  # type <- 'IGBPgen'   # 'IGBPdet'   # 'IGBPgen'
  # iMon <- 8
  # iTrn <- 2
  # path <- 'results/final_products'
  
  # Get the data needed in a dataframe
  target_file <- paste0(path, '/', iVar, '_', type, '_ext.nc')
  nc <- nc_open(target_file)
  dat_v1 <- ncvar_get(nc, varid = paste0('Delta_', iVar, ''), 
                      start = c(1, 1, iMon, iTrn), count = c(-1, -1, 1, 1))
  dat_v2 <- ncvar_get(nc, varid = paste0('Delta_', iVar, '_ext'), 
                      start = c(1, 1, iMon, iTrn), count = c(-1, -1, 1, 1))
  lat <- ncvar_get(nc, varid = 'lat')
  lon <- ncvar_get(nc, varid = 'lon')
  df_dat <- data.frame(lat = rep(lat, each = length(lon)),
                       lon = rep(lon, times = length(lat)),
                       version_1 = as.vector(dat_v1),
                       version_2 = as.vector(dat_v2))
  df_dat <- df_dat[!is.nan(df_dat$version_2),]
  
  # get valid ranges
  sqshlims <- c(-1,1) * max(abs(quantile(df_dat$version_1, 
                                         probs = c(0.05, 0.95), na.rm = T)))
  fulllims <- quantile(df_dat$version_1, probs = c(0.01, 0.99), na.rm = T)
  
  # plot the map
  gmap <- ggplot(gather(df_dat, key = Source, value = delta, 
                        c('version_1','version_2'))) +
    # geom_sf(data = wmap_df_land, fill='Grey50',colour='Grey50',size=0)+
    geom_tile(aes(x = lon, y = lat, fill = delta))+
    #geom_sf(data = wmap_df_ocean, fill='Grey20',colour='Grey20',size=0)+
    facet_wrap(~Source, nc = 1)+
    scale_fill_gradientn('', colors = rev(brewer.pal(8,'RdBu')),
                         limits = sqshlims, oob = squish) +
    coord_sf(ylim = c(-58,84), expand = F)+
    theme(panel.background = element_rect(fill='grey60'),
          panel.grid = element_blank(),
          legend.position = 'right',
          legend.key.height = unit(2.1, "cm")) +
    guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5)) + 
    ggtitle(label = paste('Change in', iVar, 'under', type, 'classification'), 
            subtitle = paste('Transition', iTrn, 'for', month.name[iMon]))
  

  # plot the scatterplot
  gsct <- ggplot(df_dat)+
    geom_point(aes(x = version_1, y = version_2), 
                   colour = 'darkgreen', alpha = 0.1, size = 0.5) +
    geom_abline() +
    # ggtitle(label = paste0(iVar, '_', type), 
    #         subtitle = paste('Confrontation between versions | RMSE =', 
    #                          round(sqrt(mean((df_dat$version_1 - df_dat$version_2)^2,
    #                                          na.rm = T)), digits = 4))) +
    coord_equal(xlim = fulllims, ylim = fulllims) +
    theme(plot.background = element_rect(fill = 'white', colour = 'grey10')) 
    
  # merge them all and return output
  plt <- gmap + gsct + 
    plot_layout(design = c(area(t = 1, l = 1, b = 5, r = 6),
                           area(t = 5, l = 1, b = 5, r = 1)))
  return(plt)
}
