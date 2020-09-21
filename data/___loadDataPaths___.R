# load paths where the data can be found

sys <- Sys.info()
if(sys['nodename'] == "Gregorys-MBP.fritz.box" & 
   sys['user'] == 'greg'){ # if it's me, at home, likely the data is on a NAS
  main_data_path <- '/Volumes/home/work/data'
}

# CRU data
dpath_cru <- paste0(main_data_path, '/external_datasets/climate/CRU')

# S4T original data (v1) before gap-filling spatial extension
dpath_s4t <-  paste0(main_data_path, '/internal_datasets/bph-lulcc___S4Tdata')
