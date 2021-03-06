# Script saying where the input data originally are and setting symbolic links
# 
# NOTE: Still unsure if this is the best way, as there is no reason to set this 
#       info in a public repo


sys <- Sys.info()
if(sys['nodename'] == "Gregorys-MBP.fritz.box" & 
   sys['user'] == 'greg'){ # if it's me, at home, likely the data is on a NAS
  
  spath_s4t <- '/Volumes/home/work/data/internal_datasets/bph-lulcc___S4Tdata/'
  spath_cru <- '/Volumes/home/work/data/external_datasets/climate/CRU/'
  spath_vct <- '/Volumes/home/work/data/external_datasets/WorldVector/'
  
}

if(sys['nodename'] == "jeodpp-terminal-151p-02" & 
   sys['user'] == 'duveigr'){ # if it's me, on JEODPP
  
  spath_s4t <- '/storage/duveigr/internal_datasets/bph-lulcc___S4Tdata/v1.0/'
  spath_cru <- '/storage/duveigr/external_datasets/climate/CRU/'
  spath_vct <- '/storage/duveigr/external_datasets/WorldVector/'
}

# prepare symlinks to data on the local setup

tpath_s4t <- 'data/bph-lulcc___S4Tdata'
if(file.exists(tpath_s4t)){
  if(Sys.readlink(tpath_s4t) != spath_s4t){
    file.remove(tpath_s4t)
    file.symlink(from = spath_s4t, to = tpath_s4t)
  }
} else {
  file.symlink(from = spath_s4t, to = tpath_s4t)
}

tpath_cru <- 'data/CRU'
if(file.exists(tpath_cru)){
  if(Sys.readlink(tpath_cru) != spath_cru){
    file.remove(tpath_cru)
    file.symlink(from = spath_cru, to = tpath_cru)
  }
} else {
  file.symlink(from = spath_cru, to = tpath_cru)
}

tpath_vct <- 'data/WorldVector'
if(file.exists(tpath_vct)){
  if(Sys.readlink(tpath_vct) != spath_vct){
    file.remove(tpath_vct)
    file.symlink(from = spath_vct, to = tpath_vct)
  }
} else {
  file.symlink(from = spath_vct, to = tpath_vct)
}
