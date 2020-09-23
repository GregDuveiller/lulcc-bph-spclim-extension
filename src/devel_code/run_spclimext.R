# Script to run an spclim job 

require(here)

#### setup run #### ---- 

sys <- Sys.info()
if(sys['nodename'] == "Gregorys-MBP.fritz.box" & 
   sys['user'] == 'greg'){ # if it's me, at home, likely the data is on a NAS
  
  spath_s4t <- '/Volumes/home/work/data/internal_datasets/bph-lulcc___S4Tdata/'
  spath_cru <- '/Volumes/home/work/data/external_datasets/climate/CRU/'
  
}

if(sys['nodename'] == "jeodpp-terminal-151p-02" & 
   sys['user'] == 'duveigr'){ # if it's me, on JEODPP

  spath_s4t <- '/storage/duveigr/internal_datasets/bph-lulcc___S4Tdata/v1.0/'
  spath_cru <- '/storage/duveigr/external_datasets/climate/CRU/'
  
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



#### prepare run ----
dir.create('scratch')
dir.create('results/final_products', recursive = T)

calc_climData <- FALSE

#### get clim data #### ---- 

if(calc_climData){source('src/devel_code/get_climdata.R')}

# load CLIM data
load('results/cleaned_input_data/climate/df4ClimateSpace_1dd.RData') # dfClim

#### run #### ---- 
source('src/devel_code/spclimext_withRF.R')
spclimext_withRF(dfClim, iVar = 'HG', type = 'IGBPgen', do_checkplot = F)
spclimext_withRF(dfClim, iVar = 'LE', type = 'IGBPgen', do_checkplot = F)

## Opt parameters... (set as default)
# nu_thr <- 10
# minTHR <- 0.05
# comprLevel <- 7
# mzVal <- NA
