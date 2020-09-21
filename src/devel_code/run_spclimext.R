# Script to run an spclim job 

require(here)

#### Initialisation #### ---- 

iVar <- 'albedo'    # 'LSTnight'  # 'LSTday'  #  
type <- 'IGBPgen'   # 'IGBPdet'   # 'IGBPgen'

## Opt parameters... (set as default)
# nu_thr <- 10
# minTHR <- 0.05
# comprLevel <- 7
# mzVal <- NA


#### load data #### ---- 

# load CLIM data
load('results/cleaned_input_data/climate/df4ClimateSpace_1dd.RData') # dfClim

#### run #### ---- 
source('src/devel_code/spclimext_withRF.R')
spclimext_withRF(dfClim, iVar, type, do_checkplot = T)
