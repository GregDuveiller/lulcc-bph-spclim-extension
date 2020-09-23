# Script to run an spclim job 

require(here)

#### Initialisation #### ---- 

iVar <- 'SWreflected' # 'HG'  # 'albedo'    # 'LSTnight'  # 'LSTday'  
type <- 'IGBPdet'   # 'IGBPdet'   # 'IGBPgen'
chkp <- F

# set symlinks to connect to the data... 
# source('data/___loadDataPaths___.R')

dir.create('scratch')
dir.create('results/final_products', recursive = T)
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
spclimext_withRF(dfClim, 'SWreflected', type, do_checkplot = chkp)
#spclimext_withRF(dfClim, 'albedo', type, do_checkplot = chkp)
