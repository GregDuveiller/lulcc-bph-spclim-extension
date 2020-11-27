# Script to run the spclim job 

require(here)

#### prepare run ----
source('code/setup_datafiles.R')

dir.create('scratch', showWarnings = F)
dir.create('data', showWarnings = F)
dir.create('results/final_products', showWarnings = F, recursive = T)

calc_climData <- FALSE

#### get clim data #### ---- 

if(calc_climData){source('code/get_climdata.R')}

# load CLIM data
load('results/cleaned_input_data/climate/df4ClimateSpace_1dd.RData') # dfClim

#### run #### ---- 
source('code/spclimext_withRF.R')

iVars <- c('HG', 'LE', 'LSTday', 'LSTnight', 
           'LWemitted', 'LWsfc', 'SWreflected', 'albedo')

for(iVar in iVars){
  for(type in c('IGBPgen', 'IGBPdet')){
    spclimext_withRF(dfClim, iVar = iVar, type = type, do_checkplot = F)
  }
}


