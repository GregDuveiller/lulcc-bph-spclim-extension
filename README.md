# lulcc-bph-spclim-extension 

This repository provides the code necessary to reproduce the update of a specific dataset from its original version (v1.0) to a new version (v2.0) that provides a spatial extension based on climatological similarities. 

## Background
The original dataset provides spatialised estimates of how a series of biophysical environmental variables (e.g. daytime land surface temperature, latent heat, etc.) would locally change following potential changes in vegetation cover. A full description of the dataset is provided in the following document:

Duveiller, G., Hooker, J. & Cescatti, A. A dataset mapping the potential biophysical effects of vegetation cover change. _Sci Data_ 5, 180014 (2018). https://doi.org/10.1038/sdata.2018.14

For an illustration of what can be done with such data, one can refer to the following study:

Duveiller, G., Hooker, J. & Cescatti, A. The mark of vegetation change on Earth’s surface energy balance. _Nat Commun_ 9, 679 (2018). https://doi.org/10.1038/s41467-017-02810-8


## Motivation for version 2
The data-driven methodology needed to generate the dataset imposed a limitation on the spatial extent over which these estimates could be made, as it required both source (what is transitioned from) and target (what is transitioned to) vegetation cover types to co-exist within a local moving window. To provide a more widespread spatial coverage as required by some users, this second version of the dataset has been expanded to places where only one of the two vegetation classes involved in a given transition occurs. This is achieved using  Random Forests, a non-linear regression tree approach, to predict the change in biophysical variables based on local climate. The method employed is described in this study:

Duveiller, G., Caporaso, L., Abad-Viñas, R., Perugini, L., Grassi, G., Arneth, A., & Cescatti, A. Local biophysical effects of land use and land cover change: towards an assessment tool for policy makers. _Land Use Policy_, 91, 104382 (2020). https://doi.org/10.1016/j.landusepol.2019.104382

It was here applied using the updated version ts4.04 of the CRU climate dataset, which is described here:

Harris, I., Osborn, T.J., Jones, P. et al. Version 4 of the CRU TS monthly high-resolution gridded multivariate climate dataset. _Sci Data_ 7, 109 (2020).

With respect to what is done in Duveiller et al (2020), here we also use an extra climatic predictor describing snow presence (as the cumulated precipitation during days with negative temperatures in degree Celsius). The other predictors are mean temperature, temperature range, cumulated precipitation and an aridity index.


## Description
This repository includes the main script ```run_spclimext.R``` that orchestrates the processing for each input variable of the original dataset. Each of these input variables are in individual files, and a corresponding output file is generated. 

The repository also includes some intermediate products to faciliate the ease of use of the code. These are stored in ```results/cleaned_input_data```.

## Input data
The base files of the original dataset (v1.0) can be found here: 
https://data.jrc.ec.europa.eu/dataset/f97e4216-e81c-4dfb-b42d-d0c19070d029

The climate data used is the updated version ts4.04 of the CRU climate dataset:
https://catalogue.ceda.ac.uk/uuid/89e1e34ec3554dc98594a5732622bce9


## New dataset
The new version (v2.0) of the dataset can be found in the JRC Data Catalogue:
https://data.jrc.ec.europa.eu/dataset/f97e4216-e81c-4dfb-b42d-d0c19070d029

The updated (v2.0) data have the same format as before, but now contains an extra variable with an *_ext* suffix, which is the spatially extended version. The original data (without the *_ext* suffix) are also available in each file. Note that over the same location both will differ as the extended version is in effect performing a smoothing operation due to the climate-based interpolation. Finally, the file sizes have been substantially reduced by using the internal compression within the netCDF files.

## Versioning
This repository is available on GitHub at https://github.com/GregDuveiller/lulcc-bph-spclim-extension/. All versions of this repository, including the last release, are archived on Zenodo: [DOI: 10.xxxx/zenodo.xxxxxxx].