---
title: "Addendum: A dataset mapping the potential biophysical effects of vegetation cover change"
author: "G. Duveiller & A. Cescatti"
bibliography: addendum.bib
---

The original dataset [@Duveiller2018] provides a spatialised estimation of how a series of biophysical environmental variables (e.g. daytime land surface temperature, latent heat, etc.) would locally change following a potential change in vegetation cover from a source to a target vegetation class.
The data-driven methodology that is employed imposes a limitation on the spatial extent over which this estimation is provided, as it requires that both target and source vegetation cover types co-exist within a local moving window.
To provide a more exhaustive spatial coverage as required by some users, this second version of the dataset was expanded to places where only one of the two vegetation classes involved in a given transition occurs.
This is realised using a concept of climatic envelope, and using a machine learning technique known as a Random Forest [@Breimann] to infer the change in biophysical variables based on local climate similarity.
The method employed is described in [@Duveiller2020], but is here applied using an updated version of the climate data (CRU ts4.04 dataset []) and an extra climatic indicator describing snow presence as the cumulated precipitation during days with negative temperatures in degree Celsius  (the other indicators are mean temperature, temperature range, cumulated precipitation and aridity index).

This new version of the dataset is available here: [https://data.jrc.ec.europa.eu/]
The data has the same format as before, but now contains an extra variable with an *_ext* suffix, which is the spatially extended version. The original data (without the *_ext* suffix) are also available in each file. Note that over common locations both will differ as the extended version is in effect  performing a smoothing operation due to the climate-based interpolation. Finally, the file size have been substantially reduced by using the internal compression of the netcdf files.

The code to reproduce these operations is available here: [URL]
