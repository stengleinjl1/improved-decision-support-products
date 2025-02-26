# Leveraging modern quantitative tools to improve wildlife decision support products for natural resources agencies

**Author:** Jennifer Stenglein
**Date:** 2025-02-26

The data files and code in this repository accompany the analyses in Stenglein, J. L., Bemowski, R., Storm, D. J., and D. MacFarland. Leveraging modern quantitative tools to improve wildlife decision support products for natural resources agencies. 

## Vignette 1: Winter severity index (WSI) using publicly-available NOAA data for use in spatial interpolation

Data files required:

- df_Temp.Rds
- df_Snow.Rds
- dmus.Rds

Vignette1.R uses daily minimum temperature and snow depth data from locations in and around Wisconsin, USA from December 2021 - April 2022. The code creates daily minimum temperature and snow depth rasters using inverse distance weighting. These daily rasters are summed to get a minimum temperature and snow depth raster for the entire winter. These 2 rasters are summed to get a WSI raster for winter 2021-2022. The raster pixels are averaged to get a WSI value for each Deer Management Unit. Those results are plotted to produce the map found in Figure 3. 



