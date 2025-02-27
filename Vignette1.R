# Leveraging modern quantitative tools to improve wildlife decision support products for natural resources agencies
# Author: Jennifer Stenglein

# Vignette 1: WSI using publicly-available data and spatial interpolation

# code to integrate publicly-available NOAA data for use in spatial interpolation 
# to estimate Winter Severity Index for each Deer Management Unit in Wisconsin
# example is for winter 2021 - 2022 and recreates map in Figure 3 

library(tidyverse)
library(sp)
library(sf)
library(raster)
library(gstat)

# read in data
# these datasets already include publicly-available data from the Climate Data Online portal (https://www.ncdc.noaa.gov/cdo-web/)
# 'source' column indicates whether record is from DNR staff observation or NOAA publicly-available data
df_temp <- readRDS("data/df_Temp.Rds") 
coordinates(df_temp) <- ~ longitude + latitude
proj4string(df_temp) <- crs("+proj=longlat") 

df_snow <- readRDS("data/df_Snow.Rds") 
coordinates(df_snow) <- ~ longitude + latitude
proj4string(df_snow) <- crs("+proj=longlat")

# read in shapefile
dmu <- readRDS("data/dmus.Rds") %>%
  st_as_sf()

# create vector of winter dates
dates <- seq(as.Date("2021-12-01"),as.Date("2022-04-30"),1)

# create an empty grid 
longitude <- seq(-93.5, -86.5, .1)
latitude <- seq(42, 47.5, .1)
grid <- expand.grid(longitude=longitude, latitude=latitude) 
coordinates(grid) <- ~longitude + latitude
proj4string(grid) <- crs("+proj=longlat")
gridded(grid) <- TRUE

# create an empty RasterStack to hold the output rasters from the below loop
r_stack_temp <- raster::stack()
r_stack_snow <- raster::stack()

# set up the loop to loop over dates
for(i in 1:length(dates)) {
  
  # filter to date
  temp_filter <- df_temp[df_temp$date==dates[i],]
  snow_filter <- df_snow[df_snow$date==dates[i],]
  
  # print date
  print(dates[i])
  
  # run the idw
  this_idw_temp <- idw(formula = Temp.Pts~1,
                       locations = temp_filter,
                       newdata = grid,
                       idp = 2.0) 
  
  this_idw_snow <- idw(formula = Snow.Pts~1,
                       locations = snow_filter,
                       newdata = grid,
                       idp = 2.0) 
  
  # convert this_idw to this_raster 
  this_raster_temp <- raster(this_idw_temp)
  this_raster_snow <- raster(this_idw_snow)
  
  # add to RasterStack
  r_stack_temp <- addLayer(r_stack_temp,this_raster_temp)
  r_stack_snow <- addLayer(r_stack_snow,this_raster_snow)
}

# rename rasters in RasterStack to dates 
names(r_stack_temp) <- gsub('[-]', '.', x=dates)
names(r_stack_snow) <- gsub('[-]', '.', x=dates)

# sum the rasters in each RasterStack
r_sum_temp <- sum(r_stack_temp)
r_sum_snow <- sum(r_stack_snow)

# sum the temp and snow RasterStacks to get WSI raster
r_stack_wsi <- raster::stack(r_sum_temp,r_sum_snow) 
r_sum_wsi <- sum(r_stack_wsi)
plot(r_sum_wsi)

# extract mean pixel values by DMU
wsi_by_dmu <- extract(r_sum_wsi, dmu, fun=mean, na.rm=TRUE, df=TRUE)
wsi_by_dmu <- cbind(dmu$DMU,wsi_by_dmu)
colnames(wsi_by_dmu) = c("DMU","ID","WSI")
wsi_by_dmu <- wsi_by_dmu %>%
  mutate(WSI = as.numeric(WSI))

# plot
dmu = dmu %>%
  left_join(wsi_by_dmu)
ggplot(dmu) +
  geom_sf(aes(fill=WSI)) +
  ggthemes::theme_map() +
  scale_fill_gradient2(limits = c(min(dmu$WSI),max(dmu$WSI)),
                       low = "lightsteelblue1",
                       mid = "royalblue3",
                       high = "maroon1",
                       midpoint = (min(dmu$WSI)+max(dmu$WSI))/2)
