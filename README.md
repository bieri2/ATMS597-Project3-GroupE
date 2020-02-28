# ATMS597-Project3-GroupE
## Group Members: Carolina Bieri, Arka Mitra, and Jeff Thayer

This Python script finds austral summertime (DJF) extreme precipitation days near Jakarta, Indonesia (6.21S, 106.85E) over the period 1996-2019 using daily data from the Global Precipitation Climatology Project (GPCP), and then creates composite maps of several NCEP Reanalysis (Kalnay et al. 1996) products to analyze the associated global circulation patterns during these extreme precipitation days. Monthly-mean NCEP Reanalysis products provide a long-term mean for DJF over the period 1981-2010, which is used to construct anomalies. 

Extreme precipitation days are identified as the 95th percentile and above of daily precipitation for 1996-2019 using the grid point nearest to Jakarta, Indonesia during DJF. A cumulative distribution function (CDF) using all precipitation values (including non-extreme values) is calculated. This makes it possible to obtain an idea of the shape of the dataset and compare the 95th percentile value to the rest of the data.

NCEP Reanalysis products (i.e. humidity, temperature, winds, precipitable water, geopotential height) are then compiled for extreme precipitation days and composited to examine typical patterns associated with these extreme days. Plots are created for these variables to facilitate analysis. Daily means as well as mean daily anomalies are generated. In order to obtain an understanding of general long-term patterns, DJF means for 1981-2010 are also plotted.

GPCP data were obtained via this website: https://www.ncei.noaa.gov/data/global-precipitation-climatology-project-gpcp-daily/access/

NCEP Reanalysis data were obtained here: https://www.esrl.noaa.gov/psd/thredds/catalog/Datasets/catalog.html

Data were downloaded using tools such as `wget` and `xarray.open_dataset`. `xarray` was also used extensively to reduce data. `cartopy` functions and methods were used to create plots. 
