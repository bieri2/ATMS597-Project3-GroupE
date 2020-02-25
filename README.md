# ATMS597-Project3-GroupE
## Group Members: Carolina Bieri, Arka Mitra, and Jeff Thayer

This python module finds austral summertime (DJF) extreme precipitation days near Jakarta, Indonesia (6.21S, 106.85E) over the period 1996-2019 using daily data from the Global Precipitation Climatology Project (GPCP), and then creates composite maps of several NCEP Reanalysis products to analyze the associated global circulation patterns during these extreme precipitation days. Monthly-mean NCEP Reanalysis products provide a long-term mean for DJF over the period 1981-2010, which is used to construct anomalies. 

Extreme precipitation days are identified as the 95th percentile and above of daily precipitation for 1996-2019 using the grid point nearest to Jakarta, Indonesia during DJF. A cumulative distribution function (CDF) using all precipitation values (including non-extreme values) is calculated. This makes it possible to obtain an idea of the shape of the dataset and compare the 95th percentile value to the rest of the data.

NCEP Reanalysis products (i.e. humidity, temperature, winds, precipitable water, geopotential height) are then compiled for extreme precipitation days and composited to examine typical global circulation patterns associated with these extreme days. 
