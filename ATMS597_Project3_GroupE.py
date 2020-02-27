# libraries and modules to be imported 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt  
import glob
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker

#Code block to read from all the daily GPCP NetCDF files...
#...and merge them into a single NetCDF file

#Files must be downloaded from online server with the same directory structure...
#...using the "wget -r -nH" command 
years = np.arange(1996,2020)
gpcp_daily_data_directory = "/content/drive/My Drive/noaa_data_project3/data/global-precipitation-climatology-project-gpcp-daily/access/"
output_nc_dir = "/content/drive/My Drive/"

#Calculate total number of days in dataset 
#This takes care of any duplicate values
count = 0
times = []
for i in range(0,len(years)): 
  nc_files = sorted(glob.glob(gpcp_daily_data_directory + str(years[i]) + "/*.nc"))
  for n in range(0,len(nc_files)):
    filename = nc_files[n].replace(gpcp_daily_data_directory+str(years[i]),'')
    date = filename.split('_')[3].replace('d','')
    date = pd.to_datetime(date, format='%Y%m%d')
    times = np.append(times, date)
    
latitude = np.arange(-90.0,90.0)
longitude = np.arange(0.0,360.0)
data = np.zeros((len(times),len(latitude),len(longitude),1))

count = 0
for i in range(0, len(years)): 
  #print(i)
  nc_files = glob.glob(gpcp_daily_data_directory + str(years[i]) + "/*.nc")
  for n in range(0, len(nc_files)):
    try:
      nc = xr.open_dataset(nc_files[n])
      ncvar = nc['precip']
      precip = ncvar.mean(axis=0)
      data[count,:,:,0] = precip
    except:
      continue
    count = count + 1
data1 = data.squeeze(axis=3)
precip_agg = xr.DataArray(data1, coords=[times, latitude, longitude], dims=['time', 'latitude', 'longitude'])
precip_agg.to_netcdf(output_nc_dir+'Aggregate_GCPC_daily_1996_2019.nc') #Output NetCDF file generated and saved

