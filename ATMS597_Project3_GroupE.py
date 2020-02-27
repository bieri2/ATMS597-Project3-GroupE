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

#Function to select the necessary months of data (DJF)
def is_djf(month):
    return (month == 12) | (month <= 2)

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

count = 0 #loop over years and store daily data into 'data' array
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

#All future analysis will now be done off the combined NetCDF dataset
precip_agg = xr.open_dataset(output_nc_dir+'Aggregate_GCPC_daily_1996_2019.nc')
#Opened and stored the dataset in the xarray dataset precip_agg
    
#The following code-block selects the data for the grid cell closest to Jakarta, Indonesia
#Further that data is subsetted for our requisite months (DJF). Only valid data is retained
#Precip > 200 mm/day are rejected. 95 %-ile rainfall is then calculated as 26 mm/day
#Dates above the 95 %-ile value are reatained in the xarray dataset 'pcp_above_dates'. 
#Precipitation data from Jakarta for these extereme precip days are stored in 'extreme'.
pcp_j = precip_agg.sel(latitude=-6.21,longitude=106.85,method='nearest')
pcp_j = pcp_j.sel(time=is_djf(pcp_j['time.month']))
pcp_j = pcp_j.where(pcp_j != -99999.0)
pcp_j = pcp_j.where(pcp_j < 200.)
pcp_quant = pcp_j.quantile(0.95)
precip = pcp_j.to_array()
pcp_quant = pcp_quant.to_array().values
pcp_above_dates = pcp_j['time'][np.where(precip[0,:]>=pcp_quant[:])]
#print(pcp_above_dates)
extreme = precip[0,:][np.where(precip[0,:]>=pcp_quant[:])]
#print(extreme)


# The following code block calculates and plots and saves the CDF for daily rainfall...
#...for DJF, over Jakarta, Indonesia. The 95%-ile rainfall is specially marked.

#Calculate CDF
nbins = 1000
counts, edges = np.histogram(precip[0,:], bins=nbins, range=(0,200), density = False)
cdf = np.cumsum(counts)/len(precip[0,:])

#Find 95%-ile rainfall, and prepare to plot it over CDF
rain_95 = pcp_quant #approximately
x = np.zeros(100)
y = np.zeros(100)
x[:] = rain_95
y[:] = 0.95

#Plot CDF and 95%-ile rainfall over Jakarta in the required months
plt.figure(figsize=(6,6))
ax = plt.gca()
plt.plot(edges[1:], cdf)
plt.scatter(x,np.arange(0,1.,0.01),color='k',s=6.)
plt.scatter(np.arange(0,70.,70/100),y,color='r',s=6.)
plt.plot(y,y,'k')
plt.ylabel('CDF', fontsize=14)
plt.xlabel('Daily Average Rainfall near Jakarta (mm)', fontsize=14)
major_ticks = np.arange(0, 101, 10)
minor_ticks = np.arange(0, 101, 5)
ax.set_xticks(major_ticks)
ax.set_xticks(minor_ticks, minor=True)
ax.set_yticks(major_ticks/100.)
ax.set_yticks(minor_ticks/100., minor=True)
plt.grid()
ax.grid(which='minor', alpha=0.2)
ax.grid(which='major', alpha=0.5)
plt.xlim(-.50,70.)
plt.ylim(0.3,1.05)
plt.savefig('JAKARTA_DJF_CDF.png', dpi=300)
#plt.show()
