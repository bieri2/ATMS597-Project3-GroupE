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

#---------------------------------------Read and Merge Daily Rainfall Data into NetCDF-------------------------------------------------
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
#---------------------------------------End of Code Block----------------------------------------------------------------------------

#---------------------------------------Extreme precip, CDF for Jakarta DJF rainfall-------------------------------------------------
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
plt.title('CDF and 95 percentile rainfall \n for Jakarta, Indonesia (DJF)', fontsize=14)
plt.savefig(output_nc_dir+'JAKARTA_DJF_CDF.png', dpi=300)
#plt.show()
#---------------------------------------End of code block------------------------------------------------------------------------------

#---------------------------------------Extreme day Mean Global Composites (1981-2019)-------------------------------------------------
# Code-block to plot the long-term global composites...
# ...for days meeting and excceding the 95%-ile daily average rainfall value for DJF at Jakarta, Indonesia 
# Read data from combined fields at multiple levels (250 hPa, 500 hPa, 850 hPa, surface) into xarray datasets  
combined_Uwind_250hPa = xr.open_dataset('/content/drive/My Drive/Combined_fields/combined_uwnd_250hPa_1996to2019_ExtremePrecipDays.nc')
combined_Vwind_250hPa = xr.open_dataset('/content/drive/My Drive/Combined_fields/combined_vwnd_250hPa_1996to2019_ExtremePrecipDays.nc')
combined_GeopHgt_500hPa = xr.open_dataset('/content/drive/My Drive/Combined_fields/combined_hgt_500hPa_1996to2019_ExtremePrecipDays.nc')
combined_Uwind_500hPa = xr.open_dataset('/content/drive/My Drive/Combined_fields/combined_uwnd_500hPa_1996to2019_ExtremePrecipDays.nc')
combined_Vwind_500hPa = xr.open_dataset('/content/drive/My Drive/Combined_fields/combined_vwnd_500hPa_1996to2019_ExtremePrecipDays.nc')
combined_Omega_500hPa = xr.open_dataset('/content/drive/My Drive/Combined_fields/combined_omega_500hPa_1996to2019_ExtremePrecipDays.nc')
combined_Uwind_850hPa = xr.open_dataset('/content/drive/My Drive/Combined_fields/combined_uwnd_850hPa_1996to2019_ExtremePrecipDays.nc')
combined_Vwind_850hPa = xr.open_dataset('/content/drive/My Drive/Combined_fields/combined_vwnd_850hPa_1996to2019_ExtremePrecipDays.nc')
combined_SpecHum_850hPa = xr.open_dataset('/content/drive/My Drive/Combined_fields/combined_shum_850hPa_1996to2019_ExtremePrecipDays.nc')
combined_AirTemp_850hPa = xr.open_dataset('/content/drive/My Drive/Combined_fields/combined_air_850hPa_1996to2019_ExtremePrecipDays.nc')
combined_Uwind_10m = xr.open_dataset('/content/drive/My Drive/Combined_fields/combined_uwnd_10m_1996to2019_ExtremePrecipDays.nc')
combined_Vwind_10m = xr.open_dataset('/content/drive/My Drive/Combined_fields/combined_vwnd_10m_1996to2019_ExtremePrecipDays.nc')
combined_SkinTemp_Sfc = xr.open_dataset('/content/drive/My Drive/Combined_fields/combined_air_1000hPa_1996to2019_ExtremePrecipDays.nc')
combined_PrecipWater = xr.open_dataset('/content/drive/My Drive/Combined_fields/combined_pr_wtr_1996to2019_ExtremePrecipDays.nc')

#Rotate longitudes cyclically by 180 degrees to make plots start at -180 degrees.
lats = combined_Uwind_250hPa['lat'][:]
lons = combined_Uwind_250hPa['lon'][:]
lons_cyclic = np.roll(lons,71)

# level = 250 hPa
# Fields --> Wind barbs and speeds
uwnd_level = combined_Uwind_250hPa
uwnd_data = np.roll(uwnd_level['uwnd'].mean('time'),72,axis=1)
vwnd_level = combined_Vwind_250hPa
vwnd_data = np.roll(vwnd_level['vwnd'].mean('time'),72,axis=1)
winds = np.sqrt(uwnd_data[:,:]**2+vwnd_data[:,:]**2)
plt.figure(figsize=(16,16))
ax = plt.axes(projection=ccrs.PlateCarree())
plt.pcolormesh(lons_cyclic, lats, winds, 
              transform=ccrs.PlateCarree(), cmap='winter')
ax.barbs(lons_cyclic[::4], lats[::2], uwnd_data[::2,::4],
          vwnd_data[::2,::4], length=4)
# ax.quiver(lons_cyclic[::4], lats[::2], uwnd_data[::2,::4],
#           vwnd_data[::2,::4])
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.5, linestyle='-')
#set where the gridlines go
gl.xlocator = mticker.FixedLocator(np.arange(-180,181,30))
gl.ylocator = mticker.FixedLocator(np.arange(-90,91,30))
gl.xlabels_top = False
gl.ylabels_right = False
plt.title('250-hPa Wind Speeds (m/s) and Wind Barbs \n for extreme precipitation days over Jakarta, Indonesia', fontsize=15)
ax.coastlines()
cb = plt.colorbar(orientation='horizontal', shrink=0.5, pad=0.07)
cb.set_label(r'Winds speeds [$m s^{-1}$] at 250 hPa', fontsize=14)
plt.clim(0,40)
plt.show()

# level = 500 hPa
# Fields --> Wind barbs and Geopotential Heights
# Fields --> Wind barbs and Vertical Vorticity
uwnd_level = combined_Uwind_500hPa
uwnd_data = np.roll(uwnd_level['uwnd'].mean('time'),72,axis=1)
vwnd_level = combined_Vwind_500hPa
vwnd_data = np.roll(vwnd_level['vwnd'].mean('time'),72,axis=1)
hgt_level = combined_GeopHgt_500hPa
hgt_data = np.roll(hgt_level['hgt'].mean('time'),72,axis=1)
omega_level = combined_Omega_500hPa
omega_data = np.roll(omega_level['omega'].mean('time'),72,axis=1)
plt.figure(figsize=(16,16))
ax = plt.axes(projection=ccrs.PlateCarree())
cs2 = plt.pcolormesh(lons_cyclic, lats, hgt_data, 
              transform=ccrs.PlateCarree(), cmap='jet')
ax.barbs(lons_cyclic[::4], lats[::2], uwnd_data[::2,::4],
          vwnd_data[::2,::4], length=4)
# ax.quiver(lons_cyclic[::4], lats[::2], uwnd_data[::2,::4],
#           vwnd_data[::2,::4])
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.5, linestyle='-')
#set where the gridlines go
gl.xlocator = mticker.FixedLocator(np.arange(-180,181,30))
gl.ylocator = mticker.FixedLocator(np.arange(-90,91,30))
gl.xlabels_top = False
gl.ylabels_right = False
plt.title('500-hPa Geopotential Heights (m),\n Wind Speeds (m/s) and Wind Barbs \n for extreme precipitation days over Jakarta, Indonesia',
          fontsize=15)
ax.coastlines()
cb = plt.colorbar(cs2, orientation='horizontal', shrink=0.5, pad=0.07)
cb.set_label(r'Geopotential heights [m] at 500 hPa', fontsize=14)
plt.clim(5000.,6000.)
plt.show()

plt.figure(figsize=(16,16))
ax = plt.axes(projection=ccrs.PlateCarree())
cs2 = plt.pcolormesh(lons_cyclic, lats, omega_data, 
              transform=ccrs.PlateCarree(), cmap='seismic')
ax.barbs(lons_cyclic[::4], lats[::2], uwnd_data[::2,::4],
          vwnd_data[::2,::4], length=4)
# ax.quiver(lons_cyclic[::4], lats[::2], uwnd_data[::2,::4],
#           vwnd_data[::2,::4])
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.5, linestyle='-')
#set where the gridlines go
gl.xlocator = mticker.FixedLocator(np.arange(-180,181,30))
gl.ylocator = mticker.FixedLocator(np.arange(-90,91,30))
gl.xlabels_top = False
gl.ylabels_right = False
plt.title('500-hPa Vertical Vorticity [1/s],\n Wind Speeds (m/s) and Wind Barbs \n for extreme precipitation days over Jakarta, Indonesia',
          fontsize=15)
ax.coastlines()
cb = plt.colorbar(cs2, orientation='horizontal', shrink=0.5, pad=0.07)
cb.set_label(r'Vertical Vorticity [$s^{-1}$] at 500 hPa', fontsize=14)
plt.show()

# level = 850 hPa
# Fields --> Wind barbs and Air Temperatures
# Fields --> Wind barbs and Specific Humidity
uwnd_level = combined_Uwind_850hPa
uwnd_data = np.roll(uwnd_level['uwnd'].mean('time'),72,axis=1)
vwnd_level = combined_Vwind_850hPa
vwnd_data = np.roll(vwnd_level['vwnd'].mean('time'),72,axis=1)
shum_level = combined_SpecHum_850hPa
shum_data = np.roll(shum_level['shum'].mean('time'),72,axis=1)
air_level = combined_AirTemp_850hPa
air_data = np.roll(air_level['air'].mean('time'),72,axis=1)
plt.figure(figsize=(16,16))
ax = plt.axes(projection=ccrs.PlateCarree())
cs2 = plt.pcolormesh(lons_cyclic, lats, air_data, 
              transform=ccrs.PlateCarree(), cmap='jet')
ax.barbs(lons_cyclic[::4], lats[::2], uwnd_data[::2,::4],
          vwnd_data[::2,::4], length=4)
# ax.quiver(lons_cyclic[::4], lats[::2], uwnd_data[::2,::4],
#           vwnd_data[::2,::4])
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.5, linestyle='-')
#set where the gridlines go
gl.xlocator = mticker.FixedLocator(np.arange(-180,181,30))
gl.ylocator = mticker.FixedLocator(np.arange(-90,91,30))
gl.xlabels_top = False
gl.ylabels_right = False
plt.title('850-hPa Temperatures [K],\n Wind Speeds (m/s) and Wind Barbs \n for extreme precipitation days over Jakarta, Indonesia',
          fontsize=15)
ax.coastlines()
cb = plt.colorbar(cs2, orientation='horizontal', shrink=0.5, pad=0.07)
cb.set_label(r'Temperatures [K] at 850 hPa', fontsize=14)
plt.show()

# All levels integrated
# Fields --> Precipitable Water in the entire atmospheric column

pwtr_level = combined_PrecipWater
pwtr_data = np.roll(pwtr_level['pr_wtr'].mean('time'),72,axis=1)
plt.figure(figsize=(16,16))
ax = plt.axes(projection=ccrs.PlateCarree())
cs2 = plt.pcolormesh(lons_cyclic, lats, air_data, 
              transform=ccrs.PlateCarree(), cmap='GnBu')
# ax.quiver(lons_cyclic[::4], lats[::2], uwnd_data[::2,::4],
#           vwnd_data[::2,::4])
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.5, linestyle='-')
#set where the gridlines go
gl.xlocator = mticker.FixedLocator(np.arange(-180,181,30))
gl.ylocator = mticker.FixedLocator(np.arange(-90,91,30))
gl.xlabels_top = False
gl.ylabels_right = False
plt.title('Total Precipitable Water [$kg m^{-2}$]\n for extreme precipitation days over Jakarta, Indonesia',
          fontsize=15)
ax.coastlines()
cb = plt.colorbar(cs2, orientation='horizontal', shrink=0.5, pad=0.07)
cb.set_label(r'Total Precipitable Water [kg $m^{-2}$]', fontsize=14)
plt.show()

# level = surface
# Fields --> Wind barbs and Surface Skin temperatures

lats = combined_Uwind_10m['lat'][:]
lons = combined_Uwind_10m['lon'][:]
lons_cyclic = np.roll(lons,96)
uwnd_level = combined_Uwind_10m
uwnd_data = np.roll(uwnd_level['uwnd'].mean('time'),96,axis=1)
vwnd_level = combined_Vwind_10m
vwnd_data = np.roll(vwnd_level['vwnd'].mean('time'),96,axis=1)
air_level = combined_SkinTemp_Sfc
air_data = np.roll(air_level['skt'].mean('time'),96,axis=1)
plt.figure(figsize=(16,16))
ax = plt.axes(projection=ccrs.PlateCarree())
cs2 = plt.pcolormesh(lons_cyclic, lats, air_data, 
              transform=ccrs.PlateCarree(), cmap='jet')
ax.barbs(lons_cyclic[::4], lats[::2], uwnd_data[::2,::4],
          vwnd_data[::2,::4], length=4)
# ax.quiver(lons_cyclic[::4], lats[::2], uwnd_data[::2,::4],
#           vwnd_data[::2,::4])
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.5, linestyle='-')
#set where the gridlines go
gl.xlocator = mticker.FixedLocator(np.arange(-180,181,30))
gl.ylocator = mticker.FixedLocator(np.arange(-90,91,30))
gl.xlabels_top = False
gl.ylabels_right = False
plt.title('Surface Temperatures [K],\n Wind Speeds (m/s) and Wind Barbs \n for extreme precipitation days over Jakarta, Indonesia',
          fontsize=15)
ax.coastlines()
cb = plt.colorbar(cs2, orientation='horizontal', shrink=0.5, pad=0.07)
cb.set_label(r'Surface Temperatures [K]', fontsize=14)
plt.show()

#---------------------------------------End of Code Block----------------------------------------------------------------------------


#---------------------------------------Long term mean plots (1981-2010)-------------------------------------------------------------
## Importing Long-term Mean NCEP Reanalysis data from 1981-2010, only selecting DJF

# Get list of dates to be used when selecting data 
months = xr.cftime_range(start = '0001-01-01', end = '0001-12-01', freq = 'MS', calendar = 'standard') # selecting long-term mean data for 1981-2010 with cftime.DatetimeGregorian format
months = months[(months.month == 12) | (months.month == 1) | (months.month == 2)] # selecting long-term mean data for DJF in 1981-2010

# Define URLs from which to download data
url_p = 'https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.derived/pressure/'
url_s = 'https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.derived/surface/'

# Extract data - winds, geopotential height, specific humidity, temperature, precipitable water at different heights 
ds_Uwind_250hPa_LTM   = xr.open_dataset(url_p + 'uwnd.mon.1981-2010.ltm.nc', engine  = 'netcdf4').sel(level = 250, time = months)
ds_Vwind_250hPa_LTM   = xr.open_dataset(url_p + 'vwnd.mon.1981-2010.ltm.nc', engine  = 'netcdf4').sel(level = 250, time = months)
ds_GeopHgt_500hPa_LTM = xr.open_dataset(url_p + 'hgt.mon.1981-2010.ltm.nc', engine   = 'netcdf4').sel(level = 500, time = months)
ds_Uwind_500hPa_LTM   = xr.open_dataset(url_p + 'uwnd.mon.1981-2010.ltm.nc', engine  = 'netcdf4').sel(level = 500, time = months)
ds_Vwind_500hPa_LTM   = xr.open_dataset(url_p + 'vwnd.mon.1981-2010.ltm.nc', engine  = 'netcdf4').sel(level = 500, time = months)
ds_Uwind_850hPa_LTM   = xr.open_dataset(url_p + 'uwnd.mon.1981-2010.ltm.nc', engine  = 'netcdf4').sel(level = 850, time = months)
ds_Vwind_850hPa_LTM   = xr.open_dataset(url_p + 'vwnd.mon.1981-2010.ltm.nc', engine  = 'netcdf4').sel(level = 850, time = months)
ds_SpecHum_850hPa_LTM = xr.open_dataset(url_p + 'shum.mon.1981-2010.ltm.nc', engine  = 'netcdf4').sel(level = 850, time = months)
ds_AirTemp_850hPa_LTM = xr.open_dataset(url_p + 'air.mon.1981-2010.ltm.nc', engine   = 'netcdf4').sel(level = 850, time = months)
ds_Uwind_sig995_LTM   = xr.open_dataset(url_s + 'uwnd.sig995.mon.1981-2010.ltm.nc', engine    = 'netcdf4').sel(time = months)
ds_Vwind_sig995_LTM   = xr.open_dataset(url_s + 'vwnd.sig995.mon.1981-2010.ltm.nc', engine    = 'netcdf4').sel(time = months)
ds_SkinTemp_Sfc_LTM   = xr.open_dataset(url_s + '_gauss/skt.sfc.mon.1981-2010.ltm.nc', engine = 'netcdf4').sel(time = months)
ds_PrecipWater_LTM    = xr.open_dataset(url_s + 'pr_wtr.eatm.mon.1981-2010.ltm.nc', engine    = 'netcdf4').sel(time = months)

#---------------------------------------Plotting 250 mb wind speed and direction------------------------------------------------------
# Plot 250 mb wind speed and direction
fig = plt.figure(figsize = (15, 8))
ax  = fig.add_subplot(1, 1, 1, projection = ccrs.PlateCarree(central_longitude = 180))
ax.coastlines()

# Calculate means over DJF
uwind_djf_ltm_250 = ds_Uwind_250hPa_LTM.mean(dim = 'time')['uwnd']
vwind_djf_ltm_250 = ds_Vwind_250hPa_LTM.mean(dim = 'time')['vwnd']

# Calculate wind speed using U and V components
wind_speed_250    = np.sqrt(uwind_djf_ltm_250**2 + vwind_djf_ltm_250**2)

# Plot filled contours
c = ax.contourf(wind_speed_250.lon, wind_speed_250.lat,
            wind_speed_250, 10, transform=ccrs.PlateCarree(), cmap='RdPu', alpha=0.8)

# Plot winds using barbs
ax.barbs(uwind_djf_ltm_250.lon.values[::4], uwind_djf_ltm_250.lat.values[::4], 
         uwind_djf_ltm_250[::4,::4], vwind_djf_ltm_250[::4,::4], 
         length = 5, sizes = dict(emptybarb = 0.25, spacing = 0.2, height = 0.5),
         linewidth = 0.95, alpha = 0.5)

# Plot color bar
cbar = fig.colorbar(c)
cbar.ax.tick_params(labelsize=14) 
cbar.set_label(label='Wind speed (m/s)',fontsize=14)

plt.title('250 mb wind speed and direction - DJF', fontsize=15)
#plt.show()

#---------------------------------------Plotting 500 mb wind speed and direction------------------------------------------------------
# Plot 500 mb geopotential height and winds
fig = plt.figure(figsize = (15, 8))
ax  = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(central_longitude = 180))
ax.coastlines()

# Calculate means for DJF
uwind_djf_ltm_500 = ds_Uwind_500hPa_LTM.mean(dim = 'time')['uwnd']
vwind_djf_ltm_500 = ds_Vwind_500hPa_LTM.mean(dim = 'time')['vwnd']
gh_djf_ltm_500    = ds_GeopHgt_500hPa_LTM.mean(dim = 'time')['hgt']

# Plot filled contours
c = ax.contourf(gh_djf_ltm_500.lon, gh_djf_ltm_500.lat,
            gh_djf_ltm_500, 10, transform=ccrs.PlateCarree(), cmap='Purples', alpha=0.8)

# Plot winds using barbs
ax.barbs(uwind_djf_ltm_500.lon.values[::4], uwind_djf_ltm_500.lat.values[::4], 
         uwind_djf_ltm_500[::4,::4], vwind_djf_ltm_500[::4,::4], 
         length = 5,sizes = dict(emptybarb = 0.25, spacing = 0.2, height = 0.5),
         linewidth = 0.95, alpha = 0.6)

# Plot color bar
cbar = fig.colorbar(c)
cbar.ax.tick_params(labelsize=14) 
cbar.set_label(label = 'Geopotential height (m)',fontsize=14)

plt.title('500 mb geopotential height and winds - DJF', fontsize=15)
#plt.show()

#---------------------------------------Plotting 850 mb wind speed, temperature, and humidity------------------------------------------

# Plot 850 mb temperature, specific humidity, and winds

# Plot winds and temperature 
fig = plt.figure(figsize = (15, 15))
ax  = fig.add_subplot(2, 1, 1, projection = ccrs.PlateCarree(central_longitude = 180))
ax.coastlines()

# Calculate mean values for DJF
uwind_djf_ltm_850 = ds_Uwind_850hPa_LTM.mean(dim = 'time')['uwnd']
vwind_djf_ltm_850 = ds_Vwind_850hPa_LTM.mean(dim = 'time')['vwnd']
t_djf_ltm_850     = ds_AirTemp_850hPa_LTM.mean(dim = 'time')['air']

# Plot filled contours 
c = ax.contourf(t_djf_ltm_850.lon, t_djf_ltm_850.lat,
            t_djf_ltm_850 + 273.15,10, transform = ccrs.PlateCarree(), cmap = 'Reds', alpha=0.8)

# Plot winds using barbs
ax.barbs(uwind_djf_ltm_850.lon.values[::4], uwind_djf_ltm_850.lat.values[::4], 
         uwind_djf_ltm_850[::4, ::4], vwind_djf_ltm_850[::4, ::4], 
         length = 5,sizes = dict(emptybarb = 0.25, spacing = 0.2, height = 0.5),
         linewidth = 0.95, alpha = 0.6)

cbar = fig.colorbar(c)
cbar.ax.tick_params(labelsize = 14) 
cbar.set_label(label = 'Air temperature (K)', fontsize = 14)

plt.title('850 mb temperature and winds - DJF', fontsize = 15)

# Plot specific humidity
ax2 = fig.add_subplot(2, 1, 2, projection=ccrs.PlateCarree(central_longitude = 180))

# Calculate mean for DJF
sh_djf_ltm_850 = ds_SpecHum_850hPa_LTM.mean(dim = 'time')['shum']
ax2.coastlines()

# Plot filled contours
c2 = ax2.contourf(sh_djf_ltm_850.lon, sh_djf_ltm_850.lat,
            sh_djf_ltm_850, 10, transform = ccrs.PlateCarree(), cmap = 'Blues', alpha = 0.8)

# Plot winds using barbs
ax2.barbs(uwind_djf_ltm_850.lon.values[::4], uwind_djf_ltm_850.lat.values[::4], 
         uwind_djf_ltm_850[::4, ::4], vwind_djf_ltm_850[::4, ::4], 
         length = 5,sizes = dict(emptybarb = 0.25, spacing = 0.2, height = 0.5),
         linewidth = 0.95, alpha = 0.6)

cbar2 = fig.colorbar(c2)
cbar2.ax.tick_params(labelsize = 14) 
cbar2.set_label(label = 'Specific humidity (g/kg)', fontsize = 14)

ax2.set_title('850 mb specific humidity and winds - DJF', fontsize = 15)

#plt.show()

#---------------------------------------Plotting surface temperature and winds---------------------------------------------------------

fig = plt.figure(figsize = (15, 8))
ax = fig.add_subplot(1, 1, 1, projection = ccrs.PlateCarree(central_longitude = 180))
ax.coastlines()

# Calculate means for DJF
uwind_djf_ltm_sfc = ds_Uwind_sig995_LTM.mean(dim = 'time')['uwnd']
vwind_djf_ltm_sfc = ds_Vwind_sig995_LTM.mean(dim = 'time')['vwnd']
t_sfc = ds_SkinTemp_Sfc_LTM.mean(dim = 'time')['skt']

# Plot filled contours
c = ax.contourf(t_sfc.lon, t_sfc.lat,
                t_sfc, 10, transform=ccrs.PlateCarree(), cmap = 'rainbow', alpha = 0.8)

# Plot winds using barbs
ax.barbs(uwind_djf_ltm_sfc.lon.values[::4], uwind_djf_ltm_sfc.lat.values[::4], 
         uwind_djf_ltm_sfc[::4, ::4], vwind_djf_ltm_sfc[::4, ::4], 
         length = 5, sizes = dict(emptybarb = 0.25, spacing = 0.2, height = 0.5),
         linewidth = 0.95, alpha = 0.5)

cbar = fig.colorbar(c)
cbar.ax.tick_params(labelsize = 14) 
cbar.set_label(label = 'Skin temperature (deg C)',fontsize = 14)

plt.title('Surface skin temperature and winds - DJF', fontsize = 15)
#plt.show()

#---------------------------------------Plotting precipitable water---------------------------------------------------------------------
fig = plt.figure(figsize = (15, 8))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(central_longitude = 180))
ax.coastlines()

pw = ds_PrecipWater_LTM.mean(dim = 'time')['pr_wtr']
pw = pw.where(pw > 0)

c = ax.contourf(pw.lon, pw.lat, pw, 10, transform = ccrs.PlateCarree(), cmap = 'Greens', alpha = 0.8)

cbar = fig.colorbar(c)
cbar.ax.tick_params(labelsize = 14) 
cbar.set_label(label = 'Precipitable water ($kg/m^{2}$)',fontsize = 14) 

plt.title('Precipitable water - DJF', fontsize = 15)
#plt.show()
