## Import NCEP Reanalysis fields directly from the NOAA THREDDS server, selecting only the extreme precipitation days

years = np.arange(1996,2020)

datasets_Uwind_250hPa = []
datasets_Vwind_250hPa = []
datasets_GeopHgt_500hPa = []
datasets_Uwind_500hPa = []
datasets_Vwind_500hPa = []
datasets_Omega_500hPa = []
datasets_Uwind_850hPa = []
datasets_Vwind_850hPa = []
datasets_SpecHum_850hPa = []
datasets_AirTemp_850hPa = []
datasets_Uwind_sig995 = []
datasets_Vwind_sig995 = []
datasets_SkinTemp_Sfc = []
datasets_PrecipWater = []

for iyr in years:

    print('working on '+str(iyr))

    if len(ExtremePrecip_dates[ExtremePrecip_dates['time.year'].values ==iyr]) > 0:
      dates_year = ExtremePrecip_dates[ExtremePrecip_dates['time.year'].values ==iyr]

      # Extract data
      ds_Uwind_250hPa = xr.open_dataset('https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.dailyavgs/pressure/uwnd.'+str(iyr)+'.nc',engine='netcdf4').sel(\
                            level=250, time=dates_year)
      ds_Vwind_250hPa = xr.open_dataset('https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.dailyavgs/pressure/vwnd.'+str(iyr)+'.nc',engine='netcdf4').sel(\
                            level=250, time=dates_year)
      ds_GeopHgt_500hPa = xr.open_dataset('https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.dailyavgs/pressure/hgt.'+str(iyr)+'.nc',engine='netcdf4').sel(\
                            level=500, time=dates_year)
      ds_Uwind_500hPa = xr.open_dataset('https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.dailyavgs/pressure/uwnd.'+str(iyr)+'.nc',engine='netcdf4').sel(\
                            level=500, time=dates_year)
      ds_Vwind_500hPa = xr.open_dataset('https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.dailyavgs/pressure/vwnd.'+str(iyr)+'.nc',engine='netcdf4').sel(\
                            level=500, time=dates_year)
      ds_Omega_500hPa = xr.open_dataset('https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.dailyavgs/pressure/omega.'+str(iyr)+'.nc',engine='netcdf4').sel(\
                            level=500, time=dates_year)
      ds_Uwind_850hPa = xr.open_dataset('https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.dailyavgs/pressure/uwnd.'+str(iyr)+'.nc',engine='netcdf4').sel(\
                            level=850, time=dates_year)
      ds_Vwind_850hPa = xr.open_dataset('https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.dailyavgs/pressure/vwnd.'+str(iyr)+'.nc',engine='netcdf4').sel(\
                            level=850, time=dates_year)
      ds_SpecHum_850hPa = xr.open_dataset('https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.dailyavgs/pressure/shum.'+str(iyr)+'.nc',engine='netcdf4').sel(\
                            level=850, time=dates_year)
      ds_AirTemp_850hPa = xr.open_dataset('https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.dailyavgs/pressure/air.'+str(iyr)+'.nc',engine='netcdf4').sel(\
                            level=850, time=dates_year)
      ds_Uwind_sig995 = xr.open_dataset('https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.dailyavgs/surface/uwnd.sig995.gauss.'+str(iyr)+'.nc',engine='netcdf4').sel(\
                            time=dates_year)
      ds_Vwind_sig995 = xr.open_dataset('https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.dailyavgs/surface/vwnd.sig995.gauss.'+str(iyr)+'.nc',engine='netcdf4').sel(\
                            time=dates_year)
      ds_SkinTemp_Sfc = xr.open_dataset('https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.dailyavgs/surface_gauss/skt.sfc.gauss.'+str(iyr)+'.nc',engine='netcdf4').sel(\
                            time=dates_year)
      ds_PrecipWater = xr.open_dataset('https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.dailyavgs/surface/pr_wtr.eatm.'+str(iyr)+'.nc',engine='netcdf4').sel(\
                            time=dates_year)

      # Append data
      datasets_Uwind_250hPa.append(ds_Uwind_250hPa)
      datasets_Vwind_250hPa.append(ds_Vwind_250hPa)
      datasets_GeopHgt_500hPa.append(ds_GeopHgt_500hPa)
      datasets_Uwind_500hPa.append(ds_Uwind_500hPa)
      datasets_Vwind_500hPa.append(ds_Vwind_500hPa)
      datasets_Omega_500hPa.append(ds_Omega_500hPa)
      datasets_Uwind_850hPa.append(ds_Uwind_850hPa)
      datasets_Vwind_850hPa.append(ds_Vwind_850hPa)
      datasets_SpecHum_850hPa.append(ds_SpecHum_850hPa)
      datasets_AirTemp_850hPa.append(ds_AirTemp_850hPa)
      datasets_Uwind_sig995.append(ds_Uwind_sig995)
      datasets_Vwind_sig995.append(ds_Vwind_sig995)
      datasets_SkinTemp_Sfc.append(ds_SkinTemp_Sfc)
      datasets_PrecipWater.append(ds_PrecipWater)


# Combine DJF extreme precipitation days from each yearly file into 1 file
combined_Uwind_250hPa = xr.concat(datasets_Uwind_250hPa, dim='time')
combined_Vwind_250hPa = xr.concat(datasets_Vwind_250hPa, dim='time')
combined_GeopHgt_500hPa = xr.concat(datasets_GeopHgt_500hPa, dim='time')
combined_Uwind_500hPa = xr.concat(datasets_Uwind_500hPa, dim='time')
combined_Vwind_500hPa = xr.concat(datasets_Vwind_500hPa, dim='time')
combined_Omega_500hPa = xr.concat(datasets_Omega_500hPa, dim='time')
combined_Uwind_850hPa = xr.concat(datasets_Uwind_850hPa, dim='time')
combined_Vwind_850hPa = xr.concat(datasets_Vwind_850hPa, dim='time')
combined_SpecHum_850hPa = xr.concat(datasets_SpecHum_850hPa, dim='time')
combined_AirTemp_850hPa = xr.concat(datasets_AirTemp_850hPa, dim='time')
combined_Uwind_sig995 = xr.concat(datasets_Uwind_sig995, dim='time')
combined_Vwind_sig995 = xr.concat(datasets_Vwind_sig995, dim='time')
combined_SkinTemp_Sfc = xr.concat(datasets_SkinTemp_Sfc, dim='time')
combined_PrecipWater = xr.concat(datasets_PrecipWater, dim='time')

# Convert Dataset files into Netcdf 
combined_Uwind_250hPa.to_netcdf('combined_Uwind_250hPa_1996to2019_ExtremePrecipDays.nc')
combined_Vwind_250hPa.to_netcdf('combined_Vwind_250hPa_1996to2019_ExtremePrecipDays.nc')
combined_GeopHgt_500hPa.to_netcdf('combined_GeopHgt_500hPa_1996to2019_ExtremePrecipDays.nc')
combined_Uwind_500hPa.to_netcdf('combined_Uwind_500hPa_1996to2019_ExtremePrecipDays.nc')
combined_Vwind_500hPa.to_netcdf('combined_Vwind_500hPa_1996to2019_ExtremePrecipDays.nc')
combined_Omega_500hPa.to_netcdf('combined_Omega_500hPa_1996to2019_ExtremePrecipDays.nc')
combined_Uwind_850hPa.to_netcdf('combined_Uwind_850hPa_1996to2019_ExtremePrecipDays.nc')
combined_Vwind_850hPa.to_netcdf('combined_Vwind_850hPa_1996to2019_ExtremePrecipDays.nc')
combined_SpecHum_850hPa.to_netcdf('combined_SpecHum_850hPa_1996to2019_ExtremePrecipDays.nc')
combined_AirTemp_850hPa.to_netcdf('combined_AirTemp_850hPa_1996to2019_ExtremePrecipDays.nc')
combined_Uwind_sig995.to_netcdf('combined_Uwind_sig995_1996to2019_ExtremePrecipDays.nc')
combined_Vwind_sig995.to_netcdf('combined_Vwind_sig995_1996to2019_ExtremePrecipDays.nc')
combined_SkinTemp_Sfc.to_netcdf('combined_SkinTemp_Sfc_1996to2019_ExtremePrecipDays.nc')
combined_PrecipWater.to_netcdf('combined_PrecipWater_1996to2019_ExtremePrecipDays.nc')

# Move newly-created Netcdf files into your Google Drive 
!mv *.nc "/content/drive/My Drive/"

##############################################################################################################################

## Importing Long-term Mean NCEP Reanalysis data from 1981-2010, only selecting DJF

months = xr.cftime_range(start='0001-01-01', end='0001-12-01', freq='MS', calendar = 'standard') # selecting long-term mean data for 1981-2010 with cftime.DatetimeGregorian format
months = months[(months.month==12)|(months.month==1)|(months.month==2)] # selecting long-term mean data for DJF in 1981-2010

# Extract data
ds_Uwind_250hPa_LTM = xr.open_dataset('https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.derived/pressure/uwnd.mon.1981-2010.ltm.nc',engine='netcdf4').sel(\
                      level=250, time = months)
ds_Vwind_250hPa_LTM = xr.open_dataset('https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.derived/pressure/vwnd.mon.1981-2010.ltm.nc',engine='netcdf4').sel(\
                      level=250, time = months)
ds_GeopHgt_500hPa_LTM = xr.open_dataset('https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.derived/pressure/hgt.mon.1981-2010.ltm.nc',engine='netcdf4').sel(\
                      level=500, time = months)
ds_Uwind_500hPa_LTM = xr.open_dataset('https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.derived/pressure/uwnd.mon.1981-2010.ltm.nc',engine='netcdf4').sel(\
                      level=500, time = months)
ds_Vwind_500hPa_LTM = xr.open_dataset('https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.derived/pressure/vwnd.mon.1981-2010.ltm.nc',engine='netcdf4').sel(\
                      level=500, time = months)
ds_Omega_500hPa_LTM = xr.open_dataset('https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.derived/pressure/omega.mon.1981-2010.ltm.nc',engine='netcdf4').sel(\
                      level=500, time = months)
ds_Uwind_850hPa_LTM = xr.open_dataset('https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.derived/pressure/uwnd.mon.1981-2010.ltm.nc',engine='netcdf4').sel(\
                      level=850, time = months)
ds_Vwind_850hPa_LTM = xr.open_dataset('https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.derived/pressure/vwnd.mon.1981-2010.ltm.nc',engine='netcdf4').sel(\
                      level=850, time = months)
ds_SpecHum_850hPa_LTM = xr.open_dataset('https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.derived/pressure/shum.mon.1981-2010.ltm.nc',engine='netcdf4').sel(\
                      level=850, time = months)
ds_AirTemp_850hPa_LTM = xr.open_dataset('https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.derived/pressure/air.mon.1981-2010.ltm.nc',engine='netcdf4').sel(\
                      level=850, time = months)
ds_Uwind_sig995_LTM = xr.open_dataset('https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.derived/surface/uwnd.sig995.mon.1981-2010.ltm.nc',engine='netcdf4').sel(\
                      time = months)
ds_Vwind_sig995_LTM = xr.open_dataset('https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.derived/surface/vwnd.sig995.mon.1981-2010.ltm.nc',engine='netcdf4').sel(\
                      time = months)
ds_SkinTemp_Sfc_LTM = xr.open_dataset('https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.derived/surface_gauss/skt.sfc.mon.1981-2010.ltm.nc',engine='netcdf4').sel(\
                      time = months)
ds_PrecipWater_LTM = xr.open_dataset('https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.derived/surface/pr_wtr.eatm.mon.1981-2010.ltm.nc',engine='netcdf4').sel(\
                      time = months)
