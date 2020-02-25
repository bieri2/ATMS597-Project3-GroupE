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
datasets_Uwind_10m = []
datasets_Vwind_10m = []
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
      ds_Uwind_10m = xr.open_dataset('https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.dailyavgs/surface_gauss/uwnd.10m.gauss.'+str(iyr)+'.nc',engine='netcdf4').sel(\
                            time=dates_year)
      ds_Vwind_10m = xr.open_dataset('https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.dailyavgs/surface_gauss/vwnd.10m.gauss.'+str(iyr)+'.nc',engine='netcdf4').sel(\
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
      datasets_Uwind_10m.append(ds_Uwind_10m)
      datasets_Vwind_10m.append(ds_Vwind_10m)
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
combined_Uwind_10m = xr.concat(datasets_Uwind_10m, dim='time')
combined_Vwind_10m = xr.concat(datasets_Vwind_10m, dim='time')
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
combined_Uwind_10m.to_netcdf('combined_Uwind_10m_1996to2019_ExtremePrecipDays.nc')
combined_Vwind_10m.to_netcdf('combined_Vwind_10m_1996to2019_ExtremePrecipDays.nc')
combined_SkinTemp_Sfc.to_netcdf('combined_SkinTemp_Sfc_1996to2019_ExtremePrecipDays.nc')
combined_PrecipWater.to_netcdf('combined_PrecipWater_1996to2019_ExtremePrecipDays.nc')

# Move newly-created Netcdf files into your Google Drive 
!mv *.nc "/content/drive/My Drive/"
