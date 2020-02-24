## Import NCEP Reanalysis fields directly from the NOAA THREDDS server, selecting only the extreme precipitation days

#days = pd.date_range(start='1996-01-01', end='2019-12-31', freq='D') # selecting daily data from 1996-2019
#days = days[(days.month==12)|(days.month==1)|(days.month==2)] # selecting only DJF from 1996-2019

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
datasets_Uwind_Sfc = []
datasets_Vwind_Sfc = []
datasets_SkinTemp_Sfc = []
datasets_PrecipWater = []

for iyr in years:
    print('working on '+str(iyr))
    dates_year = ExtremePrecip_dates[ExtremePrecip_dates['time.year'].values == iyr]

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
    ds_Uwind_Sfc = xr.open_dataset('https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.dailyavgs/surface_gauss/uwnd.10m.gauss.'+str(iyr)+'.nc',engine='netcdf4').sel(\
                          time=dates_year)
    ds_Vwind_Sfc = xr.open_dataset('https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.dailyavgs/surface_gauss/vwnd.10m.gauss.'+str(iyr)+'.nc',engine='netcdf4').sel(\
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
    datasets_Uwind_Sfc.append(ds_Uwind_Sfc)
    datasets_Vwind_Sfc.append(ds_Vwind_Sfc)
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
combined_Uwind_Sfc = xr.concat(datasets_Uwind_Sfc, dim='time')
combined_Vwind_Sfc = xr.concat(datasets_Vwind_Sfc, dim='time')
combined_SkinTemp_Sfc = xr.concat(datasets_SkinTemp_Sfc, dim='time')
combined_PrecipWater = xr.concat(datasets_PrecipWater, dim='time')
