import numpy as np
import datetime
from datetime import timedelta
import Evapotranspiration as fpet
from astral import *
from astral.sun import sun



def Grid_ETp(fname_in, lat, lon, z=0.0):
	"""
	lat: latitude
	lon: longitude
	"""
	data = Dataset(fname_in,'r')
	
	z = 1368
	
	a = LocationInfo()
	a.name = 'WalnutGulch'
	a.region = 'US'
	a.latitude = lat
	a.longitude = lon
	a.timezone = 'US/Central'
	a.elevation = z
			
	time_in = num2date(data['time'][:], units = data['time'].units, calendar = data['time'].calendar)
	
	time_in = time_in+pd.to_timedelta(Local_hour, 'H')
	
	daylight_t = np.ones(len(time_in))*0.5
		
	endi = -1#time_in.tail(1).index[0]
	
	startdate = datetime.datetime(time_in[0].year,time_in[0].month,time_in[0].day,time_in[0].hour,0)
	
	enddate = datetime.datetime(time_in[endi].year,time_in[endi].month,time_in[endi].day,time_in[endi].hour,0)
	
	idate = startdate
	
	while idate < enddate:
	
		sunrise = sun(a,idate)['sunrise'].replace(tzinfo = None)
		
		sunset = sun(a,idate)['sunset'].replace(tzinfo = None)
	
		iperioday = np.where((time_in <= sunset) & (time_in > sunrise))[0]
	
		daylight_t[iperioday] = 0.1
		
		idate = idate + pd.to_timedelta(1, 'd')
	
	# Create new netcdf
	
	dataset = Dataset(fname_out, 'w', format = 'NETCDF4_CLASSIC')
	
	dataset.createDimension('time', None)
	
	dataset.createDimension('lat', len(data['latitude']))
	
	dataset.createDimension('lon', len(data['longitude']))
	
	#dataset.createDimension('grid', len(data['latitude'])*len(data['longitude']))
	
	# Create coordinate variables for 4-dimensions
	
	time = dataset.createVariable('time', np.float64, ('time',))
	
	lon = dataset.createVariable('lon', np.float64, ('lon',))
	
	lat = dataset.createVariable('lat', np.float64, ('lat',))
	
	# Create the actual 4-d variable
	
	#pet = dataset.createVariable('pet', np.float32,('time','grid'))
	
	pet = dataset.createVariable('pet', np.float32,('time','lat','lon'))
	
	time.units = 'hours since 1980-01-01 00:00:00'
	
	time.calendar = 'gregorian'
	
	for i in range(len(time_in)):
	
		Tdew = -273.15+np.array(data['d2m'][i,:,:])#.flatten() # Dewpoint temperature in C
		
		Ta = -273.15+np.array(data['t2m'][i,:,:])#.flatten() # Temperature in C
		
		ws_u = np.abs(np.array(data['u10'][i,:,:]))#.flatten() # Wind speed in m/s in x-direction
		
		ws_v = np.abs(np.array(data['v10'][i,:,:]))#.flatten() # Wind speed in m/s in y-direction
		
		P = 0.001*np.array(data['sp'][i,:,:])#.flatten() # Pressure in kPa
		
		Rns = (1e-6)*np.array(data['ssr'][i,:,:])#.flatten() # Net solar radiation in MJ m**-2 n*-1
		
		Rnl = (1e-6)*np.array(data['str'][i,:,:])#.flatten() # Net thermal radiation in MJ m**-2 n*-1
		
		PEv = (1e-3)*np.array(data['pev'][i,:,:])#.flatten() # Potential evaporation in MJ m**-2 n*-1
		
		Ev = (1e-3)*np.array(data['e'][i,:,:])#.flatten() # Evaporation in MJ m**-2 n*-1
		
		ws = np.sqrt(np.power(ws_u,2)+np.power(ws_v,2))
		
		# Calculate net radiation
		Rn = Rns+Rnl
		# calculate wind speed
		ws = fpet.wind_2m_from_wz(ws, 10.)	
		# Calculate Saturation vapour pressure [kPa]
		svp = fpet.svp_from_t(Ta)
		# Calculate Actual vapour pressure [kPa]
		avp = fpet.svp_from_t(Tdew)
		# Calculate Slope of saturation vapour pressure curve [kPa degC-1]
		delta_svp = fpet.delta_svp(Ta)
		# Calculate psy: Psychrometric constant [kPa deg C]
		psy = 0.000665*P
		
		shf = daylight_t[i]*Rn
		
		ETo = fpet.PM_FAO_hourly(Rn, Ta, ws, svp, avp, delta_svp, psy, shf)