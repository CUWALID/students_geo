import numpy as np
import datetime
from datetime import timedelta
import Evapotranspiration as fpet
from astral import *
from astral.sun import sun



def csv_ETp(fname_in, lat, lon, z=0.0):
	"""
	lat: latitude [degree]
	lon: longitude (degree)
	z: elevation (m)
	"""
	
	a = LocationInfo()
	a.name = 'WalnutGulch'
	a.region = 'US'
	a.latitude = lat
	a.longitude = lon
	a.timezone = 'US/Central'
	a.elevation = z
			
	
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
	

def get_pet_PM(data, daylight_t=0.5):
	"""
	data: csv file containing all ERA-land data
	"""
	
	
	Tdew = -273.15+data['d2m'].values # Dewpoint temperature in C
	
	Ta = -273.15+data['t2m'].values # Temperature in C
	
	ws_u = np.abs(data['u10'].values)) # Wind speed in m/s in x-direction
	
	ws_v = np.abs(data['v10'].values # Wind speed in m/s in y-direction
	
	P = 0.001*data['sp'].values # Pressure in kPa
	
	Rns = (1e-6)*data['ssr'].values # Net solar radiation in MJ m**-2 n*-1
	
	Rnl = (1e-6)*data['str'].values # Net thermal radiation in MJ m**-2 n*-1
	
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
	
	shf = daylight_t*Rn
	
	ETo = fpet.PM_FAO_hourly(Rn, Ta, ws, svp, avp, delta_svp, psy, shf)
		
	return Eto