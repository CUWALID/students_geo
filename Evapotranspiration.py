"""
FAO Penman-Monteith equation

"""


import pandas as pd
import numpy as np
from pandas import read_csv
from pandas import datetime
from matplotlib import pyplot
import datetime
from datetime import timedelta
import math
	
def psy_const(atmos_pres):
    """
    Calculate the psychrometric constant.
    This method assumes that the air is saturated with water vapour at the
    minimum daily temperature. This assumption may not hold in arid areas.
    Based on equation 8, page 95 in Allen et al (1998).
    :param atmos_pres: Atmospheric pressure [kPa]. Can be estimated using
        ``atm_pressure()``.
    :return: Psychrometric constant [kPa degC-1].
    :rtype: float
    """
	return(0.000665*atmos_pres)
	
def Ra_hourly(lat,lon,J,t):
	# lat: latitud
	# lon: longitud
	# J: Julian day
	# t: hour
	# Lz: longitude of the centre of the local time zone
	# Sc: seasonal correction
	# J julian day
	b=2*np.pi*(J-81)/364 # Eq. 33
	Sc = 0.1645*np.sin(2*b) - 0.1255*np.cos(b) - 0.025*np.sin(b) # Eq. 32
	
	# solar time angle at midpoint
	# Lz longitude of the centre of the local time zone
	Lz=15*np.around(lon/15)#+7.5
	Lm = lon
	# Lm: longitude of the measurement site
	# w: Solar time angle at the midpoint of period
	w=(np.pi/12)*((t+0.06667*(Lz-Lm)+Sc)-12)
	
	# solar time angles at the beginning and end of the period
	t_1=1. # for one hour period
	# w 1 solar time angle at beginning of period [rad] (Equation 29),
	# w 2 solar time angle at end of period [rad] (Equation 30).
	w_1=w-np.pi*t_1/24
	w_2=w+np.pi*t_1/24
	
	# Gsc solar constant = 0.0820 MJ m-2 min-1,
	Gsc = 0.0820 # [MJ m-2 min-1]
	
	# dr inverse relative distance Earth-Sun (Equation 23),
	dr = 1+ 0.033*np.cos(2*np.pi*J/365)
	
	# d solar declination [rad] (Equation 24)
	d=0.409*np.sin(2*np.pi*J/365-1.39)
	
	# j latitude [rad] (Equation 22)
	j=(np.pi/180)*lat
	
	# Sunset hour angle (Equation 25)
	ws=np.arccos(-np.tan(d)*np.tan(j))
	
	aux_1=np.ones(len(J))
	
	aux_1[np.where(w>ws)]=0.0
	aux_1[np.where(w<-ws)]=0.0
	# Ra extraterrestrial radiation in the hour (or shorter) period [MJ m-2 hour-1],
	
	Ra=aux_1*(12*60/np.pi)*Gsc*dr*((w_2-w_1)*np.sin(d)*np.sin(j)+np.cos(d)*np.cos(j)*(np.sin(w_2)-np.sin(w_1)))
	
	return Ra,ws,w

def Ra_daily(lat,lon,J):
	# lat: latitud
	# lon: longitud
	# J: Julian day
	# seasonal correction

	# Gsc solar constant = 0.0820 MJ m-2 min-1,
	
	Gsc = 0.0820 # [MJ m-2 min-1]
	
	# dr inverse relative distance Earth-Sun (Equation 23),
	
	dr = 1+ 0.033*np.cos(2*np.pi*J/365)
	
	# d solar declination [rad] (Equation 24)
	
	d=0.409*np.sin(2*np.pi*J/365-1.39)
	
	# j latitude [rad] (Equation 22)
	
	j=(np.pi/180)*lat
	
	# Sunset hour angle (Equation 25)

	ws=np.arccos(-np.tan(d)*np.tan(j))

	# Ra extraterrestrial radiation in the hour (or shorter) period [MJ m-2 hour-1], (Equation 28)
	
	Ra=(24*60/np.pi)*Gsc*dr*(ws*np.sin(d)*np.sin(j)+np.cos(d)*np.cos(j)*np.sin(ws))
	
	return Ra

def avp_from_rh(svp, rh):
    """
    Estimate actual vapour pressure (*e*a) from saturation vapour pressure at
    daily minimum temperature and maximum relative humidity
    Based on FAO equation 18 in Allen et al (1998).
    :param svp_tmin: Saturation vapour pressure at daily minimum temperature
        [kPa]. Can be estimated using ``svp_from_t()``.
    :param rh_max: Maximum relative humidity [%]
    :return: Actual vapour pressure [kPa]
    :rtype: float
    """
	return svp * (rh / 100.0)
	
def avp_from_Tdew(tdew):
    """
    Estimate actual vapour pressure (*e*a) from saturation vapour pressure at
    daily minimum temperature and maximum relative humidity
    Based on FAO equation 18 in Allen et al (1998).
    :param svp_tmin: Saturation vapour pressure at daily minimum temperature
        [kPa]. Can be estimated using ``svp_from_t()``.
    :param rh_max: Maximum relative humidity [%]
    :return: Actual vapour pressure [kPa]
    :rtype: float
    """
	return 0.6108 * np.exp((17.27 * tdew) / (tdew + 237.3))

def svp_from_t(t):
    """
    Estimate saturation vapour pressure (*es*) from air temperature.
    Based on equations 11 and 12 in Allen et al (1998).
    :param t: Temperature [deg C]
    :return: Saturation vapour pressure [kPa]
    :rtype: float
    """
    return 0.6108 * np.exp((17.27 * t) / (t + 237.3))


def delta_svp(t):
    """
    Estimate the slope of the saturation vapour pressure curve at a given
    temperature.
    Based on equation 13 in Allen et al (1998). If using in the Penman-Monteith
    *t* should be the mean air temperature.
    :param t: Air temperature [deg C]. Use mean air temperature for use in
        Penman-Monteith.
    :return: Saturation vapour pressure [kPa degC-1]
    :rtype: float
    """
    tmp = 4098 * (0.6108 * np.exp((17.27 * t) / (t + 237.3)))
	
    return tmp / ((t + 237.3)**2)

def Rsn_from_Ra(Rs,a):
	if a is None:
		a=0.23
	return (1-a)*Rs
	
def P_T_PETp_daily(net_rad, G, psy, delta_svp):
    """
    Estimate reference evapotranspiration (ETo) from a hypothetical
    short grass reference surface using the FAO-56 Penman-Monteith equation.
    Based on equation 6 in Allen et al (1998).
    :param net_rad: Net radiation at crop surface [MJ m-2 day-1]. If
        necessary this can be estimated using ``net_rad()``.
    :param t: Air temperature at 2 m height [deg Kelvin].
    :param ws: Wind speed at 2 m height [m s-1]. If not measured at 2m,
        convert using ``wind_speed_at_2m()``.
    :param svp: Saturation vapour pressure [kPa]. Can be estimated using
        ``svp_from_t()''.
    :param avp: Actual vapour pressure [kPa]. Can be estimated using a range
        of functions with names beginning with 'avp_from'.
    :param delta_svp: Slope of saturation vapour pressure curve [kPa degC-1].
        Can be estimated using ``delta_svp()``.
    :param psy: Psychrometric constant [kPa deg C]. Can be estimatred using
        ``psy_const_of_psychrometer()`` or ``psy_const()``.
    :param shf: Soil heat flux (G) [MJ m-2 day-1] (default is 0.0, which is
        reasonable for a daily or 10-day time steps). For monthly time steps
        *shf* can be estimated using ``monthly_soil_heat_flux()`` or
        ``monthly_soil_heat_flux2()``.

	gamma: latent heat of vaporization
    :return: Reference evapotranspiration (ETo) from a hypothetical
        grass reference surface [mm day-1].
    :rtype: float
    """
	
	G=0.0 # For daily values
	
	alpha=1.74 # For arid regions
	
	gamma=2.45 # MJ kg**-1

	PET=(net_rad-G)*(1/gamma)*delta_svp/(delta_svp+psy)
	
	return alpha*PET
	

def P_T_PETp_hourly(net_rad, G, psy, delta_svp,alpha):
    """
	Priestley Taylor Method
	
    :param net_rad: Net radiation at crop surface [MJ m-2 day-1]. If
        necessary this can be estimated using ``net_rad()``.
    :param t: Air temperature at 2 m height [deg Kelvin].
    :param ws: Wind speed at 2 m height [m s-1]. If not measured at 2m,
        convert using ``wind_speed_at_2m()``.
    :param svp: Saturation vapour pressure [kPa]. Can be estimated using
        ``svp_from_t()''.
    :param avp: Actual vapour pressure [kPa]. Can be estimated using a range
        of functions with names beginning with 'avp_from'.
    :param delta_svp: Slope of saturation vapour pressure curve [kPa degC-1].
        Can be estimated using ``delta_svp()``.
    :param psy: Psychrometric constant [kPa deg C]. Can be estimatred using
        ``psy_const_of_psychrometer()`` or ``psy_const()``.
    :param shf: Soil heat flux (G) [MJ m-2 day-1] (default is 0.0, which is
        reasonable for a daily or 10-day time steps). For monthly time steps
        *shf* can be estimated using ``monthly_soil_heat_flux()`` or
        ``monthly_soil_heat_flux2()``.

	gamma: latent heat of vaporization
    :return: Reference evapotranspiration (ETo) from a hypothetical
        grass reference surface [mm day-1].
    :rtype: float
    """
	if alpha is None:
		alpha=1.74
		
	#alpha=1.74 # For arid regions
	
	gamma=2.45 # MJ kg**-1

	PET=(net_rad-G)*(1/gamma)*delta_svp/(delta_svp+psy)
	
	return alpha*PET


	
def PM_FAO_hourly(net_rad, t, ws, svp, avp, delta_svp, psy, shf):
    """
	Penman-Monteith (PM) from FAO-56 at hourly time steps

    Estimate reference evapotranspiration (ETo) from a hypothetical
    short grass reference surface using the FAO-56 Penman-Monteith equation.
    Based on equation 6 in Allen et al (1998).
    :param net_rad: Net radiation at crop surface [MJ m-2 day-1]. If
        necessary this can be estimated using ``net_rad()``.
    :param t: Air temperature at 2 m height [deg Kelvin].
    :param ws: Wind speed at 2 m height [m s-1]. If not measured at 2m,
        convert using ``wind_speed_at_2m()``.
    :param svp: Saturation vapour pressure [kPa]. Can be estimated using
        ``svp_from_t()''.
    :param avp: Actual vapour pressure [kPa]. Can be estimated using a range
        of functions with names beginning with 'avp_from'.
    :param delta_svp: Slope of saturation vapour pressure curve [kPa degC-1].
        Can be estimated using ``delta_svp()``.
    :param psy: Psychrometric constant [kPa deg C]. Can be estimatred using
        ``psy_const_of_psychrometer()`` or ``psy_const()``.
    :param shf: Soil heat flux (G) [MJ m-2 day-1] (default is 0.0, which is
        reasonable for a daily or 10-day time steps). For monthly time steps
        *shf* can be estimated using ``monthly_soil_heat_flux()`` or
        ``monthly_soil_heat_flux2()``.
    :return: Reference evapotranspiration (ETo) from a hypothetical
        grass reference surface [mm day-1].
    :rtype: float
    """

	b=(delta_svp + (psy * (1 + 0.34 * ws)))
	
	a1 = (0.408 * (net_rad - shf) * delta_svp )/b

	a2 = (37 * ws * (svp - avp) * psy/ (t+273))/b
	
	return a1 + a2	


def PM_FAO_daily(net_rad, t, ws, svp, avp, delta_svp, psy, shf):
    """
	Penman-Monteith (PM) from FAO-56 at hourly time steps

    Estimate reference evapotranspiration (ETo) from a hypothetical
    short grass reference surface using the FAO-56 Penman-Monteith equation.
    Based on equation 6 in Allen et al (1998).
    :param net_rad: Net radiation at crop surface [MJ m-2 day-1]. If
        necessary this can be estimated using ``net_rad()``.
    :param t: Air temperature at 2 m height [deg Kelvin].
    :param ws: Wind speed at 2 m height [m s-1]. If not measured at 2m,
        convert using ``wind_speed_at_2m()``.
    :param svp: Saturation vapour pressure [kPa]. Can be estimated using
        ``svp_from_t()''.
    :param avp: Actual vapour pressure [kPa]. Can be estimated using a range
        of functions with names beginning with 'avp_from'.
    :param delta_svp: Slope of saturation vapour pressure curve [kPa degC-1].
        Can be estimated using ``delta_svp()``.
    :param psy: Psychrometric constant [kPa deg C]. Can be estimatred using
        ``psy_const_of_psychrometer()`` or ``psy_const()``.
    :param shf: Soil heat flux (G) [MJ m-2 day-1] (default is 0.0, which is
        reasonable for a daily or 10-day time steps). For monthly time steps
        *shf* can be estimated using ``monthly_soil_heat_flux()`` or
        ``monthly_soil_heat_flux2()``.
    :return: Reference evapotranspiration (ETo) from a hypothetical
        grass reference surface [mm day-1].
    :rtype: float
    """

	b=(delta_svp + (psy * (1 + 0.34 * ws)))
	
	a1 = (0.408 * (net_rad - shf) * delta_svp )/b

	a2 = (900 * ws * (svp - avp) * psy/ (t+273))/b
	
	return a1 + a2	


def wind_2m_from_wz(v_z,z):
	"""
	w_z: wind at z meters from ground surface
	z : elevation from ground surface
	"""
	u_2m=v_z*4.87 / np.log(67.8*z-5.42)
	
	return u_2m