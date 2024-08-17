from sys import path
path.append('./Tools')
import numpy as np
import Planetary_data as pd

#Import functions where needed
from sv_coe import coe_from_sv
#from geo_func import *
#from Kep_Anom import *
#from plot_func import *
#from misc_maths import *

cb = pd.earth
mu = cb['mu']
re = cb['radius']
f = cb['flat']

#Observer Information
st_l = np.deg2rad(300)
H = 0
lat = np.deg2rad(60)
#//

#Observed spacecraft information
pmag = 2551 #range
Azi = np.deg2rad(90)
elev_a = np.deg2rad(30)
pmag_rate = 0   #Rads/s
Azi_rate = 1.973e-3 #Rads/s
elev_a_rate = 9.864e-4  #Rads/s
#//

#Step 1: Calculate Geocentric Position of the observer, using H, lat, and st_l
const_1 = np.sqrt(1 - (2*f - (f**2))*(np.sin(lat)**2))
const_ij = ((re/const_1) + H)*np.cos(lat)
const_k = ((re*((1-f)**2))/const_1) + H
R = np.array([const_ij*np.cos(st_l), const_ij*np.sin(st_l), const_k*np.sin(lat)])
#//

#Step 2: Calculate the topocentric declination
sin_topo_dec = np.cos(lat)*np.cos(Azi)*np.cos(elev_a) + np.sin(lat)*np.sin(elev_a)
topo_dec = np.arcsin(sin_topo_dec)
#//

#Step 3: Calculate the topocentric right ascension
cos_h_ad = (np.cos(lat)*np.sin(elev_a) - np.sin(lat)*np.cos(Azi)*np.cos(elev_a))/np.cos(topo_dec)
if 0 < Azi < np.pi:
    h_ad = 2*np.pi - np.arccos(cos_h_ad)
elif np.pi <= Azi <= 2*np.pi:
    h_ad = np.arccos(cos_h_ad)
topo_ra = st_l - h_ad
#//

#Step 4: Caculate the range direction cosine unit vector
p_hat = np.array([np.cos(topo_ra)*np.cos(topo_dec),
                  np.sin(topo_ra)*np.cos(topo_dec),
                  np.sin(topo_dec)])
#//

#Step 5: Calculate the geocentric position vector of the spacecraft
r = R + pmag*p_hat
#//

#Step 6: Calculate the inertial velocity
Omega = np.array([0, 0, cb['dphi']])
R_dot = np.cross(Omega,R)
#//

#Step 7: Calculate the declination rate
decdot_p1 = -Azi_rate*np.cos(lat)*np.sin(Azi)*np.cos(elev_a)
decdot_p2 = np.sin(lat)*np.cos(elev_a) - np.cos(lat)*np.cos(Azi)*np.sin(elev_a)
dec_dot = (1/np.cos(topo_dec))*(decdot_p1 + elev_a_rate*decdot_p2)
#//

#Step 8: Calculate the right ascension rate
radot_p1 = Azi_rate*np.cos(Azi)*np.cos(elev_a) - elev_a_rate*np.sin(Azi)*np.sin(elev_a) + dec_dot*np.sin(Azi)*np.cos(elev_a)*np.tan(topo_dec)
radot_p2 = np.cos(lat)*np.sin(elev_a) - np.sin(lat)*np.cos(Azi)*np.cos(elev_a)
ra_dot = cb['dphi'] + (radot_p1/radot_p2)
#//

#Step 9: Caculate the range rate direction cosine unit vector
p_hat_dot = np.array([-ra_dot*np.sin(topo_ra)*np.cos(topo_dec) - dec_dot*np.cos(topo_ra)*np.sin(topo_dec),
                      ra_dot*np.cos(topo_ra)*np.cos(topo_dec) - dec_dot*np.sin(topo_ra)*np.sin(topo_dec),
                      dec_dot*np.cos(topo_dec)])
#//

#Step 10: Calculate the geocentric velocity
v = R_dot + pmag*p_hat_dot + pmag_rate*p_hat
print(r,'\n',v)

#//Determine the orbital elements
a,h,hmag,i,lon_an,e,emag,pearg,ta1 = coe_from_sv(mu,r,v)
print('\nOrbital Elements\nh (km^2/s) =',hmag,'\na (km) =',a,'\ni (Deg) =',np.rad2deg(i),'\nRAAN (Deg) =',np.rad2deg(lon_an),'\ne =',emag,'\nArg of Pe (Deg) =',np.rad2deg(pearg),'\nTa (Deg) =',np.rad2deg(ta1),'\n')
#//