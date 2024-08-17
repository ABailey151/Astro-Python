from sys import path
path.append('./Tools')
import numpy as np
import Planetary_data as pd
from misc_maths import JulianD

GD = np.array([3,3,2004,4,30,0])
lon = np.deg2rad(139.8)
print('\nGreg Date(UT) =',GD,', Eastward Longitude from Greenwich =',np.rad2deg(lon),'deg')

#Step 1: Calculate to Julian day at 0h UT
J0,UT = JulianD(GD)
#Step 2: Calculate the number opf Julian Centuries since J2000
T0 = (J0 - 2451545)/36525
print('\nJ0 =',J0,'days, t0 =',T0)

#Step 3: Calculate the Greenwich Sidereal Time at 0h UT and correct to 0 and 360 degree range
st_G0 = np.deg2rad(100.4606184 + 36000.77004*T0 + 0.000387933*(T0**2) - (2.583e-8)*(T0**3))
if st_G0 > (2*np.pi):
    st_G0 = st_G0 - int(st_G0/(2*np.pi))*(2*np.pi)
if st_G0 < 0:
    st_G0 = st_G0 - np.floor(st_G0/(2*np.pi))*(2*np.pi)
print('\nGreenwich Sidereal Time at 0h UT =',np.rad2deg(st_G0),'deg')

#Step 4: find the Greenwich Sidereal Time at any other UTC
st_G = st_G0 + np.deg2rad(360.98564724)*(UT/24)
print('Greenwich Sidereal Time at UT',UT,'=',np.rad2deg(st_G),'deg')

#Step 5: Find the Local Sidereal Time and correct to 0 and 360 degree range
st_l = st_G + lon
if st_l > (2*np.pi):
    st_l = st_l - int(st_l/(2*np.pi))*(2*np.pi)
print('\nLocal Sidereal Time =',np.rad2deg(st_l),'deg, or',st_l*(24/(2*np.pi)),'h\n')
