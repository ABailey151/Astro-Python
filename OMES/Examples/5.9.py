from sys import path
path.append('./Tools')
import numpy as np
import Planetary_data as pd

#Import functions where needed
#from sv_coe import *
#from geo_func import *
#from Kep_Anom import *
#from plot_func import *
#from misc_maths import *

cb = pd.earth
re = cb['radius']
f = cb['flat']

#Given Geocentric Position Vector
r = np.array([-2032.4, 4591.2, -4544.8])

#Fine the azimuth and elevation angle in the Topocentric horizon frame at sea level
H = 0
lat = np.deg2rad(-40)
st_l = np.deg2rad(110)

#Calculate Geocentric Position of the observer
const_1 = np.sqrt(1 - (2*f - (f**2))*(np.sin(lat)**2))
const_ij = ((re/const_1) + H)*np.cos(lat)
const_k = ((re*((1-f)**2))/const_1) + H
R = np.array([const_ij*np.cos(st_l), const_ij*np.sin(st_l), const_k*np.sin(lat)])

#Therefore, find the spacecraft position relative to the observer in the geocentric frame
pX = r - R

#Define the transformation matrix to go from the geocentric frame to the topocentric horizon frame
QXx = np.array(([-np.sin(st_l),              np.cos(st_l),              0],
                [-np.sin(lat)*np.cos(st_l), -np.sin(lat)*np.sin(st_l),  np.cos(lat)],
                [np.cos(lat)*np.cos(st_l),  np.cos(lat)*np.sin(st_l),   np.sin(lat)]))

#Use Matrix to transform between frames
pxx = pX[0]*QXx[0,0] + pX[1]*QXx[0,1] + pX[2]*QXx[0,2]
pxy = pX[0]*QXx[1,0] + pX[1]*QXx[1,1] + pX[2]*QXx[1,2]
pxz = pX[0]*QXx[2,0] + pX[1]*QXx[2,1] + pX[2]*QXx[2,2]
px = np.array([pxx,pxy,pxz])
unit_px = px/np.linalg.norm(px)
print('\nTopocentric Horizon Position vector =',px,'km')
print('Unit Vector =',unit_px)

#from eq.(5.58), the elevation angle is given by the k component of the position vector in the topocentric horizon frame
elev_a = np.arcsin(unit_px[2])
print('\nElevation Angle =',np.rad2deg(elev_a),'deg')
#The azimuth is found from the i and j components of the same equation
sinA = unit_px[0]/np.cos(elev_a)
cosA = unit_px[1]/np.cos(elev_a)
azi = np.arccos(cosA)
print('\nsinA =',sinA,'cosA =',cosA)
print('Azimuth =',np.rad2deg(azi),'or',360-np.rad2deg(azi),'deg\n')