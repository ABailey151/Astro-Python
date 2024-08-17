from sys import path
path.append('./Tools')
import numpy as np
import Planetary_data as pd

#Import functions where needed
#from sv_coe import *
#from geo_func import *
#from Kep_Anom import *
from plot_func import *
from misc_maths import JulianD

def JD_Rm(JD):

    T0 = (JD - 2451545.0)/36525 #Number of Julian centuries since J2000

    ob = np.deg2rad(23.439 - (0.0130042*T0))

    LpCoefs = np.array(([0,     218.32,  481267.881, 0,     0,      0,         0.9508, 0,      0],
                        [6.29,  135.0,   477198.87,  5.13,  93.3,   483202.03, 0.0518, 135.0,  477198.87],
                        [-1.27, 259.3,  -413335.36,  0.28,  220.2,  960400.89, 0.0095, 259.3, -413335.38],
                        [0.66,  235.7,   890534.22, -0.28,  318.3,  6003.15,   0.0078, 253.7,  890534.22],
                        [0.21,  269.9,   954397.74, -0.17,  217.6, -407332.21, 0.0028, 269.9,  954397.70],
                        [-0.19, 357.5,   35999.05,  0,      0,     0,          0,      0,      0],
                        [-0.11, 106.5,   966404.03, 0,      0,     0,          0,      0,      0]))

    i = 0
    for m in range(0,7):
        i = i + (LpCoefs[m,0] * np.sin(np.deg2rad(LpCoefs[m,1] + (LpCoefs[m,2]*T0))))
    Lelon = np.deg2rad((LpCoefs[0,1] + (LpCoefs[0,2]*T0) + i)%360)
    
    j = 0
    for b in range(1,5):
        j = j + (LpCoefs[b,3] * np.sin(np.deg2rad(LpCoefs[b,4] + (LpCoefs[b,5]*T0))))
    Lelat = np.deg2rad(j%360)

    k = 0
    for f in range(1,5):
        k = k + (LpCoefs[f,6] * np.sin(np.deg2rad(LpCoefs[f,7] + (LpCoefs[f,8]*T0))))
    HP = np.deg2rad(LpCoefs[0,6] + k%360)
    
    rm = 6378.14/np.sin(HP)

    Uxyz = np.array([np.cos(Lelat)*np.cos(Lelon),
                    (np.cos(ob)*np.cos(Lelat)*np.sin(Lelon)) - (np.sin(ob)*np.sin(Lelat)),
                    (np.sin(ob)*np.cos(Lelat)*np.sin(Lelon)) + (np.cos(ob)*np.sin(Lelat))])
    
    Rm = rm * Uxyz
    return Rm

epoch = np.array([25,7,2013,8,0,0])
J0,UT = JulianD(epoch)
JD0 = J0 + (UT/24)

ttotal = 27.322*24*3600
tstep = 600
N = int(np.ceil(ttotal/tstep))
n = 0

Rms = np.zeros([N,3])
while n <= N-1:
    Rm = JD_Rm(JD0)
    Rms[n] = Rm
    JD0 = JD0 + (tstep/86400)
    n = n + 1

##Plots
fig = plt.figure(figsize=(10,8))
axes = plt.subplot(projection='3d')
#p1 = Traj(axes,np.linalg.norm(Rms),Rms,ttotal,N)
p1 = Orbit_Plot(axes,6378.14,Rms,ttotal,N,'0.4','b')
plt.subplots_adjust(left=0.05,right=0.95,bottom=0.09,top=0.95)
plt.show()