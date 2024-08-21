from sys import path
path.append('./Tools')
import Planetary_data as pd
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
from sv_coe import sv_from_coe
from plot_func import *
from sv_coe import sv_to_elems
from misc_maths import time_arr
#from Luna_3b import *

def JulianD(gregdate):
        
    UT = gregdate[3] + (gregdate[4]/60) + (gregdate[5]/3600)

    j0 = (367*gregdate[2]) - int((7*(gregdate[2]+int((gregdate[1]+9)/12)))/4) + int((275*gregdate[1])/9) + gregdate[0] + 1721013.5

    JD = j0 + (UT/24)
    return JD

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


def prop_luna_3b(bods,N,tstep,r0,v0,epoch):
    def f(t,y,JD):

        r = y[0:3]
        v = y[3:]

        am = -(bods[0]['mu']/(np.linalg.norm(r)**3))*r

        Rm = JD_Rm(JD)
        RmS = Rm - r
        pm = bods[1]['mu'] * ((RmS/(np.linalg.norm(RmS)**3)) - (Rm/(np.linalg.norm(Rm)**3)))

        a = am + pm

        dydt = np.concatenate((v,a), axis=0)
        return dydt

    y0 = np.append(r0,v0)

    JD = JulianD(epoch)
    rs = np.zeros([N,3])
    vs = np.zeros([N,3])
    Rms = np.zeros([N,3])
    rs[0] = r0
    vs[0] = v0
    Rms[0] = JD_Rm(JD)
    n = 1

    prop = sp.integrate.ode(f,jac=None).set_integrator('dopri5')
    while n < N:
        
        prop.set_initial_value(y0,0).set_f_params(JD)
        prop.integrate(tstep)

        rs[n] = prop.y[0:3]
        vs[n] = prop.y[3:6]

        Rms[n] = JD_Rm(JD)
        y0 = prop.y
        
        JD = JD + ((tstep)/86400)
        n = n + 1
        #print(n,'/',N)
    
    return rs,vs,Rms

#Orbit Initial conditions
epoch = np.array([25,3,2024,8,0,0])

mb = pd.earth
sb = pd.moon
bods = np.array([mb,sb])

ra = mb['radius'] + 305000
rp = mb['radius'] + 220000
a = (ra+rp)/2
emag = (ra-rp)/(ra+rp)
if emag < 1:
    hmag = np.sqrt((mb['mu']*a)*(1-(emag**2)))
else:
    hmag = np.sqrt((mb['mu']*(-a))*(1-(emag**2)))
i = np.deg2rad(25)
lon_an = np.deg2rad(0)
pearg = np.deg2rad(68)
ta0 = np.deg2rad(30)
T = ((2*np.pi)/np.sqrt(mb['mu']))*(a**(3/2))
print('\nOrbital Elements\nh (km^2/s) =',hmag,'\na (km) =',a,'\ni (Deg) =',np.rad2deg(i),'\nRAAN (Deg) =',np.rad2deg(lon_an),'\ne =',emag,'\nArg of Pe (Deg) =',np.rad2deg(pearg),'\nTa (Deg) =',np.rad2deg(ta0),'\n')

r0,v0 = sv_from_coe(mb['mu'],emag,hmag,i,lon_an,pearg,ta0)

#Set time period and step length/amount
ttotal = 2.5*27.321661*24*3600 #One Lunar orbit
tstep = 1200
N = int(np.ceil(ttotal/tstep) + 1)
print('\n'+str(N)+' Steps at '+str(tstep)+' Seconds\n',2.5*27.321661)

rs,vs,Rms = prop_luna_3b(bods,N,tstep,r0,v0,epoch)
elems = sv_to_elems(mb,rs,vs,N)
ts = time_arr(N,tstep)

Rm_S = np.zeros([N,3])
Rm_Smags = np.zeros([N])
b = 0
while b <= (N-1):
    Rm_S[b] = Rms[b] - rs[b]
    Rm_Smags[b] = np.linalg.norm(Rm_S[b]) - bods[1]['radius']
    b = b + 1
#'''

#Plots
elemfig, elemaxs = elem_plot(elems,ts,hours=False)
fig = plt.figure(figsize=(8,8))
ax_3d = plt.subplot(projection='3d')
tr_plot = Orbit_Plot(ax_3d,mb['radius'],rs,ttotal,N,tr_col='r',cb_col=mb['col'])
ax_3d.legend()


#Moon
Traj(ax_3d,Rms,ttotal,N,'0.4')
_u,_v = np.mgrid[0:2*np.pi:20j,0:np.pi:20j]
_x = sb['radius'] * np.cos(_u)*np.sin(_v) + Rms[-1,0]
_y = sb['radius'] * np.sin(_u)*np.sin(_v) + Rms[-1,1]
_z = sb['radius'] * np.cos(_v) + Rms[-1,2 ]
ax_3d.plot_wireframe(_x,_y,_z,color=sb['col'])

Rm_smagsfig, Rm_smagsax = plt.subplots()
Rm_smagsax.grid()
Rm_smagsax.set_title('Distance to Luna Surface')
Rm_smagsax.plot(ts/(3600*24), Rm_Smags, linewidth=2.0)
Rm_smagsax.set_xlabel('Time (Days)')

plt.show()