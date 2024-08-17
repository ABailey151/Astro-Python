from sys import path
path.append('./Tools')
import numpy as np
import scipy as sp
import Planetary_data as pd

#Import functions where needed
from sv_coe import *
#from geo_func import *
#from Kep_Anom import *
from plot_func import *
from misc_maths import time_arr

def atmosphere(z):

    #Geometric altitudes (km)
    h = np.array([0, 25, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 180, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000])

    #Corresponding densities (kg/m^3)
    r = np.array([1.225, 4.008e-2, 1.841e-2, 3.996e-3, 1.027e-3, 3.097e-4, 8.283e-5, 1.846e-5, 3.416e-6, 5.606e-7, 9.708e-8, 2.222e-8, 8.152e-9, 3.831e-9, 2.076e-9, 5.194e-10, 2.541e-10, 6.073e-11, 1.916e-11, 7.014e-12, 2.803e-12, 1.184e-12, 5.215e-13, 1.137e-13, 3.070e-14, 1.136e-14, 5.759e-15, 3.561e-15])

    #Scale heights (H)
    H = np.array([7.310, 6.427, 6.546, 7.360, 8.342, 7.583, 6.661, 5.927, 5.533, 5.703, 6.782, 9.973, 13.243, 16.322, 21.652, 27.974, 34.934, 43.342, 49.755, 54.513, 58.019, 60.980, 65.654, 76.377, 100.587, 147.203, 208.020])

    #Handle altitudes outside of the range
    if z > 1000:
        z = 1000
    elif z < 0:
        z = 0

    #Determine the interpolation interval
    for j in range(27):

        if z >= h[j] and z < h[j+1]:
            i = j

    if z == 1000:
        i = 26

    density = r[i] * np.exp(-(z - h[i])/H[i])

    return density

#Central Body
cb = pd.earth
Re = cb['radius']
we = np.array([0, 0, cb['dphi']])

#Orbit Initial conditions
ra = cb['radius'] + 250
rp = cb['radius'] + 215
a = (ra+rp)/2
e = (ra-rp)/(ra+rp)
if e < 1:
    hmag = np.sqrt((cb['mu']*a)*(1-(e**2)))
else:
    hmag = np.sqrt((cb['mu']*(-a))*(1-(e**2)))
i = np.deg2rad(65.1)
lon_an = np.deg2rad(340)
pearg = np.deg2rad(58)
ta0 = np.deg2rad(332)
T = ((2*np.pi)/np.sqrt(cb['mu']))*(a**(3/2))

r0,v0 = sv_from_coe(cb['mu'],e,hmag,i,lon_an,pearg,ta0)

y0 = np.concatenate([r0,v0])

#Spacecraft properties
Cd = 2.2
d = 1
A = np.pi * (d**2)
m = 500

def f(t,y):

    r = y[0:3]
    rmag = np.linalg.norm(r)
    v = y[3:6]

    rho = atmosphere(rmag-Re)
    vrel = v - np.cross(we,r)
    vrelmag = np.linalg.norm(vrel)

    p_atmos = -(1/2) * rho * (vrelmag*1000)**2 * ((Cd*A)/m) * (vrel/vrelmag)

    ag = -(cb['mu']/(np.linalg.norm(r)**3)) * r

    a = ag + (p_atmos/1000)

    dydt = np.concatenate([v,a])
    return dydt

#Set time period and step length/amount
ttotal = 2*24*(3600)
tstep = 100
N = int(np.ceil(ttotal/tstep))
#print(N,'Steps')
n = 1

#Initialise arrays and write initial conditions

rs = np.zeros([N,3])
vs = np.zeros([N,3])
ras = np.zeros([N])
rps = np.zeros([N])

rs = r0
vs = v0
ras = ra - Re
rps = rp - Re

prop = sp.integrate.ode(f,jac=None).set_integrator('dopri5')
prop.set_initial_value(y0,0)
while np.linalg.norm(prop.y[0:3]) > (Re + 100): #Propagate until the spacecraft falls to 100km altittude

    prop.set_initial_value(y0,0)
    prop.integrate(tstep)
    
    '''
    rs[n] = prop.y[0:3]
    vs[n] = prop.y[3:6]

    _,hmagi = ang_mom(prop.y[0:3],prop.y[3:6])
    _,emagi = eccen(cb['mu'],prop.y[0:3],prop.y[3:6])
    ras[n] = ((hmagi**2)/cb['mu'])*(1/(1-emagi)) - Re
    rps[n] = ((hmagi**2)/cb['mu'])*(1/(1+emagi)) - Re
    '''
    rs = np.vstack([rs,prop.y[0:3]])
    vs = np.vstack([vs,prop.y[3:6]])
    _,hmagi = ang_mom(prop.y[0:3],prop.y[3:6])
    _,emagi = eccen(cb['mu'],prop.y[0:3],prop.y[3:6])
    ras = np.vstack([ras,((hmagi**2)/cb['mu'])*(1/(1-emagi)) - Re])
    rps = np.vstack([rps,((hmagi**2)/cb['mu'])*(1/(1+emagi)) - Re])

    y0 = prop.y
    n = n + 1

ts = time_arr(n,tstep)

Elems = plt.figure(figsize=(8,6))
ax2d1 = Elems.add_subplot(111)
ax2d1.grid()
ax2d1.plot((ts/(24*(60**2))), ras,'r',label='Apogee')
plt.title(('Orbit decay of '+str(m)+'kg Spacecraft'))
plt.xlabel('Time (Days)')
plt.ylabel('Altitude (km)')
#ax2d2 = Elems.add_subplot(212,sharex=ax2d1)
ax2d1.plot((ts/(24*(60**2))), rps,'b',label='Perigee')
plt.subplots_adjust(left=0.1,right=0.95,bottom=0.1,top=0.94,hspace=0.25)
ax2d1.legend()

orb3d = plt.figure(figsize=(8,8))
ax3d = plt.subplot(projection='3d')
plt.subplots_adjust(left=0.05,right=0.95,bottom=0.1,top=0.94)
tr_plot = Orbit_Plot(ax3d,cb['radius'],rs,int(ts[-1]),n,tr_col='r',cb_col=cb['col'])
ax3d.set_title(str(cb['name'])+' Centered Inertial Frame')
ax3d.legend()
plt.show()
