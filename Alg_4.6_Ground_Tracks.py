from sys import path
path.append('./Tools')
import numpy as np
import matplotlib.pyplot as plt
import Planetary_data as pd

#Import functions where needed
from sv_coe import sv_from_coe
#from geo_func import *
from Kep_Anom import kepler_E
#from plot_func import *
from misc_maths import ra_and_dec_from_r
from frame_transform import int_to_fix

cb = pd.earth


#Orbit Initial conditions
ra = 10000
rp = 6700
a = 26500#(ra+rp)/2
emag = 0.75#(ra-rp)/(ra+rp)
if emag < 1:
    hmag = np.sqrt((cb['mu']*a)*(1-(emag**2)))
else:
    hmag = np.sqrt((cb['mu']*(-a))*(1-(emag**2)))
i = np.deg2rad(63.4)
rightasc = np.deg2rad(270)
pearg = np.deg2rad(270)
ta0 = np.deg2rad(230)
T = ((2*np.pi)/np.sqrt(cb['mu']))*(a**(3/2))
'''
#Orbit Initial conditions
ra = 10000
rp = 6700
a = 42164#(ra+rp)/2
emag = 0.3#(ra-rp)/(ra+rp)
if emag < 1:
    hmag = np.sqrt((cb['mu']*a)*(1-(emag**2)))
else:
    hmag = np.sqrt((cb['mu']*(-a))*(1-(emag**2)))
i = np.deg2rad(63.4)
rightasc = np.deg2rad(270)
pearg = np.deg2rad(270)
ta0 = np.deg2rad(230)
T = ((2*np.pi)/np.sqrt(cb['mu']))*(a**(3/2))
'''
dt = 1
tstop = 2*T
N = int(np.ceil(tstop/dt))

r_RAs = np.zeros([N+1])
r_decs = np.zeros([N+1])

phi0 = 0

#Step 1: Find j2 perturbations on righasc and pearg. \/both in rad/s
d_RightAsc = -((3/2) * ((np.sqrt(cb['mu'])*cb['j2']*(cb['radius']**2)) / (((1-(emag**2))**2)*(a**(7/2))))) * np.cos(i)
d_pearg = -((3/2) * ((np.sqrt(cb['mu'])*cb['j2']*(cb['radius']**2)) / (((1-(emag**2))**2)*(a**(7/2))))) * (((5/2)*(np.sin(i)**2))-2)

#Step 2: Find t0 since periapsis passage\
#a. Find Eccentric Anomaly
E = 2 * np.arctan(np.tan(ta0/2)*np.sqrt((1-emag)/(1+emag)))
#b. Find the Mean Anomaly
M = E - emag*np.sin(E)
#c. Find t0
t0 = (M*T)/(2*np.pi)
ta = ta0
t = t0
phi = phi0
n = 0

while n <= N:
    r,v = sv_from_coe(cb['mu'],emag,hmag,i,rightasc,pearg,ta)
    r_fix = int_to_fix(phi,r)
    ra,dec = ra_and_dec_from_r(r_fix)
    r_RAs[n] = ra * (180/np.pi)
    r_decs[n] = dec * (180/np.pi)

    #Step 3: At time t = t0 + dt, calculate r_righasc and r_dec
    #a. Calculate the true anomaly
    Mt = (2*np.pi*t)/T
    E = kepler_E(emag,Mt)
    ta = 2 * np.arctan(np.tan(E/2)*np.sqrt((1+emag)/(1-emag)))

    rightasc_t = rightasc + (d_RightAsc*dt)
    pearg_t = pearg + (d_pearg*dt)
    t = t + dt
    phi = phi + cb['dphi']*dt
    n = n + 1

Rm_smagsfig, Rm_smagsax = plt.subplots()
Rm_smagsax.grid()
Rm_smagsax.plot(r_RAs, r_decs, linewidth=2.0)
Rm_smagsax.set_xlim(0, 360)  # X-axis range from 0 to 5
Rm_smagsax.set_ylim(-90,90)
plt.show()