from sys import path
path.append('./Tools')
import numpy as np
import scipy as sp
import Planetary_data as pd
from spher_harm import *
from misc_maths import time_arr
from sv_coe import sv_from_coe
from plot_func import *


#Orbit Initial conditions
cb = pd.earth
ra = cb['radius'] + 3062
rp = cb['radius'] + 300
a = (ra+rp)/2
emag = (ra-rp)/(ra+rp)
if emag < 1:
    hmag = np.sqrt((cb['mu']*a)*(1-(emag**2)))
else:
    hmag = np.sqrt((cb['mu']*(-a))*(1-(emag**2)))
i = np.deg2rad(28)
lon_an = np.deg2rad(45)
pearg = np.deg2rad(30)
ta0 = np.deg2rad(40)

T = ((2*np.pi)/np.sqrt(cb['mu']))*(a**(3/2))

y0 = np.array([hmag,emag,ta0,lon_an,i,pearg])

def guass_dt(cb,rmag,hmag,emag,ta,lon_an,i,pearg):
    p = j2_rsw(cb,rmag,i,pearg,ta)

    u = pearg + ta

    dh_dt = rmag * p[1]

    de_dt = ((hmag/cb['mu']) * np.sin(ta) * p[0]) + ((1/(cb['mu']*hmag)) * ((((hmag**2) + (cb['mu']*rmag)) * np.cos(ta)) + (cb['mu']*emag*rmag)) * p[1])

    dtta_dt = (hmag/(rmag**2)) + ((1/(emag*hmag)) * ((((hmag**2)/cb['mu']) * np.cos(ta) * p[0]) - ((rmag + ((hmag**2)/cb['mu'])) * np.sin(ta) * p[1])))

    dtlonan_dt = (rmag/(hmag*np.sin(i))) * np.sin(u) * p[2]

    di_dt = (rmag/hmag) * np.cos(u) * p[2]

    dtpearg_dt = (-(1/(emag*hmag)) * ((((hmag**2)/cb['mu']) * np.cos(ta) * p[0]) - ((rmag + ((hmag**2)/cb['mu'])) * np.sin(ta) * p[1]))) - (((rmag*np.sin(u))/(hmag*np.tan(i))) * p[2])

    return np.array([dh_dt, de_dt, dtta_dt, dtlonan_dt, di_dt, dtpearg_dt])

def f(t,y):
    hmag,emag,ta,lon_an,i,pearg = y
    rmag = (hmag**2) / (cb['mu'] * (1 + (emag * np.cos(ta))))
    
    df_dt = guass_dt(cb,rmag,hmag,emag,ta,lon_an,i,pearg)
    return df_dt

#Set time period and step length/amount
ttotal = 2*24*(60**2)
tstep = ttotal/1000
N = int(np.ceil(ttotal/tstep) + 1)
print(tstep,'Seconds per',N,'Steps')
ts = time_arr(N,tstep)

hmags = np.zeros([N,1])
emags = np.zeros([N,1])
tas = np.zeros([N,1])
lon_ans = np.zeros([N,1])
incs = np.zeros([N,1])
peargs = np.zeros([N,1])
elems = np.zeros([N,6])

hmags[0] = hmag
emags[0] = emag
tas[0] = ta0
lon_ans[0] = lon_an
incs[0] = i
peargs[0] = pearg
elems[0] = y0

n = 1

prop = sp.integrate.ode(f,jac=None).set_integrator('dopri5')
while n < N:

    prop.set_initial_value(y0,0)
    prop.integrate(tstep)
    
    #hmags[n] = prop.y[0]
    #emags[n] = prop.y[1]
    #tas[n] = prop.y[2]
    lon_ans[n] = prop.y[3]
    #incs[n] = prop.y[4]
    peargs[n] = prop.y[5]
    
    elems[n] = prop.y

    y0 = prop.y
    n = n + 1

rs = np.zeros([N,3])
vs = np.zeros([N,3])
si = 0
while si < N:
    hmag = np.squeeze(elems[si,0])
    emag = np.squeeze(elems[si,1])
    ta = np.squeeze(elems[si,2])
    lon_an = np.squeeze(elems[si,3])
    i = np.squeeze(elems[si,4])
    pearg = np.squeeze(elems[si,5])
    rs[si],vs[si] = sv_from_coe(cb['mu'],emag,hmag,i,lon_an,pearg,ta)
    si = si + 1

lon_an_rate = np.rad2deg((lon_ans[-1]-lon_ans[0])/(ttotal/3600))
pearg_rate = np.rad2deg((peargs[-1]-peargs[0])/(ttotal/3600))
print('Lon_AN deg/hour =',lon_an_rate,'Cowell delta =',lon_an_rate - -0.17230583)
print('Pe_Arg deg/hour =',pearg_rate,'Cowell delta =',pearg_rate - 0.28216778)

#3d Plot
fig = plt.figure(figsize=(8,8))
ax = plt.subplot(projection='3d')
tr_plot = Orbit_Plot(ax,cb['radius'],rs,ttotal,N,tr_col='r',cb_col=cb['col'])
plt.show()