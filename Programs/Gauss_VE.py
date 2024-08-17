from sys import path
path.append('./Tools')
import numpy as np
import Planetary_data as pd
from misc_maths import time_arr
from plot_func import *
from int_props import *

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

#Set time period and step length/amount
ttotal = 2*24*(60**2)
tstep = ttotal/1000
N = int(np.ceil(ttotal/tstep) + 1)
ts = time_arr(N,tstep)

#Propagate Elements
elems,rs,vs = gauss_ve_prop(cb,N,tstep,hmag,emag,ta0,lon_an,i,pearg)

lon_an_rate = np.rad2deg((elems[-1,3]-elems[0,3])/(ttotal/3600))
pearg_rate = np.rad2deg((elems[-1,5]-elems[0,5])/(ttotal/3600))
print('Lon_AN deg/hour =',lon_an_rate)
print('Pe_Arg deg/hour =',pearg_rate)

#Plots
elemfig,elemaxs = elem_plot(elems,ts)

fig = plt.figure(figsize=(8,8))
ax = plt.subplot(projection='3d')
gauss_plot = Orbit_Plot(ax,cb['radius'],rs,ttotal,N,tr_col='r',cb_col=cb['col'])
ax.legend()

plt.show()