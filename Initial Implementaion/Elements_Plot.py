from sys import path
path.append('./Tools')
from Planetary_data import *

from int_props import *
from sv_coe import sv_from_coe,sv_to_elems
from misc_maths import time_arr
from plot_func import *

#Orbit Initial conditions
cb = pd.earth
ra = cb['radius'] + 150
rp = cb['radius'] + 80
a = (ra+rp)/2
e = (ra-rp)/(ra+rp)
if e < 1:
    hmag = np.sqrt((cb['mu']*a)*(1-(e**2)))
else:
    hmag = np.sqrt((cb['mu']*(-a))*(1-(e**2)))
i = np.deg2rad(90)
lon_an = np.deg2rad(0)
pearg = np.deg2rad(270)
ta0 = np.deg2rad(180)

T = ((2*np.pi)/np.sqrt(cb['mu']))*(a**(3/2))

r0,v0 = sv_from_coe(cb['mu'],e,hmag,i,lon_an,pearg,ta0)

#Set time period and step length/amount
ttotal = 2*24*(60**2)
tstep = ttotal/1000
N = int(np.ceil(ttotal/tstep) + 1)
ts = time_arr(N,tstep)

rs, vs = J2_Cowell(cb,N,tstep,r0,v0)
elems = sv_to_elems(cb,rs,vs,N)

#Plots
elemfig, elemaxs = elem_plot(elems,ts)
Orbfig, Orbaxs = plt.subplots(1, 1, subplot_kw=dict(projection='3d'), figsize=(8,8))
plt.subplots_adjust(left=0.05,right=0.95,bottom=0.06,top=0.930,wspace=0.2)
plteci = Orbit_Plot(Orbaxs,cb['radius'],rs,ttotal,N,tr_col='r',cb_col=cb['col'],wireframe=True)
Orbaxs.set_title((cb['name']+"-Centered Inertial Frame"))
plt.show()
