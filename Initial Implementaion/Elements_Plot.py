from sys import path
path.append('./Tools')
from Planetary_data import *

from int_props import *
from sv_coe import sv_from_coe,sv_to_elems
from misc_maths import time_arr
from plot_func import *

#Orbit Initial conditions
cb = pd.earth
ra = cb['radius'] + 30620
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

r0,v0 = sv_from_coe(cb['mu'],emag,hmag,i,lon_an,pearg,ta0)

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
