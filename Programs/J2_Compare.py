from sys import path
path.append('./Tools')
from Planetary_data import *

import numpy as np
from int_props import *
from sv_coe import sv_from_coe
from misc_maths import time_arr
from plot_func import *
from sv_coe import sv_to_elems

#Orbit Initial conditions
cb = pd.earth
ra = cb['radius'] + 3062
rp = cb['radius'] + 300
a = (ra+rp)/2
e = (ra-rp)/(ra+rp)
if e < 1:
    hmag = np.sqrt((cb['mu']*a)*(1-(e**2)))
else:
    hmag = np.sqrt((cb['mu']*(-a))*(1-(e**2)))
i = np.deg2rad(28)
lon_an = np.deg2rad(45)
pearg = np.deg2rad(30)
ta0 = np.deg2rad(40)

r0,v0 = sv_from_coe(cb['mu'],e,hmag,i,lon_an,pearg,ta0)

ttotal = 48*(3600)
tstep = ttotal/1000
#N = int(np.ceil(ttotal/tstep) + 1)

#Cowell Prop
tstep_cowell = tstep
N_cowell = int(np.ceil(ttotal/tstep_cowell) + 1)
ts_cowell = time_arr(N_cowell,tstep_cowell)
rs_cowell, vs_cowell = J2_Cowell(cb,N_cowell,tstep_cowell,r0,v0)
elems_cowell = sv_to_elems(cb,rs_cowell,vs_cowell,N_cowell)
lon_an_rate_cowell = np.rad2deg((elems_cowell[-1,3]-elems_cowell[0,3])/(ttotal/3600))
pearg_rate_cowell = np.rad2deg((elems_cowell[-1,5]-elems_cowell[0,5])/(ttotal/3600))


#Encke Prop
tstep_encke = tstep
#ddt = 1
N_encke = int(np.ceil(ttotal/tstep_encke) + 1)
ts_encke = time_arr(N_encke,tstep_encke)
rs_encke, vs_encke = J2_Encke(cb,N_encke,tstep_encke,r0,v0)
#rs_encke,vs_encke, N_encke = J2_Encke_variable_step(cb,tstep_encke,ddt,ttotal,r0,v0)
elems_encke = sv_to_elems(cb,rs_encke,vs_encke,N_encke)
lon_an_rate_encke = np.rad2deg((elems_encke[-1,3]-elems_encke[0,3])/(ttotal/3600))
pearg_rate_encke = np.rad2deg((elems_encke[-1,5]-elems_encke[0,5])/(ttotal/3600))
print('Encke Lon_AN deg/hour =',lon_an_rate_encke,', Cowell delta =',lon_an_rate_encke - lon_an_rate_cowell,',',((lon_an_rate_encke - lon_an_rate_cowell)/(lon_an_rate_cowell)*100),'%')
print('Encke Pe_Arg deg/hour =',pearg_rate_encke,', Cowell delta =',pearg_rate_encke - pearg_rate_cowell,',',((pearg_rate_encke - pearg_rate_cowell)/(pearg_rate_cowell)*100),'%')

def make_plots(axs):

    def axes_labels(ax):
        ax.set_xlabel('$km$')
        ax.set_ylabel('$km$')
        ax.set_zlabel('$km$')

    #Plot ECF
    ax1 = axs[0]
    plt_cowell = Orbit_Plot(ax1,cb['radius'],rs_cowell,ttotal,N_cowell,tr_col='r',cb_col=cb['col'],wireframe=True)
    ax1.set_title(('Cowell Propagation'))
    axes_labels(ax1)

    #Plot ECI
    ax2 = axs[1]
    plteci = Orbit_Plot(ax2,cb['radius'],rs_encke,ttotal,N_encke,tr_col='r',cb_col=cb['col'],wireframe=True)
    ax2.set_title(('Encke Propagation'))
    ax2.legend()
    axes_labels(ax2)
Orbfig, Orbaxs = plt.subplots(1, 2, subplot_kw=dict(projection='3d'), figsize=(12,6))
plt.subplots_adjust(left=0.05,right=0.95,bottom=0.06,top=0.930,wspace=0.2)
make_plots(Orbaxs)
plt.show()