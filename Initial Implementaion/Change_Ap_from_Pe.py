from sys import path
path.append('./Tools')
import numpy as np
import Planetary_data as pd

#Import functions where needed
from sv_coe import *
#from geo_func import *
#from Kep_Anom import *
from plot_func import *
#from misc_maths import *

fig = plt.figure(figsize=(8,8))
ax = plt.subplot(projection='3d')

cb = pd.earth

#Orbit Initial conditions
cb = pd.earth
ra0 = cb['radius'] + 1500
rp0 = cb['radius'] + 300
a0 = (ra0+rp0)/2
emag0 = (ra0-rp0)/(ra0+rp0)
if emag0 < 1:
    hmag0 = np.sqrt((cb['mu']*a0)*(1-(emag0**2)))
else:
    hmag0 = np.sqrt((cb['mu']*(-a0))*(1-(emag0**2)))
i0 = np.deg2rad(28)
lon_an0 = np.deg2rad(45)
pearg0 = np.deg2rad(30)
ta0 = np.deg2rad(0.01)
T = ((2*np.pi)/np.sqrt(cb['mu']))*(a0**(3/2))

r0,v0 = sv_from_coe(cb['mu'],emag0,hmag0,i0,lon_an0,pearg0,ta0) #Return state vector from initial orbit elements

#Target orbit elements
target_ra = cb['radius'] + 3500 #Change desired Apogee altitude
target_a = (target_ra+rp0)/2
target_emag = (target_ra-rp0)/(target_ra+rp0)
if target_emag < 1:
    target_hmag = np.sqrt((cb['mu']*target_a)*(1-(target_emag**2)))
else:
    target_hmag = np.sqrt((cb['mu']*(-a0))*(1-(target_emag**2)))
tar_r0,tar_v0 = sv_from_coe(cb['mu'],target_emag,target_hmag,i0,lon_an0,pearg0,ta0)
T_target = orb_period(cb,target_a)
#//


'''
v0_unit, Use the inital velocity unit vector for the direction of the impulse
correct, Initial correction multipier is tiny to avoid overshooting
dv, The correction value is the magnitude of the velocity change
'''
v0_unit = (v0/np.linalg.norm(v0))   
correct = (cb['radius']/target_ra**2)*((target_ra-ra0)/ra0)   
dv = correct*v0_unit    

N = 0
while abs(correct) > 1e-8 and N < 100:
    a,_,hmag,_,_,_,emag,_,_ = coe_from_sv(cb['mu'],r0,v0+dv)
    rp = ((hmag**2)/cb['mu'])*(1/(1+emag))
    ra = ((hmag**2)/cb['mu'])*(1/(1-emag))
    if a < 0:
        T_inter = T_target #Avoid tryin to calculate the orbit period of intermediate hyperbolic trajectories
    else:
        T_inter = orb_period(cb,a)
    print(N,ra-cb['radius'],target_ra-ra,correct,dv*1000,np.linalg.norm(dv)*1000)
    
    step_div = 1.6535911**((hmag/rp)+1.35632997)
    Kep_traj_Prop_Plot(ax,cb,T_inter/2,r0,v0+dv,'g',step_div,legend=False) #Plot each intermediate targeting orbit up to Apogee

    correct = (cb['radius']/target_ra)*((target_ra-ra)/ra)
    N = N + 1
    dv = dv + correct*v0_unit

print('\n| Step | Current Apogee altitude | Delta to target apogee | dV Correction value | dV Vector(m/s) | dV Magnitude(m/s) |\n')

Kep_Prop_Plot(ax,cb,T,r0,v0,step_div=100,legend=False,show_cb=True) #Plot Initial orbit with the Earth and axis
#Kep_traj_Prop_Plot(ax,cb,T_inter/2,r0,v0+dv,'g',step_div,legend=False)  #Plot the final targeting orbit 
Kep_traj_Prop_Plot(ax,cb,T_target,tar_r0,tar_v0,step_div=1.6535911**((target_hmag/rp0)+1.35632997),legend=False,tr_col='b') #Plot the full target orbit
ax.set_title('Apogee Change dV Targeting')
plt.show()

