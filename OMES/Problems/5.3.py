from sys import path
path.append('./Tools')
import numpy as np
import Planetary_data as pd

#Import functions where needed
from sv_coe import *
#from geo_func import *
#from Kep_Anom import *
#from plot_func import *
#from misc_maths import *
from Lambert_Functions import *
from plot_func import *

cb = pd.earth
mu = cb['mu']
re = cb['radius']

r1mag = re + 400
r2mag = re + 1000

dt = 30*60
del_ta = np.deg2rad(120)

A = np.sin(del_ta) * np.sqrt((r1mag*r2mag)/(1-np.cos(del_ta)))

#//Find z
z0 = 0
Fz0 = lambert_F_z(r1mag,r2mag,mu,dt,A,z0)
z = lambert_z(r1mag,r2mag,mu,dt,A,z0,Fz0)

#Create position vectors in an arbitrary frame with the x axis on r1
r1 = np.array([r1mag, 0, 0])
r2 = np.array([r2mag*np.cos(del_ta), r2mag*np.sin(del_ta), 0])

#Calculate the velocities
v1,v2 = lambert_vs_from_rs(r1,r2,r1mag,r2mag,z,A,mu)

h,hmag = ang_mom(r1,v1)
e,emag = eccen(mu,r1,v1)

rp = ((hmag**2)/mu)*(1/(1+emag))
ra = ((hmag**2)/mu)*(1/(1-emag))
a = 0.5*(rp+ra)
T = ((2*np.pi)/np.sqrt(cb['mu']))*(a**(3/2))
print('\nPerigee altitude =',rp-re,'km\n')

fig = plt.figure(figsize=(8,8))
ax = plt.subplot(projection='3d')
Kep_Prop_Plot(ax,cb,T,r1,v1)
plt.show()
