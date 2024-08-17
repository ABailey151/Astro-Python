from sys import path
path.append('./Tools')
import numpy as np
import Planetary_data as pd
from sv_coe import *

def rv_from_r0v0_ta(r0,v0,dta):
    #Magnitudes of r and v
    r0mag = np.linalg.norm(r0)
    v0mag = np.linalg.norm(v0)

    #radial component of v0 by prjecting onto r0
    vr0 = np.dot(r0,v0) / r0mag
    
    #Angular momentum magnitude
    hmag = r0mag * np.sqrt((v0mag**2)-(vr0**2))

    #calculate new postition magnitude
    d = 1 + ((((hmag**2)/(u*r0mag))-1)*np.cos(dta)) - (((hmag*vr0)/u)*np.sin(dta))
    rmag = ((hmag**2)/u) * (1/d)

    #Calculate lagrange coefficients
    f = 1 - (((u*rmag)/(hmag**2))*(1-np.cos(dta)))
    g = ((rmag*r0mag)/hmag)*np.sin(dta)

    b = (u/hmag) * ((1-np.cos(dta))/np.sin(dta))
    fdot = b * (((u/(hmag**2))*(1-np.cos(dta))) - (1/r0mag) - (1/rmag))
    gdot = 1 - ((u*r0mag)/(hmag**2))*(1-np.cos(dta))

    #Calculate new r and v
    R = f*r0 + g*v0
    V = fdot*r0 + gdot*v0
    return R,V

cb = pd.earth
u = cb['mu']

#A Use Algorithm 2.3
r0 = np.array([7000,0,0])
v0 = np.array([7,7,0])
r0mag = np.linalg.norm(r0)
dta = 90*((2*np.pi)/360)
evec, e = eccen(u,r0,v0)

R, V = rv_from_r0v0_ta(r0,v0,dta)
Rmag = np.linalg.norm(R)

print("\nA ---------------------------------------------------------------------------------------")
print("Position after 90 degree ta increase =",R,"km")

#B Find the initial true anomaly
ta0rad = np.arccos(np.dot(r0,evec)/(r0mag*e))
ta0 = np.rad2deg(ta0rad)

print("\nB -----------------------------------------------")
print("Initial True Anomaly =",ta0,"Degrees\n")