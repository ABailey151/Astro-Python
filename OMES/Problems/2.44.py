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
r0 = np.array([3450,-1700,7750])
v0 = np.array([5.4,-5.4,1])
r0mag = np.linalg.norm(r0)
v0mag = np.linalg.norm(v0)
evec, e = eccen(u,r0,v0)
ta0 = np.rad2deg(np.arccos((np.dot(r0,evec))/(r0mag*e)))
dta = np.deg2rad(82)

R, V = rv_from_r0v0_ta(r0,v0,dta)
Rmag = np.linalg.norm(R)
Vmag = np.linalg.norm(V)
ta = np.rad2deg(np.arccos((np.dot(R,evec))/(Rmag*e)))

print("\nInitial Distance, Speed, and True Anomaly\nr0 =",r0mag,"km\nv0 =",v0mag,"km/s\nTA0 =",ta0,"Degrees")

print("\nA ---------------------------------------------")
print("Distance, Speed, and True Anomaly after 82 degree TA increase\nR =",Rmag,"km\nV =",Vmag,"km/s\nTA =",ta,"Degrees")

#B Complare h and sme values
h0 = np.cross(r0,v0)
sme0 = (v0mag**2)/2 - u/r0mag

h = np.cross(R,V)
sme = (Vmag**2)/2 - u/Rmag

print("\nB ---------------------------------------------")
print("h0 =",h0,"km^2/s\nh(dta82) =",h,"km^2/s\n")
print("sme0 =",sme0,"km^2/s^2\nsme(dta82) =",sme,"km^2/s^2\n")
print("Therefore, both the Specific Angular Momentum and Mechanical Energy are conserved\n")