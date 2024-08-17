from sys import path
path.append('./Tools')
import numpy as np
import Planetary_data as pd
from sv_coe import *
from kep_Anom import ta_from_F

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

#A Calculate Orbit elements from state vectors
r0 = np.array([6320,0,7750])
v0 = np.array([0,11,0])
r0mag = np.linalg.norm(r0)
v0mag = np.linalg.norm(v0)
H, hmag = ang_mom(r0,v0)
evec, e = eccen(u,r0,v0)

#Unused since SC starts a perigee
#costa0 = np.dot(r0,evec)/(e*r0mag)
#ta0 = np.arccos(costa0)

#Calculate Mean anomaly from time interval
dt = 60*10
M = ((u**2)/(hmag**3))*(((e**2-1))**(3/2))*(dt)

#Use kepler_F to find Hyperbolic Eccentric Anomaly and new True Anomaly
F, dta = ta_from_F(e,M) #- ta0
dtadeg = np.rad2deg(dta)

#From observation, ta = dta since costa0 = 1 (r0=rp)
#print(costa0,ta0)

R, V = rv_from_r0v0_ta(r0,v0,dta)


print("\n-------------------------------------------------------------------------")
print("Initial Position Vector =",r0,"km, at 0 degrees True Anomaly")
print("\nA ------------------------------------------------------------------------------------")
print("The Position Vector after 10 minutes =",R,"km")
print("\nB -----------------------------------------------------")
print("The change in true anomaly =",dtadeg,"degrees\n")