from sys import path
path.append('./Tools')
import numpy as np
u = 398600.4

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
    #print(rmag,f,g,fdot,gdot)
    R = f*r0 + g*v0
    V = fdot*r0 + gdot*v0
    return R,V

#Define initial position and velocity conditions
r0 = np.array([8182.4,-6865.9,0])
v0 = np.array([0.47572,8.8116,0])
dta = 120*((2*np.pi)/360)

R, V = rv_from_r0v0_ta(r0,v0,dta)

print("R =",R,"km\nV =",V,"km/s")