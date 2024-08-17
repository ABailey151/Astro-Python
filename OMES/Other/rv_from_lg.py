from sys import path
path.append('./Tools')
from sv_coe import *
from kep_Anom import *
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
    print(rmag,f,g,fdot,gdot)
    R = f*r0 + g*v0
    V = fdot*r0 + gdot*v0
    return R,V

#Calculate Orbit elements from state vectors
r0 = np.array([8182.4, -6865.9, 0])
v0 = np.array([0.47572, 8.8116, 0])
r0mag = np.linalg.norm(r0)
v0mag = np.linalg.norm(v0)
H, hmag = ang_mom(r0,v0)

evec, e = eccen(u,r=r0,v=v0)
sme = ((v0mag**2)/2) - (u/r0mag)
a = -(u/(2*sme))
T = ((2*np.pi)/np.sqrt(u))*(a**(3/2))

#Mean Anomaly
dt = 60*60                  #This method only holds for time intervals less that 10 minutes

if e < 1:
    M = (2*np.pi)*(dt/T)
    E, dta = ta_from_E(e,M)

else:
    M = ((u**2)/(hmag**3))*(((e**2)-1)**(3/2))*(dt)
    E, dta = ta_from_F(e,M)

dtadeg = np.rad2deg(dta)
R, V = rv_from_r0v0_ta(r0,v0,dta)
print(np.linalg.norm(R))


print("\n-------------------------------------------------------------------------")
print("Initial Position Vector =",r0,"km")
print("\nA ------------------------------------------------------------------------------------")
print("The Position Vector after 10 minutes =",R,"km")
print("\nB -----------------------------------------------------")
print("The change in true anomaly =",dtadeg,"degrees\n")
