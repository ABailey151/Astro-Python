from sys import path
path.append('./Tools')
import numpy as np
import Planetary_data as pd
from sv_coe import *
from Kep_Anom import *

u = pd.earth['mu']

r0 = np.array([3231.6687, -5597.41438, 2352.46243])
r0mag = np.linalg.norm(r0)
v0 = np.array([0, 0, 10])
v0mag = np.linalg.norm(v0)
#Calculate orbital elements
evec, e = eccen(u,r0,v0)
H = np.cross(r0,v0)
h = np.linalg.norm(H)
k = np.array([0, 0, 1])
i = np.rad2deg(np.arccos(np.dot(H,k)/(h)))
if e < 1:
    a = ((h**2)/u)*(1/(1-e**2))
else:
    a = -((h**2)/u)*(1/((e**2)-1))

if a < 0:
    T = "N/A"
else:
    T = ((2*np.pi)/np.sqrt(u))*(a**(3/2))

#Calculate the initial True anomaly
cos𝜃0 = (h**2)/(e*r0mag*u) - (1/e)
if cosθ0 > 1:
    cosθ0 = 1
elif cosθ0 < -1:
    cosθ0 = -1

tat0 = np.sign(np.dot(H,np.cross(evec,r0)))
if tat0 >= 0:
    𝜃0 = np.rad2deg(np.arccos(cosθ0))
else:
    𝜃0 = 360 - np.rad2deg(np.arccos(cosθ0))

#Additional variable for Universal Anomaly algorithm
vr0 = np.dot(r0,v0)/r0mag
alpha = 2/r0mag - (v0mag**2)/u
dt = 30*60
if e < 1:
    if dt >= T:
        n = dt/T

#Use UA and fg/fdgd algorithms to find final R and V vectors
X = kepler_U(u,dt,r0mag,vr0,alpha)

'''
f, g = fg_from_ua(dt,r0mag,alpha,X)
R = f*r0 + g*v0
Rmag = np.linalg.norm(R)

fdot, gdot = fdgd_from_ua(r0mag,Rmag,alpha,X)
V = fdot*r0 + gdot*v0
Vmag = np.linalg.norm(V)
'''
R, V = rv_from_kepler_U(u,dt,r0,v0)
Rmag = np.linalg.norm(R)
Vmag = np.linalg.norm(V)

#Calculate True Anomaly after time interval
cos𝜃 = (h**2)/(e*Rmag*u) - (1/e)
if cos𝜃 > 1:
    cos𝜃 = 1
elif cosθ < -1:
    cosθ = -1

tat = np.sign(np.dot(H,np.cross(evec,R)))
if tat >= 0:
    𝜃 = np.rad2deg(np.arccos(cosθ))
else:
    𝜃 = 360 - np.rad2deg(np.arccos(cosθ))

#Print Section
print("\n---------------------------------------------------------------------------")
print("Initial Position Vector and Magnitude,\nr0 =",r0[0],"i +",r0[1],"j +",r0[2],"k\nr0mag =",r0mag,"km\n")
print("Inital Velocity Vector and Magnitude,\nv0 =",v0[0],"i +",v0[1],"j +",v0[2],"k\nv0mag =",v0mag,"km/s\n")
print("e =",e)
print("a =",a,"km")
print("i =",i,"Degrees")
if e < 1:
    print("T =",T,"Seconds,(",T/60,"Minutes)")
print("Initial True Anomaly =",np.deg2rad(𝜃0),"RADs /",𝜃0,"Degress")
print("---------------------------------------------------------------------------")
print("dt =",dt,"Seconds,(",dt/60,"Minutes)")
print("Universal Anomaly =",X,"km^1/2")
print("---------------------------------------------------------------------------")
print("New Position Vector and Magnitude,\nR =",R[0],"i +",R[1],"j +",R[2],"k\nRmag =",Rmag,"km\n")
print("New Velocity Vector and Magnitude,\nV =",V[0],"i +",V[1],"j +",V[2],"k\nVmag =",Vmag,"km\n")
print("New True Anomaly =",np.deg2rad(𝜃),"RADs /",𝜃,"Degress")
if e < 1:
    if dt >= T:
        print("The Spacecraft completed",n,"revolution(s)")
print("---------------------------------------------------------------------------\n")

